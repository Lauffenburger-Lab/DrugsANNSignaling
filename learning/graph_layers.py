"""
Taken from https://github.com/NickMeim/deepSIBA_pytorch.
"""
from __future__ import absolute_import,division, print_function
import numpy as np
from numpy import inf, ndarray
import pandas as pd
import torch
import torch.nn as nn
from torch.functional import F
from tqdm.auto import tqdm
import os
import random
import sklearn
from sklearn.model_selection import KFold
from sklearn import metrics
import re
from tempfile import TemporaryFile
from rdkit import Chem
from functools import partial
from multiprocessing import cpu_count, Pool
from copy import deepcopy
from functools import reduce

def get_default_device():
    if torch.cuda.is_available():
        #print('cuda mode')
        return torch.device('cuda')
    else:
        #print('cpu mode')
        return torch.device('cpu')

def temporal_padding(x, padding=(1, 1), padvalue=0):
    """Pads the middle dimension of a 3D tensor.
    Arguments:
        x: Tensor or variable.
        padding: Tuple of 2 integers, how many zeros to
            add at the start and end of dim 1.
    Returns:
        A padded 3D tensor.
    """
    assert len(padding) == 2
    pattern = (0, 0, padding[0], padding[1], 0, 0)
    return F.pad(x, pattern, "constant",padvalue)


def neighbour_lookup(atoms, edges,maskvalue=0,include_self=False):
    ''' Looks up the features of an all atoms neighbours, for a batch of molecules.
    # Arguments:
        atoms (K.tensor): of shape (batch_n, max_atoms, num_atom_features)
        edges (K.tensor): of shape (batch_n, max_atoms, max_degree) with neighbour
            indices and -1 as padding value
        maskvalue (numerical): the maskingvalue that should be used for empty atoms
            or atoms that have no neighbours (does not affect the input maskvalue
            which should always be -1!)
        include_self (bool): if True, the featurevector of each atom will be added
            to the list feature vectors of its neighbours
    # Returns:
        neigbour_features (K.tensor): of shape (batch_n, max_atoms(+1), max_degree,
            num_atom_features) depending on the value of include_self
    # Todo:
        - make this function compatible with Tensorflow, it should be quite trivial
            because there is an equivalent of `T.arange` in tensorflow.
    '''

    device = get_default_device()

    # The lookup masking trick: We add 1 to all indices, converting the
    #   masking value of -1 to a valid 0 index.
    masked_edges = edges + 1
    masked_edges = masked_edges.detach()
    # We then add a padding vector at index 0 by padding to the left of the
    #   lookup matrix with the value that the new mask should get
    masked_atoms = temporal_padding(atoms, (1, 0), padvalue=maskvalue)

    # Import dimensions
    atoms_shape = masked_atoms.shape
    batch_n = atoms_shape[0]
    lookup_size = atoms_shape[1]
    num_atom_features = atoms_shape[2]

    edges_shape = masked_edges.shape
    max_atoms = edges_shape[1]
    max_degree = edges_shape[2]

    # create broadcastable offset
    offset_shape = (batch_n, 1, 1)
    offset = torch.arange(start=0,end=batch_n).reshape(offset_shape)
    offset *= lookup_size
    offset = offset.to(device)
    offset = offset.detach()

    # apply offset to account for the fact that after reshape, all individual
    #   batch_n indices will be combined into a single big index
    flattened_atoms = torch.reshape(masked_atoms,(-1,num_atom_features))
    flattened_edges = torch.reshape(torch.add(masked_edges, offset), (batch_n, -1)).detach()

    # Gather flattened
    flattened_result=flattened_atoms[flattened_edges]

    # Unflatten result
    output_shape = (batch_n, max_atoms, max_degree, num_atom_features)
    output = torch.reshape(flattened_result, output_shape)

    if include_self:
        return torch.cat([torch.unsqueeze(atoms, dim=2), output], dim=2)
    return output


class NeuralGraphHidden(nn.Module):
    ''' Hidden Convolutional layer in a Neural Graph (as in Duvenaud et. al.,
    2015). This layer takes a graph as an input. The graph is represented as by
    three tensors.

    - The atoms tensor represents the features of the nodes.
    - The bonds tensor represents the features of the edges.
    - The edges tensor represents the connectivity (which atoms are connected to
        which)

    It returns the convolved features tensor, which is very similar to the atoms
    tensor. Instead of each node being represented by a num_atom_features-sized
    vector, each node now is represented by a convolved feature vector of size
    conv_width.

    # Example
        Define the input:
        ```python
            atoms0 = Input(name='atom_inputs', shape=(max_atoms, num_atom_features))
            bonds = Input(name='bond_inputs', shape=(max_atoms, max_degree, num_bond_features))
            edges = Input(name='edge_inputs', shape=(max_atoms, max_degree), dtype='int32')
        ```

        The `NeuralGraphHidden` can be initialised in three ways:
        1. Using an integer `conv_width` and possible kwags (`Dense` layer is used)
            ```python
            atoms1 = NeuralGraphHidden(conv_width, activation='relu', bias=False)([atoms0, bonds, edges])
            ```
        2. Using an initialised `Dense` layer
            ```python
            atoms1 = NeuralGraphHidden(Dense(conv_width, activation='relu', bias=False))([atoms0, bonds, edges])
            ```
        3. Using a function that returns an initialised `Dense` layer
            ```python
            atoms1 = NeuralGraphHidden(lambda: Dense(conv_width, activation='relu', bias=False))([atoms0, bonds, edges])
            ```

        Use `NeuralGraphOutput` to convert atom layer to fingerprint

    # Arguments
        inner_layer_arg: Either:
            1. an int defining the `conv_width`, with optional kwargs for the
                inner Dense layer
            2. An initialised but not build (`Dense`) keras layer (like a wrapper)
            3. A function that returns an initialised keras layer.
        kwargs: For initialisation 1. you can pass `Dense` layer kwargs

    # Input shape
        List of Atom and edge tensors of shape:
        `[(samples, max_atoms, atom_features), (samples, max_atoms, max_degrees,
          bond_features), (samples, max_atoms, max_degrees)]`
        where degrees referes to number of neighbours

    # Output shape
        New atom featuers of shape
        `(samples, max_atoms, conv_width)`

    # References
        - [Convolutional Networks on Graphs for Learning Molecular Fingerprints](https://arxiv.org/abs/1509.09292)

    '''

    def __init__(self, num_at_feats,num_bond_feats, hidden_size,max_degree,activ, bias , **kwargs):
        # Initialise based on one of the three initialisation methods
        super(NeuralGraphHidden, self).__init__()

        # inner_layer_arg is conv_width
        if (isinstance(num_at_feats,int) and isinstance(num_bond_feats,int) and isinstance(hidden_size, int) and isinstance(max_degree,int)):
          self.input_size = int(num_at_feats+num_bond_feats)
          self.hidden_size = hidden_size
          self.activ=activ
          #self.create_inner_layer_fn = nn.Linear(self.input_size, self.hidden_size, bias=bias)
          if activ is not None:
            self.relu= nn.ReLU()
        else:
            raise ValueError('NeuralGraphHidden has to be initialised with 4 integers for: num_at_feats, num_bond_feats,hidden size,and max_degree')

        self.max_degree = max_degree
        # Add the dense layers (that contain trainable params)
        #   (for each degree we convolve with a different weight matrix)
        self.inner_3D_layers = nn.ModuleList([nn.Linear(self.input_size, self.hidden_size, bias=bias) for degree in range(self.max_degree)])

        self.init_weights(NeuralGraphHidden)

    def init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.xavier_normal_(m.weight)
            nn.init.zeros_(m.bias)

    def forward(self, atoms, bonds, edges,mask=None):
        #atoms, bonds, edges = x

        # Import dimensions
        num_samples = list(atoms.size())[0]
        max_atoms = list(atoms.size())[1]
        num_atom_features = list(atoms.size())[-1]
        num_bond_features = list(bonds.size())[-1]

        # Create a matrix that stores for each atom, the degree it is
        #atom_degrees = torch.sum((edges.eq(1)), dim=-1, keepdim=True) - 1
        atom_degrees = torch.sum((~edges.eq(-1)), dim=-1, keepdim=True)
        atom_degrees = atom_degrees.detach()

        # For each atom, look up the features of it's neighbour
        neighbour_atom_features = neighbour_lookup(atoms, edges, include_self=True)

        # Sum along degree axis to get summed neighbour features
        summed_atom_features = torch.sum(neighbour_atom_features, dim=-2).type(torch.cuda.FloatTensor)
        #summed_atom_features = torch.bmm(edges,atoms)

        # Sum the edge features for each atom
        summed_bond_features = torch.sum(bonds, dim=-2).type(torch.cuda.FloatTensor)
        #summed_bond_features = torch.bmm(edges,torch.sum(bonds, dim=-2))

        # Concatenate the summed atom and bond features
        summed_features = torch.cat([summed_atom_features, summed_bond_features], dim=-1)

        # For each degree we convolve with a different weight matrix
        new_features_by_degree = []
        for degree in range(self.max_degree):

            # Create mask for this degree
            atom_masks_this_degree = atom_degrees.eq(degree+1)

            # Multiply with hidden merge layer
            #   (use time Distributed because we are dealing with 2D input/3D for batches)
            #print(summed_features.shape)
            #print(self.inner_3D_layers[degree].weight.shape)
            new_unmasked_features = self.inner_3D_layers[degree](summed_features)
            if self.activ is not None:
              new_unmasked_features=self.relu(new_unmasked_features)

            # Do explicit masking because TimeDistributed does not support masking
            new_masked_features = new_unmasked_features * atom_masks_this_degree

            new_features_by_degree.append(new_masked_features)

        # Finally sum the features of all atoms
        new_features = torch.stack(new_features_by_degree,dim=0).sum(dim=0)
        #new_features = reduce(torch.add, new_features_by_degree)

        return new_features
