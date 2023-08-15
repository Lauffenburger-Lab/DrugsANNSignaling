"""
Taken from https://github.com/BioSysLab/deepSIBA/tree/master.
"""
from __future__ import division, print_function
import numpy as np
from numpy import inf, ndarray
import pandas as pd
import os
import random
import re
from tempfile import TemporaryFile
from rdkit import Chem
from functools import partial
from multiprocessing import cpu_count, Pool
from copy import deepcopy

def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return list(map(lambda s: x == s, allowable_set))

def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))

def atom_features(atom):
    return np.array(one_of_k_encoding_unk(atom.GetSymbol(),
                                      ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na',
                                       'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb',
                                       'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H',    # H?
                                       'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr',
                                       'Cr', 'Pt', 'Hg', 'Pb', 'Unknown']) +
                    one_of_k_encoding_unk(atom.GetDegree(), [0, 1, 2, 3, 4, 5]) +
                    one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4]) +
                    one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5]) +
                    [atom.GetIsAromatic()])

def bond_features(bond):
    bt = bond.GetBondType()
    return np.array([bt == Chem.rdchem.BondType.SINGLE,
                     bt == Chem.rdchem.BondType.DOUBLE,
                     bt == Chem.rdchem.BondType.TRIPLE,
                     bt == Chem.rdchem.BondType.AROMATIC,
                     bond.GetIsConjugated(),
                     bond.IsInRing()])

def num_atom_features():
    # Return length of feature vector using a very simple molecule.
    m = Chem.MolFromSmiles('CC')
    alist = m.GetAtoms()
    a = alist[0]
    return len(atom_features(a))

def num_bond_features():
    # Return length of feature vector using a very simple molecule.
    simple_mol = Chem.MolFromSmiles('CC')
    Chem.SanitizeMol(simple_mol)
    return len(bond_features(simple_mol.GetBonds()[0]))



def padaxis(array, new_size, axis, pad_value=0, pad_right=True):
    ''' Padds one axis of an array to a new size
    This is just a wrapper for np.pad, more usefull when only padding a single axis
    # Arguments:
        array: the array to pad
        new_size: the new size of the specified axis
        axis: axis along which to pad
        pad_value: pad value,
        pad_right: boolean, pad on the right or left side
    # Returns:
        padded_array: np.array
    '''
    add_size = new_size - array.shape[axis]
    assert add_size >= 0, 'Cannot pad dimension {0} of size {1} to smaller size {2}'.format(axis, array.shape[axis], new_size)
    pad_width = [(0,0)]*len(array.shape)

    #pad after if int is provided
    if pad_right:
        pad_width[axis] = (0, add_size)
    else:
        pad_width[axis] = (add_size, 0)

    return np.pad(array, pad_width=pad_width, mode='constant', constant_values=pad_value)

def tensorise_smiles(smiles, max_degree=5, max_atoms=None):
    '''Takes a list of smiles and turns the graphs in tensor representation.
    # Arguments:
        smiles: a list (or iterable) of smiles representations
        max_atoms: the maximum number of atoms per molecule (to which all
            molecules will be padded), use `None` for auto
        max_degree: max_atoms: the maximum number of neigbour per atom that each
            molecule can have (to which all molecules will be padded), use `None`
            for auto
        **NOTE**: It is not recommended to set max_degree to `None`/auto when
            using `NeuralGraph` layers. Max_degree determines the number of
            trainable parameters and is essentially a hyperparameter.
            While models can be rebuilt using different `max_atoms`, they cannot
            be rebuild for different values of `max_degree`, as the architecture
            will be different.
            For organic molecules `max_degree=5` is a good value (Duvenaud et. al, 2015)
    # Returns:
        atoms: np.array, An atom feature np.array of size `(molecules, max_atoms, atom_features)`
        bonds: np.array, A bonds np.array of size `(molecules, max_atoms, max_neighbours)`
        edges: np.array, A connectivity array of size `(molecules, max_atoms, max_neighbours, bond_features)`
    TODO:
        * Arguments for sparse vector encoding
    '''

    # import sizes
    n = len(smiles)
    n_atom_features = num_atom_features()
    n_bond_features = num_bond_features()

    # preallocate atom tensor with 0's and bond tensor with -1 (because of 0 index)
    # If max_degree or max_atoms is set to None (auto), initialise dim as small
    #   as possible (1)
    atom_tensor = np.zeros((n, max_atoms or 1, n_atom_features),dtype='float32')
    bond_tensor = np.zeros((n, max_atoms or 1, max_degree or 1, n_bond_features),dtype='float32')
    edge_tensor = -np.ones((n, max_atoms or 1, max_degree or 1), dtype=int)
    #edge_tensor = np.zeros((n, max_atoms or 1, max_atoms or 1),dtype='float32')
    #di=np.diag_indices(max_atoms)
    #for i in range(n):
    #    x=edge_tensor[i]
    #    x[di]=1
    #    edge_tensor[i]=x

    for mol_ix, s in enumerate(smiles):

        #load mol, atoms and bonds
        mol = Chem.MolFromSmiles(s)
        assert mol is not None, 'Could not parse smiles {}'.format(s)
        atoms = mol.GetAtoms()
        bonds = mol.GetBonds()

        # If max_atoms is exceeded, resize if max_atoms=None (auto), else raise
        if len(atoms) > atom_tensor.shape[1]:
            assert max_atoms is None, 'too many atoms ({0}) in molecule: {1}'.format(len(atoms), s)
            atom_tensor = padaxis(atom_tensor, len(atoms), axis=1)
            bond_tensor = padaxis(bond_tensor, len(atoms), axis=1)
            edge_tensor = padaxis(edge_tensor, len(atoms), axis=1, pad_value=-1)

        rdkit_ix_lookup = {}
        connectivity_mat = {}

        for atom_ix, atom in enumerate(atoms):
            # write atom features
            atom_tensor[mol_ix, atom_ix, : n_atom_features] = atom_features(atom)

            # store entry in idx
            rdkit_ix_lookup[atom.GetIdx()] = atom_ix

        # preallocate array with neighbour lists (indexed by atom)
        connectivity_mat = [ [] for _ in atoms]

        for bond in bonds:
            # lookup atom ids
            a1_ix = rdkit_ix_lookup[bond.GetBeginAtom().GetIdx()]
            a2_ix = rdkit_ix_lookup[bond.GetEndAtom().GetIdx()]

            # lookup how many neighbours are encoded yet
            a1_neigh = len(connectivity_mat[a1_ix])
            a2_neigh = len(connectivity_mat[a2_ix])

            # If max_degree is exceeded, resize if max_degree=None (auto), else raise
            new_degree = max(a1_neigh, a2_neigh) + 1
            if new_degree > bond_tensor.shape[2]:
                assert max_degree is None, 'too many neighours ({0}) in molecule: {1}'.format(new_degree, s)
                bond_tensor = padaxis(bond_tensor, new_degree, axis=2)
                edge_tensor = padaxis(edge_tensor, new_degree, axis=2, pad_value=-1)

            # store bond features
            bond_features_var = np.array(bond_features(bond), dtype=int)
            bond_tensor[mol_ix, a1_ix, a1_neigh, :] = bond_features_var
            bond_tensor[mol_ix, a2_ix, a2_neigh, :] = bond_features_var

            #add to connectivity matrix
            connectivity_mat[a1_ix].append(a2_ix)
            connectivity_mat[a2_ix].append(a1_ix)

        #store connectivity matrix
        for a1_ix, neighbours in enumerate(connectivity_mat):
            degree = len(neighbours)
            #edge_tensor[mol_ix, a1_ix, neighbours] = 1.
            edge_tensor[mol_ix, a1_ix, : degree] = neighbours

    return atom_tensor, bond_tensor, edge_tensor