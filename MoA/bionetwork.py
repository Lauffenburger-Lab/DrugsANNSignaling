'''
Adapted from https://github.com/Lauffenburger-Lab/LEMBAS
'''
import torch
import torch.nn as nn
import scipy.sparse
import numpy.random
import pandas
import numpy
from scipy.sparse.linalg import eigs
from scipy.linalg import eig

########Activation function
def activation(x, leak):
    #x[x<0] = leak*x[x<0]
    #x[x>=0.5] = 0.5 * (1 + (1./(0.5/(x[x>=0.5]-0.5) + 1)))
    x = numpy.where(x < 0, x * leak, x)
    x = numpy.where(x > 0.5, 1-0.25/x, x) #Pyhton will display division by zero warning since it evaluates both before selecting
    return x

def deltaActivation(x, leak):
#    middleIndex = numpy.logical_and(x<0.5, x>0)
#    x[x>=0.5] = 0.25/(((x[x>=0.5]-0.5) + 0.5)**2)
#    x[middleIndex] = 1
#    x[x<0] = leak
    # middleIndex = numpy.logical_and(x<=0.5, x>0)
    # x = numpy.where(x > 0.5, 0.25/(((x-0.5) + 0.5)**2), x)
    # x[middleIndex] = 1
    # x = numpy.where(x < 0, leak, x)
    y = numpy.ones(x.shape) #derivative = 1 if nothing else is stated
    y = numpy.where(x <= 0, leak, y)  #let derivative be 0.01 at x=0
    #y = numpy.where(x > 0.5, 0.25/(((x-0.5) + 0.5)**2), y)
    y = numpy.where(x > 0.5, 0.25/(x**2), y)
    return y



def oneStepDeltaActivationFactor(yhatFull, leak):
    y = torch.ones(yhatFull.shape, dtype=yhatFull.dtype)
    piece1 = yhatFull<=0
    piece3 = yhatFull>0.5
    y[piece1] = torch.tensor(leak, dtype=yhatFull.dtype) #there is a bug in torch that sets this to 0 if piece1 all true, will probably never happen
    y[piece3] = 0.25/((-0.25/(yhatFull[piece3]-1))**2)
    return y

# def oneStepActivationFactor(yhatFull, leak):
#     y = torch.ones(yhatFull.shape, dtype=yhatFull.dtype)
#     piece1 = yhatFull<=0
#     piece3 = yhatFull>0.5
#     y[piece1] = torch.tensor(leak, dtype=yhatFull.dtype) #there is a bug in torch that sets this to 0 if piece1 all true, will probably never hapen
#     y[piece3] = 4 * (yhatFull[piece3] - yhatFull[piece3]**2)
#     return y

# def activationFactor(x, leak):
#     #x[x<0] = leak*x[x<0]
#     #x[x>=0.5] = 0.5 * (1 + (1./(0.5/(x[x>=0.5]-0.5) + 1)))
#     #y = numpy.ones(x.shape)
#     y = numpy.where(x <= 0, leak, 1)
#     y = numpy.where(x > 0.5, 0.5 * (1 + (1/(0.5/(x-0.5) + 1)))/x, y) #Pyhton will display division by zero warning since it evaluates both before selecting
#     return y

def invActivation(x, leak):
    if leak>0:
        x = numpy.where(x < 0, x/leak, x)
    else:
        x = numpy.where(x < 0, 0, x)
    x = numpy.where(x > 0.5, -0.25/(x-1), x) #Pyhton will display division by zero warning since it evaluates both before selecting
    return x

# def activation(x, leak):
#     x = numpy.where(x <= 0, x * leak, x)
#     x = numpy.where(x > 0, 1/((1/x) + 1), x) #Pyhton will display division by zero warning since it evaluates both before selecting
#     return x

# def activationFactor(x, leak):
#     y = numpy.where(x <= 0, leak, 1)
#     y = numpy.where(x > 0, (1/((1/x) + 1))/x, y) #Pyhton will display division by zero warning since it evaluates both before selecting
#     return y

# def deltaActivation(x, leak):
#     y = numpy.ones(x.shape) #derivative = 1 if nothing else is stated
#     y = numpy.where(x <= 0, leak, y)  #let derivative be 0.01 at x=0
#     y = numpy.where(x > 0, 1/((x + 1)**2), y)
#     return y

# def invActivation(x, leak):
#     if leak>0:
#         x = numpy.where(x < 0, x/leak, x)
#     else:
#         x = numpy.where(x < 0, 0, x)
#     x = numpy.where(x > 0, -(1*x)/(x - 1), x)
#     return x



def gradCliping(grad, n):
    clipingFilter = grad<-n
    grad[clipingFilter] = numpy.tanh(grad[clipingFilter]+n) - n
    clipingFilter = grad>n
    grad[clipingFilter] = numpy.tanh(grad[clipingFilter]-n) + n
    return grad




##########Spectral radius

class spectralRadius(torch.autograd.Function):
    @staticmethod
    def forward(ctx, weights, A, networkList):
        ctx.dim = weights.shape
        ctx.tol = 10**-6
        weights = weights.detach().numpy().flatten()

        ctx.networkList = networkList
        ctx.weights = weights
        ctx.A = A
        ctx.A.data = ctx.weights


        try:
            e, v = eigs(ctx.A, k=1, which='LM', ncv=100, tol=ctx.tol)
            v = v[:,0]
            e = e[0]
        except  (KeyboardInterrupt, SystemExit):
            raise
        except:
            print('Forward fail (did not find any eigenvalue with eigs)')
            tmpA = ctx.A.toarray()
            e, v, w = lreig(tmpA) #fall back to solving full eig problem

        spectralRadius = numpy.abs(e)
        ctx.e = e
        ctx.v = v
        ctx.w = numpy.empty(0)

        return torch.from_numpy(numpy.asarray(spectralRadius))

    @staticmethod
    def backward(ctx, grad_output):
        v = ctx.v
        e = ctx.e
        w = ctx.w
        networkList = ctx.networkList
        tmpA = ctx.A
        tmpA.data = ctx.weights
        tmpA = tmpA.T  #tmpA.T.toarray()

        if w.shape[0]==0:
            try:
                eT = e
                if numpy.isreal(eT): #does for some reason not converge if imag = 0
                    eT = eT.real
                e2, w = eigs(tmpA, k=1, sigma=eT, OPpart='r', tol=ctx.tol)
                selected = 0 #numpy.argmin(numpy.abs(e2-eT))
                w = w[:,selected]
                e2 = e2[selected]
                #Check if same eigenvalue
                if abs(e-e2)>(ctx.tol*10):
                    print('Backward fail (eigs left returned different eigenvalue)')
                    w = numpy.empty(0)
                    #e, v, w = lreig(tmpA) #fall back to solving whole eig problem
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                print('Backward fail (did not find any eigenvalue with eigs)')
                #e, v, w = lreig(tmpA) #fall back to solving full eig problem
                delta = numpy.zeros(ctx.weights.shape)


        if w.shape[0] != 0:
            divisor = w.T.dot(v).flatten()
            if abs(divisor) == 0:
                delta = numpy.zeros(ctx.weights.shape)
                print('Empty eig')
            else:
                delta = numpy.multiply(w[networkList[0]], v[networkList[1]])/divisor
                direction = e/numpy.abs(e)
                delta = (delta/direction).real
        else:
            #print('Empty eig')
            delta = numpy.zeros(ctx.weights.shape)

        #deltaFilter = numpy.not_equal(numpy.sign(delta), numpy.sign(ctx.weights))
        #delta[deltaFilter] = 0

        delta = torch.tensor(delta, dtype = grad_output.dtype)

        constrainNorm = True
        if constrainNorm:
            norm = torch.norm(delta, 2)
            if norm>10:
                delta = delta/norm #typical seems to be ~0.36
            #delta = delta * numpy.abs(ctx.weights)
            #delta = delta/norm(delta)


        dW = grad_output * delta

        return dW, None, None, None



class bionetworkFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, weights, bias, A, networkList, parameters):
        # Load into memory
        ctx.weights = weights.detach().numpy()
        A.data = ctx.weights
        ctx.networkList = networkList
        ctx.A = A

        bIn = x.transpose(0, 1).detach().numpy() + bias.detach().numpy()

        xhat = numpy.zeros(bIn.shape, dtype=bIn.dtype)
        # xhat = numpy.random.rand(bIn.shape[0], bIn.shape[1])

        for i in range(parameters['iterations']):
            if i > 40:  # normally takes around 40 iterations to reach steady state
                if i > 41:
                    if numpy.sum(numpy.abs(xhat - xhatBefore)) < 1e-6:
                        break
                xhatBefore = xhat.copy()
            xhat = A.dot(xhat)
            xhat += bIn
            xhat = activation(xhat, parameters['leak'])

        output = torch.from_numpy(xhat)
        output.requires_grad_()
        output = output.transpose(0, 1)

        # Pass to backward
        ctx.xRaw = A.dot(xhat) + bIn  # When converged this is the same as taking inv(activation(xhat))
        ctx.x = xhat
        ctx.parameters = parameters
        return output

    @staticmethod
    def backward(ctx, grad_output):

        AT = ctx.A
        AT.data = ctx.weights.data
        AT = AT.T

        gradIn = grad_output.transpose(0, 1).detach().numpy()
        grad = numpy.zeros(gradIn.shape)
        # grad = numpy.random.rand(gradIn.shape[0], gradIn.shape[1])
        deltaX = ctx.xRaw.copy()
        deltaX = deltaActivation(deltaX, ctx.parameters['leak'])
        gradBefore = grad.copy()

        for i in range(ctx.parameters['iterations']):
            if i > 20:  # normally takes around 30 iterations to reach steady state
                if i > 21:
                    if numpy.sum(numpy.abs(grad - gradBefore)) < 1e-6:
                        break
                gradBefore = grad.copy()
            grad = deltaX * (AT.dot(grad) + gradIn)
            # as a precaution clipping with tanh for |gradients| > clipping val
            grad = gradCliping(grad, ctx.parameters['clipping'])
            # norm = numpy.linalg.norm(grad, ord='fro')
            # normFactor = numpy.where(norm>ctx.parameters['clipping'], ctx.parameters['clipping']/norm, 1)
            # grad = grad * normFactor

        output = torch.from_numpy(grad).transpose(0, 1)

        # Construct gradients
        grad_weight = torch.from_numpy(
            numpy.sum(numpy.multiply(ctx.x[ctx.networkList[1], :], grad[ctx.networkList[0], :]), axis=1))
        grad_bias = torch.from_numpy(numpy.sum(grad, axis=1)).unsqueeze(1)

        return output, grad_weight, grad_bias, None, None, None

class model(torch.nn.Module): ### CHANGED IT TO ADD THE DRUG LAYER: ADDED druginName and targetNames and targetList
    def __init__(self, networkList, nodeNames, modeOfAction, inputAmplitude, projectionFactor,inName, outName, bionetParams, valType):

        super(model, self).__init__()

        # if (drugSim is None) and (batch_size is None):
        #     raise Exception("One of the following arguments `drugSim` and `batch_size` needs to have a value. Got None for both.")

        #self.batch_size = batch_size
        #self.drugSim = drugSim
        #self.drugLayer = DrugstoSignal(drugTargets,len(druginName),len(targetNames),drugSim=drugSim,druginName=druginName,useMask=useMask,useDrugs=useDrugs)
        self.inputLayer = projectInput(nodeNames, inName, inputAmplitude, valType)
        self.network = bionet(networkList, len(nodeNames), modeOfAction, bionetParams,nodeNames, inName, valType)
        self.projectionLayer = projectOutput(nodeNames, outName, projectionFactor, valType)

    def forward(self, X):
        #Xin = self.drugLayer(X)
        fullX = self.inputLayer(X)
        fullY = self.network(fullX)
        Yhat = self.projectionLayer(fullY)
        return Yhat, fullY

def spectralLoss(model, YhatFull, weights, expFactor = 20, lb=0.5):
    bionetParams = model.network.parameters

    randomIndex = numpy.random.randint(YhatFull.shape[0])
    activationFactor = oneStepDeltaActivationFactor(YhatFull[randomIndex,:], bionetParams['leak']).detach()
    weightFactor = activationFactor[model.network.networkList[0]]
    multipliedWeightFactor = weights * weightFactor
    spectralRadius = model.network.getSpectralRadius(multipliedWeightFactor)
    #spectralClampFactor = 1/torch.max(spectralRadius.detach()/bionetParams['spectralLimit'], torch.tensor(1.0).double()) #Prevents infinte penalty
    #spectralRadiusLoss =  (1/(1 - spectralClampFactor*spectralRadius) - 1)

    # if spectralRadius>bionetParams['spectralTarget']:
    #     spectralRadiusLoss = torch.abs(spectralRadius)
    # else:
    #     spectralRadiusLoss = torch.tensor(0.0)

    scaleFactor = 1/numpy.exp(expFactor *  bionetParams['spectralTarget'])

    if spectralRadius>lb:
        spectralRadiusLoss = scaleFactor * (torch.exp(expFactor*spectralRadius)-1)
    else:
        spectralRadiusLoss = torch.tensor(0.0)

    return spectralRadiusLoss, spectralRadius


def uniformLoss(curState, dataIndex, YhatFull, targetMin = 0, targetMax = 0.99, maxConstraintFactor = 10):
    data = curState.detach().clone()
    data[dataIndex, :] = YhatFull

    targetMean = (targetMax-targetMin)/2
    targetVar= (targetMax-targetMin)**2/12

    factor = 1
    meanFactor = factor
    varFactor = factor
    minFactor = factor
    maxFactor = factor
    maxConstraintFactor = factor * maxConstraintFactor

    nodeMean = torch.mean(data, dim=0)
    nodeVar = torch.mean(torch.square(data-nodeMean), dim=0)
    maxVal, _ = torch.max(data, dim=0)
    minVal, _ = torch.min(data, dim=0)

    meanLoss = meanFactor * torch.sum(torch.square(nodeMean - targetMean))
    varLoss =  varFactor * torch.sum(torch.square(nodeVar - targetVar))
    maxLoss = maxFactor * torch.sum(torch.square(maxVal - targetMax))
    minloss = minFactor * torch.sum(torch.square(minVal- targetMin))
    maxConstraint = -maxConstraintFactor * torch.sum(maxVal[maxVal.detach()<=0]) #max value should never be negative

    loss = meanLoss + varLoss + minloss + maxLoss + maxConstraint
    return loss

def generateRandomInput(model, N, simultaniousInput):
    sizeX = model.inputLayer.weights.shape[0]
    X = torch.zeros((N, sizeX), dtype=torch.double)
    for i in range(1, N): #leave first block blank
        selected = numpy.random.randint(sizeX, size=simultaniousInput)
        X[i, selected] = torch.rand(simultaniousInput, dtype=torch.double)
    return X

def reduceSpectralRadius(model, spectralTarget, localX, maxIter = 100):
    N = localX.shape[0]
    leak = model.network.parameters['leak']
    networkList = model.network.networkList

    localY, localFull = model(localX)
    localY = localY.clone().detach()
    criterion = torch.nn.MSELoss(reduction='mean')
    maxSr = numpy.inf
    optimizer = torch.optim.Adam(model.parameters(), lr = 0.0005, weight_decay=0)
    noise = 1e-8
    srFactor =  1e-5

    print('Reducing spectral radius for all input (target ', spectralTarget, '):')
    for k in range(maxIter):
        optimizer.zero_grad()
        Yhat, YhatFull = model(localX)
        fitLoss = criterion(Yhat, localY)
        signConstraint = torch.sum(torch.abs(model.network.weights[model.network.getViolations()]))
        spectralRadius = torch.zeros(N, dtype=torch.double)
        spectralRegulation = torch.zeros(N, dtype=torch.double)
        for i in range(N):
            activationFactor = oneStepDeltaActivationFactor(localFull[i,:].flatten(), leak).detach()
            weightFactor = activationFactor[networkList[0]]
            multipliedWeightFactor = model.network.weights * weightFactor
            sr = model.network.getSpectralRadius(multipliedWeightFactor)
            spectralRadius[i] = sr.item()
            if spectralRadius[i]>spectralTarget:
                spectralRegulation[i] = srFactor * torch.abs(sr)
                #spectralClampFactor = 1/torch.max(spectralRadius[i]/bionetParams['spectralLimit'], torch.tensor(1.0).double()) #Prevents infinte penalty
                #spectralRegulation[i] = (1/(1 - spectralClampFactor*sr) - 1)
            else:
                spectralRegulation[i] = torch.tensor(0, dtype=torch.double, requires_grad=True)

        loss = 1e-5 * fitLoss + torch.sum(spectralRegulation) + signConstraint
        loss.backward()
        optimizer.step()
        model.network.weights.data = model.network.weights.data + torch.randn(model.network.weights.shape) * noise
        maxSr = numpy.max(spectralRadius.detach().numpy())
        print(k, maxSr)
        if maxSr<spectralTarget:
            break

    return model

def oneCycle(e, maxIter, maxHeight = 2e-3, minHeight = 1e-8, peak = 1000):
    phaseLength = 0.95 * maxIter

    if e<=peak:
        effectiveE = e/peak
        lr = (maxHeight-minHeight) * 0.5 * (numpy.cos(numpy.pi*(effectiveE+1))+1) + minHeight
    elif e<=phaseLength:
        effectiveE = (e-peak)/(phaseLength-peak)
        lr = (maxHeight-minHeight) * 0.5 * (numpy.cos(numpy.pi*(effectiveE+2))+1) + minHeight
    else:
        lr = minHeight

    return lr

def getSamples(N, batchSize):
    order = numpy.random.permutation(N)
    outList = []
    while len(order)>0:
        outList.append(order[0:batchSize])
        order = order[batchSize:]
    return outList

def getAllSpectralRadius(model, YhatFull):
    leak = model.network.parameters['leak']
    sr = numpy.zeros(YhatFull.shape[0])
    activationFactor = oneStepDeltaActivationFactor(YhatFull, leak)
    for i in range(len(sr)):
        weightFactor = activationFactor[i, model.network.networkList[0]]
        sr[i] = model.network.getSpectralRadius(model.network.weights * weightFactor).item()
    return sr

class bionet(nn.Module):
    def __init__(self, networkList, size, modeOfAction, parameters,nodeNames, inName, dtype):
        super().__init__()
        self.parameters = parameters

        self.size_in = size
        self.size_out = size
        self.networkList = networkList
        self.modeOfAction = torch.tensor(modeOfAction)
        self.type = dtype
        self.nodeNames = nodeNames
        self.inName = inName

        # initialize weights and biases
        weights, bias = self.initializeWeights()
        #weights = self.initializeWeights()

        #for this to work as intended network list must be sorted on index 0
        self.A = scipy.sparse.csr_matrix((weights.detach().numpy(), networkList), shape=(size, size), dtype='float64')

        self.weights = nn.Parameter(weights)
        self.bias = nn.Parameter(bias)



    def forward(self, x):
        return bionetworkFunction.apply(x, self.weights, self.bias, self.A, self.networkList, self.parameters) # deleted self.bias

    def getWeight(self, nodeNames, source, target):
        self.A.data = self.weights.detach().numpy()
        locationSource = numpy.argwhere(numpy.isin(nodeNames, source))[0]
        locationTarget = numpy.argwhere(numpy.isin(nodeNames, target))[0]
        weight = self.A[locationTarget, locationSource][0]
        return weight

    def getViolations(self, weights = None):
        if weights == None:
            weights = self.weights.detach()
        wrongSignActivation = torch.logical_and(weights<0, self.modeOfAction[0] == True)#.type(torch.int)
        wrongSignInhibition = torch.logical_and(weights>0, self.modeOfAction[1] == True)#.type(torch.int)
        return torch.logical_or(wrongSignActivation, wrongSignInhibition)

    # def getSpectralLoss(self):
    #     self.A.data = self.weights.detach().numpy()
    #     e, v, w = lreigs(self.A)
    #     spectralRadius = numpy.abs(e)
    #     divisor = w.T.dot(v)

    #     if abs(divisor)>10**-4:
    #         w = w/divisor #Ensures wT*v=1
    #         delta = numpy.multiply(w[self.networkList[0]], v[self.networkList[1]])
    #         delta = numpy.squeeze(delta)
    #         direction = spectralRadius/e
    #         delta = (direction * delta).real
    #         deltaFilter = numpy.not_equal(numpy.sign(delta), numpy.sign(self.weights.detach().numpy()))
    #         delta[deltaFilter] = 0
    #     else:
    #         print('missmatch using L2 instead')
    #         delta = numpy.sign(self.weights.detach().numpy())  #Could occure if degenerate eigenvalues
    #     delta = delta/norm(delta, 1)
    #     delta = torch.tensor(delta).reshape([-1,1])
    #     return delta, spectralRadius
    def getSpectralRadius(self, weights):
        return spectralRadius.apply(weights, self.A, self.networkList)

    def getRevSpectralRadius(self, weights):
        return spectralRadius.apply(weights, self.A.T, self.networkList)

    def initializeWeights(self):
        # The idea is to scale input to a source so that well connected nodes have lower weights
        # weights = 0.5 + 0.1 * (torch.rand(self.networkList.shape[1])-0.5)
        # bias = 0.01 + 0.001 * (torch.rand(self.size_in,1)-0.5)

        weights = 0.1 + 0.1 * torch.rand(self.networkList.shape[1], dtype=self.type)
        weights[self.modeOfAction[1, :]] = -weights[self.modeOfAction[1, :]]
        bias = 1e-3 * torch.ones((self.size_in, 1), dtype=self.type) # was 1e-3 instead of 0.5
        dictionary = dict(zip(self.nodeNames, list(range(len(self.nodeNames)))))
        nodeOrder = numpy.array([dictionary[x] for x in self.inName])
        bias[nodeOrder] = bias[nodeOrder]*500. #put 0.5 as input in targets
        values, counts = numpy.unique(self.networkList[0,:], return_counts=True)

        for i in range(self.size_in):
            affectedIn = self.networkList[0, :] == i
            if numpy.any(affectedIn):
                if torch.all(weights[affectedIn] < 0):
                    bias.data[i] = 1  # only affected by inhibition, relies on bias for signal

        # for i in range(self.size_in):
        #     affectedIn = self.networkList[0,:] == i
        #     fanIn = max(sum(affectedIn), 1)
        #     affectedOut = self.networkList[0,:] == i
        #     fanOut = max(sum(affectedOut), 1)
        #     weights.data[affectedIn] = weights.data[affectedIn] * numpy.sqrt(2.0/numpy.sqrt(fanIn * fanOut))

        return weights, bias


    def preScaleWeights(self, targetRadius = 0.8):
        spectralRadius = self.getSpectralRadius(self.weights)
        factor = targetRadius/spectralRadius.item()
        self.weights.data = self.weights.data * factor
        # print('Pre-scaling eig')
        # optimizer = torch.optim.Adam([self.weights], lr=0.001)
        # weightFactor = self.weights
        # for i in range(1000):
        #     optimizer.zero_grad()
        #     spectralRadius = self.getSpectralRadius(weightFactor)
        #     if i % 20 == 0:
        #         print('i={:.0f}, e={:.4f}'.format(i, spectralRadius.item()))

        #     if spectralRadius.item()>targetRadius:
        #         spectralRadius.backward()
        #         optimizer.step()
        #     else:
        #         break


    # def getSpectralR(self, Xfull):
    #     rawX =


# def ismember(a, b):
#     bind = {}
#     for i, elt in enumerate(b):
#         if elt not in bind:
#             bind[elt] = i
#     return [bind.get(itm, None) for itm in a]

class projectInput(nn.Module):
    def __init__(self, nodeList, inputNames, amplitude, type):
        super().__init__()

        self.size_in = len(inputNames)
        self.size_out = len(nodeList)
        self.type = type
        dictionary = dict(zip(nodeList, list(range(len(nodeList)))))
        self.nodeOrder = numpy.array([dictionary[x] for x in inputNames])
        weights = amplitude * torch.ones(len(inputNames), dtype=type)
        self.weights = nn.Parameter(weights)

    def forward(self, x):
        curIn = torch.zeros([x.shape[0],  self.size_out], dtype=self.type)
        curIn[:, self.nodeOrder] = self.weights * x
        return curIn

class projectOutput(nn.Module):
    def __init__(self, nodeList, outputNames, scaleFactor, type):
        super().__init__()

        self.size_in = len(nodeList)
        self.size_out = len(outputNames)
        self.type = type

        dictionary = dict(zip(nodeList, list(range(len(nodeList)))))
        self.nodeOrder = numpy.array([dictionary[x] for x in outputNames])

        #bias = torch.zeros(len(outputNames), dtype=type)
        weights = scaleFactor * torch.ones(len(outputNames), dtype=type)
        self.weights = nn.Parameter(weights)

    def forward(self, x):
        curOut = self.weights * x[:, self.nodeOrder]
        return curOut


def getEigenvalue(model):
    try:
        eigenValue = abs(eigs(model.A, k=1)[0])
        eigenvalue = torch.from_numpy(eigenValue)
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        eigenvalue = numpy.nan
    return eigenvalue

def lreig(A):
    #fall back if eigs fails
    e, w, v = eig(A, left = True)
    selected = numpy.argmax(numpy.abs(e))
    eValue = e[selected]
    # selected = (e == eValue)

    # if numpy.sum(selected) == 1:
    w = w[:,selected]
    v = v[:,selected]
    # else:
    #     w = numpy.sum(w[:,selected], axis=1, keepdims=True)
    #     v = numpy.sum(v[:,selected], axis=1, keepdims=True)
    #     w = w/norm(w)
    #     v = v/norm(v)
    return eValue, v, w

def getRandomNet(networkSize, sparsity):
    network = scipy.sparse.random(networkSize, networkSize, sparsity)
    scipy.sparse.lil_matrix.setdiag(network, 0)
    networkList = scipy.sparse.find(network)
    networkList = numpy.array((networkList[1], networkList[0])) #we flip the network for rowise ordering
    nodeNames = [str(x+1) for x in range(networkSize)]
    #weights = torch.from_numpy(networkList[2])
    return networkList, nodeNames

def loadNetwork(filename, banList = []):
    net = pandas.read_csv(filename, sep='\t', index_col=False)
    net = net[~ net["source"].isin(banList)]
    net = net[~ net["target"].isin(banList)]

    sources = list(net["source"])
    targets = list(net["target"])
    stimulation = numpy.array(net["stimulation"])
    inhibition = numpy.array(net["inhibition"])
    modeOfAction = 0.1 * numpy.ones(len(sources))
    modeOfAction[stimulation==1] = 1
    modeOfAction[inhibition==1] = -1

    networkList, nodeNames, weights = makeNetworkList(sources, targets, modeOfAction)  #0 == Target 1 == Source due to numpy sparse matrix structure
    modeOfAction = numpy.array([[weights==1],[weights==-1]]).squeeze()

    return networkList, nodeNames, modeOfAction

def makeNetworkList(sources, targets, weights):
    nodeNames = list(numpy.unique(sources + targets))
    dictionary = dict(zip(nodeNames, list(range(len(nodeNames)))))
    sourceNr = numpy.array([dictionary[x] for x in sources]) #colums
    targetNr = numpy.array([dictionary[x] for x in targets]) #rows
    size = len(nodeNames)
    A = scipy.sparse.csr_matrix((weights, (sourceNr, targetNr)), shape=(size, size))
    networkList = scipy.sparse.find(A)
    weights = networkList[2]
    networkList = numpy.array((networkList[1], networkList[0]))  #0 == Target 1 == Source due to numpy sparse matrix structure
    return networkList, nodeNames, weights


def trainingParameters(**attributes):
    #set defaults
    params = {'iterations': 150, 'leak': 0.01, 'clipping': 1,  'targetPrecision': 1e-4}

    for curKey in params.keys():
        if curKey in attributes.keys():
            params[curKey] = attributes[curKey]

    if 'spectralTarget' in attributes.keys():
        params[curKey] = attributes[curKey]
    else:
        params['spectralTarget'] = numpy.exp(numpy.log(params['targetPrecision'])/params['iterations'])

    return params



def saveParam(model, nodeList, fileName):
    nodeList = numpy.array(nodeList).reshape([-1, 1])

    #Weights
    networkList = model.network.networkList
    sources = nodeList[networkList[1]]
    targets = nodeList[networkList[0]]
    paramType = numpy.array(['Weight'] * len(sources)).reshape([-1, 1])
    values = model.network.weights.detach().numpy().reshape([-1, 1])
    data1 = numpy.concatenate((sources, targets, paramType, values), axis=1)

    #Bias
    sources = nodeList
    targets = numpy.array([''] * len(sources)).reshape([-1, 1])
    paramType = numpy.array(['Bias'] * len(sources)).reshape([-1, 1])
    values = model.network.bias.detach().numpy().reshape([-1, 1])
    data2 = numpy.concatenate((sources, targets, paramType, values), axis=1)

    #Projection
    projectionList = model.projectionLayer.nodeOrder
    sources = nodeList[projectionList]
    targets = numpy.array([''] * len(sources)).reshape([-1, 1])
    paramType = numpy.array(['Projection'] * len(sources)).reshape([-1, 1])
    values = model.projectionLayer.weights.detach().numpy().reshape([-1, 1])
    data3 = numpy.concatenate((sources, targets, paramType, values), axis=1)

    #Input projection
    projectionList = model.inputLayer.nodeOrder
    sources = nodeList[projectionList]
    targets = numpy.array([''] * len(sources)).reshape([-1, 1])
    paramType = numpy.array(['Input'] * len(sources)).reshape([-1, 1])
    values = model.inputLayer.weights.detach().numpy().reshape([-1, 1])
    data4 = numpy.concatenate((sources, targets, paramType, values), axis=1)

    data = numpy.concatenate((data1, data2, data3, data4))

    pd = pandas.DataFrame(data)
    pd.columns = ['Source', 'Target', 'Type', 'Value']
    pd.to_csv(fileName, sep = '\t', quoting = None, index = False)

def loadParam(fileName, model, nodeNames):
    dictionary = dict(zip(nodeNames, list(range(len(nodeNames)))))
    data = pandas.read_csv(fileName, delimiter = '\t')

    #Reset model to zero
    model.inputLayer.weights.data = torch.zeros(model.inputLayer.weights.shape)
    model.network.weights.data = torch.zeros(model.network.weights.shape)
    model.network.bias.data = torch.zeros(model.network.bias.shape)
    model.projectionLayer.weights.data = torch.zeros(model.projectionLayer.weights.shape)


    inputLookup = model.inputLayer.nodeOrder
    networkLookup = model.network.networkList #model.network.A.nonzero()
    projectionLookup = model.projectionLayer.nodeOrder

    for i in range(data.shape[0]):
        curRow = data.iloc[i,:]
        source = dictionary[curRow['Source']]
        value = curRow['Value']
        if curRow['Type'] == 'Weight':
            target = dictionary[curRow['Target']]
            weightNr = numpy.argwhere(numpy.logical_and(networkLookup[1,:] == source, networkLookup[0,:] == target))
            model.network.weights.data[weightNr] = value
        elif curRow['Type'] == 'Bias':
            model.network.bias.data[source] = value
        elif curRow['Type'] == 'Projection':
            model.projectionLayer.weights.data[projectionLookup == source] = value
        elif curRow['Type'] == 'Input':
            model.inputLayer.weights.data[inputLookup == source] = value
    return model


def getMeanLoss(criterion, Y):
    averagePerOutput = torch.mean(Y, dim=0)
    medianPerOutput = torch.from_numpy(numpy.median(Y, axis=0))
    errorFromPredictingTheMean = criterion(Y, averagePerOutput)
    errorFromPredictingTheMedian = criterion(Y, medianPerOutput)
    return (errorFromPredictingTheMean.item(), errorFromPredictingTheMedian.item())


def sensitivityAnalysis(model, nodeNames, X, conditionName, fileName):
    downValue = -10
    upValue = 10
    nodeNames = numpy.array(nodeNames, dtype=object)
    n = len(nodeNames)

    Xfull = model.inputLayer(X)
    ctrlMatrix = numpy.zeros((2, n))
    downMatrix = downValue * numpy.identity(n)
    upMatrix = upValue * numpy.identity(n)

    joinedMatrix = numpy.concatenate((ctrlMatrix, downMatrix, upMatrix))
    upNames = nodeNames + '_u'
    downNames = nodeNames + '_d'
    ctrlName = numpy.array(['null', 'ctrl'])
    joinedNames = numpy.concatenate((ctrlName, downNames, upNames))
    joinedMatrix = torch.tensor(joinedMatrix)
    joinedMatrix = joinedMatrix + Xfull
    joinedMatrix[0,:] = 0

    sensitivityFull = model.network(joinedMatrix)
    df = pandas.DataFrame(sensitivityFull.detach().numpy(), index=joinedNames, columns=nodeNames)
    df = df.round(decimals=3)
    df.to_csv(fileName, sep='\t')

def generateConditionNames(X, inName):
    inName = numpy.array(inName)
    X = X.detach().numpy()
    names = numpy.empty(X.shape[0], dtype=object)
    for i in range(len(names)):
        curSelection = X[i,:]>0
        if sum(curSelection) == 0:
            names[i] = '(none)'
        else:
            curLigands = list(inName[curSelection])
            names[i] = '_'.join(curLigands)
    return names


#For reference:
class bionetworkAutoGrad(nn.Module):
    def __init__(self, networkList, size, reps=150):
        super().__init__()
        self.size_in = size
        self.size_out = size

        #requires_grad=False
        bias = torch.Tensor(size, 1)
        self.bias = nn.Parameter(bias, requires_grad = True)
        self.reps = reps
        self.leak = 0.01

        # initialize weights and biases
        weights = torch.Tensor(len(networkList[0]))
        nn.init.uniform_(weights, -0.1, 0.1) # weight init
        nn.init.uniform_(self.bias, -0.1, 0.1)  # bias init

        self.A = torch.sparse.FloatTensor(torch.from_numpy(networkList).long(), weights, torch.Size([size, size])).double()
        self.A = nn.Parameter(self.A, requires_grad = True)
        self.A.values = weights

    def forward(self, x):
        #Load into memory
        bIn = x.transpose(0, 1)
        curB = torch.add(bIn, self.bias)
        xhat = torch.zeros(bIn.shape, dtype=bIn.dtype)

        #print(self.A)
        for i in range(self.reps):
            xhat = torch.sparse.mm(self.A, xhat)
            xhat = torch.add(xhat, curB)
            xhat[xhat<0] = self.leak*xhat[xhat<0]
            xhat[xhat>0.5] = 0.5 * (1 + (1/(0.5/(xhat[xhat>0.5]-0.5) + 1)))
            #xhat[xhat>0] = 1/(0.5/(xhat[xhat>0]) + 1)

        output = xhat.transpose(0, 1)
        return output