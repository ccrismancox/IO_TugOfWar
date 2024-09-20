"""
Functions for the estimation. Altered for the tit-for-tat model
Likelihood, constraint, equilibrium, others
"""
import numpy as np
from usaParam_t4t import *
from adolc import adouble 
from adolc import erf 
from attackProbs import *  

def constraint(x, given):
    betaH = x[0:np.int_(given['nBeta']/2)]
    betaF = x[len(betaH):(given['nBeta'])]
    kappaH = x[(given['nBeta']):(given['nBeta'] + np.int_(given['nKappa']/2))]
    kappaF = x[(given['nBeta'] + np.int_(given['nKappa']/2)):(given['nBeta'] +given['nKappa'])]
    
    U = usaParam(betaH, betaF, kappaH, kappaF,given)
    
    V0 = x[given['nReal']:(given['nAll']+1)].reshape((-1,1), order='F')
    V = PsiQRE(V0, U, given['P'], given['delta'])
    dV = V- np.ravel(V0)
    
    return dV
    
  
  
def logLik(x, given):
  V = x[given['nReal']:(given['nAll']+1)].reshape((-1,1), order='F')
  S = given['P'].shape[1]

  VH = np.reshape(V[:(2*S)], (1,-1), order='F')
  VF = np.reshape(V[(2*S):], (1,-1), order='F')
  
  vHm = VH.reshape((-1, S), order='F')
  vFm = VF.reshape((-1, S), order='F')
  cSAD1 = given['stateActionsInData'][:(S*2),:].reshape((-1,1), order='F')    
  cSAD2 = given['stateActionsInData'][(S*2):,:].reshape((-1,1), order='F')    
  
  LL = VH.dot(cSAD1) -  np.sum(np.log(np.sum(np.exp(vHm), axis=0).T) * given['statesInData'])
  LL = LL + VF.dot(cSAD2)  - np.sum(np.log(np.sum(np.exp(vFm), axis=0).T) * given['statesInData'])
  
  return -LL




def PsiQRE(v, U, P, delta):
  #####################################################################################
  # INFO :: Value function
  # INPUTS ::
  # # # v : a vector of actor-state-action values
  # # # # # # v = c(v1,v2), vi state-action values for i
  # # # U : a matrix of state-action-profile utilities, ncol=2
  # # # P : a transition matrix, see model['P']
  # # # delta : a discount factor
  # OUTPUTS :: v, return values
  #####################################################################################

  S= P.shape[1]
  
  # Values, vi[j,k] := action j, state k
  vH = np.reshape(v[:(2*S)], (-1,S), order='F')
  vF = np.reshape(v[(2*S):], (-1,S), order='F')
  

  evHn = np.exp(vH)
  evFn = np.exp(vF)
  
  #Conditional Choice Probabilites
  PH = (evHn/np.sum(evHn, axis=0)).reshape((-1,1), order='F')
  PF = (evFn/np.sum(evFn, axis=0)).reshape((-1,1), order='F')


  #U_i(a|s)
  UHsa = (U[:,0] + delta * (P.dot(np.log(np.sum(evHn, axis=0)) ))).reshape((-1,1))
  UFsa = (U[:,1] + delta * (P.dot(np.log(np.sum(evFn, axis=0)) ))).reshape((-1,1))
  
  VH=np.kron(PF.reshape((-1,S), order='F'), np.ones((1,2))).reshape((-1,1), order='F') * UHsa
  VF=np.kron(PH.reshape((-1,S), order='F'), np.ones((2,1))).reshape((-1,1),order='F') * UFsa
  
  VH = np.sum(VH.reshape((-1,2)), axis=1)
  VF = np.sum(np.array(np.split(VF, S)).reshape(S, 2,2), axis=1)

  return np.concatenate((VH, np.ravel(VF, order='C')))



def logLik_pointwise(v, Y, given):
  V = v.reshape((-1,1), order='F')
  VH, VF = attackProbs(V)
  LL = np.empty(Y.shape[0], dtype="O")
  N = Y.shape[0]  
  
  for i in range(N):
      prAH = VH[np.where(given['states']==Y[i,0])[0]]
      pH = prAH*Y[i,1] + (1-prAH)*(1-Y[i,1])
      prAF = VF[np.where(given['states']==Y[i,0])[0]]
      pF = prAF*Y[i,2] + (1-prAF)*(1-Y[i,2])
      LL[i] = (np.log( pH) + np.log(pF))[0]
      
  return -LL
