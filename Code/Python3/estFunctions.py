"""
Functions for the estimation
Likelihood, constraint, equilibrium, others
"""
import numpy as np
from usaParam import *
from adolc import adouble 
from adolc import erf 
from attackProbs import *  
from adolc import condassign

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



def norm_cdf(x):
  return 0.5 + 0.5*(erf(x/np.sqrt(2.0)))

def gamma2trans(gamma,sigma, given):
  nrows = given['G'].shape[0]
  G = given['G']
  states = given['states']
  Z= np.hstack((np.ones((nrows, 1)),G[:,np.array([1,2,0])], np.reshape(G[:,0] * G[:,1], (-1,1)), np.reshape(G[:,0] * G[:,2], (-1,1))))
  d = (states[3]-states[2])/2.0
  mu = Z.dot(gamma)
  pr = adouble(np.zeros((nrows,nrows//4)))
  for i in range(nrows):
    for j in range(nrows//4):
        hi = norm_cdf((states[j] + d - mu[i])/sigma)
        lo = norm_cdf((states[j] - d - mu[i])/sigma)
        pr[i,j] =  hi*1*(j==0) + (hi-lo)*1*(j>0)*(j<(nrows/4-1)) + (1-lo)*1*(j==(nrows/4-1))
    
  return pr

def constraint_parameter(beta, kappa, gamma,v, given):
    betaH = beta[:np.int_(given['nBeta']/2)]
    betaF = beta[len(betaH):(given['nBeta'])]
    kappaH = kappa[:np.int_(given['nKappa']/2)]
    kappaF = kappa[len(kappaH):(given['nKappa'])]

    betaH = np.array([beta[0]])
    betaF = np.array([beta[1]])
    kappaH = np.array([kappa[0]])
    kappaF = np.array([kappa[1]])
    sigma = gamma[-1]
    gamma1 = gamma[:-1]
    P = gamma2trans(gamma1, sigma, given)
    U = usaParam(betaH, betaF, kappaH, kappaF,given)

    
    V0 = v.reshape((-1,1), order='F')
    V = PsiQRE(V0, U,P, given['delta'])
    dV = V- np.ravel(V0)
    
    return dV
    


def logLik_parameter(beta, kappa, gamma,v, Y, given):
  V = v.reshape((-1,1), order='F')
  VH, VF = attackProbs(V)
  LL = np.empty(Y.shape[0], dtype="O")
  N = Y.shape[0]
  states = given['states']
  Z= np.hstack((np.ones((N-1, 1)),Y[:(N-1),np.array([1,2,0])], np.reshape(Y[:(N-1),0] * Y[:(N-1),1], (-1,1)), np.reshape(Y[:(N-1),0] * Y[:(N-1),2], (-1,1))))
  d = (states[3]-states[2])/2.0

  mu = Z.dot(gamma[:-1])
  sigma = gamma[-1]
  pr = adouble(np.ones((N,)))
  pr0= adouble(np.ones((N,)))
  pr1= adouble(np.ones((N,)))
  for i in range(N-1):
      hi = norm_cdf((Y[1:,0][i] + d - mu[i])/sigma)
      lo = norm_cdf((Y[1:,0][i] - d - mu[i])/sigma)
      pr[i+1] = condassign(pr[i+1], (hi-lo)-1e-5, hi-lo, hi-lo+1e-5)
      pr0[i+1] = condassign(pr0[i+1], (hi)-1e-5, hi, hi+1e-5)
      pr1[i+1] = condassign(pr1[i+1], (1-lo)-1e-5, 1-lo, 1-lo+1e-5)
      

  lowest = 1*(Y[:,0]==np.min(given['stateActions'][:,0]))
  highest = 1*(Y[:,0]==np.max(given['stateActions'][:,0]))
  trans = pr*(1-lowest)*(1-highest) + pr0*(lowest)*(1-highest) + pr1*(1-lowest)*(highest)
  
  
  
  for i in range(N):
      prAH = VH[np.where(given['states']==Y[i,0])[0]]
      pH = prAH*Y[i,1] + (1-prAH)*(1-Y[i,1])
      prAF = VF[np.where(given['states']==Y[i,0])[0]]
      pF = prAF*Y[i,2] + (1-prAF)*(1-Y[i,2])
      LL[i] = (np.log( pH) + np.log(pF))[0]
      
  LL = LL + np.log(trans)     
  return -LL



def logLik_parameter2(beta, kappa, gamma,v, Y, given):
  V = v.reshape((-1,1), order='F')
  VH, VF = attackProbs(V)
  LL = np.empty((Y.shape[0],3), dtype="O")
  N = Y.shape[0]
  states = given['states']
  Z= np.hstack((np.ones((N-1, 1)),Y[:(N-1),np.array([1,2,0])], np.reshape(Y[:(N-1),0] * Y[:(N-1),1], (-1,1)), np.reshape(Y[:(N-1),0] * Y[:(N-1),2], (-1,1))))
  d = (states[3]-states[2])/2.0

  mu = Z.dot(gamma[:-1])
  sigma = gamma[-1]
  pr = adouble(np.ones((N,)))
  pr0= adouble(np.ones((N,)))
  pr1= adouble(np.ones((N,)))
  for i in range(N-1):
      hi = norm_cdf((Y[1:,0][i] + d - mu[i])/sigma)
      lo = norm_cdf((Y[1:,0][i] - d - mu[i])/sigma)
      pr[i+1] = condassign(pr[i+1], (hi-lo)-1e-5, hi-lo, hi-lo+1e-5)
      pr0[i+1] = condassign(pr0[i+1], (hi)-1e-5, hi, hi+1e-5)
      pr1[i+1] = condassign(pr1[i+1], (1-lo)-1e-5, 1-lo, 1-lo+1e-5)
      

  lowest = 1*(Y[:,0]==np.min(given['stateActions'][:,0]))
  highest = 1*(Y[:,0]==np.max(given['stateActions'][:,0]))
  trans = pr*(1-lowest)*(1-highest) + pr0*(lowest)*(1-highest) + pr1*(1-lowest)*(highest)
  
  
  
  for i in range(N):
      prAH = VH[np.where(given['states']==Y[i,0])[0]]
      pH = prAH*Y[i,1] + (1-prAH)*(1-Y[i,1])
      prAF = VF[np.where(given['states']==Y[i,0])[0]]
      pF = prAF*Y[i,2] + (1-prAF)*(1-Y[i,2])
      LL[i,0] = (np.log( pH))[0]
      LL[i,1] = (np.log( pF))[0]
      LL[i,2] = np.log(trans[i])
      
  LL = np.ravel(LL)
  return -LL



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

