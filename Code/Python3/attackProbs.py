"""
Function for computing the choice probabilities
for each actor
"""
import numpy as np
def attackProbs(V):
  S = len(V)//4
  vH = np.reshape(V[:(2*S)], (-1,S), order='F')
  vF = np.reshape(V[(2*S):], (-1,S), order='F')
  
  mvH = np.amax(vH, axis=0)
  mvF = np.amax(vF, axis=0)
  evHn = np.exp(vH - mvH)
  evFn = np.exp(vF - mvF)
  
  #Conditional Choice Probabilites
  PH = (evHn/np.sum(evHn, axis=0)).reshape((-1,1), order='F')
  PF = (evFn/np.sum(evFn, axis=0)).reshape((-1,1), order='F')
  
  prAH = PH.reshape(2,S, order="F")[1,:]
  prAF = PF.reshape(2,S, order="F")[1,:]
  return (prAH, prAF)
