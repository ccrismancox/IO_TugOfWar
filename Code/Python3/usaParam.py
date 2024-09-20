import numpy as np

def usaParam(betaH, betaF, kappaH, kappaF,given):
  G = given['G']
  
 
  Xh = np.vstack((G[:,0], G[:, 1]))
  Xf = np.vstack((G[:,0], G[:, 2]))
  Uh = Xh.T.dot(np.concatenate((betaH, kappaH))).reshape(-1,1)
  Uf = Xf.T.dot(np.concatenate((betaF,kappaF))).reshape(-1,1)

  usa = np.hstack((Uh, Uf))
  return usa


