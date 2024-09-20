import numpy as np

def usaParam(betaH, betaF, kappaH, kappaF,given):
  G = given['G']
   
  GH = np.bitwise_or(G[:,0]==2, G[:,0]==4) #Did Fatah attack last period?
  GF = np.bitwise_or(G[:,0]==3, G[:,0]==4) #Did Hammas attack last period?

  Xh = np.vstack((GH*G[:, 1], G[:, 1]))
  Xf = np.vstack((GF*G[:, 2], G[:, 2]))
  Uh = Xh.T.dot(np.concatenate((betaH, kappaH))).reshape(-1,1)
  Uf = Xf.T.dot(np.concatenate((betaF,kappaF))).reshape(-1,1)

  usa = np.hstack((Uh, Uf))
  return usa



