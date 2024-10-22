"""
Code for fitting the main model
"""
from genGiven import *
from estFunctions import *
from attackProbs import *
import numpy as np
import os
import shutil
from time import time
import scipy.sparse as sps
import adolc 
import pyipopt
from numpy.linalg import inv
import scipy.stats as stats
import argparse
import pandas as pd

P = np.genfromtxt("../R4.4.1/ipoptTEMP/trans.csv", delimiter=",")[1:,:]
Xdata = np.genfromtxt("../R4.4.1/ipoptTEMP/regData.csv", delimiter=",")[1:,:]
'''
columns in X are
states(continuous); Hattacks; Fattacks; lag.states; states(discrete)
'''
Y = Xdata[:,np.array([4,1,2])]


states =  np.genfromtxt("../R4.4.1/ipoptTEMP/statespace.csv", delimiter=",")[1:]


'''parameters'''
parameters = pd.read_csv("../R4.4.1//ipoptTEMP/params.csv")
delta = parameters['delta'][0]
nK = int(parameters['nkappa'])
nB = int(parameters['nbeta'])
try:
  twostep = parameters['twostep'][0]
except:
  twostep = True
  
try:
  maxit = parameters['maxit'][0]
except:
  maxit = 5000
  
try:
  tol = parameters['tol'][0]
except:
  tol = 1e-5



given = genGiven(states, Y, P, nBeta=nB, nKappa=nK, delta=delta)
  
start =np.genfromtxt("../R4.4.1/ipoptTEMP/start.csv", delimiter=",")[1:]
x0 = start[:given['nAll']]
xL = np.zeros(given['nAux'])

    
def LL(x):
  return  logLik(x, given)


def const(x):
  """ constraint function """
  return constraint(x, given)

adolc.trace_on(1)
ax = adolc.adouble(x0)
adolc.independent(ax)
ay = LL(ax)
adolc.dependent(ay)
adolc.trace_off()

# trace constraint function
adolc.trace_on(2)
ax = adolc.adouble(x0)
adolc.independent(ax)
ay = const(ax)
adolc.dependent(ay)
adolc.trace_off()


def LLadolc(x):
    return adolc.function(1,x)

def grLLadolc(x):
    return adolc.gradient(1,x)

def const_adolc(x):
    return adolc.function(2,x)

class jac_c_adolc:
    
    def __init__(self, x):
        #options = np.array([1,1,0,0],dtype=int)
        options = None
        result = adolc.colpack.sparse_jac_no_repeat(2,x,options)
        
        self.nnz  = result[0]     
        self.rind = np.asarray(result[1],dtype=int)
        self.cind = np.asarray(result[2],dtype=int)
        self.values = np.asarray(result[3],dtype=float)
        
    def __call__(self, x, flag, user_data=None):
        if flag:
            return (self.rind, self.cind)
        else:
            result = adolc.colpack.sparse_jac_repeat(2, x, self.nnz, self.rind,
                self.cind, self.values)
            return result[3]

Jac_c_adolc = jac_c_adolc(x0)

    
    
nvar = len(x0)
x_L = np.array([-np.inf]*nvar, dtype=float)
x_U = np.array([np.inf]*nvar, dtype=float)
 
# define the inequality constraints
ncon = const(x0).shape[0]
g_L = np.array([0]*ncon, dtype=float)
g_U = np.array([0]*ncon, dtype=float)


idx = Jac_c_adolc(x0, True)
val = Jac_c_adolc(x0, False)
Jac = sps.coo_matrix( (val, (idx[0], idx[1])), shape= (given['nAll']-given['nReal'], given['nAll']))
# define the number of nonzeros in the jacobian 
nnzj = Jac.nnz
print("Jacobian has %i non-zero elements (%f%% dense)"%(nnzj, 100*nnzj/((given['nAll']-given['nReal'])*given['nAll'])))
   

             
def lagrangian(x, lagrange, obj_factor):
    return  obj_factor*LL(x) + np.dot(lagrange, const(x))   
  
  
# trace lagrangian function
adolc.trace_on(3)
ax = adolc.adouble(x0)
adolc.independent(ax)
ay = lagrangian(ax, xL, 1.0)
adolc.dependent(ay)
adolc.trace_off()


class hessLag_adolc_sp:
  def __init__(self, x, given):
      #options = np.array([0,1],dtype=int)
      options=None
      result = adolc.colpack.sparse_hess_no_repeat(3,x,options)
          
      self.cind = np.asarray(result[2],dtype=int)        
      #self.mask = np.where(self.cind < given['nAll'])    
      self.rind = np.asarray(result[1],dtype=int)
      self.cind = self.cind
      self.values = np.asarray(result[3],dtype=float)
      self.mask = np.where(self.cind < given['nAll'])    
  def __call__(self, x, lagrange,obj_factor,flag, user_data=None):
      if flag:
          return (self.rind[self.mask], self.cind[self.mask])
      else:
       #   x = np.hstack([x,lagrange,obj_factor])
          result = adolc.colpack.sparse_hess_repeat(3, x, self.rind,
                                                      self.cind, self.values)
      return result[3][self.mask]
      
      
      
hessLag_adolc = hessLag_adolc_sp(x0, given)


t0 = time()
H2 = hessLag_adolc(x0, xL, 1.0, False)
t1 = time() - t0
print("Hessian time: %d seconds"%t1)
nnzh = len(H2)
print("Hessian has %i non-zero elements (%f%% dense)"%(nnzh, 100*nnzh/(given['nAll']**2)))
    
          
          
# create the nonlinear programming model
nlp2 = pyipopt.create(
nvar,
x_L,
x_U,
ncon,
g_L,
g_U,
nnzj,
nnzh,
LLadolc,
grLLadolc,
const_adolc,
Jac_c_adolc,
hessLag_adolc
)
        

nlp2.num_option('expect_infeasible_problem_ctol', 1e-15)
nlp2.num_option('max_cpu_time', 6e5)
nlp2.num_option('tol', tol)
nlp2.str_option('warm_start_init_point', 'yes')
nlp2.int_option('max_iter', maxit)

   
t0 = time()
results = nlp2.solve(x0)
t1 = time()-t0
print("Estimation time: %d seconds" %t1)

# free the model
nlp2.close()
out = {'results': results}
    


Btheta = adolc.hessian(1, results[0])


idx = Jac_c_adolc(results[0], True)
val = Jac_c_adolc(results[0], False)
Htheta = sps.coo_matrix( (val, (idx[0], idx[1])), shape= (given['nAll']-given['nReal'], given['nAll']))
Htheta = Htheta.T

W = Btheta + Htheta.dot(Htheta.T)
top = sps.hstack((W, -Htheta))
bottom = sps.hstack((-Htheta.T, np.zeros((Htheta.shape[1], Htheta.shape[1]))))
BorderHessian = sps.vstack((top, bottom)).todense()
VCOV = inv(BorderHessian)
VCOV = VCOV #correct for the weighted likelihood
SE = np.sqrt(np.array(VCOV.diagonal())[0,:given['nReal']]).reshape((-1,1))[:len(results[0][:given['nReal']])][:,0]
tstat = results[0][:given['nReal']]/SE
pval = stats.norm.sf(np.abs(tstat))*2
output = np.concatenate((results[0][:given['nReal']], SE,
                    tstat, pval)).reshape((-1, 4), order="F")

Vout = VCOV[:given['nReal'], :given['nReal']]

output0 = output
np.savetxt("../R4.4.1/ipoptTEMP/regtable.csv", output0, delimiter=",")
np.savetxt("../R4.4.1/ipoptTEMP/convergence.csv", np.array([results[5], results[4]]), delimiter=",")
np.savetxt("../R4.4.1/ipoptTEMP/v.csv", np.array(results[0][given['nReal']:]), delimiter=",")
np.savetxt("../R4.4.1/ipoptTEMP/VCOV1.csv", Vout, delimiter=",")






if twostep:
  V2=VCOV
  x0=results[0]
  agamma = adolc.adouble(np.genfromtxt("../R4.4.1/ipoptTEMP/gamma.csv", delimiter=",")[1:])
  gamma = np.genfromtxt("../R4.4.1/ipoptTEMP/gamma.csv", delimiter=",")[1:]
  V1 = np.genfromtxt("../R4.4.1/ipoptTEMP/V1.csv", delimiter=",")[1:,:]
  given['states'] = states
  def LLgamma(gamma, x, Y, given):
    beta = x[:given['nBeta']]
    kappa = x[given['nBeta']:given['nReal']]
    v= x[given['nReal']:]    
    return logLik_parameter(beta, kappa, gamma, v, Y, given)

  adolc.trace_on(3)
  adolc.independent(agamma)
  ay = LLgamma(agamma, x0, Y, given)
  adolc.dependent(ay)
  adolc.trace_off()

  def const_gamma(gamma, x, given):
    beta = x[:given['nBeta']]
    kappa = x[given['nBeta']:given['nReal']]
    v= x[given['nReal']:]    
    return constraint_parameter(beta, kappa, gamma, v,given)
  
  adolc.trace_on(4)
  adolc.independent(agamma)
  ay=const_gamma(agamma, x0, given)
  adolc.dependent(ay)
  adolc.trace_off()
  
  
  def LLtheta(x, gamma, Y, given):
    beta = x[:given['nBeta']]
    kappa = x[given['nBeta']:given['nReal']]
    v= x[given['nReal']:]    
    return logLik_parameter(beta, kappa, gamma, v, Y, given)
  
  
  adolc.trace_on(5)
  ax = adolc.adouble(x0)
  adolc.independent(ax)
  ay = LLtheta(ax, agamma, Y, given)
  adolc.dependent(ay)
  adolc.trace_off()
  
  
  class jac_theta_adolc:
    
    def __init__(self, x):
        #options = np.array([1,1,0,0],dtype=int)
        options = None
        result = adolc.colpack.sparse_jac_no_repeat(5,x,options)
        
        self.nnz  = result[0]     
        self.rind = np.asarray(result[1],dtype=int)
        self.cind = np.asarray(result[2],dtype=int)
        self.values = np.asarray(result[3],dtype=float)
        
    def __call__(self, x, flag, user_data=None):
        if flag:
            return (self.rind, self.cind)
        else:
            result = adolc.colpack.sparse_jac_repeat(5, x, self.nnz, self.rind,
                self.cind, self.values)
            return result[3]

  Jac_theta_adolc = jac_theta_adolc(x0)

    



  idx = Jac_theta_adolc(x0, True)
  val = Jac_theta_adolc(x0, False)
  JLL_theta = sps.coo_matrix( (val, (idx[0], idx[1])), shape= (300,  given['nAll']))

  JLL_gamma = adolc.jacobian(3, gamma)
  JV_gamma  = adolc.jacobian(4, gamma)
  
  Bgamma = JLL_theta.T.dot(JLL_gamma)
  
  Hgamma = Htheta.dot(JV_gamma)
  W = Bgamma + Hgamma
  R = np.diag(const(x0)).dot(JV_gamma) #zeros
  C = np.vstack((W, R))
  


  meat = C.dot(V1).dot(C.T)
  twostep=V2 + V2.dot(meat).dot(V2)
  np.set_printoptions(precision=3, suppress=True)
  
  SE = np.sqrt(np.array(twostep.diagonal())[0,:given['nReal']]).reshape((-1,1))[:len(results[0][:given['nReal']])][:,0]
  tstat = results[0][:given['nReal']]/SE
  pval = stats.norm.sf(np.abs(tstat))
  output = np.concatenate((results[0][:given['nReal']], SE, tstat, pval)).reshape((-1, 4), order="F")
  
  

  np.savetxt("../R4.4.1/ipoptTEMP/regtable2.csv", output, delimiter=",")
  Vout = V2[:given['nReal'], :given['nReal']]
  np.savetxt("../R4.4.1/ipoptTEMP/VCOV2.csv", twostep, delimiter=",")

