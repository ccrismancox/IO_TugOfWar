'''
Code to fit the model with changing deltas
'''
from genModel import *
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


P = np.genfromtxt("ipoptTEMP_delta/trans.csv", delimiter=",")[1:,:]
Xdata = np.genfromtxt("ipoptTEMP_delta/regData.csv", delimiter=",")[1:,:]
'''
columns in X are
states(continuous); Hattacks; Fattacks; lag.states; states(discrete)
'''
Y = Xdata[:,np.array([4,1,2])]


states =  np.genfromtxt("ipoptTEMP_delta/statespace.csv", delimiter=",")[1:]


'''parameters'''
parameters = np.genfromtxt("ipoptTEMP_delta/params.csv", delimiter=",")[1:,:]
delta = parameters[0][0]
nK = int(parameters[0][1])
nB = int(parameters[0][2])


given = genGiven(states, Y, P, nBeta=nB, nKappa=nK, delta=delta)
  
start =np.genfromtxt("ipoptTEMP_delta/start.csv", delimiter=",")[1:]
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


def  hessLagBuild(x, lagrange, obj_factor):
    idx = hessLag_adolc(x, lagrange, obj_factor, True)
    val = hessLag_adolc(x, lagrange, obj_factor, False)
    return sps.coo_matrix( (val, (idx[0], idx[1])), shape=(x0.shape*2))


t0 = time()
H2 = hessLag_adolc(x0, xL, 1.0, False)
t1 = time() - t0
print("Hessian time: %d seconds" %t1)
nnzh = len(H2)
    
    
          
          
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
nlp2.num_option('max_cpu_time', 1e5)
nlp2.num_option('tol', 1e-4)
nlp2.str_option('warm_start_init_point', 'yes')
if delta < .999:
  nlp2.int_option('max_iter', 5000)
else:
  nlp2.int_option('max_iter', 50)

   
t0 = time()
results = nlp2.solve(x0)
t1 = time()-t0
print("Estimation time: %d seconds" %t1)

# free the model
nlp2.close()
print(results[0][:given['nReal']])
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
print(np.round(output,2))


np.savetxt("ipoptTEMP_delta/regtable.csv", output, delimiter=",")
np.savetxt("ipoptTEMP_delta/convergence.csv", np.array([results[5], results[4]]), delimiter=",")
np.savetxt("ipoptTEMP_delta/v.csv", np.array(results[0][given['nReal']:]), delimiter=",")








