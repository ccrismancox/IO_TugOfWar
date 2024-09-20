import numpy as np
import itertools as it


def genGiven(states, ts, P, nBeta, nKappa, delta, agent=None):
    given = {}
    given['G'] = np.array(list(it.product(list(states),[0,1] , [0,1])))
    given['P'] = P
    given['delta'] = delta
    given['nBeta'] = nBeta
    given['nKappa'] = nKappa
    given['agent']= agent
    
    S = P.shape[1]
    given['nReal'] = given['nBeta'] + given['nKappa']
    given['nAux'] = S*2*2
    given['nAll'] = given['nReal'] + given['nAux']
    
    
    counts = np.hstack((states.reshape(-1,1), np.zeros((S,1))))
    counts2 = np.unique(ts[:,0], return_counts=True)
    counts[np.array([counts[i,0] in counts2[0] for i in range(S)]) ,1] = counts2[1]
    given['statesInData'] = counts[:,1]
    stateActions = np.array(list(it.product(list(states), [0,1])))

    stateActionsInData = np.zeros((S*2*2,1))
    for i in range(stateActions.shape[0]):
      addon = stateActions.shape[0]
      stateActionsInData[i,:] = np.sum(np.logical_and(ts[:,0] == stateActions[i,0], ts[:,1]==stateActions[i,1]))
      stateActionsInData[i+addon,:] = np.sum(np.logical_and(ts[:,0] == stateActions[i,0], ts[:,2]==stateActions[i,1]))
        
      
    given['stateActionsInData'] = stateActionsInData
    given['stateActions'] = stateActions
    return given
