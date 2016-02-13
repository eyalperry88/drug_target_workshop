""" Propagting information on the graph """

import time
import numpy as np
import random

# Debug structure used to time the propagation
times = {'prior_time' : 0,
'copy_time' : 0,
'propogate_time' : 0,
'converge_time' : 0}

def propagate(g, prior, EPSILON = 0.0001, ALPHA = 0.9, MAX_ITERATIONS=40, RANDOM_PRIORS=0, KNOCKOUT_IDX = []):
    """ 
    Propagate information on the graph.
    
    If prior == 'MT' - propagation will start from the nodes marked as mutated.
    If prior == 'DE' - propagation will start from the nodes marked as differentially expressed.
    If prior == 'RAND_MT' - propagation will start from RANDOM_PRIORS random nodes.
    
    We stop propagating when we reached MAX_ITERATIONS or convergence is less than EPSILON
    
    The function return a tuple (F, iterations), where F is a vector representing
    the signal on the network post-propagation. iterations is the number of iterations it took.
    """
    global times
    
    n = len(g.nodes)
    
    Y = np.zeros(n)
    count_prop_start = 0
    if prior == 'RAND_MT':
        random_indices = random.sample(range(n), RANDOM_PRIORS)
        for idx in random_indices:
            Y[idx] = (1 - ALPHA)
    else:
        for node in g.nodes:
            if prior == 'GE':
                if g.nodes[node].expression_level != None:
                    Y[g.gene2index[node]] = (1 - ALPHA)
                    count_prop_start += 1
            elif prior == 'MT':
                if g.nodes[node].mutation_type != None:
                    Y[g.gene2index[node]] = (1 - ALPHA)
                    count_prop_start += 1
        print(str(count_prop_start) + ' initial genes out of ' + str(n))

    F = np.zeros(n)
    
    iterations = 1
    while True:
        newF = g.W.dot(ALPHA * F) + Y
        for knock in KNOCKOUT_IDX:
            newF[knock] = 0
        
        summ = np.sum((newF - F)**2)
        F = newF
        
        if (summ < EPSILON or iterations >= MAX_ITERATIONS):
            # print('Converged after ' + str(iterations) + ' iterations')
            break
        iterations += 1

    return F, iterations
