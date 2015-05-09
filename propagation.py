import time
import numpy as np
import random

times = {'prior_time' : 0,
'copy_time' : 0,
'propogate_time' : 0,
'converge_time' : 0}

def propagate(g, prior, EPSILON = 0.0001, ALPHA = 0.9, MAX_ITERATIONS=40, RANDOM_PRIORS=0, KNOCKOUT_IDX = -1):
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
        if (KNOCKOUT_IDX >= 0):
            newF[KNOCKOUT_IDX] = 0
        
        summ = np.sum((newF - F)**2)
        F = newF
        
        if (summ < EPSILON or iterations >= MAX_ITERATIONS):
            print('Converged after ' + str(iterations) + ' iterations')
            break
        iterations += 1

    return F, iterations



def propagate_old(g, prior, EPSILON = 0.0001, ALPHA = 0.75, MAX_ITERATIONS=40):
    global times
    
    prior_knowledge = {}
    
    # fill prior knowledge by GE
    # differentially expressed -> Y = 1
    start = time.clock()
    count_prop_start = 0
    for node in g.nodes:
        if (g.nodes[node].expression_level if prior == 'GE' else g.nodes[node].mutation_type) != None:
            prior_knowledge[node] = (1 - ALPHA)
            count_prop_start += 1
        else:
            prior_knowledge[node] = 0
    times['prior_time'] += (time.clock() - start)

    print(str(count_prop_start) + ' initial genes out of ' + str(len(g.nodes)))

    propagation_scores = {}
    for node in g.nodes:
        propagation_scores[node] = 0

    iterations = 1
    while True:
        start = time.clock()
        new_propagation_scores = propagation_scores.copy()
        times['copy_time'] += (time.clock() - start)
        
        start = time.clock()
        for node in propagation_scores:
            summ = 0
            for neighbor, w in g.nodes[node].neighbors.items():
                summ += (propagation_scores[neighbor] * w)
            new_propagation_scores[node] = ALPHA * summ + prior_knowledge[node]
        times['propogate_time'] += (time.clock() - start)
        
        start = time.clock()
        summ = 0
        for node in propagation_scores:
            summ += (propagation_scores[node] - new_propagation_scores[node]) ** 2
        times['converge_time'] += (time.clock() - start)
        
        start = time.clock()
        propagation_scores = new_propagation_scores.copy()
        times['copy_time'] += (time.clock() - start)
        
        # print('Dist^2 after ' + str(iterations) + ' is ' + str(summ))
        if (summ < EPSILON or iterations >= MAX_ITERATIONS):
            print('Converged after ' + str(iterations) + ' iterations')
            break
        iterations += 1

    return propagation_scores, iterations
