

def propagate(g, prior, ITERATIONS = 10, ALPHA = 0.75):
    prior_knowledge = {}
    
    # fill prior knowledge by GE
    # differentially expressed -> Y = 1
    count_prop_start = 0
    for node in g.nodes:
        if (g.nodes[node].expression_level if prior == 'GE' else g.nodes[node].mutation_type) != None:
            prior_knowledge[node] = 1
            count_prop_start += 1
        else:
            prior_knowledge[node] = 0

    print(str(count_prop_start) + ' initial genes out of ' + str(len(g.nodes)))

    propagation_scores = {}
    for node in g.nodes:
        propagation_scores[node] = 0

    iterations = 1
    while True:
        new_propagation_scores = propagation_scores.copy()
        for node in propagation_scores:
            summ = 0
            for neighbor in g.nodes[node].neighbors:
                summ += (propagation_scores[neighbor] * g.getEdgeWeight(node, neighbor))
            new_propagation_scores[node] = ALPHA * summ + (1 - ALPHA) * prior_knowledge[node]
        
        summ = 0
        for node in propagation_scores:
            summ += (propagation_scores[node] - new_propagation_scores[node]) ** 2
        propagation_scores = new_propagation_scores.copy()
        
        # print('Dist^2 after ' + str(iterations) + ' is ' + str(summ))
        if (summ < 0.0001):
            print('Converged after ' + str(iterations) + ' iterations')
            break
        iterations += 1

    return propagation_scores
