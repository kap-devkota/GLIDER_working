import numpy as np


def individual_k_acc(k_ratios,
                     prots,
                     glide_net,
                     prot_go_map,
                     verbose = False):
    """
    Functions that chooses the best k among the choices listed in `k_ratios`. Accuracy 
    is the preferred metric of choice

    k_choices   - [int] of k neighbors
    prots       - [proteins]
    glide_net   - network formed from GLIDE scores in the algorithm described in the paper. A networkx object
    prot_go_map - map prot -> [GO]
    """
    def log(strng):
        if verbose:
            print(strng)
    
    # Compute the average k-value, and use it to generate the actual k-scores.
    avg_k    = np.average(list(glide_net.degree()))
    k_choices = [int(k * avg_k) for k in k_ratios]

    log(f" Choices of k: {k_choices}")
    
    """
    Get neighbors 
    """
    prot_neighbors = {}
    for p in prots:
        # Get neighbors
        prot_neighbors[p] = sorted(glide_net[p].items(),
                                   reverse = True,
                                   key = lambda edge: edge[1]["weight"])        
    
    def get_acc(prot, k):
        """
        Given a protein given by prot and a k-value, return the accuracy
        """

        labels_dict = {}

        """
        Find the k nearest neighbors.
        """
        
        neighbors   = prot_neighbors[prot]
        if len(neighbors) > k:
            neighbors = neighbors[: k]
        for n, _ in neighbors:
            # Ignore neighbors with no go association available
            if n not in prot_go_map:
                continue
            for go in prot_go_map[n]:
                labels_dict[go] = 1 if go not in prot_go_map else labels_dict[go] + 1

        # If the neighbors did not produce any GO association, return 0
        if len(labels_dict) == 0:
            return 0
        
        best_label = max(labels_dict, key = labels_dict.get)
        """
        For a protein, the accuracy is 1 the best GO label is selected, else 0
        """
        return 1 if best_label in prot_go_map[prot] else 0
    
    k_acc  = []
    counts = 0 
    for k in k_choices:
        acc = 0.0
        for prot in prots:
            # If no go association available, ignore
            if prot not in prot_go_map:
                continue
            counts += 1
            acc += get_acc(prot, k)
        acc /= counts
        k_acc.append(acc)

    """
    Return accuracy of the individual k-values.
    """
    return k_choices, k_acc
