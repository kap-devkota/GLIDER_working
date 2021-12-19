import numpy as np


def individual_k_acc(k_choices,
                     prots,
                     glide_mat,
                     prot_go_map):
    """
    Functions that chooses the best k among the choices listed in `k_choices`. Accuracy 
    is the preferred metric of choice

    k_choices   - [int] of k neighbors
    prots       - [proteins]
    glide_mat   - A_nxn matrix of glide similarity
    prot_go_map - map prot -> [GO]
    """

    sorted_mat = np.argsort(-glide_mat, axis = 1)

    def get_acc(prot, k):
        """
        Given protein, given by `prot` and the number of `k`, returns the 
        accuracy. 
        """
        neighbors   = sorted_mat[prot][:k]
        labels_dict = {}
        for n in neighbors:
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
    return k_acc
