import numpy as np

def glide_generate_graph(k,
                         glide_mat,
                         id_symbol_map):
    """
    Generate a new graph with each node having degree `k` using the glide scores represented in the 
    `glide_mat` matrix.
    """

    sorted_mat = np.argsort(-glide_mat, axis = 1)

    edges      = []
    n_nodes    = glide_mat.shape[0]

    """
    For each protein id, return a list of neighbors and add it 
    to the edgelist
    """
    for i in range(n_nodes):
        neighbors = sorted_mat[i][: k]
        edges  += [(id_symbol_map[i], id_symbol_map[n], glide_mat[i, n]) for n in neighbors]
    return edges
    
