from creation_operator import (
    constraint_lowering1d,
    one_body_bound_state_creation1d
)
from scipy.linalg import expm

def vacuum_state1d(n_site: int, state_vecs: "np.array"):
    """prepare the vacuum state |φ>=1/sqrt(2) * (|Z2>+(-1)^N|Z2'>) for 1d pxp model
    Args:
      n_site (int): the number of sites (for spin 1/2)
      state_vecs ("np.array"): np.arrays representing states under consideration
    Returns:
      the vacuum state which is np.array
    """
    if n_site % 2 == 1:
        raise TypeError("only even number of sites is allowed!")
    else:
        all_l = "01" * int(n_site / 2)
        all_r = "10" * int(n_site / 2)
        idx_l = state_vecs.index(all_l)
        idx_r = state_vecs.index(all_r)
        vacuum = np.zeros(len(state_vecs), dtype="float32")
        vacuum[idx_l] = 1 / np.sqrt(2)
        vacuum[idx_r] = (-1) ** int(n_site / 2) /np.sqrt(2)
        return vacuum


def condensate_ansatz1d(n_site: int, state_vecs: "np.array"):
    """prepare the (zeroth order) ansatz |ψ^0_scar> for the zero energy scarred state for 1d pxp model
    Return:
      np.array representing |ψ^0_scar>
    """
    vacuum = vacuum_state1d(n_site, state_vecs)
    condensate = np.zeros((len(state_vecs), len(state_vecs)), dtype="float32")
    for k in range(int(n_site / 2)):
        for l in range(int(n_site / 2)):
            if k > l:
                condensate += constraint_lowering1d(n_site, state_vecs, k, l)
    condensate = - 1 / (int(n_site / 2) - 1) * condensate
    condensate_ansatz = expm(condensate).dot(vacuum)
    return condensate_ansatz
  

def condensate_with_bound_state1d(n_site: int, state_vecs: "np.array"):
    """prepare the (first order) ansatz |ψ^1_scar> by adding bound state operators
    Return:
      np.array representing |ψ^1_scar>  
    """
    
