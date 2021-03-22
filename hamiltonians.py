def pxp_1d(n_site: int, state_vecs: "np.array"):
  """create PXP Hamiltonian in 1 dimensional lattice.
  
  Args:
    n_site:
      The number of sites. For 1 dimensional case, this is the same as the number of spins.
    state_vecs (np.array):
      A certain subspace of the whole Hamiltonian. For example, you can restrict the Hilbert space by discarding "dark states".
      Each element of it should be a spin configuration such as "110101010101".
      
  Return:
    pxp (scipy.sparse.csc_matrix of (# of state_vecs, # of state_vecs)):
      pxp[i][j] corresponds to a matrix element of the PXP model, i.e., <i|H|j> where |i> = state_vecs[i]      
  """
  pxp = np.zeros((len(state_vecs), len(state_vecs)))
  
  if n_site % 2 == 1:
      raise TypeError("the number of particles must be even.")
    
  else:
      for i in range(len(state_vecs)):
          tmp_str = state_vecs[i]
          for k in range(n_site):
              cond_str = "".join([tmp_str[(k - 1) % n_site], tmp_str[(k + 1) % n_site]])
              if cond_str == "11":
                  tmp_str_list = list(tmp_str)
                  if tmp_str_list[k] == "0":
                      tmp_str_list[k] = "1"
                  else:
                      tmp_str_list[k] = "0"
                  act_str = "".join(tmp_str_list)
                  
                  idx = state_vecs.index(act_str)
                  pxp[idx][i] += 1
      return csr_matrix(pxp)  
    
def perturbation_1d(
    n_site: int,
    state_vecs: "np.array",
    cpl_const: float
) -> "scipy.sparse.csr_matrix":
    """create perturbations of 1d PXP model
    Args:
      n_site (int): the number of sites.
      state_vecs (np.array)
      cpl_const (float): the coupling constant 
    Returns:
      perturb (scipy.sparse.csr_matrix): H_pert = cpl_const * Σ (PXPP + PPXP)
    """
    perturb = np.zeros((len(state_vecs), len(state_vecs)), dtype=float32)
    for state_idx, state_str in enumerate(state_vecs):
        for k in range(len(state_str)):
            cond_str = "".join([state_str[(k - i) % n_site] for i in range(-2, 2, 1)])
            if cond_str == "1111":
                act_str1_list = list(state_str)
                act_str2_list = list(state_str)
                act_str1_list[(k - 1) % n_site] = "0"
                act_str2_list[k] = "0"
                act_str1 = "".join(act_str1_list)
                act_str2 = "".join(act_str2_list)
                
                idx1 = state_vecs.index(act_str1)
                idx2 = state_vecs.index(act_str2)
                
                perturb[idx1][state_idx] += 1
                perturb[idx2][state_idx] += 1
            elif cond_str == "1011":
                act_str_list = list(state_str)
                act_str_list[(k - 1) % n_site] = "1"
                act_str = "".join(act_str_list)
                
                idx = state_vecs.index(act_str)
                
                perturb[idx][state_idx] += 1
            elif cond_str == "1101":
                act_str_list = list(state_str)
                act_str_list[k] = "1"
                act_str = "".join(act_str_list)
                
                idx = state_vecs.index(act_str)
                
                perturb[idx][state_idx] += 1
                
    return csr_matrix(cpl_const * perturb)
        
