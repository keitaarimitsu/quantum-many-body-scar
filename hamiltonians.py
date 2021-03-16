def pxp_1d(n_site: int, state_vecs):
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
