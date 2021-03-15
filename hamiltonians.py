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
              if k == 0:
                  if (tmp_str[k + 1] == "1") and (tmp_str[int(n_site / 2 - 1)] == "1"):
                      if tmp_str[k] == "0":
                          aft_str = tmp_str[k + 1:]
                          act_str = "1" + aft_str
                      else:
                          aft_str = tmp_str[k + 1:]
                          act_str = "0" + aft_str
                          
                      idx = state_vecs.index(act_str)
                      pxp[idx][i] += 1
              elif k == int(n_site / 2 -1):
                  if (tmp_str[k - 1] == "1") and (tmp_str[0] == "1"):
                      if tmp_str[k] == "0":
                          bef_str = tmp_str[0:k]
                          act_str = bef_str + "1"
                      else:
                          bef_str = tmp_str[0:k]
                          act_str = bef_str + "0"
                          
                      idx = state_vecs.index(act_str)
                      pxp[idx][i] += 1
              else:
                  if (tmp_str[k - 1] == "1") and (tmp_str[k + 1] == "1"):
                      if tmp_str[k] == "0":
                          bef_str = tmp_str[0:k]
                          aft_str = tmp_str[k + 1:]
                          act_str = bef_str + "1" + aft_str
                      else:
                          bef_str = tmp_str[0:k]
                          aft_str = tmp_str[k + 1:]
                          act_str = bef_str + "0" + aft_str
                              
                      idx = state_vec.index(act_str)
                      pxp[idx][i] += 1
      return pxp  
