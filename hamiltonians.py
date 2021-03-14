def pxp_1d(n_site, state_vecs):
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
  return pxp
