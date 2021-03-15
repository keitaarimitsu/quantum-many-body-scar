def constraint_lowering_ordered1d(
    n_site: int,
    state_vecs: "np.array",
    pos1: int,
    pos2: int,
) -> "scipy.sparse.csc_matrix":
    """construct a pair of islands creation operator, i.e., a_[pos1]^dag a_[pos2]^dag
    Args: 
        n_site (int):
            the number of spins
        state_vecs (np.array):
            1d array which represents states under consideration, 
            i.e., state_vecs[i] = "10101110101011..." (ith configuration)
        pos1 (int):
            position of 1st island
        pos2 (int):
            position of 2nd island, which must be bigger than pos1.
    Return:
        operator (scipy.sparse.csc_matrix): a matrix corresponding a_[pos1]^dag a_[pos2]^dag
    """
    operator = np.zeros((len(state_vecs), len(state_vecs)), dtype="int32")
    if pos1 => pos2:
        raise TypeError("pos1 must be less than pos2")
    else:
        for i in range(len(state_vecs)):
            tmp_str = state_vecs[i]
            if pos1 == 0:
                if pos2 == 2:
                    if (tmp_str[n_site - 1] == "1") and (tmp_str[pos1 + 1] == "1") and (tmp_str[pos2 + 1] == "1"):
                        if (tmp_str[pos1] == "0") and (tmp_str[pos2] == "0"):
                            aft_str = tmp_str[pos2 + 1:]
                            act_str = "111" + aft_str
                            idx = state_vecs.index(act_str)
                            operator[idx][i] += 1
                elif pos2 == n_site - 2:
                    if (tmp_str[n_site - 3] == "1") and (tmp_str[n_site - 1] == "1") and (tmp_str[pos1 + 1] == "1"):
                        if (tmp_str[pos1] == "0") and (tmp_str[pos2] == "0"):
                            aft_str = tmp_str[pos1 + 1: pos2]
                            act_str = "1" + aft_str + "11"
                            idx = state_vecs.index(act_str)
                            operator[idx][i] += 1
                else:
                    if (tmp_str[n_site - 1] == "1") and (tmp_str[pos1 + 1] == "1") and (tmp_str[pos2 - 1] == "1") and (tmp_str[pos2 + 1] == "1"):
                        if (tmp_str[pos1] == "0") and (tmp_str[pos2] == "0"):
                            aft
        
