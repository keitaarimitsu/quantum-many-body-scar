def constraint_lowering1d(
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
            position of 2nd island
    Return:
        operator (scipy.sparse.csc_matrix): a matrix corresponding a_[pos1]^dag a_[pos2]^dag
    """
    operator = np.zeros((len(state_vecs), len(state_vecs)), dtype="int32")
    
    for i in range(len(state_vecs)):
        tmp_str = state_vecs[i]
        cond_str_pos1 = "".join([tmp_str[(pos1 - 1) % n_site], tmp_str[(pos1 + 1) % n_site]])
        cond_str_pos2 = "".join([tmp_str[(pos2 - 1) % n_site], tmp_str[(pos2 + 1) % n_site]])
        cond_str = cond_str_pos1 + cond_str_pos2
        if (cond_str == "1111") and (tmp_str[pos1] == "0") and (tmp_str[pos2] == "0"):
            tmp_str_list = list(tmp_str)
            tmp_str_list[pos1] = "1"
            tmp_str_list[pos2] = "1"
            act_str = "".join(tmp_str_list)
            idx = state_vecs.index(act_str)
            operator[idx][i] += 1
    return csc_matrix(operator)
        

def constraint_bound_creation1d(
    n_site: int,
    state_vecs: "np.array",
    pos1: int,
    pos2: int
) -> "scipy.sparse.csc_matrix":
     """construct a pair of an island @ pos1 and the bound state @ pos2, i.e., a_[pos1]^dag γ_[pos2]^dag
    Args: 
        n_site (int):
            the number of spins
        state_vecs (np.array):
            1d array which represents states under consideration, 
            i.e., state_vecs[i] = "10101110101011..." (ith configuration)
        pos1 (int):
            position of an island
        pos2 (int):
            position of a bound state
    Return:
        operator (scipy.sparse.csc_matrix): a matrix corresponding a_[pos1]^dag γ_[pos2]^dag
    """
    operator = np.zeros((len(state_vecs), len(state_vecs)), dtype="int32")
    for i in range(len(state_vecs)):
        tmp_str = state_vecs[i]
        cond_str_pos1 = "".join([tmp_str[(pos1 + k) % n_site] for k in range(-1, 2, 1)])
        cond_str_pos2 = "".join([tmp_str[(pos2 + k) % n_site] for k in range(-2, 3, 1)])
        cond_str = cond_str_pos1 + cond_str_pos2
        if cond_str == "10111111":
            tmp_str_list = list(tmp_str)
            tmp_str_list[pos1] = "1"
            tmp_str_list[pos2] = "0"
            act_str = "".join(tmp_str_list)
            idx = state_vecs.index(act_str)
            operator[idx][i] += 1
    return operator


def one_body_island_creation1d(
    n_site: int,
    state_vecs: "np.array",
    pos: int
) -> "scipy.sparse.csc_matrix":
    operator = np.zeros((len(state_vecs), len(state_vecs)), dtype="int32")
    for i in range(len(state_vecs)):
        tmp_str = state_vecs[i]
        cond_str = "".join([tmp_str[(pos + k) % n_site] for k in range(-1, 2, 1)])
        if cond_str == "101":
            tmp_str_list  = list(tmp_str)
            tmp_str_list[pos] = "1"
            act_str = "".join(tmp_str_list)
            idx = state_vecs.index(act_str)
            operator[idx][i] += 1
    return operator


def one_body_bound_state_creation1d(
    n_site: int,
    state_vecs: "np.array",
    pos: int,
) -> "scipy.sparse.csc_matrix":
    """construct a bound state in 1d pxp chain, i.e., P_[pos-2]P_[pos-1]s^+_[pos]P_[pos+1]P_[pos+2]
    Args:
        n_site (int): 
            the number of spins
        state_vecs (np.array):
            1d array represents states
        pos (int):
            position of the bound state
    Return:
        operator (scipy.sparse.csc_matrix):
            P_[pos-2]P_[pos-1]s^+_[pos]P_[pos+1]P_[pos+2]
    """
    operator = np.zeros((len(state_vecs), len(state_vecs)), dtype="int32")
    for i in range(len(state_vecs)):
        tmp_str = state_vecs[i]
        cond_str = "".join([tmp_str[(pos + k) % n_site] for k in range(-2, 3, 1)]) 
        if cond_str == "11111":
            tmp_str_list = list(tmp_str)
            tmp_str_list[pos] = "0"
            act_str = "".join(tmp_str_list)
            idx = state_vecs.index(act_str)
            operator[idx][i] += 1
    return csc_matrix(operator)       


def spectrum_generating_algebra1d(
    n_site: int,
    state_vecs: "np.array"
) -> "scipy.sparse.csc_matrix":
    operator = np.zeros((len(state_vecs), len(state_vecs)), dtype="float32")
    for idx, state in enumerate(state_vecs):
        for k in range(len(state)):
            cond_str = "".join([state[(k + l) % n_site] for l in range(-1, 2, 1)])
            if cond_str == "111":
                act_str_list = list(state)
                act_str_list[k] = "0"
                act_str = "".join(act_str_list)
                act_idx = state_vecs.index(act_str)
                operator[idx][idx] += (-1) ** k * (-1)
                operator[act_idx][idx] += (-1) ** k * (-1) / np.sqrt(2)
            elif cond_str == "101":
                act_str_list = list(state)
                act_str_list[k] = "1"
                act_str = "".join(act_str_list)
                act_idx = state_vecs.index(act_str)
                operator[idx][idx] += (-1) ** k
                operator[act_idx][idx] += (-1) ** k / np.sqrt(2)
    return csc_matrix(operator)
