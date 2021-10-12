from sage.modules.free_module_integer import IntegerLattice

def help():
    print("---Tool list---")
    print("lattice_attack(A, B, C, p, lb, ub) -> vector : (Simple explanation) solve A_1 * x_1 + ... + A_n * x_n + B = C * y mod p")


# Directly taken from rbtree's LLL repository
# From https://oddcoder.com/LOL-34c3/, https://hackmd.io/@hakatashi/B1OM7HFVI
def Babai_CVP(mat, target):
    M = IntegerLattice(mat, lll_reduce=True).reduced_basis
    G = M.gram_schmidt()[0]
    diff = target
    for i in reversed(range(G.nrows())):
        diff -=  M[i] * ((diff * G[i]) / (G[i] * G[i])).round()
    return target - diff

def solve(mat, lb, ub, weight = None):
    num_var  = mat.nrows()
    num_ineq = mat.ncols()

    max_element = 0 
    for i in range(num_var):
        for j in range(num_ineq):
            max_element = max(max_element, abs(mat[i, j]))

    if weight == None:
        weight = num_ineq * max_element

    # sanity checker
    if len(lb) != num_ineq:
        print("Fail: len(lb) != num_ineq")
        return

    if len(ub) != num_ineq:
        print("Fail: len(ub) != num_ineq")
        return

    for i in range(num_ineq):
        if lb[i] > ub[i]:
            print("Fail: lb[i] > ub[i] at index", i)
            return

        # heuristic for number of solutions
    DET = 0

    if num_var == num_ineq:
        DET = abs(mat.det())
        num_sol = 1
        for i in range(num_ineq):
            num_sol *= (ub[i] - lb[i])
        if DET == 0:
            print("Zero Determinant")
        else:
            num_sol //= DET
            # + 1 added in for the sake of not making it zero...
            print("Expected Number of Solutions : ", num_sol + 1)
            
    # scaling process begins
    max_diff = max([ub[i] - lb[i] for i in range(num_ineq)])
    applied_weights = []

    for i in range(num_ineq):
        ineq_weight = weight if lb[i] == ub[i] else max_diff // (ub[i] - lb[i])
        applied_weights.append(ineq_weight)
        for j in range(num_var):
            mat[j, i] *= ineq_weight
        lb[i] *= ineq_weight
        ub[i] *= ineq_weight

    # Solve CVP
    target = vector([(lb[i] + ub[i]) // 2 for i in range(num_ineq)])
    result = Babai_CVP(mat, target)

    for i in range(num_ineq):
        if (lb[i] <= result[i] <= ub[i]) == False:
            print("Fail : inequality does not hold after solving")
            break
    
        # recover x
    fin = None

    if DET != 0:
        mat = mat.transpose()
        fin = mat.solve_right(result)
    
    ## recover your result
    return result, applied_weights, fin

# n coefficients, m individual formula
# len(A) = m, len(A[0]) = n
def lattice_attack(A, B, C, p, lb, ub):
    assert len(A) == len(B) and len(B) == len(C)
    m, n = len(A), len(A[0])
    assert len(lb) == n + 1
    mtrx = matrix(ZZ, (n+1)*m+2, (n+1)*m+2)
    for i in range(n*m):
        mtrx[i, i] = 1
    for i in range(m):
        for j in range(n):
            mtrx[i*n+j, n*m+i] = A[i][j]
        mtrx[n*m+i, n*m+i] = -p
        mtrx[(n+1)*m, n*m+i] = -C[i]
        mtrx[(n+1)*m+1, n*m+i] = B[i]
    mtrx[(n+1)*m, (n+1)*m] = 1
    mtrx[(n+1)*m+1, (n+1)*m+1] = 1
    _lb = lb[:-1] * m + [0] * m + [lb[-1]] + [1]
    _ub = ub[:-1] * m + [0] * m + [ub[-1]] + [1]
    res = solve(mtrx, _lb, _ub)
    return res[2]