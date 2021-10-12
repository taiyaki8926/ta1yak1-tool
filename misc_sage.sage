from sage.crypto.util import carmichael_lambda
import sys

def help():
    print("---Tool list---")
    print("pohlig_mul(g, h, p) -> x : DLP solver with multiplicative cyclic group")
    
def pohlig_mul(g, h, p):
    def bsgs_mul(g, h, n, p):
        assert pow(g, n, p) == 1
        m = int(sqrt(n))
        if pow(m, 2) != n:
            m += 1
        ls = [(pow(g, j, p)) for j in range(m)]
        gamma = h
        for i in range(m):
            if gamma in ls:
                ret = i * m + ls.index(gamma)
                assert pow(g, ret, p) == h
                return ret
            else:
                gamma *= pow(g, -m, p)
        print('[+] Not Found')
        assert 1 == 2
        return -1

    def pohlig_prime_mul(g, h, _p, e, p):
        n = pow(_p, e)
        assert pow(g, n, p) == 1
        x = 0
        gamma = pow(g, pow(_p, e-1), p)
        for k in range(e):
            hk = pow(pow(g, -x, p) * h % p, pow(_p, e-1-k), p)
            dk = bsgs_mul(gamma, hk, _p, p)
            x += pow(_p, k) * dk
        assert pow(g, x, p) == h
        print('[+] Found discrete log {} in modulo {}^{}'.format(x, _p, e))    
        return x

    def calc_order(g, p):
        n = carmichael_lambda(p)
        assert pow(g, n, p) == 1
        f = list(factor(n))
        div = 1
        for _p, e in f:
            for k in range(1, e+1):
                if pow(g, n//(pow(_p, k)), p) == 1:
                    div *= _p
                else:
                    break
        return n // div

    n = calc_order(g, p)
    assert pow(g, n, p) == 1
    xi_ls = []
    mod_ls = []
    f = factor(n)
    print('[+] factorize :',f)
    for _p, e in list(f):
        gi = pow(g, n//pow(_p, e), p)
        hi = pow(h, n//pow(_p, e), p)
        xi_ls.append(pohlig_prime_mul(gi, hi, _p, e, p))
        mod_ls.append(pow(_p, e))
        ret = crt(xi_ls, mod_ls)
        if pow(g, ret, p) == h:
            print('[+] Final result:', ret)
            return ret
    print('[+] Not Found')
    sys.exit()
    return -1