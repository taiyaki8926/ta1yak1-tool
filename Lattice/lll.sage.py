

# This file was *autogenerated* from the file ./lll.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_0p75 = RealNumber('0.75'); _sage_const_0xffffffffffffffffffffffffffffffff7fffffff = Integer(0xffffffffffffffffffffffffffffffff7fffffff); _sage_const_0xffffffffffffffffffffffffffffffff7ffffffc = Integer(0xffffffffffffffffffffffffffffffff7ffffffc); _sage_const_0x1c97befc54bd7a8b65acf89f81d4d4adc565fa45 = Integer(0x1c97befc54bd7a8b65acf89f81d4d4adc565fa45); _sage_const_0x4a96b5688ef573284664698968c38bb913cbfc82 = Integer(0x4a96b5688ef573284664698968c38bb913cbfc82); _sage_const_0x23a628553168947d59dcc912042351377ac5fb32 = Integer(0x23a628553168947d59dcc912042351377ac5fb32); _sage_const_0x0100000000000000000001f4c8f927aed3ca752257 = Integer(0x0100000000000000000001f4c8f927aed3ca752257); _sage_const_0x01 = Integer(0x01); _sage_const_16 = Integer(16); _sage_const_68570986338040913473862653949389814867104146286 = Integer(68570986338040913473862653949389814867104146286); _sage_const_32 = Integer(32); _sage_const_2359071867217229033899952554549558255586216332644160293624620302437641015683071171760925782397 = Integer(2359071867217229033899952554549558255586216332644160293624620302437641015683071171760925782397); _sage_const_5 = Integer(5); _sage_const_350 = Integer(350); _sage_const_512 = Integer(512); _sage_const_1461501637330902918203684832716283019655932542983 = Integer(1461501637330902918203684832716283019655932542983); _sage_const_10 = Integer(10); _sage_const_22 = Integer(22); _sage_const_151 = Integer(151)# From Furutsuki (https://furutsuki.hatenablog.com/entry/2021/03/21/101039)
# A_i + x_i = y * B_i (mod p)
# l := x_i.bit_length()
def lll_solve_1(a, b, p, l, debug=False):
    print("[+] ----------LLL SOLVE 1----------")
    N = len(a)
    nonce_max = _sage_const_2 **l
    assert(len(a) == len(b))
    B = matrix(QQ, N + _sage_const_2 , N + _sage_const_2 )
    B.set_block(_sage_const_0   , _sage_const_0 , matrix.identity(N) * p)
    B.set_block(N  , _sage_const_0 , matrix(QQ, _sage_const_1 , N, b))
    B.set_block(N+_sage_const_1 , _sage_const_0 , matrix(QQ, _sage_const_1 , N, a))
    if not debug:
        B[N,N] = nonce_max / p
        B[N+_sage_const_1 ,N+_sage_const_1 ] = nonce_max
        L = B.LLL()

        flag = False
        for row in list(L):
            x = [int(abs(i)) for i in row]
            #print("[+] candidate: ", x[0])
            if _sage_const_0  < x[_sage_const_0 ] < nonce_max:
                y_ls = [(a[i] + x[i]) * pow(b[i], -_sage_const_1 , p) % p for i in range(len(a))]
                y = (a[_sage_const_0 ] + x[_sage_const_0 ]) * pow(b[_sage_const_0 ], -_sage_const_1 , p) % p
                _y = (a[_sage_const_1 ] + x[_sage_const_1 ]) * pow(b[_sage_const_1 ], -_sage_const_1 , p) % p
                if y == _y:
                    #print("[+] x_i: ", x)
                    print("[+] y  : ", y_ls)
                    #return x, y
                    flag = True
        if not flag:
            print("[+] Not found")
        return _sage_const_0 , _sage_const_0 
    else:
        for i in range(_sage_const_2 *l):
            # set various B[N, N] parameter
            B[N,N] = (_sage_const_2 **i) / p
            B[N+_sage_const_1 ,N+_sage_const_1 ] = nonce_max
            L = B.LLL()
            for row in list(L):
                x = [int(abs(i)) for i in row]
                if _sage_const_0  < x[_sage_const_0 ] < nonce_max:
                    y = (a[_sage_const_0 ] + x[_sage_const_0 ]) * pow(b[_sage_const_0 ], -_sage_const_1 , p) % p
                    _y = (a[_sage_const_1 ] + x[_sage_const_1 ]) * pow(b[_sage_const_1 ], -_sage_const_1 , p) % p
                    if y == _y:
                        print("[+] parameter, y: ", i, y)
        print('[+] Debug done')
        return _sage_const_0 , _sage_const_0 

# From Yoshiking (https://yoshiking.hatenablog.jp/entry/2020/12/27/100446)
def lll_solve_2(a, b, p, l, verbose=False):
    print("[+] ----------LLL SOLVE 2----------")

    def solve_cvp(B, t, verbose=False):
        t_ = t - B.stack(t).gram_schmidt()[_sage_const_0 ].row(-_sage_const_1 )
        if verbose:
            print("Target vector projection:")
            print(numerical_approx(t_, digits=_sage_const_4 ))
        B_ = B.LLL()
        if verbose:
            print("\nLLL-reduced basis:")
            print(numerical_approx(B_, digits=_sage_const_4 ))
        c = B_.solve_left(t_)
        c_ = vector(map(round, c))
        if verbose:
            print("\nRound-off errors:")
            print(numerical_approx(vector(map(abs, c - c_)), digits=_sage_const_4 ))
        return c_ * B_
    
    N = len(a)
    assert(len(a) == len(b))
    B = matrix(QQ, N + _sage_const_1 , N + _sage_const_1 )
    B.set_block(_sage_const_0 , _sage_const_0 , matrix(QQ, _sage_const_1 , N, b))
    B.set_block(_sage_const_1 , _sage_const_0 , matrix.identity(N) * p)
    B[_sage_const_0 , N] = (_sage_const_2 **l) / p
    t = vector(QQ, a + [_sage_const_0 ])
    ans = solve_cvp(B, t, verbose)
    #print("[+] x_i:", ans[:-1])
    print("[+] y  : ", (ans[-_sage_const_1 ] * p / (_sage_const_2 **l)) % p)
    return ans[:-_sage_const_1 ], ans[-_sage_const_1 ]

# From CTF-Wiki (https://ctf-wiki.org/en/crypto/asymmetric/lattice/cvp/)
# Unlike lll_solve_1, the definition of l is "the number of bits that already known"
# i.e., p.bit_length() - x_i.bit_length()

def lll_solve_3(a, b, p, l):
    print("[+] ----------LLL SOLVE 3----------")

    # http://www.isg.rhul.ac.uk/~sdg/igor-slides.pdf
    def babai(A, w):
        A = A.LLL(delta=_sage_const_0p75 )
        G = A.gram_schmidt()[_sage_const_0 ]
        t = w
        for i in reversed(range(A.nrows())):
            c = ((t * G[i]) / (G[i] * G[i])).round()
            t -= A[i] * c
        return w - t

    N = len(a)
    assert(len(a) == len(b))
    k = len(bin(p))-_sage_const_2 -l
    B = Matrix(QQ, N+_sage_const_1 , N+_sage_const_1 )
    B.set_block(_sage_const_0 , _sage_const_0 , matrix.identity(N) * p)
    B.set_block(N, _sage_const_0 , matrix(QQ, _sage_const_1 , N, b))
    B[N, N] = _sage_const_1  / (_sage_const_2  ** (k + _sage_const_1 ))
    closest = babai(B, vector(QQ, a + [_sage_const_0 ]))
    print("[+] y  : ", (closest[-_sage_const_1 ] * (_sage_const_2  ** (k + _sage_const_1 ))) % p)
    return (closest[-_sage_const_1 ] * (_sage_const_2  ** (k + _sage_const_1 ))) % p

# From Furutsuki (https://furutsuki.hatenablog.com/entry/2021/03/21/101039)
# LINE CTF 2021 babycrypto4
def test_1():
    # From https://neuromancer.sk/std/secg/secp160r1
    p = _sage_const_0xffffffffffffffffffffffffffffffff7fffffff 
    K = GF(p)
    a = K(_sage_const_0xffffffffffffffffffffffffffffffff7ffffffc )
    b = K(_sage_const_0x1c97befc54bd7a8b65acf89f81d4d4adc565fa45 )
    E = EllipticCurve(K, (a, b))
    G = E(_sage_const_0x4a96b5688ef573284664698968c38bb913cbfc82 , _sage_const_0x23a628553168947d59dcc912042351377ac5fb32 )
    E.set_order(_sage_const_0x0100000000000000000001f4c8f927aed3ca752257  * _sage_const_0x01 )

    n = int(E.order())

    dataset = []
    with open("line_ctf_babycrypto4_output.txt") as f:
        lines = f.read().strip().split("\n")
        for l in lines:
            dataset.append( [int(x, _sage_const_16 ) for x in l.split(" ")] )
    # Given R_i, S_i, H_i and n such that
    # -inv(S_i) * H_i + k_i = inv(S_i) * R_i * d (mod n)
    # R, S, H : 160 bit
    # n       : 161 bit
    # k_i     : 32  bit
    # d       : 161 bit
    r, s, k, h = zip(*dataset)
    N = len(r)
    a = [-pow(s[i], -_sage_const_1 , n) * h[i] % n for i in range(N)]
    b = [pow(s[i], -_sage_const_1 , n) * r[i] % n for i in range(N)]

    # find this
    d = _sage_const_68570986338040913473862653949389814867104146286 
    print("[+] ans: ", d)
    return a, b, n, _sage_const_32 

# From Furutsuki (https://furutsuki.hatenablog.com/entry/2020/12/27/160652#Crypto-Wilhelmina-Says)
# Harekaze mini CTF 2020 Wilhelmina Says
def test_2():
    from Crypto.Util.number import getStrongPrime
    import random

    # Given X_i, Z_s, p such that
    # Z_i + k_i = X_i * flag (mod p)
    # Z, X, n : 512 bit
    # k_i     : 350 bit
    # flag    : 311 bit
    flag = _sage_const_2359071867217229033899952554549558255586216332644160293624620302437641015683071171760925782397 
    n = _sage_const_5 
    k = _sage_const_350 
    p = getStrongPrime(_sage_const_512 )
    xs = [random.randint(_sage_const_2 , p-_sage_const_1 ) for _ in range(n)]
    ys = [x * flag % p for x in xs]
    zs = [(y >> k) << k for y in ys]

    # find this
    print("[+] ans: ", flag)
    return zs, xs, p, k

# From CTF-Wiki (https://ctf-wiki.org/en/crypto/asymmetric/lattice/cvp/)
# BCTF 2018 guess number
def test_3():
    import random
    
    def msb(k, x, p):
        delta = p >> (k + _sage_const_1 )
        ui = random.randint(x - delta, x + delta)
        return ui

    # p := gmpy2.next_prime(2**160)
    p = _sage_const_1461501637330902918203684832716283019655932542983 

    # Given T_i, U_s, p such that
    # U_i + k_i = T_i * alpha (mod p)
    # U, T, p : 161 bit
    # k_i     : 150 bit
    # alpha   : 161 bit
    alpha = random.randint(_sage_const_1 , p - _sage_const_1 )
    t = []
    u = []
    k = _sage_const_10 
    for i in range(_sage_const_22 ):
        t.append(random.randint(_sage_const_1 , p - _sage_const_1 ))
        u.append(msb(k, alpha * t[i] % p, p))

    # find this
    print("[+] ans: ", alpha)
    return u, t, p, _sage_const_151 


def main():
    for _ in range(_sage_const_2 ):
        a, b, p, l = test_1()
        lll_solve_1(a, b, p, l)
        lll_solve_2(a, b, p, l)
        lll_solve_3(a, b, p, l)

if __name__ == "__main__":
    main()


