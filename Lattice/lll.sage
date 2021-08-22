
# From Yoshiking (https://yoshiking.hatenablog.jp/entry/2020/12/27/100446)
def lll_solve(a, b, p, l, verbose=False):

    def solve_cvp(B, t, verbose=False):
        t_ = t - B.stack(t).gram_schmidt()[0].row(-1)
        if verbose:
            print("Target vector projection:")
            print(numerical_approx(t_, digits=4))
        B_ = B.LLL()
        if verbose:
            print("\nLLL-reduced basis:")
            print(numerical_approx(B_, digits=4))
        c = B_.solve_left(t_)
        c_ = vector(map(round, c))
        if verbose:
            print("\nRound-off errors:")
            print(numerical_approx(vector(map(abs, c - c_)), digits=4))
        return c_ * B_
    
    N = len(a)
    assert(len(a) == len(b))
    B = matrix(QQ, N + 1, N + 1)
    B.set_block(0, 0, matrix(QQ, 1, N, b))
    B.set_block(1, 0, matrix.identity(N) * p)
    B[0, N] = (2^l) / p
    t = vector(QQ, a + [0])
    ans = solve_cvp(B, t, verbose)
    #print("[+] x_i:", ans[:-1])
    print("[+] y  : ", (ans[-1] * p / (2^l)) % p)
    return (ans[-1] * p / (2^l)) % p)

# From Furutsuki (https://furutsuki.hatenablog.com/entry/2021/03/21/101039)
# A_i + x_i = y * B_i (mod p)
# l := x_i.bit_length()
def lll_solve_2(a, b, p, l, debug=False):
    print("[+] ----------LLL SOLVE 2----------")
    N = len(a)
    nonce_max = 2^l
    assert(len(a) == len(b))
    B = matrix(QQ, N + 2, N + 2)
    B.set_block(0  , 0, matrix.identity(N) * p)
    B.set_block(N  , 0, matrix(QQ, 1, N, b))
    B.set_block(N+1, 0, matrix(QQ, 1, N, a))
    if not debug:
        B[N,N] = nonce_max / p
        B[N+1,N+1] = nonce_max
        L = B.LLL()

        flag = False
        for row in list(L):
            x = [int(abs(i)) for i in row]
            #print("[+] candidate: ", x[0])
            if 0 < x[0] < nonce_max:
                y_ls = [(a[i] + x[i]) * pow(b[i], -1, p) % p for i in range(len(a))]
                y = (a[0] + x[0]) * pow(b[0], -1, p) % p
                _y = (a[1] + x[1]) * pow(b[1], -1, p) % p
                if y == _y:
                    #print("[+] x_i: ", x)
                    print("[+] y  : ", y_ls)
                    #return x, y
                    flag = True
        if not flag:
            print("[+] Not found")
        return 0, 0
    else:
        for i in range(2*l):
            # set various B[N, N] parameter
            B[N,N] = (2^i) / p
            B[N+1,N+1] = nonce_max
            L = B.LLL()
            for row in list(L):
                x = [int(abs(i)) for i in row]
                if 0 < x[0] < nonce_max:
                    y = (a[0] + x[0]) * pow(b[0], -1, p) % p
                    _y = (a[1] + x[1]) * pow(b[1], -1, p) % p
                    if y == _y:
                        print("[+] parameter, y: ", i, y)
        print('[+] Debug done')
        return 0, 0

# From CTF-Wiki (https://ctf-wiki.org/en/crypto/asymmetric/lattice/cvp/)
# Unlike lll_solve_1, the definition of l is "the number of bits that already known"
# i.e., p.bit_length() - x_i.bit_length()

def lll_solve_3(a, b, p, l):
    print("[+] ----------LLL SOLVE 3----------")

    # http://www.isg.rhul.ac.uk/~sdg/igor-slides.pdf
    def babai(A, w):
        A = A.LLL(delta=0.75)
        G = A.gram_schmidt()[0]
        t = w
        for i in reversed(range(A.nrows())):
            c = ((t * G[i]) / (G[i] * G[i])).round()
            t -= A[i] * c
        return w - t

    N = len(a)
    assert(len(a) == len(b))
    k = len(bin(p))-2-l
    B = Matrix(QQ, N+1, N+1)
    B.set_block(0, 0, matrix.identity(N) * p)
    B.set_block(N, 0, matrix(QQ, 1, N, b))
    B[N, N] = 1 / (2 ** (k + 1))
    closest = babai(B, vector(QQ, a + [0]))
    print("[+] y  : ", (closest[-1] * (2 ** (k + 1))) % p)
    return (closest[-1] * (2 ** (k + 1))) % p

# From Furutsuki (https://furutsuki.hatenablog.com/entry/2021/03/21/101039)
# LINE CTF 2021 babycrypto4
def test_1():
    # From https://neuromancer.sk/std/secg/secp160r1
    p = 0xffffffffffffffffffffffffffffffff7fffffff
    K = GF(p)
    a = K(0xffffffffffffffffffffffffffffffff7ffffffc)
    b = K(0x1c97befc54bd7a8b65acf89f81d4d4adc565fa45)
    E = EllipticCurve(K, (a, b))
    G = E(0x4a96b5688ef573284664698968c38bb913cbfc82, 0x23a628553168947d59dcc912042351377ac5fb32)
    E.set_order(0x0100000000000000000001f4c8f927aed3ca752257 * 0x01)

    n = int(E.order())

    dataset = []
    with open("line_ctf_babycrypto4_output.txt") as f:
        lines = f.read().strip().split("\n")
        for l in lines:
            dataset.append( [int(x, 16) for x in l.split(" ")] )
    # Given R_i, S_i, H_i and n such that
    # -inv(S_i) * H_i + k_i = inv(S_i) * R_i * d (mod n)
    # R, S, H : 160 bit
    # n       : 161 bit
    # k_i     : 32  bit
    # d       : 161 bit
    r, s, k, h = zip(*dataset)
    N = len(r)
    a = [-pow(s[i], -1, n) * h[i] % n for i in range(N)]
    b = [pow(s[i], -1, n) * r[i] % n for i in range(N)]

    # find this
    d = 68570986338040913473862653949389814867104146286
    print("[+] ans: ", d)
    return a, b, n, 32

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
    flag = 2359071867217229033899952554549558255586216332644160293624620302437641015683071171760925782397
    n = 5
    k = 350
    p = getStrongPrime(512)
    xs = [random.randint(2, p-1) for _ in range(n)]
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
        delta = p >> (k + 1)
        ui = random.randint(x - delta, x + delta)
        return ui

    # p := gmpy2.next_prime(2**160)
    p = 1461501637330902918203684832716283019655932542983

    # Given T_i, U_s, p such that
    # U_i + k_i = T_i * alpha (mod p)
    # U, T, p : 161 bit
    # k_i     : 150 bit
    # alpha   : 161 bit
    alpha = random.randint(1, p - 1)
    t = []
    u = []
    k = 10
    for i in range(22):
        t.append(random.randint(1, p - 1))
        u.append(msb(k, alpha * t[i] % p, p))

    # find this
    print("[+] ans: ", alpha)
    return u, t, p, 151


def main():
    a, b, p, l = test_1()
    # a, b, p, l = test_2()
    # a, b, p, l = test_3()

    lll_solve(a, b, p, l)
    # lll_solve_2(a, b, p, l)
    # lll_solve_3(a, b, p, l)

if __name__ == "__main__":
    main()

