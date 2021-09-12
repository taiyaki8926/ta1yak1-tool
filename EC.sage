
def help():
    print("---Tool list---")
    print("smart(p, A, B, P, Q) -> k : Calculation of k such that Q = k*P. It is efficient when EC.order == p")
    print("pohlig(p, A, B, P, Q, bound=0) -> k : It is efficient when we can get the result of factorization of EC.order()")
    print("calc_A_B(p, P, Q) -> [A,B] : Calculation of coefficient A and B such that y^2 = x^3+Ax + B mod p, from two point P and Q")
    print("singular(p, A, P, Q)-> k : It is efficient when the elliptic curve is singular")


# from: https://hxp.io/blog/25/SharifCTF-2016-crypto350-British-Elevator-writeup/
def smart(p, A, B, P:list, Q:list):
    x_P, y_P = P
    x_Q, y_Q = Q
    E = EllipticCurve(GF(p), [0, 0, 0, A, B])
    assert E.order() == p
    Q_p = Qp(p, 2)
    E_p = EllipticCurve(Q_p, [0, 0, 0, A, B])
    
    yP_p = sqrt(Q_p(x_P) ** 3 + Q_p(A) * Q_p(x_P) + Q_p(B))
    P_p = E_p(Q_p(x_P), (-yP_p, yP_p)[yP_p[0] == y_P])
    
    yQ_p = sqrt(Q_p(x_Q) ** 3 + Q_p(A) * Q_p(x_Q) + Q_p(B))
    Q_p = E_p(Q_p(x_Q), (-yQ_p, yQ_p)[yQ_p[0] == y_Q])
    
    l_Q = E_p.formal_group().log()(- (p * Q_p)[0] // (p * Q_p)[1]) / p
    l_P = E_p.formal_group().log()(- (p * P_p)[0] // (p * P_p)[1]) / p
    e = l_Q / l_P
    assert e[0] * E(x_P, y_P) == E(x_Q, y_Q)
    return e[0]

# from: https://furutsuki.hatenablog.com/entry/2020/05/05/112207#ptr-yudai%E3%81%AE%E7%99%BA%E8%A1%A8
def pohlig(p, A, B, P:list, Q:list, bound=0):
    EC = EllipticCurve(GF(p), [A, B])
    _P = EC((P[0], P[1]))
    _Q = EC((Q[0], Q[1]))
    order = EC.order()
    factors = list(factor(order))
    mods = []
    zs = []
    cur = 1
    if bound == 0:
        bound = order-1

    for f in factors:
        p, e = f
        pe = p**e
        Pi = (order // pe) * _P
        Qi = (order // pe) * _Q
        zi = Pi.discrete_log(Qi)
        zs.append(zi)
        mods.append(pe)
        cur *= pe
        if cur > bound:
            break
    
    ret = CRT(zs, mods)
    assert(ret * _P == _Q)
    return ret


def calc_A_B(p, P, Q):
    x1, y1 = P
    x2, y2 = Q
    a = inverse_mod(x1-x2, p) * (y1^2 - y2^2 - x1^3 + x2^3) % p
    b = (y1^2 - x1^3 - a*x1) % p
    try:
        EC = EllipticCurve(GF(p), [a,b])
        assert(EC.is_on_curve(x1, y1) and EC.is_on_curve(x2, y2))
    except:
        print('This may be a singular elliptic curve. Use "singular(p, A, P, Q)" function, setting A = [0, 0, 0, a, b]')
    return [a, b]

# from: https://gitlab.com/n0tsobad/ctf-writeups/tree/master/2019-02-01-nullcon/Singular
def singular(p, A:list, P, Q):
    def calc_trans(p, A):
        a1, a2, a3, a4, a6 = A
        b2 = a1^2 + 4*a2
        b4 = 2*a4 + a1*a3
        b6 = a3^2 + 4*a6
        b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2
        Di = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
        assert(Di % p == 0)
        
        x = 413400541209677581972773119133520959089878607131
        r = mod((4*a2 + a1^2)^2 - 24 * (2*a4 + a1*a3), p).sqrt()
        if r == 0:
            x = -(4*a2 + a1**2) * pow(12, -1, p) % p
            y = (-a1 * x - a3) * pow(2, -1, p) % p
            return (x, y), True
        else:
            x1 = (-(4*a2 + a1**2) + r) * pow(12, -1, p) % p
            y1 = (-a1 * x1 - a3) * pow(2, -1, p) % p
            x2 = (-(4*a2 + a1**2) - r) * pow(12, -1, p) % p
            y2 = (-a1 * x2 - a3) * pow(2, -1, p) % p
            return ((x1, y1), (x2, y2)), False
        
    def trans(p, A, ret, _P, _Q):
        a1, a2, a3, a4, a6 = A
        assert((a1*ret[1] - 3*ret[0]^2 - 2*a2*ret[0] - a4) % p == 0)
        assert((2*ret[1] + a1*ret[0] + a3) % p == 0)
        P.<x> = GF(p)[]
        y = ret[1]
        f = x^3 + a2*x^2 + a4*x + a6 - y^2 - a1*x*y - a3*y
        f_ = f.subs(x = x + ret[0])
        if f_.subs(x = 0) != 0:
            print('[+] Not this ...')
            return 0
        else:
            r = f_.roots()
            if len(r) == 1:
                assert((0,3) in r)
                print('[+] This is cusp type')
                P_ = (_P[0] - ret[0], _P[1] - ret[1])
                Q_ = (_Q[0] - ret[0], _Q[1] - ret[1])
                u = P_[0] * pow(P_[1], -1, p) % p
                v = Q_[0] * pow(Q_[1], -1, p) % p
                return v * pow(u,-1,p) % p
            else:
                print('[+] This is node type')
                assert((0,2) in r)
                r.remove((0,2))
                P_ = (_P[0] - ret[0], _P[1] - ret[1])
                Q_ = (_Q[0] - ret[0], _Q[1] - ret[1])
                if pow(p-r[0][0], (p-1)//2, p) != 1:
                    print('[+] cannot find square root ...')
                    return 0
                t = GF(p)(p-r[0][0]).square_root()
                u = (P_[1] + t*P_[0]) * pow(P_[1] - t*P_[0], -1, p) % p
                v = (Q_[1] + t*Q_[0]) * pow(Q_[1] - t*Q_[0], -1, p) % p
                try:
                    ret = discrete_log(v, u)
                    return ret
                except:
                    print('[+] cannot find discrete log ...')
                    return 0
    
    def mul_singular_EC(p, A, k, P):
        def add_EC(p, A, P, Q):
            xP, yP = P
            xQ, yQ = Q
            if P != Q:
                s = (yP - yQ) * pow(xP - xQ, -1, p) % p
                xR = (s^2 - xP - xQ) % p
                return [xR, -(yP + s * (xR - xP)) % p]
            else:
                s = (3 * xP^2 + A) * pow(2 * yP, -1, p) % p
                xR = (s^2 - 2 * xP) % p
                return [xR, -(yP + s * (xR - xP)) % p]
        _P = P
        l = int(k).bit_length()
        ls = []
        for i in range(l):
            ls.append(_P)
            _P = add_EC(p, A, _P, _P)
        init = True
        for i, _ in enumerate(bin(k)[2:][::-1]):
            if _ == '1':
                if init:
                    ans = ls[i]
                    init = False
                else:
                    ans = add_EC(p, A, ans, ls[i])
        return ans
                
    ret, flag = calc_trans(p, A)
    if flag:
        print('[+] check ', ret)
        ret = trans(p, A, ret, P, Q)
        return ret
    else:
        r1, r2 = ret
        print('[+] check ', r1)
        ret = trans(p, A, r1, P, Q)
        if ret != 0 :
            if mul_singular_EC(p, A[3], ret, P) == Q:
                print('[+] find!')
                return ret
        print('[+] check ', r2)
        ret = trans(p, A, r2, P, Q)
        if ret != 0:
            if mul_singular_EC(p, A[3], ret, P) == Q:
                print('[+] find!')
                return ret
    return 0