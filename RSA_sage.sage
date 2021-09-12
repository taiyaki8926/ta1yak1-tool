import sys

def help():
    print("---Tool list---")
    print("coppersmith_given_d(c, e, n, d0) -> m : It is efficient when we know tha partial LSBs of d")


def coppersmith_given_d(c, e, n, d0):
    def partial_p(p0, kbits, n):
        PR.<x> = PolynomialRing(Zmod(n))
        nbits = n.nbits()

        f = 2^kbits*x + p0
        f = f.monic()
        roots = f.small_roots(X=2^(nbits//2-kbits), beta=0.3)  # find root < 2^(nbits//2-kbits) with factor >= n^0.3
        if roots:
            x0 = roots[0]
            p = gcd(2^kbits*x0 + p0, n)
            return ZZ(p)

    def find_p(d0, kbits, e, n):
        X = var('X')
        for k in range(1, e+1):
            results = solve_mod([e*d0*X - k*X*(n-X+1) + k*n == X], 2^kbits)
            for x in results:
                p0 = ZZ(x[0])
                try:
                    p = partial_p(p0, kbits, n)
                    if p:
                        return p
                except:
                    pass

    p = find_p(d0, d0.nbits(), e, n)
    q = n//p
    assert(n % p == 0)
    d = pow(e, -1, (p-1)*(q-1))
    m = pow(c, d, n)
    return m



if __name__ == '__main__':
    # sample
    n = 91228313843871422767839961657666790992268686181318288926107142629455536678341935900541382795206832349055585535906809836634229957830715848801103620493532865352996568386538805904907414117294419594023947589369412777921689835783244330332449834233642171173739935465571238598520641743485122408455438394294985480639
    e = 17
    d0 = 11431345612482449038094888498936316909461350091333222118655485374375384706307876941389738498046341068089191534780236606397204076298227974869985785520393937
    c = 3514670141058599920940442500404554889150141069826897009811842212403284442746185817636482379495201627390855088639977010471092729461125359933776264573183725753359659913176329404374125862300034280900990690348106060362523466667237379073824116366226768297382408622204390184563632005893478202279135233583985544305
    