import gmpy2
import owiener # if there is no module, then "python3 -m pip install owiener"
from sympy import *


def help():
    print("---Tool list---")
    print("calc_prime_from_d(d, e, n) -> (p, q) : Calculation of primes p and q from n, e and d")
    print("lowpub(c, e, n) -> m : It is efficient when e is very small")
    print("wiener(c, e, n) -> m : It is efficient when e is very large")


def calc_prime_from_d(d, e, n):
    # check that p is larger than 3 or not
    if n % 2 == 0:
        return 2, n // 2
    elif n % 3 == 0:
        return 3, n // 3
    else:
        k_min = (e * d - 1) // n
        for k in range(k_min, 2 * k_min + 1):
            if (e * d - 1) % k == 0:
                phi = (e * d - 1) // k
                _sum = n - phi  + 1
                x = Symbol('x')
                p, q = solve(x ** 2 - _sum * x + n)
                if p.is_integer :
                    return p, q

def lowpub(c, e, n):
    m = gmpy2.iroot(c, e)[0]
    return m

def wiener(c, e, n):
    d = owiener.attack(e, n)
    return pow(c, d, n)