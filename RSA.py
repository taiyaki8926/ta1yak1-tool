import gmpy2
import owiener # if there is no module, then "python3 -m pip install owiener"
from sympy import *
import subprocess


def help():
    print("---Tool list---")
    print("calc_prime_from_d(d, e, n) -> (p, q) : Calculation of primes p and q from n, e and d")
    print("lowpub(c, e, n) -> m : It is efficient when e is very small")
    print("wiener(c, e, n) -> m : It is efficient when e is very large")
    print("view_lsb_decryption_oracle(): It is efficient when there exists an oracle that returns the LSB of the arbitrary plaintext")
    print("    Note : you must copy and paste this function and to your own solver, and implement 'get_oracle' function. ")
    print("           Cannot import because of the lack of my knowledge :(")


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


def view_get_oracle():
    print('''# this is the sample of get_oracle function
def get_oracle(_c:int, conn):
    conn.sendline(str(_c).encode())
    ret = eval(conn.recvline())
    assert(ret == 1 or ret == 0)
    return ret
    ''')


# from: https://kataware.hatenablog.jp/entry/2018/07/04/225426
def view_lsb_decryption_oracle():
    print('''
from fractions import Fraction
from math import ceil
import sys

def lsb_decryption_oracle(c, e, n, conn)
decision = [0,Fraction(n)]
i = 0
while decision[1] - decision[0] >= 1:
    i += 1
    print('[+] {}/{}'.format(i, n.bit_length()))
    cc = (c * pow(2,e,n)) % n
    try:
        oracle = get_oracle(cc, conn)
    except:
        print('Need to implement get_oracle function')
        print('You can view sample to use "view_get_oracle()"')
        sys.exit()
    if oracle == 1:
        decision[0] = Fraction(sum(decision)/2)
    else:
        decision[1] = Fraction(sum(decision)/2)
    c = cc
return ceil(decision[0])
    ''')
