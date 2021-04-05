from pwn import *
from sympy import var, poly, core

def help():
    print("---Tool list---")
    print("interpolation(i, o, p) -> [coeff] : Lagrange interpolation")


def interpolation(_input, _output, p):
    var('x')

    def lagrange_coeff(offset, _input, p):
        _coeff = 1
        for i in range(len(_input)):
            if i != offset:
                _coeff *= (x - _input[i]) * pow((_input[offset] - _input[i]), -1, p)
        return _coeff

    _f = 0
    for i in range(len(_input)):
        _f += _output[i] * lagrange_coeff(i, _input, p)
    
    if type(_f) == core.add.Add:
        _f = poly(_f)
        _flag = _f.all_coeffs()
        return [_f % p for _f in _flag]
    else:
        return [_f % p]
    