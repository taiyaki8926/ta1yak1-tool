# ta1yak1-tool

Various CTF crypto tools

Made by [taiyaki](https://twitter.com/ta1yak1_8926)

# Index

* RSA.py

    * `p` and `q` calculation from public key `e, n` and secret key `d`

    * Low public exponent attack : It is efficient when `e` is very small 
    
    * Wiener's attack : It is efficient when `e` is very large"
    
    * LSB decryption oracle Attack : It is efficient when there exists an oracle that returns the LSB of the arbitrary plaintext

* RSA_sage.sage

    * Partial key exposure attack : It is efficient when we know tha partial LSBs of `d`

* EC.sage
    * Smart's attack : It is efficient when EC.order() == p

    * Pohlig-Hellman attack : It is efficient when we can get the result of factorization of EC.order()

    * Coefficients `A` and `B` (such that `y^2 = x^3 + Ax + B mod p`) calculation from two point `P, Q` and modulo `p`

    * Attack for singular elliptic curve : It is efficient when the elliptic curve is singular. The definition of singular is discriminant is zero, or cannot be defined as `EllipticCurve(GF(p), A:coefficient list)` in Sage

# Usage

```Python
from ta1yak1tool import RSA
RSA.help()
```

# How to import sage file (for developer)

```
$ sage -preparse *.sage;mv *.sage.py *.py
```


# References

https://www.slideshare.net/sonickun/rsa-n-ssmjp (Only Japanese)