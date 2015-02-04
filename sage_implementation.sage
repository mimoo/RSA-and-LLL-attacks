def small_roots(pol, X=None, beta=1.0, epsilon=None, **kwds):
    r"""
    Let `N` be the characteristic of the base ring this polynomial
    is defined over: ``N = self.base_ring().characteristic()``.
    This method returns small roots of this polynomial modulo some
    factor `b` of `N` with the constraint that `b >= N^\beta`.
    Small in this context means that if `x` is a root of `f`
    modulo `b` then `|x| < X`. This `X` is either provided by the
    user or the maximum `X` is chosen such that this algorithm
    terminates in polynomial time. If `X` is chosen automatically
    it is `X = ceil(1/2 N^{\beta^2/\delta - \epsilon})`.
    The algorithm may also return some roots which are larger than `X`.
    'This algorithm' in this context means Coppersmith's algorithm for finding
    small roots using the LLL algorithm. The implementation of this algorithm
    follows Alexander May's PhD thesis referenced below.

    INPUT:

    - ``X`` -- an absolute bound for the root (default: see above)
    - ``beta`` -- compute a root mod `b` where `b` is a factor of `N` and `b
      \ge N^\beta`. (Default: 1.0, so `b = N`.)
    - ``epsilon`` -- the parameter `\epsilon` described above. (Default: `\beta/8`)
    - ``**kwds`` -- passed through to method :meth:`Matrix_integer_dense.LLL() <sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL>`.

    EXAMPLES:

    First consider a small example::

        sage: N = 10001
        sage: K = Zmod(10001)
        sage: P.<x> = PolynomialRing(K, implementation='NTL')
        sage: f = x^3 + 10*x^2 + 5000*x - 222

    This polynomial has no roots without modular reduction (i.e. over `\ZZ`)::

        sage: f.change_ring(ZZ).roots()
        []

    To compute its roots we need to factor the modulus `N` and use the Chinese
    remainder theorem::

        sage: p,q = N.prime_divisors()
        sage: f.change_ring(GF(p)).roots()
        [(4, 1)]
        sage: f.change_ring(GF(q)).roots()
        [(4, 1)]

        sage: crt(4, 4, p, q)
        4

    This root is quite small compared to `N`, so we can attempt to
    recover it without factoring `N` using Coppersmith's small root
    method::

        sage: f.small_roots()
        [4]

    An application of this method is to consider RSA. We are using 512-bit RSA
    with public exponent `e=3` to encrypt a 56-bit DES key. Because it would be
    easy to attack this setting if no padding was used we pad the key `K` with
    1s to get a large number::

        sage: Nbits, Kbits = 512, 56
        sage: e = 3

    We choose two primes of size 256-bit each::

        sage: p = 2^256 + 2^8 + 2^5 + 2^3 + 1
        sage: q = 2^256 + 2^8 + 2^5 + 2^3 + 2^2 + 1
        sage: N = p*q
        sage: ZmodN = Zmod( N )

    We choose a random key::

        sage: K = ZZ.random_element(0, 2^Kbits)

    and pad it with 512-56=456 1s::

        sage: Kdigits = K.digits(2)
        sage: M = [0]*Kbits + [1]*(Nbits-Kbits)
        sage: for i in range(len(Kdigits)): M[i] = Kdigits[i]

        sage: M = ZZ(M, 2)

    Now we encrypt the resulting message::

        sage: C = ZmodN(M)^e

    To recover `K` we consider the following polynomial modulo `N`::

        sage: P.<x> = PolynomialRing(ZmodN, implementation='NTL')
        sage: f = (2^Nbits - 2^Kbits + x)^e - C

    and recover its small roots::

        sage: Kbar = f.small_roots()[0]
        sage: K == Kbar
        True

    The same algorithm can be used to factor `N = pq` if partial
    knowledge about `q` is available. This example is from the
    Magma handbook:

    First, we set up `p`, `q` and `N`::

        sage: length = 512
        sage: hidden = 110
        sage: p = next_prime(2^int(round(length/2)))
        sage: q = next_prime( round(pi.n()*p) )
        sage: N = p*q

    Now we disturb the low 110 bits of `q`::

        sage: qbar = q + ZZ.random_element(0,2^hidden-1)

    And try to recover `q` from it::

        sage: F.<x> = PolynomialRing(Zmod(N), implementation='NTL')
        sage: f = x - qbar

    We know that the error is `\le 2^{\text{hidden}}-1` and that the modulus
    we are looking for is `\ge \sqrt{N}`::

        sage: set_verbose(2)
        sage: d = f.small_roots(X=2^hidden-1, beta=0.5)[0] # time random
        verbose 2 (<module>) m = 4
        verbose 2 (<module>) t = 4
        verbose 2 (<module>) X = 1298074214633706907132624082305023
        verbose 1 (<module>) LLL of 8x8 matrix (algorithm fpLLL:wrapper)
        verbose 1 (<module>) LLL finished (time = 0.006998)
        sage: q == qbar - d
        True

    REFERENCES:

    Don Coppersmith. *Finding a small root of a univariate modular equation.*
    In Advances in Cryptology, EuroCrypt 1996, volume 1070 of Lecture
    Notes in Computer Science, p. 155--165. Springer, 1996.
    http://cr.yp.to/bib/2001/coppersmith.pdf

    Alexander May. *New RSA Vulnerabilities Using Lattice Reduction Methods.*
    PhD thesis, University of Paderborn, 2003.
    http://www.cs.uni-paderborn.de/uploads/tx_sibibtex/bp.pdf
    """
    from sage.misc.misc import verbose
    from sage.matrix.constructor import Matrix
    from sage.rings.all import RR

    N = pol.parent().characteristic()

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    beta = RR(beta)
    if beta <= 0.0 or beta > 1.0:
        raise ValueError("0.0 < beta <= 1.0 not satisfied.")

    f = pol.change_ring(ZZ)

    P,(x,) = f.parent().objgens()

    delta = f.degree()

    if epsilon is None:
        epsilon = beta/8
    verbose("epsilon = %d"%epsilon, level=2)

    m = max(beta**2/(delta * epsilon), 7*beta/delta).ceil()
    verbose("m = %d"%m, level=2)

    t = int( ( delta*m*(1/beta -1) ).floor() )
    verbose("t = %d"%t, level=2)

    if X is None:
        X = (0.5 * N**(beta**2/delta - epsilon)).ceil()
    verbose("X = %s"%X, level=2)

    # we could do this much faster, but this is a cheap step
    # compared to LLL
    g  = [x**j * N**(m-i) * f**i for i in range(m) for j in range(delta) ]
    g.extend([x**i * f**m for i in range(t)]) # h

    B = Matrix(ZZ, len(g), delta*m + max(delta,t) )
    for i in range(B.nrows()):
        for j in range( g[i].degree()+1 ):
            B[i,j] = g[i][j]*X**j

    B =  B.LLL(**kwds)

    f = sum([ZZ(B[0,i]//X**i)*x**i for i in range(B.ncols())])
    R = f.roots()

    ZmodN = pol.base_ring()
    roots = set([ZmodN(r) for r,m in R if abs(r) <= X])
    Nbeta = N**beta
    return [root for root in roots if N.gcd(ZZ(pol(root))) >= Nbeta]

# test
Nbits, Kbits = 512, 56
e = 3
p = 2^256 + 2^8 + 2^5 + 2^3 + 1
q = 2^256 + 2^8 + 2^5 + 2^3 + 2^2 + 1
N = p*q
ZmodN = Zmod( N )
K = ZZ.random_element(0, 2^Kbits)
Kdigits = K.digits(2)
M = [0]*Kbits + [1]*(Nbits-Kbits)
for i in range(len(Kdigits)): M[i] = Kdigits[i]

M = ZZ(M, 2)
C = ZmodN(M)^e
P.<x> = PolynomialRing(ZmodN, implementation='NTL')
f = (2^Nbits - 2^Kbits + x)^e - C
print(K, small_roots(f))
