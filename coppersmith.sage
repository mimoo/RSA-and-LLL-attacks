def coppersmith_univariate(pol, bb, beta):
    """Howgrave-Graham revisited method
    using with epsilon
    """
    # init
    dd = pol.degree()
    NN = pol.parent().characteristic()

    # checks
    if not 0 < beta <= 1:
        raise ValueError("beta should belongs in (0, 1]")

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")
    """ in reality it doesn't have to be, we could play around this
    """
    
    # choose epsilon, m and t
    """ epsilon can be anything?
    if epsilon is <= 1/7 * beta
    then we can use m = ceil(beta^2/delta epsilon)
    otherwise m >= max{ beta^2/delta epsilon, 7beta/delta }
    """
    epsilon = beta / 7
    mm = ceil(beta**2 / (dd * epsilon))
    tt = floor(dd * mm * ((1/beta) - 1)) # t = 0 if beta = 1, rly?

    # change ring of pol and x
    polZ = pol.change_ring(ZZ) # shouldnt it be bb^mm ?
    x = polZ.parent().gen()
    
    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append(x**jj * bb**(mm - ii) * polZ**ii)
    hh = [] # beta=1 => t=0 => no h_i polynomials
    for ii in range(tt):
        hh.append(x**ii * polZ**mm)
    
    # compute bound X
    XX = ceil(N**((beta**2/dd) - epsilon))
    
    # construct lattice B
    nn = dd * mm + tt
    BB = Matrix(ZZ, nn) # why not use gen_lattice?
    """here sage's implementation uses rectangular matrix
    why???
    """
    for ii in range(nn):
        for jj in range(ii+1):
            # fill gg
            if ii < dd*mm:
                BB[ii, jj] = gg[ii][jj] * XX**jj
            # fill hh
            else:
                BB[ii, jj] = 0#hh[ii][jj]

    # LLL
    BB = BB.LLL()
    
    # Find shortest vector in new basis
    """ Apparently Sage doesn't sort after LLL
    """
    normn = norm(BB[0])
    norm_index = 0

    for ii in range(1, nn):
        if norm(BB[ii]) < normn:
            normn = norm(BB[ii])
            norm_index = ii

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[norm_index, ii] / XX**ii
    
    # factor polynomial
    potential_roots = new_pol.roots() # doesn't find anything...
    print(potential_roots)
    # test roots on original pol
    roots = []
    for root in potential_roots:
        result = ZZ(pol(ZZ(root[0])))

        if gcd(NN, result) >= NN**beta:
            roots.append(root[0])

    # no roots found
    return roots
    
# Test 1
# (from http://www.jscoron.fr/cours/mics3crypto/tpcop.pdf)
"""
N = 2122840968903324034467344329510307845524745715398875789936591447337206598081
C = 1792963459690600192400355988468130271248171381827462749870651408943993480816

K = Zmod(N)
R.<x> = PolynomialRing(K)
pol = (2**500 + x)**3 - C
M = coppersmith_univariate(pol, N, 1)
"""
# pol.small_roots() doesn't compute anything either

# Test 2
# (from http://www.sagemath.org/doc/reference/polynomial_rings/sage/rings/polynomial/polynomial_modn_dense_ntl.html#sage.rings.polynomial.polynomial_modn_dense_ntl.small_roots)

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

print("solution a trouver:", K)
roots = coppersmith_univariate(f, N, 1)
print(roots)
