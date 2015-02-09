def coppersmith_univariate(pol, modulus, beta, mm, tt, XX):
    """
    Howgrave-Graham revisited method of Coppersmith following theorem:
    
    modulus integer of unknown factorization,
    0 < beta <= 1
    0 < epsilon <= beta / 7
    b | modulus and b >= modulus^beta,
    pol is a monic polynomial of degree dd,
    
    THEN
    
    we can find roots of pol(x) = 0 mod b with
    |root| <= (1/2) * modulus^((beta^2/dd) - epsilon)
    """
    #
    # init
    #
    dd = pol.degree()
    nn = dd * mm + tt

    #
    # checks
    #
    if not 0 < beta <= 1:
        raise ValueError("beta should belongs in (0, 1]")

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    #
    # calculate bounds and display them
    #
    """
    we want to find g(x) such that ||g(xX)|| <= b^m / sqroot(n)
    we know LLL will give us a short vector v such that:
    ||v|| <= 2^((n - 1)/4) * det(L)^(1/n)

    so we want to satisfy:
    2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqroot(n)
    so we can obtain ||v|| < N^(beta*m) / sqroot(n) <= b^m / sqroot(n)
    (it's important to play with N because we might not know b)

    
    """
    
    # t optimized?
    print "\n# Optimized t?\n"
    print "we want X^(n-1) < N^(beta*m) so that each vector is helpful"
    cond1 = RR(XX^(nn-1))
    print "* X^(n-1) = ", cond1
    cond2 = pow(modulus, beta*mm)
    print "* N^(beta*m) = ", cond2
    print "* X^(n-1) < N^(beta*m) # GOOD" if cond1 < cond2 else "* X^(n-1) >= N^(beta*m) # NOT GOOD"
    
    # bound for X
    print "\n# X bound respected?\n"
    print "we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M"
    print "* X =", XX
    cond2 = RR(modulus^(((2*beta*mm)/(nn-1)) - ((dd*mm*(mm+1))/(nn*(nn-1)))) / 2)
    print "* M =", cond2
    print "* X <= M # GOOD" if XX <= cond2 else "* X > M # NOT GOOD"

    # solution possible?
    print "\n# Solutions possible?\n"
    detL = RR(modulus^(dd * mm * (mm + 1) / 2) * XX^(nn * (nn - 1) / 2))
    print "we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqroot(n)"
    cond1 = RR(2^((nn - 1)/4) * detL^(1/nn))
    print "* 2^((n - 1)/4) * det(L)^(1/n) = ", cond1
    cond2 = RR(modulus^(beta*mm) / sqrt(nn))
    print "* N^(beta*m) / sqroot(n) = ", cond2
    print "* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqroot(n) # SOLUTION WILL BE FOUND" if cond1 < cond2 else "* 2^((n - 1)/4) * det(L)^(1/n) >= N^(beta*m) / sqroot(n) #NO SOLUTIONS MIGHT BE FOUND"
    
    #
    # Coppersmith revisited algo
    #

    # change ring of pol and x
    polZ = pol.change_ring(ZZ)
    x = polZ.parent().gen()

    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append(x**jj * modulus**(mm - ii) * polZ**ii)
    for ii in range(tt):
        gg.append(x**ii * polZ**mm)
    
    # construct lattice B
    BB = Matrix(ZZ, nn)
    """here sage's implementation uses rectangular matrix
    why???
    """
    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj] * XX**jj

    # LLL
    BB = BB.LLL()

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[0, ii] / XX**ii

    # factor polynomial
    potential_roots = new_pol.roots()

    # test roots on original pol
    roots = []
    for root in potential_roots:
        result = ZZ(pol(ZZ(root[0])))

        if gcd(modulus, result) >= modulus**beta:
            roots.append(root[0])

    # no roots found
    return roots

############################################
# Test on Stereotyped Messages
##########################################    
# (from http://www.sagemath.org/doc/reference/polynomial_rings/sage/rings/polynomial/polynomial_modn_dense_ntl.html#sage.rings.polynomial.polynomial_modn_dense_ntl.small_roots)

print "//////////////////////////////////"
print "// TEST 1"
print "////////////////////////////////"


Nbits, Kbits = 512, 56; e = 3; p = 2^256 + 2^8 + 2^5 + 2^3 + 1; q = 2^256 + 2^8 + 2^5 + 2^3 + 2^2 + 1; N = p*q; ZmodN = Zmod( N ); K = ZZ.random_element(0, 2^Kbits); Kdigits = K.digits(2); M = [0]*Kbits + [1]*(Nbits-Kbits); 
for i in range(len(Kdigits)): M[i] = Kdigits[i]; 
M = ZZ(M, 2); C = ZmodN(M)^e; P.<x> = PolynomialRing(ZmodN, implementation='NTL'); f = (2^Nbits - 2^Kbits + x)^e - C


# PLAY WITH THOSE
""" epsilon can be anything?
    if epsilon is <= 1/7 * beta
    then we can use m = ceil(beta^2/delta epsilon)
    otherwise m >= max{ beta^2/delta epsilon, 7beta/delta }
"""

dd = f.degree()
beta = 1
epsilon = beta / 7
mm = ceil(beta**2 / (dd * epsilon))
tt = floor(dd * mm * ((1/beta) - 1))
XX = ceil(N**((beta**2/dd) - epsilon))
roots = coppersmith_univariate(f, N, beta, mm, tt, XX)

# output
print "\n# Solutions"
print "we want to find:",str(K)
print "we found:", str(roots)
print "\n"

############################################
# Test on Factoring with High Bits Known
##########################################
print "//////////////////////////////////"
print "// TEST 2"
print "////////////////////////////////"

length = 512; hidden = 110; p = next_prime(2^int(round(length/2))); q = next_prime( round(pi.n()*p) ); N = p*q; qbar = q + ZZ.random_element(0,2^hidden-1); F.<x> = PolynomialRing(Zmod(N), implementation='NTL'); f = x - qbar;

# PLAY WITH THOSE:
beta = 0.5
dd = f.degree()
epsilon = beta / 7
mm = ceil(beta**2 / (dd * epsilon))
tt = floor(dd * mm * ((1/beta) - 1))
XX = ceil(N**((beta**2/dd) - epsilon))
XX += 100000000000000000000000000000000
roots = coppersmith_univariate(f, N, beta, mm, tt, XX)

# output
#d = f.small_roots(X=2^hidden-1, beta=0.5)[0]; print("we found:", qbar - d)
print "\n# Solutions"
print "we want to find:", qbar - q
print "we found:", roots
