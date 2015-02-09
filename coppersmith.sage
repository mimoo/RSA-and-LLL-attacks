def coppersmith_howgrave_univariate(pol, modulus, beta, mm, tt, XX):
    """
    Howgrave-Graham revisited method of Coppersmith following theorem:
    
    modulus integer of unknown factorization,
    0 < beta <= 1
    0 < epsilon <= beta / 7
    b | modulus and b >= modulus^beta,
    pol is a monic polynomial of degree dd,
    
    THEN
    
    we can find roots of pol(x) = 0 mod b with
    |root| <= |X|
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
    * we want to find g(x) such that ||g(xX)|| <= b^m / sqrt(n)

    * we know LLL will give us a short vector v such that:
    ||v|| <= 2^((n - 1)/4) * det(L)^(1/n)

    * we will use that vector as a coefficient vecotr for our g(x)
    
    * so we want to satisfy:
    2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)
    
    so we can obtain ||v|| < N^(beta*m) / sqrt(n) <= b^m / sqrt(n)
    (it's important to play with N because we might not know b)
    """
    
    # t optimized?
    print "\n# Optimized t?\n"
    print "we want X^(n-1) < N^(beta*m) so that each vector is helpful"
    cond1 = RR(XX^(nn-1))
    print "* X^(n-1) = ", cond1
    cond2 = pow(modulus, beta*mm)
    print "* N^(beta*m) = ", cond2
    print "* X^(n-1) < N^(beta*m) \n-> GOOD" if cond1 < cond2 else "* X^(n-1) >= N^(beta*m) \n-> NOT GOOD"
    
    # bound for X
    print "\n# X bound respected?\n"
    print "we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M"
    print "* X =", XX
    cond2 = RR(modulus^(((2*beta*mm)/(nn-1)) - ((dd*mm*(mm+1))/(nn*(nn-1)))) / 2)
    print "* M =", cond2
    print "* X <= M \n-> GOOD" if XX <= cond2 else "* X > M \n-> NOT GOOD"

    # solution possible?
    print "\n# Solutions possible?\n"
    detL = RR(modulus^(dd * mm * (mm + 1) / 2) * XX^(nn * (nn - 1) / 2))
    print "we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)"
    cond1 = RR(2^((nn - 1)/4) * detL^(1/nn))
    print "* 2^((n - 1)/4) * det(L)^(1/n) = ", cond1
    cond2 = RR(modulus^(beta*mm) / sqrt(nn))
    print "* N^(beta*m) / sqrt(n) = ", cond2
    print "* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n) \n-> SOLUTION WILL BE FOUND" if cond1 < cond2 else "* 2^((n - 1)/4) * det(L)^(1/n) >= N^(beta*m) / sqroot(n) \n-> NO SOLUTIONS MIGHT BE FOUND (but we never know)"

    # warning about X
    print "\n# Note that no solutions will be found _for sure_ if you don't respect:\n* |root| < X \n* b >= modulus^beta\n"
    
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
            gg.append((x * XX)**jj * modulus**(mm - ii) * polZ(x * XX)**ii)
    for ii in range(tt):
        gg.append((x * XX)**ii * polZ(x * XX)**mm)
    
    # construct lattice B
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj in range(ii+1):
            BB[ii, jj] = gg[ii][jj]

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

    # 
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
roots = coppersmith_howgrave_univariate(f, N, beta, mm, tt, XX)

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

# RSA gen
length = 512;
p = next_prime(2^int(round(length/2)));
q = next_prime( round(pi.n()*p) );
N = p*q;

# qbar is q + [hidden_size_random]
hidden = 110;
diff = ZZ.random_element(0, 2^hidden-1)
qbar = q + diff; 

F.<x> = PolynomialRing(Zmod(N), implementation='NTL'); 
f = x - qbar;

# PLAY WITH THOSE:
beta = 0.5 # we should have q >= N^beta
dd = f.degree()
epsilon = beta / 7
mm = ceil(beta**2 / (dd * epsilon))
tt = floor(dd * mm * ((1/beta) - 1))
XX = ceil(N**((beta**2/dd) - epsilon)) + 1000000000000000000000000000000000 # we should have |diff| < X
roots = coppersmith_howgrave_univariate(f, N, beta, mm, tt, XX)

# output
print "\n# Solutions"
print "we want to find:", qbar - q
print "we found:", roots
