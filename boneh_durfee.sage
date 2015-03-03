def boneh_durfee(pol, modulus, mm, tt, XX, YY):
    """
    Boneh and Durfee revisited by Herrmann and May
    finds a solution if:
    * |x| < e^delta
    * |y| < e^0.5
    whenever delta < 1 - sqrt(2)/2 ~ 0.292
    """


    #
    # calculate bounds and display them
    #
    '''to do'''
    #
    # Algorithm
    #

    # change ring of pol and x
    #polZ = pol.change_ring(ZZ)
    #x, y = polZ.parent().gens()
    '''useless?'''

    # substitution (Herrman and May)
    PR.<x, y, u> = PolynomialRing(ZZ)
    Q = PR.quotient(x*y + 1 - u) # u = x*y + 1
    polZ = Q(pol).lift()

    UU = XX*YY + 1

    # x-shifts
    gg = []
    for kk in range(mm + 1):
        for ii in range(mm - kk + 1):
            gg.append((x * XX)^ii * modulus^(mm - kk) * polZ(x * XX, y * YY, u * UU)^kk)

    # y-shifts (selected by Herrman and May)
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            yshift = (y * YY)^jj * polZ(x * XX, y * YY, u * UU)^kk * modulus^(mm - kk)
            gg.append(Q(yshift).lift()) # substitution

    # unravelled linerization (Herrman and May)
    monomials = []

    # x-shift
    for kk in range(mm + 1):
        for ii in range(kk + 1):
            monomials.append(u^(ii) * x^(kk-ii))

    # y-shift
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            monomials.append(u^kk * y^jj)

    # construct lattice B
    nn = len(monomials)
    return monomials,0
    BB = Matrix(ZZ, nn)

    for ii in range(nn):

        BB[ii, 0] = gg[ii](0, 0)

        for jj in range(1, ii + 1):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj])

    #
    # DET
    #
    det = 1
    for ii in range(nn):
        det *= BB[ii, ii]

    bound = modulus^(mm * (nn - 1)) / (nn * 2^nn)^((nn - 1)/2)

    print "det=", det
    print "bound=", bound
    print "det < bound?", det < bound

    return 0,0
    # LLL
    BB = BB.LLL()

    # transform shortest vectors in polynomials  
    pol1 = pol2 = 0

    for ii in range(nn):
        pol1 += monomials[ii] * BB[0, ii] / monomials[ii](X,Y)
        pol2 += monomials[ii] * BB[1, ii] / monomials[ii](X,Y)

    # resultant
    polx = pol1.resultant(pol2, y)
    poly = pol1.resultant(pol2, x)
    '''doesn't work!'''

    solx = polx.roots()
    soly = poly.roots()

    return solx, soly


############################################
# Test 
##########################################

# RSA gen
length = 512;
p = next_prime(2^int(round(length/2)));
q = next_prime( round(pi.n()*p) );
N = p*q;
phi = (p-1)*(q-1)

d = 3
while gcd(d, phi) != 1:
    d += 2
e = d.inverse_mod((p-1)*(q-1))

print "d:", d

# Problem put in equation
P.<x,y> = PolynomialRing(Zmod(e))
pol = 1 + x * (N + 1 + y)
m = 20
t = 1
"""
how to choose m and t?
"""
X = floor(e^0.292)
Y = floor(e^0.5)
"""why those values?
"""

# boneh_durfee
solx, soly = boneh_durfee(pol, e, m, t, X, Y)