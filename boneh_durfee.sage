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

    # substitution (Herrman and May)
    PR.<u, x, y> = PolynomialRing(ZZ)
    Q = PR.quotient(x*y + 1 - u) # u = x*y + 1
    polZ = Q(pol).lift()

    UU = XX*YY + 1

    # x-shifts
    gg = []

    for kk in range(mm + 1):
        for ii in range(mm - kk + 1):
            xshift = (x * XX)^ii * modulus^(mm - kk) * polZ(u * UU, x * XX, y * YY)^kk
            gg.append(xshift)
    gg.sort()

    # x-shifts monomials
    monomials = []
    for polynomial in gg:
        for monomial in polynomial.monomials():
            if monomial not in monomials:
                monomials.append(monomial)
    monomials.sort()

    # y-shifts (selected by Herrman and May)
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            yshift = (y * YY)^jj * polZ(u * UU, x * XX, y * YY)^kk * modulus^(mm - kk)
            gg.append(Q(yshift).lift()) # substitution

    # y-shifts monomials
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            monomials.append(u^kk * y^jj)

    #
    # DEBUG
    # 
    
    # what are the roots?
    #e * d = 1 + x (N + 1 -p - q)
    y_temp = -p -q
    x_temp = (e * d - 1) / (N + 1 + y_temp)
    for ii in range(len(gg) - 1):
        print gg[ii](x_temp*y_temp,x_temp,y_temp) % e
    return 0,0

    #
    # END DEBUG
    #

    # construct lattice B
    nn = len(monomials)

    BB = Matrix(ZZ, nn)

    for ii in range(nn):

        BB[ii, 0] = gg[ii](0, 0, 0)

        for jj in range(1, ii + 1):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj])
    
    #
    # DET CHECK (OPTIONAL)
    #
    
    det = 1
    for ii in range(nn):
        det *= BB[ii, ii]

    bound = modulus^(mm * (nn - 1)) / (nn * 2^nn)^((nn - 1)/2)
    bound = int(bound)

    if det >= bound:
        print "we don't have det < bound"
        #print "det - bound = ", abs(det - bound)
    
    # LLL
    BB = BB.LLL()

    # shortest vectors to polynomials  
    pol1 = pol2 = 0

    for ii in range(nn):
        pol1 += monomials[ii] * BB[0, ii] / monomials[ii](UU,XX,YY)
        pol2 += monomials[ii] * BB[1, ii] / monomials[ii](UU,XX,YY)

    # revert substitution
    u, x, y = pol1.parent().gens() #dunno why I have to do this
    pol1 = pol1.subs({u:x*y + 1})
    pol2 = pol2.subs({u:x*y + 1})

    # resultant
    polx = pol1.resultant(pol2)

    # DOESN'T WORK HERE
    print polx

    return polx, gg


############################################
# Test 
##########################################

# RSA gen
length = 512;
p = next_prime(2^int(round(length/2)));
q = next_prime( round(pi.n()*p) );
N = p*q;
phi = (p-1)*(q-1)

d = 3 # short d
while gcd(d, phi) != 1:
    d += 2
e = d.inverse_mod((p-1)*(q-1))

print "d:", d

# Problem put in equation
P.<x,y> = PolynomialRing(Zmod(e))
pol = 1 + x * (N + 1 + y)
delta = (2 - sqrt(2)) / 2
tho = (1 - 2 * delta)
m = 5
t = int(tho * m)
X = floor(e^0.292)
Y = floor(e^0.5)

# boneh_durfee
solx, soly = boneh_durfee(pol, e, m, t, X, Y)