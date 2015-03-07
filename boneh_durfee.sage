def boneh_durfee(pol, modulus, mm, tt, XX, YY):
    """
    Boneh and Durfee revisited by Herrmann and May
    finds a solution if:
    * |x| < e^delta
    * |y| < e^0.5
    whenever delta < 1 - sqrt(2)/2 ~ 0.292
    """

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
            # debug
            yshifttest = y^jj * polZ(u, x, y)^kk * modulus^(mm - kk)
            yshifttest = Q(yshifttest).lift()
            if yshifttest == 0 or yshifttest(uu,xx,yy) % e^mm != 0:
                print "OUIE"
                print jj, kk
                return gg, yshifttest
            # end debug

    # y-shifts monomials
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            monomials.append(u^kk * y^jj)

    '''debug'''
    for ii in range(36,48):
        new = gg[ii](u/UU,x/XX,y/YY)
        print ii, new(uu,xx,yy) % e^mm == 0
    return new,gg
    '''end debug'''


    # construct lattice B
    nn = len(monomials)

    BB = Matrix(ZZ, nn)

    for ii in range(nn):

        BB[ii, 0] = gg[ii](0, 0, 0)

        for jj in range(1, ii + 1):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj])

    '''debug'''
    # for ii in range(nn):
    #     '''check if vector is helpful'''
    #     if BB[ii,ii] > modulus^mm:
    #         print "vector "+str(ii)+" not helpful", BB[ii, ii]

    #     '''check triangular matrix'''  
    #     for jj in range(ii + 1, nn):
    #         if BB[ii,jj] != 0:
    #             print "ugggg", ii, jj

    '''debug det'''
    # from the lattice
    '''
    det = 1
    for ii in range(nn):
        det *= BB[ii, ii]

    # from the formula
    sx = mm^3 / 6
    sy = tho^2 * sx
    su = ((1/6) + (tho/3))*mm^3
    se = su
    detformula = XX^sx * YY^sy * UU^su * modulus^se

    print det - int(detformula)

    bound = modulus^(mm*nn)

    if det >= bound:
        print "we don't have det < bound"
        print "det - bound = ", abs(det - bound)
    else:
        print "det < bound"
    '''


    ''' debug: test if roots work on every vectors'''
    #PR.<x,y> = PolynomialRing(ZZ) # useful?
    for jj in range(nn):
        print jj
        poltest = 0
        for ii in range(nn):
            poltest += monomials[ii](u,x,y) * BB[jj, ii] / monomials[ii](UU,XX,YY)
        if poltest == 0 or poltest(uu,xx,yy) % e^mm != 0:
            print "AIE AIE AIE number 1"
            print jj
            #print poltest
            print poltest == 0
            print poltest(uu,xx,yy) != 0 % e^mm
            print poltest % e^mm
           # return poltest, BB

    # LLL
    BB = BB.LLL()

    # shortest vectors to polynomials
    pol1 = pol2 = 0

    PR.<x,y> = PolynomialRing(ZZ) # useful?

    for ii in range(nn):
        pol1 += monomials[ii](x*y+1,x,y) * BB[0, ii] / monomials[ii](UU,XX,YY)
        pol2 += monomials[ii](x*y+1,x,y) * BB[1, ii] / monomials[ii](UU,XX,YY)

    ''' debug: test if roots work on every vectors'''
    for jj in range(2,nn):
        print jj
        poltest = 0
        for ii in range(nn):
            poltest += monomials[ii](x*y+1,x,y) * BB[jj, ii] / monomials[ii](UU,XX,YY)
        if poltest == 0 or poltest(xx,yy) % e^mm != 0:
            print "AIE AIE AIE"
            print jj
            print poltest
            print poltest == 0
            print poltest(xx,yy) != 0 % e^mm
            print poltest % e^mm
            return poltest, BB
    print "done, roots work for every vectors"

    '''
    # revert substitution
    u, x, y = pol1.parent().gens() #dunno why I have to do this
    pol1 = pol1.subs({u:x*y + 1})
    pol2 = pol2.subs({u:x*y + 1})
    '''
    print pol1(xx,yy)
    print pol2(xx,yy)

    # resultant
    polx = pol1.resultant(pol2)

    # DOESN'T WORK HERE
    print polx

    return pol1, pol2


############################################
# Test 
##########################################

# RSA gen
length = 1024
p = next_prime(2^int(round(length/2)));
q = next_prime( round(pi.n()*p) );
N = p*q;
phi = (p-1)*(q-1)

d = int(2^(log(N^(0.25)))) # short d
if d % 2 == 0: d += 1
print "d starts at", d
while gcd(d, phi) != 1:
    d += 2
e = d.inverse_mod((p-1)*(q-1))

print "d:", d

# Problem put in equation
P.<x,y> = PolynomialRing(Zmod(e))
A = int((N+1)/2)
pol = 1 + x * (A + y)
delta = (2 - sqrt(2)) / 2
tho = (1 - 2 * delta)
m = 7
t = int(tho * m)
X = floor(e^0.292)
Y = 2*floor(e^0.5)

# hard debug
m = 7
t = 3

#
# debug
# 
yy = (-p -q)/2
xx = (e * d - 1) / (A + yy)
uu = xx*yy + 1

print "|y| < Y:", abs(yy) < Y
print "|x| < X:", abs(xx) < X
print "d < N^0.292", d < N^(0.292)

# boneh_durfee
solx, soly = boneh_durfee(pol, e, m, t, X, Y)