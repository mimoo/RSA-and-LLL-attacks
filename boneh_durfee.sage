debug = True

# test list of polynomials for a root
def test_polynomials(BB, monomials, x, y, UU, XX, YY, xx, yy, modulus):
    polynomials = []
    for jj in range(len(monomials)):
        poll.append(0)
        for ii in range(len(monomials)):
            poll[-1] += monomials[ii](x*y+1,x,y) * BB[jj, ii] / monomials[ii](UU,XX,YY)

    for polynomial in polynomials:
        if polynomial(xx,yy) % modulus != 0:
            print "root not working on polynomial"
            return false
    return true

# test if matrix is triangular
def matrix_is_triangular(BB):
    for ii in range(BB.dimensions()[0]):
        if BB[ii, ii] == 0:
            print "zero detected", ii
        for jj in range(ii + 1, BB.dimensions()[0]):
            if BB[ii,jj] != 0:
                print "not triangular", ii, jj

# test if matrix has full rank
def matrix_full_rank(BB):
    a = BB.dimensions()
    if a[0] != a[1] or a[0] != BB.rank():
        print "matrix is not full rank"

# display stats on helpful vectors
def helpful_vectors(BB, modulus):
    nothelpful = 0
    for ii in range(BB.dimensions()[0]):
        if BB[ii,ii] >= modulus:
            nothelpful += 1

    print nothelpful, "/", BB.dimensions()[0], " vectors are not helpful"

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
            xshift = x^ii * modulus^(mm - kk) * polZ(u, x, y)^kk
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
            yshift = y^jj * polZ(u, x, y)^kk * modulus^(mm - kk)
            yshift = Q(yshift).lift()
            gg.append(yshift) # substitution
    
    # y-shifts monomials
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            monomials.append(u^kk * y^jj)

    '''
    # IF WE USE ABOVE WE GET BAD BOUND ON DET
    # IF WE USE BELOW WE GET ONLY HELPFUL VECTORS

    # y-shifts (selected by Herrman and May)
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            yshift = y^jj * polZ(u, x, y)^kk * modulus^(mm - kk)
            yshift = Q(yshift).lift()
            yshift = yshift(u*UU, x*XX, y*YY)
            #gg.append(yshift) # substitution
            # test if helpful
            if YY^jj * UU^kk * modulus^(mm-kk) < modulus^mm:
                print "yshift(",jj,kk,") helpful"
                gg.append(yshift) # substitution
                monomials.append(y*jj * u^kk)
    '''
    # construct lattice B
    nn = len(monomials)

    BB = Matrix(ZZ, nn)

    for ii in range(nn):

        BB[ii, 0] = gg[ii](0, 0, 0)

        for jj in range(1, ii + 1):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](UU,XX,YY)

    # debug
    if debug:
        matrix_is_triangular(BB)
        matrix_full_rank(BB)
        helpful_vectors(BB, modulus^mm)
    
    # check on determinant's bound
    if debug:
        det = BB.det()
        bound = modulus^(mm*(nn-1))
        bound = bound / (nn*(2^nn))^((nn-1)/2)
        if det >= bound:
            print "we don't have det < bound"
            print "size of det - bound = ", int(log(abs(det) / abs(bound))/log(2))
        else:
            print "det < bound"

    # now we get rid of u
    PR.<x,y> = PolynomialRing(ZZ)

    # debug
    if debug:
        test_polynomials(BB, monomials, x, y, UU, XX, YY, xx, yy, modulus^mm)

    # LLL
    BB = BB.LLL()

    # debug
    if debug:
        test_polynomials(BB, monomials, x, y, UU, XX, YY, xx, yy, modulus^mm)

    '''
    # shortest vectors to polynomials
    # this approach doesn't work... we need more vectors
    pol1 = pol2 = 0

    for ii in range(nn):
        pol1 += monomials[ii](x*y+1,x,y) * BB[0, ii] / monomials[ii](UU,XX,YY)
        pol2 += monomials[ii](x*y+1,x,y) * BB[1, ii] / monomials[ii](UU,XX,YY)
    '''

    # find two vectors we can work with
    pols = []
    for ii in range(nn):
        pols.append(0)
        for jj in range(nn):
            pols[-1] += monomials[jj](x*y+1,x,y) * BB[ii, jj] / monomials[jj](UU,XX,YY)
        if pols[-1](xx,yy) != 0:
            pols.pop()
            break

    pol1 = pol2 = 0

    for ii, pol in enumerate(pols):
        for jj in range(ii + 1, len(pols)):
            if gcd(pol, pols[jj]) == 1:
                pol1 = pol
                pol2 = pols[jj]
                break

    # failure
    if pol1 == pol2 == 0:
        print "failure"
        return 0, 0

    # resultant
    PR.<x> = PolynomialRing(ZZ)
    rr = pol1.resultant(pol2)
    rr = rr(x, x)
    soly = rr.roots()[0][0]

    ss = pol1(x, soly)
    solx = ss.roots()[0][0]

    return solx, soly


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
if d % 2 == 0: d += 1 # in case d even
while gcd(d, phi) != 1:
    d += 2
e = d.inverse_mod((p-1)*(q-1))

print "d=", d

# Problem put in equation
P.<x,y> = PolynomialRing(Zmod(e))
A = int((N+1)/2)
pol = 1 + x * (A + y)
delta = (2 - sqrt(2)) / 2
tho = (1 - 2 * delta)
m = 7
t = int(tho * m) # we must have m >= t !
X = floor(e^0.292)
Y = 2*floor(e^0.5)

# hard debug
m = 4
t = 2
X = xx + 10
Y = abs(yy) + 10

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