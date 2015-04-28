"""
This code tries reducing the lattice before feeding it to LLL
It is experimental and not documented

It automatically decides if vectors should be remove with the function
`remove_unhelpful`. You can also automatically remove vectors 
(although not recommanded) with the function `manually_remove`

Most of this work has been merged to `boneh_durfee.sage`
"""

import time

debug = True
helpful_only = True

# display stats on helpful vectors
def helpful_vectors(BB, modulus):
    nothelpful = 0
    for ii in range(BB.dimensions()[0]):
        if BB[ii,ii] >= modulus:
            nothelpful += 1

    print nothelpful, "/", BB.dimensions()[0], " vectors are not helpful"

# display matrix picture with 0 and X
def matrix_overview(BB, bound):
    for ii in range(BB.dimensions()[0]):
        a = ('%02d ' % ii)
        for jj in range(BB.dimensions()[1]):
            a += '0' if BB[ii,jj] == 0 else 'X'
            if BB.dimensions()[0] < 60:
                a += ' '
        if BB[ii, ii] >= bound:
            a += '~'
        print a

# tries to remove unhelpful vectors
# we start at current = n-1 (last vector)
def remove_unhelpful(BB, monomials, bound, current):
    # end of our recursive function
    if current == -1:
        return BB

    # we start by checking from the end
    for ii in range(current, -1, -1):
        # if it is unhelpful:
        if BB[ii, ii] >= bound:
            affected_vectors = 0
            affected_vector_index = 0
            # let's check if it affects other vectors
            for jj in range(ii + 1, BB.dimensions()[0]):
                # if another vector is affected:
                # we increase the count
                if BB[jj, ii] != 0:
                    affected_vectors += 1
                    affected_vector_index = jj

            # level:0
            # if no other vectors end up affected
            # we remove it
            if affected_vectors == 0:
                print "* removing unhelpful vector", ii
                BB = BB.delete_columns([ii])
                BB = BB.delete_rows([ii])
                monomials.pop(ii)
                BB = remove_unhelpful(BB, monomials, bound, ii-1)
                return BB

            # level:1
            # if just one was affected we check
            # if it is affecting someone else
            elif affected_vectors == 1:
                affected_deeper = True
                for kk in range(affected_vector_index + 1, BB.dimensions()[0]):
                    # if it is affecting even one vector
                    # we give up on this one
                    if BB[kk, affected_vector_index] != 0:
                        affected_deeper = False
                # remove both it if no other vector was affected and
                # this helpful vector is not helpful enough
                # compared to our unhelpful one
                if affected_deeper and abs(bound - BB[affected_vector_index, affected_vector_index]) < abs(bound - BB[ii, ii]):
                    print "* removing unhelpful vectors", ii, "and", affected_vector_index
                    BB = BB.delete_columns([affected_vector_index, ii])
                    BB = BB.delete_rows([affected_vector_index, ii])
                    monomials.pop(affected_vector_index)
                    monomials.pop(ii)
                    BB = remove_unhelpful(BB, monomials, bound, ii-1)
                    return BB
    # nothing happened
    return BB

def manually_remove(BB, monomials, vectors):
    # nothing to do?
    if len(vectors) == 0:
        return BB
    print "* manually removing vectors:", vectors
    # start from last
    vectors.sort()
    vectors.reverse()
    # removing in basis matrix
    BB = BB.delete_columns(vectors)
    BB = BB.delete_rows(vectors)
    # removing in list of monomials
    for ii in vectors:
        monomials.pop(ii)
    #
    return BB

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
    
    # y-shifts (selected by David Wong)
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            yshift = y^jj * polZ(u, x, y)^kk * modulus^(mm - kk)
            yshift = Q(yshift).lift()
            gg.append(yshift) # substitution
            monomials.append(u^kk * y^jj)

    # construct lattice B
    nn = len(monomials)

    BB = Matrix(ZZ, nn)

    for ii in range(nn):

        BB[ii, 0] = gg[ii](0, 0, 0)

        for jj in range(1, ii + 1):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](UU,XX,YY)
    
    # ERASING ROWS : PROTOTYPE TO GET BETTER BOUNDS
    if helpful_only:
        # automatically remove
        BB = remove_unhelpful(BB, monomials, modulus^mm, nn-1)
        # manually remove
        #BB = manually_remove(BB, monomials, [25,26,27,28,29,30,31,32])
        # reset dimension
        nn = BB.dimensions()[0]
        if nn == 0:
            print "failure"
            return 0,0


    # check if vectors are helpful
    if debug:
        helpful_vectors(BB, modulus^mm)
    
    # check if determinant is correctly bounded
    if debug:
        det = BB.det()
        bound = modulus^(mm*nn)
        if det >= bound:
            print "We do not have det < bound. Solutions might not be found."
            diff = (log(det) - log(bound)) / log(2)
            print "size det(L) - size e^(m*n) = ", floor(diff)
            return 0,0
        else:
            print "det(L) < e^(m*n)"

    # debug: display matrix
    if debug:
        print monomials
        matrix_overview(BB, modulus^mm)

    # LLL
    BB = BB.LLL()

    # vector 1 & 2 -> polynomials 1 & 2
    PR.<w,z> = PolynomialRing(ZZ)

    pol1 = pol2 = 0
    for jj in range(nn):
        pol1 += monomials[jj](w*z+1,w,z) * BB[0, jj] / monomials[jj](UU,XX,YY)
        pol2 += monomials[jj](w*z+1,w,z) * BB[1, jj] / monomials[jj](UU,XX,YY)

    # resultant
    PR.<q> = PolynomialRing(ZZ)
    rr = pol1.resultant(pol2)

    if rr.is_zero() or rr.monomials() == [1]:
        print "failure"
        return pol1, pol2
    
    rr = rr(q, q)

    # solutions
    soly = rr.roots()[0][0]
    print "found for y_0:", soly

    ss = pol1(q, soly)
    solx = ss.roots()[0][0]
    print "found for x_0:", solx

    #
    return solx, soly


############################################
# Test 
##########################################

# RSA gen options (tweakable)
length_N = 1024
length_d = 0.28

# RSA gen (for the demo)
p = next_prime(2^int(round(length_N/2)))
q = next_prime(round(pi.n()*p))
N = p*q
phi = (p-1)*(q-1)
d = int(N^length_d) 
if d % 2 == 0: d += 1
while gcd(d, phi) != 1:
    d += 2
e = d.inverse_mod((p-1)*(q-1))

# Problem put in equation (default)
P.<x,y> = PolynomialRing(Zmod(e))
A = int((N+1)/2)
pol = 1 + x * (A + y)

# and the solutions to be found (for the demo)
yy = (-p -q)/2
xx = (e * d - 1) / (A + yy)

#
# Default values 
# you should tweak delta and m. X should be OK as well
# 
delta = 0.28              # < 0.292 (Boneh & Durfee's bound)
X = 2*floor(N^delta)      # this _might_ be too much
Y = floor(N^(1/2))        # correct if p, q are ~ same size
m = 12                     # bigger is better (but takes longer)
t = int((1-2*delta) * m)  # optimization from Herrmann and May
# Checking bounds (for the demo)
print "=== checking values ==="
print "* |y| < Y:", abs(yy) < Y
print "* |x| < X:", abs(xx) < X
print "* d < N^0.292", d < N^(0.292)
print "* size of d:", int(log(d)/log(2))
print "* size of N:", int(log(N)/log(2))
print "* delta:", delta

# boneh_durfee
print "=== running algorithm ==="
start_time = time.time()
solx, soly = boneh_durfee(pol, e, m, t, X, Y)

# Checking solutions (for the demo)
if xx == solx and yy == soly:
    print "\n=== the solutions are correct ==="
else:
    print "=== FAIL ==="

# Stats
print("=== %s seconds ===" % (time.time() - start_time))
