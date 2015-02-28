def boneh_durfee(pol, modulus, delta, mm, tt, XX, YY):
    """
    Boneh and Durfee revisited by Herrmann and May
    d = N^delta
    """

    # change ring of pol and x
    polZ = pol.change_ring(ZZ)
    x, y = polZ.parent().gens()

    # compute polynomials
    # x-shifts
    gg = []
    for kk in range(mm + 1):
        for ii in range(mm - kk + 1):
            #print kk, ii
            gg.append((x * XX)^ii * modulus^(mm - kk) * polZ(x * XX, y * YY)^kk)

    # y-shifts (selected by Herrman and May)
    ''' DOESNT WORK'''
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            gg.append((y * YY)^jj * polZ(x * XX, y * YY)^kk * modulus^(mm - kk))
    gg.sort()

    # unravelled linerization (Herrman and May)
    monomials = [] # jusqu'à nn

    # x-shift
    for ii in range(mm + 1):
        for jj in range(ii + 1):
            monomials.append(x^ii * y^jj)

    # y-shift
    for jj in range(1, tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            monomials.append((x*y)^kk * y^jj)

    # construct lattice B
    nn = len(monomials)
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj, monomial in enumerate(monomials):
            if jj == 0:
                BB[ii, jj] = gg[ii].coefficient({x:0,y:0})
            else:
                BB[ii, jj] = gg[ii].monomial_coefficient(monomial)
    print BB
    return 0

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
# Test 
##########################################


e = 2
P.<x,y> = PolynomialRing(Zmod(e))
pol = x + y
delta = 4
m = 2
t = 1
X = 3
Y = 5

result = boneh_durfee(pol, e, delta, m, t, X, Y)