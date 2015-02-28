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
    for kk in range(mm):
        for ii in range(mm - kk):
            gg.append((x * XX)^ii * modulus^(mm - kk) * polZ(x * XX, y * YY)^kk)
    gg.sort()

    # y-shifts (selected by Herrman and May)
    hh = []
    for jj in range(tt):
        for kk in range(floor(mm/tt) * jj, mm):
            hh.append((y * YY)^jj * polZ(x * XX, y * YY)^kk * modulus^(mm - kk))
    hh.sort()

    # unravelled linerization (Herrman and May)
    nn = (mm + 1) * (mm + 2) / 2 +  tt * (mm + 1)
    monomials = [] # jusqu'à nn

    # monomials
    # x-shift
    for ii in range(mm + 1):
        for jj in range(ii + 1):
            monomials.append(x^ii * y^jj)

    # y-shift
    ''' NOT WORKING '''
    for jj in range(tt + 1):
        for kk in range(floor(mm/tt) * jj, mm + 1):
            print kk, floor(mm/tt)*jj, mm+1
            monomials.append((x*y)^kk * y^jj)

    print monomials
    return 0

    # construct lattice B
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj, monomial in enumerate(monomials):
            if jj == 0:
                BB[ii, jj] = gg[ii].coefficient({x:0,y:0})
            else:
                BB[ii, jj] = gg[ii].monomial_coefficient(monomial)

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


e = 5
P.<x,y> = PolynomialRing(Zmod(e))
pol = x + y
delta = 4
m = 2
t = 1
X = 5
Y = 5

result = boneh_durfee(pol, e, delta, m, t, X, Y)