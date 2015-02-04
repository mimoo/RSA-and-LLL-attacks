def coppersmith_univariate(pol, bb, beta):
    """Howgrave-Graham revisited method
    using with epsilon
    """
    # init
    dd = pol.degree()
    NN = pol.parent().characteristic()
    
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
    BB = Matrix(ZZ, nn)

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
    """ ok this step seems unecessary, delete when sure
    """
    normn = norm(BB[0])
    norm_index = 0

    for ii in range(1, nn):
        if norm(BB[ii]) < normn:
            normn = norm(BB[ii])
            norm_index = ii
    print(norm_index)

    # transform shortest vector in polynomial    
    new_pol = 0
    for ii in range(nn):
        new_pol += x**ii * BB[norm_index, ii]
    
    # factor polynomial
    roots = new_pol.roots() # doesn't find anything...
    print("roots found", roots)
    # test roots on original pol
    """
    in thesis it says to check root like this:
    gcd(NN, pol(root)) >= NN**beta
    """
    for root in roots:
        if pol(root) == 0:
            return root

    # throw exception?
    return 0
    
# TESTS
# (from http://www.jscoron.fr/cours/mics3crypto/tpcop.pdf)
N = 122840968903324034467344329510307845524745715398875789936591447337206598081
C = 1792963459690600192400355988468130271248171381827462749870651408943993480816

K = Zmod(N)
R.<x> = PolynomialRing(K)
pol = (2**250 + x)**3 - C
M = coppersmith_univariate(pol, N, 1)

