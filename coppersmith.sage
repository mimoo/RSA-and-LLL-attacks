# INIT
R.<x> = ZZ[]

"""pol is of degree d
NN modulus a multiple of b
lower bound beta, b >= NN^beta
epsilon <= beta / 7
"""
def coppersmith_univariate(pol, NN, bb, beta, epsilon):
    # init
    dd = pol.degree()
    # choose m and t
    mm = ceil(beta**2 / dd * epsilon)
    tt = floor(dd * mm * ((1/beta) - 1))
    # compute polynomials
    gg = {}
    for ii in range(mm):
        gg[ii] = {}
        for jj in range(dd):
            g[ii][jj] = x
    #
    return mm

# TESTS

m = 010101010101 # we know part of the message
ct = 01010101010 # we know the ciphertext
pol = (m - x)**3 - ct

print(coppersmith_univariate(pol, NN, NN, 1, 1))


