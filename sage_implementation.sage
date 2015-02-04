from sage.misc.misc import verbose 
from sage.matrix.constructor import Matrix 
from sage.rings.all import RR 

N = self.parent().characteristic() 

f = self.change_ring(ZZ) 
P,(x,) = f.parent().objgens() 

beta = RR(beta) 

delta = f.degree() 
epsilon = beta/8 

# our choice of epsilon simplifies 
#  m = ceil( max(beta^2/(delta * epsilon), 7*beta/delta) ) 
# to 
m = int( ( beta**2/(delta * epsilon) ).ceil() ) 
verbose("m = %d"%m, level=2) 

t = int( ( delta*m*(1/beta -1) ).floor() ) 
verbose("t = %d"%t, level=2) 

g  = [x**j * N**(m-i) * f**i for i in range(m) for j in range(delta) ] 
g.extend([x**i * f**m for i in range(t)]) # h 

if X is None: 
    X = int( (N**(beta**2/delta - epsilon)).ceil() ) 
verbose("X = %s"%X, level=2) 

B = Matrix(ZZ, len(g), delta*m + max(delta,t) ) 
for i in range(B.nrows()): 
    for j in range( g[i].degree()+1 ): 
        B[i,j] = g[i][j]*X**j 

B =  B.LLL(**kwds) 

f = sum([ZZ(B[0,i]//X**i)*x**i for i in range(B.ncols())]) 
R = f.roots() 

return [r for r,m in R if r<=X] 
