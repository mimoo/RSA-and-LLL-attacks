
# INTRO

## reminder on RSA:

it does:
- encrypt 
- sign

it comes from:
- generate p, q primes; N = p * q; (e, N) and (d, N)
- encryption: m^e = c [N]
- decryption: c^d = m [N]

how it works:
- recall, euler theorem: a^phi(N) = 1 [N]
- so we want to find e d = k phi(N) + 1 because:
	a^(e d) = a^phi(N) + 1 = a [N]
- we can rewrite that as e d = 1 [phi(N)]
- we need to find one invertible e in (Zphi(N))* multiplicativ group
- then invert it (extended euclidian algorithm)

## attacks
* what sorts of attack we do today:
- implementation, side-channel
- mathematical (it's hard)

* mathematical is : 
- factorize N (=pq) which allows us to recover d
- recover d directly
- recover m when we have e^m (mod N)

* so what we could do is mathematical but on relaxed RSA

# PLAN OF THE TALK:

* Hastad Broadcast Attack

* Coppersmith small roots (x < N^1/e)

* Boneh-Durfee bound (d < N^0.292)

# Hastad Broadcast attack

* explanation

* how to extend that?

# COPPERSMITH

* formulate the problem for RSA

* explanation of the idea of coppersmith attack

* high view of coppersmith technique

* lattices

* low level of coppersmith

# Boneh and Durfee