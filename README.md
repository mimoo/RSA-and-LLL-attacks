# Mathematical attacks on RSA

This repo will host implementations of different RSA attacks using lattice reduction techniques (in particular **LLL**).

First, we'll see how **Coppersmith** found out that you could use lattice reduction techniques to attack a relaxed RSA (we know parts of the message, or we know parts of one of the prime).

Second we'll see how **Boneh and Durfee** used coppersmith to factor the RSA modulus under certain assumptions.

# Coppersmith.sage

I've implemented the work of **Coppersmith** (to be correct the reformulation of his attack by **Howgrave-Graham**) in [coppersmith.sage](coppersmith.sage).

I've used it in two examples in the code:

## Stereotyped messages

For example if **you know the most significant bits of the message**. You can **find the rest of the message** with this method.

The usual RSA model is this one: you have a ciphertext `c` a modulus `N` and a public exponent `e`. Find `m` such that `m^e = c mod N`.

Now, this is the **relaxed model** we can solve: you have `c = (m + x)^e`, you know a part of the message, `m`, but you don't know `x`.
For example the message is always something like "*the password today is: [password]*".
Coppersmith says that if you are looking for `N^1/e` of the message it is then a `small root` and you should be able to find it pretty quickly.

let our polynomial be `f(x) = (m + x)^e - c` which has a root we want to find `modulo N`. Here's how to do it with my implementation:

```
dd = f.degree()
beta = 1
epsilon = beta / 7
mm = ceil(beta**2 / (dd * epsilon))
tt = floor(dd * mm * ((1/beta) - 1))
XX = ceil(N**((beta**2/dd) - epsilon))
roots = coppersmith_howgrave_univariate(f, N, beta, mm, tt, XX)
```

You can play with the values until it finds the root. The default values should be a good start. If you want to tweak:
* beta is always 1 in this case.
* `XX` is your upper bound on the root. **The bigger is the unknown, the bigger XX should be**. And the bigger it is... the more time it takes.

## Factoring with high bits known

Another case is factoring `N` knowing high bits of `q`.

The Factorization problem normally is: give `N = pq`, find `q`. In our **relaxed** model we know an approximation `q'` of `q`.

Here's how to do it with my implementation:

let `f(x) = x - q'` which has a root modulo q

```
beta = 0.5
dd = f.degree()
epsilon = beta / 7
mm = ceil(beta**2 / (dd * epsilon))
tt = floor(dd * mm * ((1/beta) - 1))
XX = ceil(N**((beta**2/dd) - epsilon)) + 1000000000000000000000000000000000
roots = coppersmith_howgrave_univariate(f, N, beta, mm, tt, XX)
```

What is important here if you want to find a solution:

*  we should have `q >= N^beta`
* as usual `XX` is the upper bound of the root, so the difference should be: |diff| < X

note: diff = | q - q` |