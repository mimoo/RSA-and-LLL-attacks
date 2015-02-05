def LLL(self, float delta=LLL_DEF_DELTA, float eta=LLL_DEF_ETA,
        method=None, float_type=None, int precision=0,
        verbose=False, siegel=False, early_red=False):
    r"""
    `(\delta,\eta)`-LLL reduce this lattice.

    INPUT:

    - ``delta`` -- (default: ``0.99``) parameter `0.25 < \delta < 1.0`
    - ``eta `` -- (default: ``0.51``) parameter `0.5 \leq \eta <
      \sqrt{\delta}`
    - ``method`` -- (default: ``None``) can be one of the following:

      * ``'wrapper'`` (``None``)
      * ``'proved'``
      * ``'fast'``
      * ``'heuristic'``

    - ``float_type`` -- (default: ``None``) can be one of the following:

      * ``None`` - for automatic choice
      * ``'double'``
      * ``'long double'``
      * ``'dpe'``
      * ``'mpfr'``

    - ``precision`` -- (default: ``0`` for automatic choice) precision
      to use
    - ``verbose`` -- (default: ``False``) be verbose
    - ``siegel`` -- (default: ``False``) use Siegel conditioning
    - ``early_red`` -- (default: ``False``) use early reduction

    OUTPUT:

    Nothing is returned but the internal state is modified.

    EXAMPLES::

        sage: from sage.libs.fplll.fplll import FP_LLL
        sage: A = random_matrix(ZZ,10,10); A
        [   -8     2     0     0     1    -1     2     1   -95    -1]
        [   -2   -12     0     0     1    -1     1    -1    -2    -1]
        [    4    -4    -6     5     0     0    -2     0     1    -4]
        [   -6     1    -1     1     1    -1     1    -1    -3     1]
        [    1     0     0    -3     2    -2     0    -2     1     0]
        [   -1     1     0     0     1    -1     4    -1     1    -1]
        [   14     1    -5     4    -1     0     2     4     1     1]
        [   -2    -1     0     4    -3     1    -5     0    -2    -1]
        [   -9    -1    -1     3     2     1    -1     1    -2     1]
        [   -1     2    -7     1     0     2     3 -1955   -22    -1]

        sage: F = FP_LLL(A)
        sage: F.LLL(method="wrapper")
        sage: L = F._sage_(); L
        [   1    0    0   -3    2   -2    0   -2    1    0]
        [  -1    1    0    0    1   -1    4   -1    1   -1]
        [  -2    0    0    1    0   -2   -1   -3    0   -2]
        [  -2   -2    0   -1    3    0   -2    0    2    0]
        [   1    1    1    2    3   -2   -2    0    3    1]
        [  -4    1   -1    0    1    1    2    2   -3    3]
        [   1   -3   -7    2    3   -1    0    0   -1   -1]
        [   1   -9    1    3    1   -3    1   -1   -1    0]
        [   8    5   19    3   27    6   -3    8  -25  -22]
        [ 172  -25   57  248  261  793   76 -839  -41  376]
        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

        sage: set_random_seed(0)
        sage: A = random_matrix(ZZ,10,10); A
        [   -8     2     0     0     1    -1     2     1   -95    -1]
        [   -2   -12     0     0     1    -1     1    -1    -2    -1]
        [    4    -4    -6     5     0     0    -2     0     1    -4]
        [   -6     1    -1     1     1    -1     1    -1    -3     1]
        [    1     0     0    -3     2    -2     0    -2     1     0]
        [   -1     1     0     0     1    -1     4    -1     1    -1]
        [   14     1    -5     4    -1     0     2     4     1     1]
        [   -2    -1     0     4    -3     1    -5     0    -2    -1]
        [   -9    -1    -1     3     2     1    -1     1    -2     1]
        [   -1     2    -7     1     0     2     3 -1955   -22    -1]

        sage: F = FP_LLL(A)
        sage: F.LLL(method="proved")
        sage: L = F._sage_(); L
        [   1    0    0   -3    2   -2    0   -2    1    0]
        [  -1    1    0    0    1   -1    4   -1    1   -1]
        [  -2    0    0    1    0   -2   -1   -3    0   -2]
        [  -2   -2    0   -1    3    0   -2    0    2    0]
        [   1    1    1    2    3   -2   -2    0    3    1]
        [  -4    1   -1    0    1    1    2    2   -3    3]
        [   1   -3   -7    2    3   -1    0    0   -1   -1]
        [   1   -9    1    3    1   -3    1   -1   -1    0]
        [   8    5   19    3   27    6   -3    8  -25  -22]
        [ 172  -25   57  248  261  793   76 -839  -41  376]

        sage: L.is_LLL_reduced()
        True
        sage: L.hermite_form() == A.hermite_form()
        True

        sage: A = random_matrix(ZZ,10,10,x=-(10^5),y=10^5)
        sage: f = FP_LLL(A)
        sage: f.LLL(method="fast")
        sage: L = f._sage_()
        sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
        True
        sage: L.hermite_form() == A.hermite_form()
        True

        sage: set_random_seed(0)
        sage: A = random_matrix(ZZ,10,10); A
        [   -8     2     0     0     1    -1     2     1   -95    -1]
        [   -2   -12     0     0     1    -1     1    -1    -2    -1]
        [    4    -4    -6     5     0     0    -2     0     1    -4]
        [   -6     1    -1     1     1    -1     1    -1    -3     1]
        [    1     0     0    -3     2    -2     0    -2     1     0]
        [   -1     1     0     0     1    -1     4    -1     1    -1]
        [   14     1    -5     4    -1     0     2     4     1     1]
        [   -2    -1     0     4    -3     1    -5     0    -2    -1]
        [   -9    -1    -1     3     2     1    -1     1    -2     1]
        [   -1     2    -7     1     0     2     3 -1955   -22    -1]

        sage: F = FP_LLL(A)
        sage: F.LLL(method="fast", early_red=True)
        sage: L = F._sage_(); L
        [   1    0    0   -3    2   -2    0   -2    1    0]
        [  -1    1    0    0    1   -1    4   -1    1   -1]
        [  -2    0    0    1    0   -2   -1   -3    0   -2]
        [  -2   -2    0   -1    3    0   -2    0    2    0]
        [   1    1    1    2    3   -2   -2    0    3    1]
        [  -4    1   -1    0    1    1    2    2   -3    3]
        [   1   -3   -7    2    3   -1    0    0   -1   -1]
        [   1   -9    1    3    1   -3    1   -1   -1    0]
        [   8    5   19    3   27    6   -3    8  -25  -22]
        [ 172  -25   57  248  261  793   76 -839  -41  376]

        sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
        True
        sage: L.hermite_form() == A.hermite_form()
        True

        sage: set_random_seed(0)
        sage: A = random_matrix(ZZ,10,10); A
        [   -8     2     0     0     1    -1     2     1   -95    -1]
        [   -2   -12     0     0     1    -1     1    -1    -2    -1]
        [    4    -4    -6     5     0     0    -2     0     1    -4]
        [   -6     1    -1     1     1    -1     1    -1    -3     1]
        [    1     0     0    -3     2    -2     0    -2     1     0]
        [   -1     1     0     0     1    -1     4    -1     1    -1]
        [   14     1    -5     4    -1     0     2     4     1     1]
        [   -2    -1     0     4    -3     1    -5     0    -2    -1]
        [   -9    -1    -1     3     2     1    -1     1    -2     1]
        [   -1     2    -7     1     0     2     3 -1955   -22    -1]

        sage: F = FP_LLL(A)
        sage: F.LLL(method="heuristic")
        sage: L = F._sage_(); L
        [   1    0    0   -3    2   -2    0   -2    1    0]
        [  -1    1    0    0    1   -1    4   -1    1   -1]
        [  -2    0    0    1    0   -2   -1   -3    0   -2]
        [  -2   -2    0   -1    3    0   -2    0    2    0]
        [   1    1    1    2    3   -2   -2    0    3    1]
        [  -4    1   -1    0    1    1    2    2   -3    3]
        [   1   -3   -7    2    3   -1    0    0   -1   -1]
        [   1   -9    1    3    1   -3    1   -1   -1    0]
        [   8    5   19    3   27    6   -3    8  -25  -22]
        [ 172  -25   57  248  261  793   76 -839  -41  376]

        sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
        True
        sage: L.hermite_form() == A.hermite_form()
        True

        sage: set_random_seed(0)
        sage: A = random_matrix(ZZ,10,10); A
        [   -8     2     0     0     1    -1     2     1   -95    -1]
        [   -2   -12     0     0     1    -1     1    -1    -2    -1]
        [    4    -4    -6     5     0     0    -2     0     1    -4]
        [   -6     1    -1     1     1    -1     1    -1    -3     1]
        [    1     0     0    -3     2    -2     0    -2     1     0]
        [   -1     1     0     0     1    -1     4    -1     1    -1]
        [   14     1    -5     4    -1     0     2     4     1     1]
        [   -2    -1     0     4    -3     1    -5     0    -2    -1]
        [   -9    -1    -1     3     2     1    -1     1    -2     1]
        [   -1     2    -7     1     0     2     3 -1955   -22    -1]

        sage: F = FP_LLL(A)
        sage: F.LLL(method="heuristic", early_red=True)
        sage: L = F._sage_(); L
        [   1    0    0   -3    2   -2    0   -2    1    0]
        [  -1    1    0    0    1   -1    4   -1    1   -1]
        [  -2    0    0    1    0   -2   -1   -3    0   -2]
        [  -2   -2    0   -1    3    0   -2    0    2    0]
        [   1    1    1    2    3   -2   -2    0    3    1]
        [  -4    1   -1    0    1    1    2    2   -3    3]
        [   1   -3   -7    2    3   -1    0    0   -1   -1]
        [   1   -9    1    3    1   -3    1   -1   -1    0]
        [   8    5   19    3   27    6   -3    8  -25  -22]
        [ 172  -25   57  248  261  793   76 -839  -41  376]

        sage: L.is_LLL_reduced(eta=0.51,delta=0.99)
        True
        sage: L.hermite_form() == A.hermite_form()
        True
    """
    _check_delta(delta)
    _check_eta(eta)
    _check_precision(precision)

    cdef LLLMethod method_
    if method == "wrapper" or method is None:
        method_ = LM_WRAPPER
    elif method == "proved":
        method_ = LM_PROVED
    elif method == "heuristic":
        method_ = LM_HEURISTIC
    elif method == "fast":
        method_ = LM_FAST
    else:
        raise ValueError("method '{}' unknown".format(method))

    cdef int flags = LLL_DEFAULT

    if verbose:
        flags |= LLL_VERBOSE
    if early_red:
        flags |= LLL_EARLY_RED
    if siegel:
        flags |= LLL_SIEGEL

    if float_type is None and method_ == LM_FAST:
        float_type = "double"

    if method_ == LM_WRAPPER and check_float_type(float_type) != FT_DEFAULT:
        raise ValueError("fpLLL's LLL wrapper function requires "
                         "float type None")
    if method_ == LM_FAST and \
            check_float_type(float_type) not in (FT_DOUBLE, FT_LONG_DOUBLE):
        raise ValueError("fpLLL's LLL fast function requires "
                         "float type 'double' or 'long double'")

    sig_on()
    cdef int r = lllReduction(self._lattice[0], delta, eta, method_,
                              check_float_type(float_type), precision, flags)
    sig_off()
    if r:
        raise RuntimeError( str(getRedStatusStr(r)) )
