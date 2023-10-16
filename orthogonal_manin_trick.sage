

# def dec_into_elem_matrices(a,c):
#     A = Matrix(2, 2, [a, 1, c, 0])
#     a, b, c, d = A.list()
#     dec = []
#     while c != 0:
#         q, r = quo_rem_in_gaussian_integers(a,c)
#         dec.append(Matrix(2,2,[1, -q, 0, 1]).inverse())
#         dec.append(Matrix(2,2,[0, -1, 1, 0]))
#         a, b, c, d = -c, -d, r, b-q*d
#     dec.append(Matrix(2, 2, [a, b, c, d]))
#     assert prod(dec) == A
#     return dec

def quo_rem_in_gaussian_integers(a,b):
    if b == 0:
        raise ZeroDivisionError
    K = a.parent()
    q = ((a/b)[0]).round() + ((a/b)[1]).round()*K.gen()
    r = a - q*b
    assert r.norm() < b.norm()
    return (q, r)

# computes the continued fraction of a/c, where a and c belong to Z[i]
def continued_fraction(a, c):
    a_over_c = a/c
    cf = []
    while c != 0:
        q, r = quo_rem_in_gaussian_integers(a,c)
        cf.append(q)
        a = c
        c = r
    assert a_over_c == eval_cf(cf)
    return cf

# given a continued fraction, computes its value as an element of Q(i); only for debugging
def eval_cf(cf):
    if len(cf) == 1:
        return cf[0]
    if len(cf) == 2:
        return cf[0] + 1/cf[1]
    new_cf = [cf[i+1] for i in range(len(cf)-1)]
    return eval_cf([cf[0],eval_cf(new_cf)])

# computes the convergents of a continued fraction, following Hardy-Wright theorem 149
def compute_convergents(cf):
    N = len(cf)
    if N == 1:
        return([[cf[0]], [1]])
    p = [cf[0], cf[1] * cf[0] +1]
    q = [1, cf[1]]
    for n in range(2, N):
        p.append(cf[n] * p[n-1] + p[n-2])
        q.append(cf[n] * q[n-1] + q[n-2])
    assert all([p[n]*q[n-1] - p[n-1]*q[n] == (-1)^(n-1) for n in range(1,N)])
    #assert p[-1]/q[-1] == eval_cf(cf)
    return(p, q)

# given the convergents [p_i] and [q_j] of a continued fraction for alpha, computes the matrices that Mi such that {Infinity, alpha} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def compute_Ms(p, q):
    Ms = []
    # Ms.append(Matrix(2,2,[1,0,0,1])) # j = -1, not necessary since we do (Infinity, alpha) instead of (0, alpha) 
    Ms.append(Matrix(2,2,[-p[0], 1, -q[0], 0])) # j = 0
    for j in range(1, len(p)):
        Ms.append(Matrix(2,2,[(-1)^(j-1)*p[j],p[j-1],(-1)^(j-1)*q[j],q[j-1]]))
    assert all([A.det() == 1 for A in Ms])
    return Ms

# given a, c in Z[i] returns the matrices that Mi such that {Infinity, a/c} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def matrices_for_unimodular_path(a,c):
    cf = continued_fraction(a, c)
    p, q = compute_convergents(cf)
    return compute_Ms(p, q)


# checking
def act_mobius(a,b,c,d,tau):
    return act_matrix(matrix([[a,b],[c,d]]), tau)

def act_matrix(gamma, tau):
    if isinstance(tau, list):
        return [act_matrix(gamma, t) for t in tau]
    a, b, c, d = gamma.list()
    K = tau.parent()
    if tau == 0:
        return K(b / d) if d else Infinity
    elif tau == Infinity:
        return K(a / c) if c else Infinity
    else:
        return (K(a)*tau + K(b)) / (K(c)*tau + K(d))


# K.<i> = QuadraticField(-1)
# while True:
#     a = K.random_element()
#     c = K.random_element()
#     if c == 0:
#         continue
#     break

# Ms = matrices_for_unimodular_path(a,c)

# alpha = a/c
# print('alpha = ', alpha)
# for M in Ms:
#     print('( ',act_matrix(M,0),', ', act_matrix(M,Infinity),')')
