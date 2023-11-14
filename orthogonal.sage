# from tqdm import tqdm
from sage.rings.padics.precision_error import PrecisionError
from multiprocessing import Process, Manager, Pool
from concurrent import futures


# Related to Manin Trick

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

### Define functions
def all_elements_of_norm(F, n):
    if n == 0:
        return [F(0)]
    else:
        return [u * beta for beta in F.elements_of_norm(n) for u in units]

def vectors_in_lattice(n):
    resp = []
    resm = []
    for mac in range(1,n+1):
        RHS = n - mac
        for mc in divisors(mac):
            a = mac // mc
            for beta in all_elements_of_norm(F, RHS):
                if a % p == 0 and \
                   mc % p == 0 and \
                   (beta / p).is_integral():
                    continue
                else:
                    new_mat = Matrix(F,2,2,[beta.conjugate(), mc, a, -beta])
                    if a % 4 == 1:
                        resp.append(new_mat)
                    elif a % 4 == 3:
                        resm.append(-new_mat)
    return resp, resm, len(resp) - len(resm)

def reduce2_mod(f):
    t = f.parent().gen()
    return sum(reduce_mod(o) * t**i for o, i in zip(f.coefficients(), f.exponents()))

def reduce_mod(f):
    t = f.parent().gen()
    return sum((o % (p**M)) * t**i for o, i in zip(f.coefficients(), f.exponents()))

def inv(FF):
    R = FF.parent()
    a0inv = ~Rp(FF(0)(0))
    pw0 = 1-a0inv * FF
    pw = 1
    ans = 0
    for n in range(M):
        ans += pw
        pw *= pw0
    ans *= a0inv
    return ans

@parallel(ncpus=ncpus)
def compute_level1_contribution(A, Ap, sgn, Rp):
    Rpol = PolynomialRing(F,2,names='u,v')
    t1, t2 = Rpol.gens()
    pol = vector(Rpol, [-1,t1]) * A * vector(Rpol, [t2,1])
    polp = pol.change_ring(phi).change_ring(GF(p))
    R1 = PolynomialRing(GF(p),names='x')
    x = R1.gen()
    to_x = pol.parent().hom([x,x], codomain=R1,check=False)
    d1, d2 = polp.degrees()
    if d1 == 0 and d2 == 0:
        j1, j2, tau1, tau2 = p, p, t1, t2
    if d1 == 1 and d2 == 0:
        j1, j2 = polp.change_ring(to_x).roots()[0][0].lift(), p
        tau1, tau2 = j1 + 1 / t1, t2
    if d1 == 0 and d2 == 1:
        j1, j2 = p, polp.change_ring(to_x).roots()[0][0].lift()
        tau1, tau2 = t1, j2 + 1 / t2
    if d1 == 1 and d2 == 1:
        ff = polp.factor()
        assert len(ff) == 2,'Does not factor'
        f1 = ff[1][0]
        f2 = ff[0][0]
        assert f1.degrees() == (1,0) and f2.degrees() == (0,1)
        j1 = f1.change_ring(to_x).roots()[0][0].lift()
        j2 = f2.change_ring(to_x).roots()[0][0].lift()
        tau1, tau2 = j1 + 1 / t1, j2 + 1 / t2
    ans = (pol.parent()(pol(tau1,tau2) * t1**d1 * t2**d2))
    ans = ans.change_ring(phi)
    while j1 < 0:
        j1 += p
    while j2 < 0:
        j2 += p
    ans = ans.change_ring(Rp)
    ans = map_poly(ans)
    return j1, j2, ans, sgn

def inv_par(kyFF):
    ky, FF = kyFF
    R = FF.parent()
    a0inv = ~Rp(FF(0)(0))
    pw0 = 1-a0inv * FF
    pw = 1
    ans = 0
    for n in range(M):
        ans += pw
        pw *= pw0
    ans *= a0inv
    return ky, ans

def Level1(V):
    res = {(i,j) : [Ruv(1), Ruv(1)] for i in range(p+1) for j in range(p+1)}
    dg = {(i,j) : 0 for i in range(p+1) for j in range(p+1)}
    pM = p**M
    input_vec = []
    for (A, B) in zip(*V):
        Ap = A.apply_map(lambda x : Rp(phi(x).lift()))
        Bp = B.apply_map(lambda x : Rp(phi(x).lift()))
        input_vec.extend([(A,Ap,1,Rp), (B,Bp,-1,Rp)])
    for A, Ap, sgn, _ in input_vec:
        i, j, FF, sgn = compute_level1_contribution(A, Ap, sgn, Rp)
        if sgn == 1:
            res[i,j][0] *= FF
        else:
            res[i,j][1] *= FF
        dg[i,j] += sgn
    with futures.ProcessPoolExecutor() as executor:
        # Calculate inverses
        print('Calculating inverses...')
        future = {executor.submit(inv_par, (ky, val[1])) : ky for ky, val in res.items()}
        for fut in futures.as_completed(future):
            ky, val = fut.result()
            res[ky] = res[ky][0] * val
    return res, dg

@cached_function
def ApplySingle(A, i, z, check=True):
    a, b, c, d = A.change_ring(ZZ).list()
    if check:
        assert A.determinant().valuation() == 0, "Error in embeddings?"
    if i == p:
        u = GF(p)(a)
        v = GF(p)(c)
    else:
        u = GF(p)(a*i+b)
        v = GF(p)(c*i+d)
    if not check and (u == 0 and v == 0):
        return 1, None
    ii = p if v == 0 else (u / v).lift()
    assert 0 <= ii <= p
    if i == p:
        if ii == p:
            r = Rp(c/a) * z
            substx = ~Rp(a) * (d*z - b)
        else:
            r = Rp((-c*ii + a)/c) * z
            substx = -~Rp(c) * (d + (d*ii-b)*z)
    else:
        if ii == p:
            r = Rp((c*i+d)/(a*i+b)) * z
            substx = -Rp(~(a*i+b)) * (a-c*z)
        else:
            r = Rp(-(ii - (a*i+b)/(c*i+d))) * z
            substx = Rp(~(c*i+d)) * (-c + (a-c*ii)*z)
    if r != 0:
        substx *= sum(r**k for k in range(M+1))
    return ii, substx

def Transform(outky):
    global gF, input_list
    res = Ruv(1)
    resinv = Ruv(1)
    for inky, s1, s2, sgn in input_list[outky]:
        f = gF[inky]
        newres = sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] for fi, i in zip(f.coefficients(), f.exponents()))
        if sgn == 1:
            res *= newres
        else:
            resinv *= newres
    return res * inv(resinv)

def Next(F):
    global gF
    gF = F
    res = {(i,j) : 1 for i in range(p+1) for j in range(p+1)}
    if parallelize:
        with futures.ProcessPoolExecutor() as executor:
            # Iteration
            future_dict = {executor.submit(Transform, outky) : outky for outky, inps in input_list.items()}
            for fut in futures.as_completed(future_dict):
                outky = future_dict[fut]
                res[outky] = fut.result()
    else:
        for outky, inps in input_list.items():
            res[outky] = Transform(outky)
    return res

def mul_dict(x,y):
    return {ky : x[ky] * y[ky] for ky in x}

def RMC(F):
    FF = F
    res = [F]
    for j in range(1,M):
        print(f'Iteration {j}')
        t = walltime()
        FF = Next(Next(FF))
        res.append(FF)
        print(f'..done in {walltime(t)} seconds.')
    print(f'Now computing product...')
    t = walltime()
    if parallelize:
        with futures.ProcessPoolExecutor() as executor:
            while len(res) > 1:
                future_dict = {executor.submit(mul_dict, res[i], res[i+1]) : i for i in range(0,len(res)-1, 2) }
                res = [] if len(res) % 2 == 0 else [res[-1]]
                for fut in futures.as_completed(future_dict):
                    res.append(fut.result())
        ans = res[0]
    else:
        ans = {ky : prod(r[ky] for r in res) for ky in res[0] }
    print(f'..done in {walltime(t)} seconds.')
    R0 = Ruv.base_ring()
    phi = R0.hom([ZZ['u'].gen()],base_map=lambda x:x.lift(),check=False)
    ans = {ky : f.change_ring(phi) for ky, f in ans.items()}
    return ans

def solve_quadratic(f, K = None, return_all = False):
    a,b,c = f[2], f[1], f[0]
    if f.degree() != 2:
        raise ValueError
    disc = b*b-4*a*c
    discrt = our_sqrt(disc,K,False)
    if return_all:
        return [(-b + discrt)/(2*a), (-b - discrt)/(2*a)]
    else:
        return (-b + discrt)/(2*a)

def our_sqrt(xx,K = None,return_all = False):
    if K is None:
        K = xx.parent()
    else:
        xx = K(xx)
    if xx == 0:
        if return_all:
            return [xx]
        else:
            return xx
    p=K.base_ring().prime()
    prec = K.precision_cap()
    valp = xx.valuation()
    valpi = xx.ordp()
    try:
        eK = K.ramification_index()
    except AttributeError:
        eK = 1
    if valp * eK % 2 != 0:
        if return_all:
            return []
        else:
            raise ValueError('Not a square')
    x = K.uniformizer()**(-valp) * xx
    try:
        z = K.unramified_generator()
    except AttributeError:
        z = K.gen()
    deg = K.residue_class_field().degree()
    found = False
    if p != 2:
        ppow = p if p != 2 else 8
        minval = 1 if p != 2 else 3
        for avec in cartesian_product_iterator([srange(ppow) for _ in range(deg)]):
            y0 = K(avec[0])
            for a in avec[1:]:
                y0 = y0*z + K(a)
            if (y0**2-x).valuation() >= minval:
                found = True
                break
        if found == False:
            if return_all:
                return []
            else:
                raise ValueError('Not a square: %s'%x)
        y, y1 = 0, y0
        while y != y1:
            y, y1 = y1, (y1+x/y1)/2
    else:
        ppow = 8
        minval = 1
        for avec in cartesian_product_iterator([srange(ppow) for _ in range(deg)]):
            y0 = K(avec[0])
            for a in avec[1:]:
                y0 = y0*z + K(a)
            if (y0**2-y0 + (1-x)/4).valuation() >= minval:
                found = True
                break
        if found == False:
            if return_all:
                return []
            else:
                raise ValueError('Not a square: %s'%x)
        y, y1 = 0, y0
        while y != y1:
            y, y1 = y1, (y1**2 - (1-x)/4)/(2*y1 - 1) #(y1+x/y1)/2
        y = 2 * y - 1
    ans = K.uniformizer()**(ZZ(valp/2)) * y
    assert ans**2 == xx
    if return_all:
        ans = [ans, -ans]
    return ans

def fixed_point(g, phi):
    a, b, c, d = g.list()
    K = g.parent()
    Kp = phi.codomain()
    x = Kp['x'].gen()
    f = phi(c)*x**2 + phi(d-a)*x - phi(b)
    p = Kp.prime()
    L = Qq(p**2, Kp.precision_cap(),names='b')
    return solve_quadratic(f.change_ring(phi).change_ring(L), L, return_all=True)

def Eval0(L0, tau):
    t0, t1 = tau
    vv = vector([-1, t0])
    ww = vector([t1, 1])
    KK = vv.parent().base_ring()
    ans = prod((vv * num.apply_morphism(phi).apply_map(lambda o : KK(o)) * ww) / (vv * den.apply_morphism(phi).apply_map(lambda o : KK(o)) * ww) for num, den in zip(*L0))
    return ans

def EvalPS(f, x0, x1, prec=M):
    ans = 1
    x0s = [1, x0]
    x1s = [1, x1]
    for _ in range(prec):
        x0s.append(x0 * x0s[-1])
        x1s.append(x1 * x1s[-1])
    return sum(sum(o * x0s[j] for o, j in zip(a.coefficients(), a.exponents()) if j < prec) * x1s[i] for a, i in zip(f.coefficients(), f.exponents()) if i < prec)

def Eval(tau, prec=M):
    global J
    t0, t1 = tau
    ans = 1
    for ky in J:
        x0 = t0 if ky[0] == p else 1/(t0 - ky[0])
        x1 = t1 if ky[1] == p else 1/(t1 - ky[1])
        val = J[ky]
        ans *= EvalPS(val, x0, x1, prec)
    return ans

def good_matrices(m):
    a, b, c, d = m.list()
    F = m.parent().base_ring()
    r = F.gen()
    I = F.ideal(2)
    if F == QQ or F == ZZ:
        raise NotImplementedError
    if c in I:
        return [ (m,-1) ]
    elif d in I:
        m = matrix([[b, a], [d, c]])
        return [ (m, 1) ]
    else:
        for d1, d2 in cartesian_product_iterator([[1,-1,r,-r] for _ in range(2)]):
            delta = -(m.determinant())
            A, B, C, D = m.list()
            y = matrix([[-C, d1], [D, d2]]).determinant()/delta
            if y in I:
                x = matrix([[d1, A], [d2, -B]]).determinant()/delta
                m1 = matrix([[A,x],[C,y]])
                m2 = matrix([[x,B],[y,D]])
                m1.rescale_col(1, -d1**-1)
                m2.rescale_col(0, d2**-1)
                m1.swap_columns(0,1)
                assert m1.determinant() == 1
                assert m2.determinant() == 1
                return [(m1,1),(m2,-1)]
    raise RuntimeError('No matrix found!')

def smallCMcycle(D):
    A = compute_gamma_tau(D).change_ring(F)
    t0, t1 = fixed_point(A, phi)
    return A, [t0, t1]

def smallRMcycle(D):
    A = compute_gamma_tau(D).change_ring(F)
    t0, t1 = fixed_point(A, phi)
    return A, [t0, t0] # not a typo!

def RMCEval(D, cycle_type = 'smallCM', prec=M):
    global L0, J
    if cycle_type == 'smallCM':
        A, tau0 = smallCMcycle(D)
    else:
        if cycle_type != 'smallRM':
            raise NotImplementedError('Cycle type should be either "smallCM" or "smallRM"')
        A, tau0 = smallRMcycle(D)
    mlist0 = matrices_for_unimodular_path(A[0,0], A[1,0])
    mlist = sum((good_matrices(m) for m in mlist0),[])
    res0 = prod(Eval0(L0, act_matrix(m.apply_morphism(phi).adjugate(), tau0))**sgn for m,sgn in mlist)
    if parallelize:
        res1 = 1
        with futures.ProcessPoolExecutor() as executor:
            future_dict = {executor.submit(Eval, act_matrix(m.apply_morphism(phi).adjugate(), tau0), prec) : sgn for m, sgn in mlist}
            for fut in futures.as_completed(future_dict):
                res1 *= fut.result()**(future_dict[fut])
    else:
        res1 = prod(Eval(act_matrix(m.apply_morphism(phi).adjugate(), tau0))**sgn for m,sgn in mlist)
    return (res0 * res1).add_bigoh(prec)

def list_powers(x, M,j):
    if x is None:
        return j, x
    plist = [x.parent()(1)]
    for i in range(M):
        plist.append(x*plist[-1])
    return j, plist

def calculate_Tp_matrices(P):
    global input_list
    Tplist = [matrix(2,2,[P, a, 0, 1]) for a in range(P.norm()) ] + [matrix(2,2,[1,0,0,P])]
    MS = []
    for m in Tplist:
        assert m[1,0] == 0
        mlist0 = matrices_for_unimodular_path(-m[0,1], m[0,0])
        mlist = sum((good_matrices(o) for o in mlist0),[])
        if not all((m * o).apply_map(lambda x : x.trace() - x).apply_morphism(phi).determinant().valuation() == 0 for o,sgn in mlist):
            raise RuntimeError('Problem with embedding, try the other one.')

        MS.extend((m * o, sgn) for o, sgn in mlist)
    apply_single_dict = {}
    S = Ruv
    R = S.base_ring()
    z1 = R.gen()
    z2 = S.gen()
    future = {}
    if False:
        with futures.ProcessPoolExecutor() as executor:
            for m, sgn in MS:
                m.set_immutable()
                mconj = m.apply_map(lambda x : x.trace() - x)
                A = m.apply_morphism(phi)
                Aconj = mconj.apply_morphism(phi)
                for i in range(p+1):
                    ii, subst1 = ApplySingle(A, i, z1, check=False)
                    jj, subst2 = ApplySingle(Aconj, i, z2, check=True)
                    future[executor.submit(list_powers, subst1, M, ii)] = (m,i,0)
                    future[executor.submit(list_powers, subst2, M, jj)] = (m,i,1)
            print(f'All {len(future)} jobs submitted, now computing in parallel')
            for fut in futures.as_completed(future):
                apply_single_dict[future[fut]] = fut.result()
            print('Done')
    else:
        cnt = 0
        for cnt, (m, sgn) in enumerate(MS):
            update_progress(float(cnt)/len(MS))
            m.set_immutable()
            mconj = m.apply_map(lambda x : x.trace() - x)
            A = m.apply_morphism(phi)
            Aconj = mconj.apply_morphism(phi)
            for i in range(p+1):
                ii, subst1 = ApplySingle(A, i, z1, check=False)
                jj, subst2 = ApplySingle(Aconj, i, z2, check=True)
                apply_single_dict[m,i,0] = list_powers(subst1,M,ii)
                apply_single_dict[m,i,1] = list_powers(subst2,M,jj)
        print('Done')

    input_list = {(i,j) : list() for i in range(p+1) for j in range(p+1)}
    for m, sgn in MS:
        for inky0 in range(p+1):
            outky0, s1 = apply_single_dict[(m, inky0, 0)]
            if s1 is None:
                continue
            for inky1 in range(p+1):
                inky = (inky0, inky1)
                outky1, s2 = apply_single_dict[(m, inky1, 1)]
                if s2 is None:
                    continue
                outky = (outky0, outky1)
                input_list[outky].append((inky, s1, s2, sgn))
    return MS


# given a non-necessary fundamental discriminant D, computes the matrix gamma_tau associated to an optimal embedding of the ring of discriminant D to M_0(2)
def compute_gamma_tau(D):
    Dsqf = D.squarefree_part()
    try:
        c = ZZ((D / fundamental_discriminant(D)).sqrt())
    except TypeError:
        raise ValueError('D is not the discrimininant of a quadratic order')
    F.<r> = QuadraticField(Dsqf)
    w = F.maximal_order().ring_generators()[0]
    # O_c has Z-basis <1, cw>, we use this basis to embed O_c into M_2(Z)
    coords = (c*w).coordinates_in_terms_of_powers()
    D0 = Matrix([[0,1], coords((c*w)^2)]).transpose() # the image of cw in M_2(Z)
    aa, bb, cc, dd = D0.list()
    # if D0 is not in M_0(2), we conjugate the embedding to force this condition
    A = Matrix(2,2,[0,1,-1,0])
    B = Matrix(2,2,[2,1,1,1])
    if cc % 2 == 0:
        D1 = D0
    elif bb % 2 == 0:
        D1 = A * D0 * A.inverse()
    elif (aa - dd) % 2 == 0:
        D1 = B * D0 * B.inverse()
    else:
        raise RuntimeError('It seems that there is no optimal embedding of O_c into Gamma_0(2)')
    assert D1.minpoly() == (c*w).minpoly()
    assert D1[1][0] % 2 == 0
    # now we find the fundamental unit of O_c
    eps = F.units()[0]
    n = 1
    while True:
        cs = coords(eps**n)
        if all([o.is_integer() for o in cs]):
            gamma_tau = cs[0] + cs[1] * D1
            assert gamma_tau.minpoly() == (eps^n).minpoly()
            assert gamma_tau[1][0] % 2 == 0
            return gamma_tau
        n +=1

def compute_gammas(max_D=1000):
    gamma_list = []
    for D in srange(2, max_D):
        if D.is_square():
            continue
        try:
            gamma_tau = compute_gamma_tau(D)
        except (TypeError,RuntimeError,NotImplementedError):
            continue
        gamma_list.append((D,gamma_tau))
    return gamma_list

def vanishing_functional(N = 10):
    g = ModularForms(20).cuspidal_subspace().gens()[0].q_expansion(N)
    E = lambda n : sum([ o for o in ZZ(n).divisors() if o % 4 != 0])
    A = Matrix([[ZZ(g[o]) for o in range(1,N)], [E(o) for o in range(1,N)]])
    # vectors in the kernel correspond to functionals that vanish on g and on the Eisenstein series
    # the first position of the vector in the kernel corresponds to a_1, and so on
    L = IntegralLattice(N-1).sublattice(A.right_kernel().basis())
    length_cap = 8
    all_vectors = sum(L.short_vectors(length_cap),[])
    ans = []
    W = L.submodule(ans)
    i = 0
    while W.rank() < L.rank():
        print(W.rank(),L.rank())
        try:
            while L.submodule(ans + [all_vectors[i]]).rank() == W.rank():
                i += 1
        except IndexError:
            length_cap *= 2
            all_vectors = sum(L.short_vectors(length_cap),[])
        v = all_vectors[i]
        print(f'Adding {v = }')
        ans.append(v)
        W = L.submodule(ans)
    return ans

def label(func):
    pos = ''
    neg = ''
    for i, o in enumerate(func):
        if o > 0:
            pos = pos + o * str(i+1)
        else:
            neg = neg + (-o) * str(i+1)
    return pos + '_' + neg

def initial_seed(v, p):
    L0 = [[], []]
    V = [[], []]
    for (im1, vi) in enumerate(v):
        i = im1 + 1
        if vi > 0:
            Li = vectors_in_lattice(i)[:2]
            Vi = vectors_in_lattice(p * i)[:2]
            L0[0].extend(vi * Li[0])
            L0[1].extend(vi * Li[1])
            V[0].extend(vi * Vi[0])
            V[1].extend(vi * Vi[1])
        elif vi < 0:
            Li = vectors_in_lattice(i)[:2]
            Vi = vectors_in_lattice(p * i)[:2]
            L0[1].extend((-vi) * Li[0])
            L0[0].extend((-vi) * Li[1])
            V[1].extend((-vi) * Vi[0])
            V[0].extend((-vi) * Vi[1])
    return L0, V


def find_value_one(maxD, cycle_type='smallCM'):
    ones = []
    not_ones = []
    errors = []
    for D in srange(2, maxD):
        print(D)
        if D.is_square():
            continue
        try:
            c2 = ZZ(D / fundamental_discriminant(D))
            c2 = c2.sqrt()
        except TypeError:
            continue
        try:
            x1 = RMCEval(D,cycle_type=cycle_type,prec=10)
            if x1 == 1:
                print(f'{D = } yields one')
                ones.append(D)
            else:
                print(f'{D = } yields {x1}')
                not_ones.append(D)
        except (TypeError, ValueError, RuntimeError, PrecisionError) as e:
            errors.append((D,str(e)))
    return ones, not_ones, errors
