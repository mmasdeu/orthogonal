# from tqdm import tqdm
from sage.rings.padics.precision_error import PrecisionError
load('orthogonal_manin_trick.sage')
from multiprocessing import Process, Manager, Pool
from concurrent import futures

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
    if sgn == 1:
        return j1, j2, ans, sgn
    else:
        return j1, j2, inv(ans), sgn

def Level1(V):
    res = {(i,j) : Ruv(1) for i in range(p+1) for j in range(p+1)}
    dg = {(i,j) : 0 for i in range(p+1) for j in range(p+1)}
    pM = p**M
    input_vec = []
    for (A, B) in zip(*V):
        Ap = A.apply_map(lambda x : Rp(phi(x).lift()))
        Bp = B.apply_map(lambda x : Rp(phi(x).lift()))
        input_vec.extend([(A,Ap,1,Rp), (B,Bp,-1,Rp)])
    if False: # parallelize:
        I = compute_level1_contribution(input_vec)
        for _, (i, j, FF, sgn) in I:
            res[i,j] *= FF
            dg[i,j] += sgn
    else:
        for A, Ap, sgn, _ in input_vec:
            i, j, FF, sgn = compute_level1_contribution(A, Ap, sgn, Rp)
            res[i,j] *= FF
            dg[i,j] += sgn
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

def NewTransform(Finpsinpky):
    F, Finv, inps, inpky = Finpsinpky
    res = 1
    for ky, s1, s2, sgn in inps:
        f = F[ky] if sgn == 1 else Finv[ky]
        S = f.parent()
        R = S.base_ring()
        if s1 is None:
            continue
        res *= sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] for fi, i in zip(f.coefficients(), f.exponents()))
    return inpky, res

def Next(F):
    res = {}
    if parallelize:
        with futures.ProcessPoolExecutor() as executor:
            # Calculate inverses
            print('Calculating inverses...')
            it = [(ky, val) for ky, val in F.items()]
            future = {executor.submit(inv_par, (ky,val)) : ky for ky, val in it}
            Finv = {}
            print('Looking at as_completed')
            for fut in futures.as_completed(future):
                ky, val = fut.result()
                print('ky = ',ky)
                Finv[ky] = val
                
            # Finv = {ky : val for ky, val in executor.map(inv_par, it,chunksize=len(it))}

            # Iteration
            print('Main iteration')
            it = [(deepcopy(F), deepcopy(Finv), deepcopy(inps), ky) for inps, ky in input_list]
            res = {}
            future = {executor.submit(NewTransform, i) : i[-1] for i in it}
            print('All submitted now')
            for fut in futures.as_completed(future):
                ky, val = fut.result()
                print('Computed ',ky, future[fut])
                res[ky] = val
            # for ky, val in executor.map(NewTransform, it, chunksize=len(it)):
            #     print('Computed ', ky)
            #     res[ky] = val
            print('Done')
            # res = {ky : val for ky, val in executor.map(NewTransform, it)}
    else:
        Finv = {ky : inv(val) for ky, val in F.items()}
        res = {ky : NewTransform((F, Finv, inps, None))[1] for inps, ky in input_list}
    return res

def RMC(F):
    FF = F
    res = F
    for j in range(1,M):
        d2 = max(val.degree() for ky, val in F.items())
        d1 = max(o.degree() for ky, val in F.items() for o in val.coefficients())
        print(f'Iteration {j}')
        t = walltime()
        FF = Next(FF)
        FF = Next(FF)
        for ky in res:
            res[ky] *= FF[ky]
        print(f'..done in {walltime(t)} seconds.')
    return res

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
    fphi = f.map_coefficients(phi)
    return solve_quadratic(f.change_ring(phi).change_ring(L), L, return_all=True)

def Eval0(L0, tau):
    t0, t1 = tau
    vv = vector([-1, t0])
    ww = vector([t1, 1])
    KK = vv.parent().base_ring()
    ans = prod((vv * num.apply_morphism(phi).apply_map(lambda o : KK(o)) * ww) / (vv * den.apply_morphism(phi).apply_map(lambda o : KK(o)) * ww) for num, den in zip(*L0))
    return ans

def Eval(J, tau):
    t0, t1 = tau
    ans = 1
    for ky, val in J.items():
        x0 = t0 if ky[0] == p else 1/(t0 - ky[0])
        x1 = t1 if ky[1] == p else 1/(t1 - ky[1])
        ans *= sum(sum(o.lift() * x0**j for o, j in zip(a.coefficients(), a.exponents())) * x1**i for a, i in zip(J[ky].coefficients(), J[ky].exponents()))
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

def RMCEval(L0, J, A):
    A = A.change_ring(F)
    tau0 = fixed_point(A, phi)
    mlist0 = matrices_for_unimodular_path(A[0,0], A[1,0])
    mlist = sum((good_matrices(m) for m in mlist0),[])
    res0 = prod(Eval0(L0, act_matrix(m.apply_morphism(phi).adjugate(), tau0))**sgn for m,sgn in mlist)
    res1 = prod(Eval(J, act_matrix(m.apply_morphism(phi).adjugate(), tau0))**sgn for m,sgn in mlist)
    res = res0 * res1
    return res, res0, res1

def calculate_Tp_matrices(P):
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
    for i in range(p+1):
        for m, sgn in MS:
            m.set_immutable()
            mconj = m.apply_map(lambda x : x.trace() - x)
            A = m.apply_morphism(phi)
            Aconj = mconj.apply_morphism(phi)
            ii, subst1 = ApplySingle(A, i, z1, check=False)
            try:
                subst1pow = [subst1.parent()(1)]
                for k in range(M):
                    subst1pow.append(subst1pow[-1] * subst1)
            except AttributeError:
                subst1pow = None
            apply_single_dict[(m, i, 0)] = (ii, subst1pow)
            jj, subst2 = ApplySingle(Aconj, i, z2, check=True)
            try:
                subst2pow = [subst2.parent()(1)]
                for k in range(M):
                    subst2pow.append(subst2pow[-1] * subst2)
            except AttributeError:
                subst2pow = None
            apply_single_dict[(m, i, 1)] = (jj, subst2pow)
    aux_dict = {(i,j) : list() for i in range(p+1) for j in range(p+1)}
    for m, sgn in MS:
        for ky0 in range(p+1):
            for ky1 in range(p+1):
                ii, s1pow = apply_single_dict[(m, ky0, 0)]
                jj, s2pow = apply_single_dict[(m, ky1, 1)]
                # if s1pow is None:
                #     s12pow = None
                # else:
                #     s12pow = {(r, s) : s1pow[r] * s2pow[s] for r in range(M+1) for s in range(M+1)}
                aux_dict[(ii,jj)].append(((ky0,ky1), s1pow, s2pow, sgn))
    return MS, [(inps, ky) for ky, inps in aux_dict.items()]

# given a non-necessary fundamental discriminant D, computes the matrix gamma_tau associated to an optimal embedding of the ring of discriminant D to M_0(2)
def compute_gamma_tau(D):
    Dsqf = D.squarefree_part()
    c = ZZ((D / fundamental_discriminant(D)).sqrt())
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


