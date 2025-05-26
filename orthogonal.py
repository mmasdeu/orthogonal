# from tqdm import tqdm
from sage.rings.padics.precision_error import PrecisionError
from multiprocessing import cpu_count
from concurrent import futures
from sage.misc.timing import cputime
from sage.all import *
from util import *

x = QQ['x'].gen()
valid_ATR_fields = dict([(41, x**4 + 4*x**3 + x**2 - 6*x - 8), (113, x**4 - 28*x**2 - 256), (137, x**4 + 4*x**3 + 105*x**2 + 202*x - 224), (313, x**4 + 4*x**3 - 111*x**2 - 230*x - 3032), (337, x**4 + 4*x**3 - 3*x**2 - 14*x - 72), (409, x**4 - 12*x**2 - 1600), (457, x**4 + 4*x**3 - 15*x**2 - 38*x - 24), (521, x**4 - 44*x**2 - 1600), (569, x**4 + 4*x**3 - 111*x**2 - 230*x - 8216), (593, x**4 - 92*x**2 - 256), (809, x**4 + 4*x**3 + x**2 - 6*x - 200), (857, x**4 + 4*x**3 - 255*x**2 - 518*x - 584), (881, x**4 + 4*x**3 - 219*x**2 - 446*x - 5408), (1153, x**4 + 4*x**3 - 27*x**2 - 62*x - 48), (1201, x**4 + 4*x**3 - 219*x**2 - 446*x - 11888), (1217, x**4 - 124*x**2 - 1024), (1249, x**4 - 60*x**2 - 4096), (1321, x**4 + 4*x**3 - 39*x**2 - 86*x - 26288), (1553, x**4 - 92*x**2 - 4096), (1657, x**4 - 68*x**3 - 1276*x**2 + 135712*x + 1478656), (1777, x**4 - 156*x**2 - 1024), (1889, x**4 + 4*x**3 - 11*x**2 - 30*x - 416), (1993, x**4 - 172*x**2 - 576), (2113, x**4 + 4*x**3 - 291*x**2 - 590*x - 21032), (2129, x**4 - 92*x**2 - 6400), (2273, x**4 + 4*x**3 + 53*x**2 + 98*x + 32), (2377, x**4 + 4*x**3 - 183*x**2 - 374*x - 39392), (2521, x**4 - 140*x**2 - 5184), (2593, x**4 + 4*x**3 - 147*x**2 - 302*x - 46808), (2617, x**4 + 4*x**3 + 465*x**2 + 922*x + 136), (2657, x**4 + 4*x**3 - 435*x**2 - 878*x - 5624), (2689, x**4 + 4*x**3 - 27*x**2 - 62*x - 432), (2729, x**4 + 4*x**3 - 39*x**2 - 86*x - 54800), (2833, x**4 - 92*x**2 - 9216), (2953, x**4 + 4*x**3 - 47*x**2 - 102*x - 88), (3001, x**4 - 204*x**2 - 1600), (3089, x**4 + 4*x**3 + 501*x**2 + 994*x - 800), (3217, x**4 + 4*x**3 - 3*x**2 - 14*x - 792), (3361, x**4 - 60*x**2 - 12544), (3433, x**4 - 100*x**3 - 2684*x**2 + 369056*x + 6718464), (3593, x**4 + 4*x**3 - 47*x**2 - 102*x - 248), (3761, x**4 + 4*x**3 - 219*x**2 - 446*x - 63728), (3881, x**4 - 236*x**2 - 1600), (3929, x**4 - 140*x**2 - 10816), (4049, x**4 - 220*x**2 - 4096), (4073, x**4 + 4*x**3 - 327*x**2 - 662*x - 55088), (4177, x**4 + 4*x**3 - 75*x**2 - 158*x - 83024), (4273, x**4 + 4*x**3 - 507*x**2 - 1022*x - 21248), (4289, x**4 + 4*x**3 - 59*x**2 - 126*x - 80), (4513, x**4 - 188*x**2 - 9216), (4657, x**4 - 148*x**3 - 476*x**2 + 589472*x + 8856576), (4721, x**4 + 4*x**3 - 219*x**2 - 446*x - 83168), (4793, x**4 + 4*x**3 - 111*x**2 - 230*x - 93752), (4801, x**4 + 4*x**3 - 59*x**2 - 126*x - 208), (4817, x**4 + 4*x**3 - 363*x**2 - 734*x - 63872), (4969, x**4 + 4*x**3 - 327*x**2 - 662*x - 73232), (4993, x**4 - 244*x**3 + 13348*x**2 + 347168*x + 589824), (5233, x**4 - 28*x**2 - 20736), (5393, x**4 + 4*x**3 - 651*x**2 - 1310*x - 1952), (5641, x**4 + 4*x**3 + 81*x**2 + 154*x + 72), (5657, x**4 + 4*x**3 - 543*x**2 - 1094*x - 39752), (5801, x**4 + 4*x**3 + x**2 - 6*x - 1448), (5897, x**4 - 44*x**2 - 23104), (6073, x**4 + 4*x**3 - 71*x**2 - 150*x - 112), (6217, x**4 + 4*x**3 - 183*x**2 - 374*x - 117152), (6329, x**4 + 4*x**3 - 71*x**2 - 150*x - 176), (6353, x**4 + 4*x**3 - 651*x**2 - 1310*x - 21392), (6449, x**4 - 28*x**2 - 25600), (6473, x**4 - 164*x**3 - 2172*x**2 + 936608*x + 19784704), (6529, x**4 + 4*x**3 - 579*x**2 - 1166*x - 47240), (6689, x**4 + 4*x**3 - 11*x**2 - 30*x - 1616), (7001, x**4 - 140*x**2 - 23104), (7121, x**4 - 220*x**2 - 16384), (7193, x**4 - 260*x**3 + 12036*x**2 + 862496*x + 5914624), (7321, x**4 + 4*x**3 - 543*x**2 - 1094*x - 73448), (7369, x**4 + 4*x**3 - 759*x**2 - 1526*x - 3680), (7393, x**4 - 188*x**2 - 20736), (7417, x**4 - 76*x**2 - 28224), (7489, x**4 + 4*x**3 - 291*x**2 - 590*x - 129896), (7793, x**4 - 28*x**2 - 30976), (7841, x**4 - 316*x**2 - 6400), (8009, x**4 + 4*x**3 - 79*x**2 - 166*x - 280), (8089, x**4 - 268*x**2 - 14400), (8161, x**4 + 4*x**3 - 75*x**2 - 158*x - 480), (8209, x**4 - 220*x**2 - 20736), (8273, x**4 - 84*x**3 - 13532*x**2 + 907168*x + 58491904), (8297, x**4 + 4*x**3 + 825*x**2 + 1642*x + 496), (8329, x**4 - 300*x**2 - 10816), (8369, x**4 + 4*x**3 - 19*x**2 - 46*x - 1960), (8377, x**4 - 204*x**2 - 23104), (8521, x**4 + 4*x**3 - 759*x**2 - 1526*x - 27008), (8609, x**4 - 188*x**2 - 25600), (8681, x**4 - 364*x**2 - 1600), (9137, x**4 - 284*x**2 - 16384), (9257, x**4 - 228*x**3 + 1924*x**2 + 1558432*x + 30647296), (9377, x**4 - 316*x**2 - 12544), (9433, x**4 + 4*x**3 - 831*x**2 - 1670*x - 16712), (9473, x**4 + 4*x**3 - 867*x**2 - 1742*x - 2168), (9497, x**4 + 4*x**3 - 55*x**2 - 118*x - 1504), (9601, x**4 - 380*x**2 - 2304), (9689, x**4 - 132*x**3 - 12284*x**2 + 1408288*x + 69222400), (9697, x**4 + 4*x**3 - 75*x**2 - 158*x - 864), (9817, x**4 + 4*x**3 + 105*x**2 + 202*x + 96), (9929, x**4 + 4*x**3 - 759*x**2 - 1526*x - 55520), (10009, x**4 - 12*x**2 - 40000), (10169, x**4 + 4*x**3 - 111*x**2 - 230*x - 202616), (10177, x**4 - 124*x**2 - 36864), (10337, x**4 - 316*x**2 - 16384), (10369, x**4 - 252*x**2 - 25600), (10433, x**4 + 4*x**3 - 91*x**2 - 190*x - 352), (10513, x**4 + 4*x**3 - 651*x**2 - 1310*x - 105632)])

def get_F_and_alpha(D):
    f = valid_ATR_fields[D]
    F = NumberField(x*x - D, names ='r')
    r = F.gen()
    E = NumberField(f, names = 'gE')
    gE = E.gen()
    EF = E.relativize(F.embeddings(E)[0], 'gEF')
    r0 = f.roots(EF)[0][0]
    alpha = r0.minpoly().discriminant()
    assert alpha.relative_norm().squarefree_part() == -1
    return F, alpha

# load('find_E_and_L.sage')

def enumerate_matrices(r, bound=10):
    F = r.parent()
    for a in F.elements_of_bounded_height(bound=bound):
        for c in F.elements_of_bounded_height(bound=bound):
            gg, d, mb = xgcd_quadratic(a, c)
            if gg.norm() != 1:
                continue
            M = Matrix(2,2,[a,-mb / gg,c,d / gg])
            assert M.determinant() == 1
            M.set_immutable()
            yield M

def find_all_cosets(r, level=2):
    F = r.parent()
    index = F.ideal(level).norm()
    for P, _ in F.ideal(level).factor():
        index *= (1 + 1/P.norm())
    id = Matrix(F,2,2,1)
    id.set_immutable()
    ans = [id]
    I = enumerate_matrices(r)
    while len(ans) < index:
        m = next(I)
        if find_coset_rep(m, ans, level=level) is None:
            ans.append(m)
    return ans

def find_coset_rep(m, coset_list, level=2):
    global F
    minv = m.adjugate()
    I = F.ideal(level)
    for c in coset_list:
        if (c * minv)[1,0] in I:
            return c
    return None


class Cocycle(SageObject):
    def __init__(self, F, p, prec, level=2, initialize=True):
        self.F = F
        self.p = p
        self.level = level
        self.prec = prec
        r = F.gen()
        self.coset_list = find_all_cosets(r, level)
        self.L0 = None
        self.data = None
        if initialize:
            self._initialize_data()

    # The same functionality as initial_seed, but takes the list
    # of coset reps for Gamma \ SL2(Z[i]) and outputs a dictionary
    def _initialize_data(self, v):
        if isinstance(v, str):
            v = functional_from_label(v)
        v = tuple(v)
        reps = self.coset_list
        L0 = {{g : []} for g in reps}
        L1 = {{g : []} for g in reps}
        for i, vi in enumerate(v, start=1):
            if vi == 0:
                continue
            for g in reps:
                Li = vectors_in_lattice(self.F, self.p, i, g)
                Lpi = vectors_in_lattice(self.F, self.p, self.p * i, g)
                L0[g].extend([(o, vi * sgn) for o, sgn in Li])
                L1[g].extend([(o, vi * sgn) for o, sgn in Lpi])
        assert all(sum(e for o, e in L0[g]) == 0 for g in reps)
        assert all(sum(e for o, e in L1[g]) == 0 for g in reps)
        self.L0 = L0
        self.data = {}        
        for g in self.coset_list:
            self.data[g] = Level1(self.p, L1[g])
        return

    def __getitem__(self, ky):
        r, s = ky
        # Write r->s as sum g0 -> goo using Manin trick
        mlist = manin_trick(r,s)
        return prod(self._eval_translate_zero_inf(m)**sgn for m, sgn in mlist)

    # Evaluate the cocycle at m0 -> moo
    def _eval_translate_zero_inf(self, g):
        global Ruv, phi
        gi = find_coset_rep(g, self.coset_list, level=self.level)
        m = g * gi.adjugate()
        if m == 1:
            return self.data[gi]

        S = Ruv
        R = S.base_ring()
        z1 = R.gen()
        z2 = S.gen()
        m.set_immutable()

        # Evaluate m * data[gi]
        mconj = m.apply_map(lambda x : x.trace() - x)
        A = m.apply_morphism(phi)
        Aconj = mconj.apply_morphism(phi)
        tmp_dict = {}
        for i in range(self.p+1):
            j1, subst1 = ApplySingle(A, i, z1, self.prec, check=False)
            j2, subst2 = ApplySingle(Aconj, i, z2, self.prec, check=True)
            tmp_dict[i,0] = (j1, list_powers(subst1,self.prec))
            tmp_dict[i,1] = (j2, list_powers(subst2,self.prec))

        transformation_list = {}
        for inky0 in range(self.p+1):
            outky0, s1 = tmp_dict[inky0, 0]
            if s1 is None:
                continue
            for inky1 in range(self.p+1):
                inky = (inky0, inky1)
                outky1, s2 = tmp_dict[inky1, 1]
                if s2 is None:
                    continue
                outky = (outky0, outky1)
                transformation_list[outky].append((inky, s1, s2, sgn))
        gF = self.data[gi]
        ans = {}
        for outky in transformation_list:
            res = Ruv(1)
            resinv = Ruv(1)
            for inky, s1, s2, sgn in transformation_list[outky]:
                f = gF[inky]
                newres = sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] \
                    for fi, i in zip(f.coefficients(), f.exponents()))
                if sgn == 1:
                    res *= newres
                else:
                    resinv *= newres
            ans[outky] = res * inv(resinv)
        return ans
    
    def next(self, **kwargs):
        timing = kwargs.get('timing', False)
        coset_list = kwargs.get('coset_list', None)
        if coset_list is None:
            coset_list = self.coset_list
        p = kwargs.get('p', None)
        if p is None:
            p = self.p
        parallelize = kwargs.get('parallelize', True)
        ncpus = kwargs.get('ncpus', cpu_count())
        res = {(g,i,j) : 1 for g in coset_list for i in range(p+1) for j in range(p+1)}
        times = vector([0, 0, 0])
        if parallelize:
            with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
                # Iteration
                future_dict = {executor.submit(Transform, outky) : outky for outky in input_list}
                for fut in futures.as_completed(future_dict):
                    outky = future_dict[fut]
                    ans, (t1p, t2p, t3p) = fut.result()
                    times += vector([t1p, t2p, t3p])
                    res[outky] = ans
        else:
            for outky in input_list:
                res[outky], (t1p, t2p, t3p) = Transform(outky)
                times += vector([t1p, t2p, t3p])
        if timing:
            return res, tuple(times)
        else:
            return res

    def RMC(self, **kwargs):
        global ncpus
        coset_list = kwargs.get('coset_list', self.coset_list)
        p = kwargs.get('p', self.p)        
        try:
            ncpus = ncpus
        except NameError:
            ncpus = cpu_count()
        FF = F
        res = [F]
        d = (-1,-1)
        j = 0
        while d != (0,0):
            j += 1
            print(f'Iteration {j}')
            t = walltime()
            FF0, t0 = self.next(timing=True, **kwargs)
            FF, t1 = self.next(timing=True, **kwargs)
            t1, t2, t3 = vector(t0) + vector(t1)
            res.append(FF)
            d = degrees(FF)
            print(f'..done in {walltime(t)} seconds. {t1 = :.2f}, {t2 = :.2f}, {t3 = :.2f}. Degrees: {d}')
        print(f'Now computing product...')
        t = walltime()
        if parallelize:
            with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
                while len(res) > 1:
                    future_dict = {executor.submit(lambda x, y: {ky : x[ky] * y[ky] for ky in x}, res[i], res[i+1]) : i for i in range(0,len(res)-1, 2) }
                    res = [] if len(res) % 2 == 0 else [res[-1]]
                    for fut in futures.as_completed(future_dict):
                        res.append(fut.result())
            ans = res[0]
        else:
            ans = {ky : prod(r[ky] for r in res) for ky in res[0] }
        print(f'..done in {walltime(t)} seconds.')
        R0 = Ruv.base_ring()
        psi = R0.hom([ZZ['u'].gen()],base_map=lambda x:x.lift(),check=False)
        ans = {ky : f.change_ring(psi) for ky, f in ans.items()}
        return ans

def get_predicted_field_and_prime_list(F, D, n, typ, char, names='z', prime_bound=600):
    r'''
    If smallRM:
    - If char == triv: return F
    - If char == conj: return HCF of Q(sqrt(-D))

    If smallCM:
    - If char == triv: return HCF of Q(sqrt(-D))
    - If char == conj: return Q
    '''
    if typ in ['smallRM','smallCM'] and n != 1:
        return None, None
    if char not in ['triv', 'conj']:
        raise ValueError('Parameter "char" must be either "triv" or "conj"')
    x = QQ['x'].gen()
    if typ == 'smallCM':
        D = -D
    if typ == 'bigRM':
        # alpha = ATR_alpha(D,n)
        alpha = get_F_and_alpha(D)[1]
        L0 = alpha.parent()
        L1 = tau_ATR_field(alpha)
        L = L1.absolute_field(names='t')
    else:
        L = QuadraticField(D, names='a')
        M = (L.composite_fields(F, names='t0')[0]).absolute_field(names='t1')
    if typ == 'smallCM' and char == 'conj':
        H = QQ
    else:
        try:
            H = (L.hilbert_class_field(names='tt').composite_fields(F)[0]).absolute_field(names=names)
        except AttributeError:
            H = NumberField(magma(L).HilbertClassField().AbsoluteField().DefiningPolynomial()._sage_(),names='tt').composite_fields(F)[0].absolute_field(names=names)
    prime_list = [p for p in prime_range(prime_bound) if len(L.ideal(p).factor()) < L.degree()]
    if typ != 'bigRM':
        assert all(len(M.ideal(p).factor()) < M.degree() for p in prime_list)
    return H, prime_list

# Related to Manin Trick
def quo_rem_quadratic(a,b):
    if b == 0:
        raise ZeroDivisionError
    K = a.parent()
    q0 = ((a/b)[0]).floor() + ((a/b)[1]).floor()*K.gen()
    for e0, e1 in product([0,1], repeat=2):
        q = q0 + e0 + e1*K.gen()
        r = a - q*b
        if r.norm() < b.norm():
            return (q, r)
    raise RuntimeError

def xgcd_quadratic(a, b):
    v1 = vector([1,0])
    v2 = vector([0,1])
    if b == 0:
        return (a, v1[0], v1[1])
    if a == 0:
        return (b, v2[0], v2[1])
    a0 = a
    b0 = b
    q, r = quo_rem_quadratic(a0,b0)
    while True:
        d = a0 - q*b0
        if d == 0:
            assert b0 == v2[0]*a + v2[1]*b
            return b0, v2[0], v2[1]
        v3 = v1 - q* v2
        v1 = vector([v2[0], v2[1]])
        v2 = vector([v3[0], v3[1]])
        a0 = b0
        b0 = d
        q, r = quo_rem_quadratic(a0,b0)
    return

# computes the continued fraction of a/c, where a and c belong to Z[i]
def continued_fraction(a, c):
    a_over_c = a/c
    cf = []
    while c != 0:
        q, r = quo_rem_quadratic(a,c)
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
    assert all([p[n]*q[n-1] - p[n-1]*q[n] == (-1)**(n-1) for n in range(1,N)])
    #assert p[-1]/q[-1] == eval_cf(cf)
    return p, q

# given the convergents [p_i] and [q_j] of a continued fraction for alpha, computes the matrices that Mi such that {Infinity, alpha} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def compute_Ms(p, q):
    Ms = []
    Ms.append(Matrix(2,2,[-p[0], 1, -q[0], 0])) # j = 0
    for j in range(1, len(p)):
        Ms.append(Matrix(2,2,[(-1)**(j-1)*p[j],p[j-1],(-1)**(j-1)*q[j],q[j-1]]))
    assert all([A.det() == 1 for A in Ms])
    return Ms

def manin_trick(x, y):
    ylist = matrices_for_unimodular_path(y.numerator(), y.denominator())
    xlist = matrices_for_unimodular_path(x.numerator(), x.denominator())
    return [(o, 1) for o in ylist] + [(o, -1) for o in xlist]

# given a, c in Z[i] returns the matrices that Mi such that {Infinity, a/c} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def matrices_for_unimodular_path(a,c):
    return compute_Ms(*compute_convergents(continued_fraction(a,c)))

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
        eps = F.unit_group().gens()[0]
        units = [F(eps**i) for i in range(eps.multiplicative_order())]
        return [u * beta for beta in F.elements_of_norm(n) for u in units]

# Given a prime p, returns two lists corresponding to
# p-primitive elements of norm n
# This corresponds to the matrices determining a region
# M which intersects
# the path g·0 -> g·∞
def vectors_in_lattice(F, p, n, g = 1):
    res = []
    for mac in range(1,n+1):
        RHS = n - mac
        betas = all_elements_of_norm(F, RHS)
        for mc in divisors(mac):
            a = mac // mc
            condition = a % p == 0 and mc % p == 0
            for beta in betas:
                if condition and (beta / p).is_integral():
                    continue
                else:
                    new_mat = Matrix([[beta.conjugate(), mc], [a, -beta]])
                    if a % 4 == 1:
                        res.append((g * new_mat, 1))
                    elif a % 4 == 3:
                        res.append((- g * new_mat, -1))
    return res

def inv(FF):
    # global Rp
    a0inv = ~(FF(0)(0))
    y1 = 0
    y = a0inv
    while y != y1:
        y1 = y
        y = y1 * (2 - y1 * FF)
    return y


@parallel(ncpus=cpu_count())
def compute_level1_contribution(F, p, A, exponent, phi, map_poly, Rp):
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
        j1, j2 = to_x(polp).roots()[0][0].lift(), p
        tau1, tau2 = j1 + 1 / t1, t2
    if d1 == 0 and d2 == 1:
        j1, j2 = p, to_x(polp).roots()[0][0].lift()
        tau1, tau2 = t1, j2 + 1 / t2
    if d1 == 1 and d2 == 1:
        ff = polp.factor()
        assert len(ff) == 2,'Does not factor'
        f1 = ff[1][0]
        f2 = ff[0][0]
        assert f1.degrees() == (1,0) and f2.degrees() == (0,1)
        # print(f1.change_ring(to_x))
        j1 = to_x(f1).roots()[0][0].lift()
        j2 = to_x(f2).roots()[0][0].lift()
        tau1, tau2 = j1 + 1 / t1, j2 + 1 / t2
    ans = (pol.parent()(pol(tau1,tau2) * t1**d1 * t2**d2)).change_ring(phi)
    while j1 < 0:
        j1 += p
    while j2 < 0:
        j2 += p
    ans = ans.change_ring(Rp)
    ans = map_poly(ans)
    return j1, j2, ans, exponent

def Level1(p, V, phi, map_poly, Rp):
    global ncpus
    try:
        ncpus = ncpus
    except NameError:
        ncpus = cpu_count()
    res = {(i,j) : [1, 1] for i in range(p+1) for j in range(p+1)}
    dg = {(i,j) : 0 for i in range(p+1) for j in range(p+1)}
    input_vec = []
    for A, exponent in V:
        Ap = A.apply_map(lambda x : Rp(phi(x).lift()))
        input_vec.append((A, Ap,exponent,Rp))
    for A, Ap, sgn, _ in input_vec:
        F = A.parent().base_ring()
        i, j, FF, exponent = compute_level1_contribution(F, p, A, exponent, phi, map_poly, Rp)
        if exponent > 0:
            res[i,j][0] *= FF**exponent
        else:
            res[i,j][1] *= FF**(-exponent)
        dg[i,j] += sgn
    with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
        # Calculate inverses
        print('Calculating inverses...')
        future = {executor.submit(inv, val[1]) : ky for ky, val in res.items()}
        for fut in futures.as_completed(future):
            ky = future[fut]
            val = fut.result()
            res[ky] = res[ky][0] * val
    return res, dg

@cached_function
def ApplySingle(A, i, z, M, check=True):
    Rp = z.parent().base_ring()
    if not Rp.is_finite():
        Rp = Rp.base_ring()
    p = ZZ(len(Rp)).factor()[0][0]
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
    t1 = 0
    t2 = 0
    for inky, s1, s2, sgn in input_list[outky]:
        f = gF[inky]
        t = cputime()
        newres = sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] for fi, i in zip(f.coefficients(), f.exponents()))
        t1 += cputime(t)
        t = cputime()
        if sgn == 1:
            res *= newres
        else:
            resinv *= newres
        t2 += cputime(t)
    t = cputime()
    ans = res * inv(resinv)
    t3 = cputime(t)
    return ans, (t1,t2,t3)


def degrees(res):
    dumax = max(ff.degree() for ff in res.values())
    dvmax = max(o.degree() for ff in res.values() for o in ff.coefficients())
    return dumax, dvmax


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
    vv = vector([-1, tau[0]])
    ww = vector([tau[1], 1])
    KK = vv.parent().base_ring()
    return prod((vv * num.apply_morphism(phi).apply_map(lambda o : KK(o)) * ww)**ZZ(exponent) for num, exponent in L0)


def Eval(tau, prec):
    global J
    t0, t1 = tau
    p = t0.parent().prime()
    ans = 1
    x0s = {[1] for _ in range(p+1)}
    x1s = {[1] for _ in range(p+1)}
    for i in range(p+1):
        x0 = t0 if i == p else 1/(t0 - i)
        x1 = t1 if i == p else 1/(t1 - i)
        for _ in range(prec+1):
            x0s[i].append(x0 * x0s[i][-1])
            x1s[i].append(x1 * x1s[i][-1])
    for ky in J:
        f = J[ky]
        x00 = x0s[ky[0]]
        x11 = x1s[ky[1]]
        ans *= sum(sum(o * x00[j] for o, j in zip(a.coefficients(), a.exponents()) if j < prec) * x11[i] \
            for a, i in zip(f.coefficients(), f.exponents()) if i < prec)
    return ans


def excellent_matrices(m, coset_list):
    a, b, c, d = m.list()
    F = m.parent().base_ring()
    one = coset_list[0]
    assert one == 1
    r = F.gen()
    I = F.ideal(2)
    if F == QQ or F == ZZ:
        raise NotImplementedError
    if c in I:
        return [ (m,-1, one) ]
    m = matrix([[b, a], [d, c]])        
    g = find_coset_rep(m, coset_list)
    if g is not None:
        return [ (m, 1, g) ]
    else:
        delta = -(m.determinant())
        A, B, C, D = m.list()
        for d1, d2 in cartesian_product_iterator([[1,-1,r,-r] for _ in range(2)]):
            y = matrix([[-C, d1], [D, d2]]).determinant() / delta
            if y in I:
                x = matrix([[d1, A], [d2, -B]]).determinant() / delta
                m1 = matrix([[x,A],[y,C]])
                m2 = matrix([[x,B],[y,D]])
                m1.rescale_col(0, -d1**-1)
                m2.rescale_col(0, d2**-1)
                assert m1.determinant() == 1
                assert m2.determinant() == 1
                return [(m1,1, one),(m2,-1, one)]
    raise RuntimeError(f'No good matrix found for m = {m.list()}!')

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
        delta = -(m.determinant())
        A, B, C, D = m.list()
        for d1, d2 in cartesian_product_iterator([[1,-1,r,-r] for _ in range(2)]):
            y = matrix([[-C, d1], [D, d2]]).determinant() / delta
            if y in I:
                x = matrix([[d1, A], [d2, -B]]).determinant() / delta
                m1 = matrix([[A,x],[C,y]])
                m2 = matrix([[x,B],[y,D]])
                m1.rescale_col(1, -d1**-1)
                m2.rescale_col(0, d2**-1)
                m1.swap_columns(0,1)
                assert m1.determinant() == 1
                assert m2.determinant() == 1
                return [(m1,1),(m2,-1)]
    raise RuntimeError(f'No good matrix found for m = {m.list()}!')

def smallCMcycle(D):
    A = compute_gamma_tau(D).change_ring(F)
    t0, t1 = fixed_point(A, phi)
    hE = QuadraticField(-D,'w').class_number()
    return A, [t0, t1], hE

def smallRMcycle(D):
    A = compute_gamma_tau(D).change_ring(F)
    t0, t1 = fixed_point(A, phi)
    hE = QuadraticField(D,'w').class_number()
    return A, [t0, t0], hE # not a typo!

def bigRMcycle_old(D, alpha=None, n=1):
    if alpha is None:
        alpha = ATR_alpha(D, n)
    else:
        assert n == 1, 'n would not be used'
    A, tau0, tau1 = compute_gamma_tau_ATR(phi, alpha)
    E = tau_ATR_field(alpha)
    return A, [tau0, tau1], E.class_number()

def bigRMcycle(D, alpha=None, n=1):
    assert alpha is None and n == 1
    try:
        alpha = get_F_and_alpha(D)[1]
    except KeyError:
        raise ValueError(f'No bigRM cycle computed for discriminant {D}')
    A, tau0, tau1 = compute_gamma_tau_ATR(phi, alpha)
    E = tau_ATR_field(alpha)
    return A, [tau0, tau1], E.class_number()

def good_matrices_random(m, limit=10**5):
    r = m[0,0].parent().gen()
    for it, c in enumerate(random_candidate_matrices(r, limit=limit)):
        try:
            return good_matrices(c * m), c
        except RuntimeError:
            continue
    raise RuntimeError('No good matrix found')

def RMCEval(D, cycle_type, prec, alpha=None, n=1, return_class_number=False):
    global L0, J, ncpus
    try:
        ncpus = ncpus
    except NameError:
        ncpus = cpu_count()
    if cycle_type == 'smallCM':
        A, tau0, hE = smallCMcycle(D)
    elif cycle_type == 'smallRM':
        A, tau0, hE = smallRMcycle(D)
    else:
        if cycle_type != 'bigRM':
            raise NotImplementedError('Cycle type should be either "smallCM" or "smallRM" or "bigRM"')
        A, tau0, hE = bigRMcycle(D, alpha=alpha, n=n)

    if any(t.trace() == 2 * t for t in tau0):
        raise ValueError(f'Tuple {(D, alpha) = } is not admissible for {cycle_type} cycle: the resulting tau is not in Hp.')

    found = False 
    r = A.parent().base_ring().gen()
    mlist0 = matrices_for_unimodular_path(A[0,0], A[1,0])
    try:
        mlist = []
        for m in mlist0:
            tmp, c = good_matrices_random(m,limit=100)
            mlist.extend(tmp)
    except RuntimeError:
        print('gamma_tau = ', A.list())
        raise RuntimeError

    print('Now computing...')
    res0 = prod(Eval0(L0, act_matrix(m.apply_morphism(phi).adjugate(), tau0))**sgn for m,sgn in mlist)
    if parallelize:
        res1 = 1
        with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
            future_dict = {executor.submit(Eval, act_matrix(m.apply_morphism(phi).adjugate(), tau0), prec) : sgn for m, sgn in mlist}
            for fut in futures.as_completed(future_dict):
                sgn = future_dict[fut]
                res1 *= fut.result()**sgn
    else:
        res1 = prod(Eval(act_matrix(m.apply_morphism(phi).adjugate(), tau0), prec)**sgn for m, sgn in mlist)
    ans = res0 * res1
    ans = ans.add_bigoh(prec + ans.valuation())
    if return_class_number:
        return ans, hE
    else:
        return ans


def change_sign_in_matrices(ms):
    return [(m, -s) for m, s in ms]

def random_candidate_matrices(r, limit=-1):
    found_matrices = 0
    yield matrix(2,2,1)
    FF = r.parent()
    i = 0
    while i != limit:
        i += 1
        a = FF.random_element()
        co2 = FF.random_element()
        c = 2*co2
        x,d,b = xgcd_quadratic(a, c)
        if not x in [1, -1, r, -r]:
            continue
        d = d/x
        b = b/x
        m = matrix(2,2,[a,-b,c,d])
        yield m
        assert m.det() == 1
        m2 = matrix(2,2,[r*a, -b, c, -r*d])
        assert m2.det() == 1
        yield m2
    return


# given a non-necessary fundamental discriminant D, computes the matrix gamma_tau associated to an optimal embedding of the ring of discriminant D to M_0(2)
def compute_gamma_tau(D):
    Dsqf = ZZ(D).squarefree_part()
    try:
        c = ZZ((ZZ(D) / fundamental_discriminant(D)).sqrt())
    except TypeError:
        raise ValueError('D is not the discrimininant of a quadratic order')
    F = QuadraticField(Dsqf, names='r')
    w = F.maximal_order().ring_generators()[0]
    # O_c has Z-basis <1, cw>, we use this basis to embed O_c into M_2(Z)
    coords = (c*w).coordinates_in_terms_of_powers()
    D0 = Matrix([[0,1], coords((c*w)**2)]).transpose() # the image of cw in M_2(Z)
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
    if eps.norm() == -1:
        eps = eps**2 # DEBUG: is this really needed?
    n = 1
    while True:
        u = eps**n
        cs = coords(u)
        if all([o.is_integer() for o in cs]):
            gamma_tau = cs[0] + cs[1] * D1
            assert gamma_tau.minpoly() == u.minpoly()
            assert gamma_tau[1][0] % 2 == 0
            return gamma_tau
        n +=1


def ATR_alpha(D, n=1):
    try:
        return QuadraticField(D, names='t').elements_of_norm(-1)[0] * n
    except IndexError:
        raise ValueError(f'Discriminant (={D}) not admissible: no elements of norm -1 in QQ(sqrt(D))')


def tau_ATR_field(alpha, names='w'):
    FF = alpha.parent()
    y = FF['y'].gen()
    return NumberField(y*y - alpha, names=names)

# it accepts a real quadratic field FF and an element alpha in FF of norm -1; then E = FF(sqrt(alpha) is the ATR extension and K = Q(i) is contained in the galois closure of E
def compute_gamma_tau_ATR(phi, alpha):
    Kp = phi.codomain()
    p = Kp.prime()
    # FF = alpha.parent()
    # x = QQ['x'].gen()
    # R.<y> = PolynomialRing(FF)
    E = tau_ATR_field(alpha,'w')
    w = E.gen()
    K = F
    i = K.gen()
    MM = E.galois_closure(names = 'gMM')
    if MM.disc() % p == 0:
        raise NotImplementedError('tau lives in a ramified extension') # p ramifies in MM so we would need to work with ramified extensions
    # Now we redefine E, because we want to view it as a subfield of MM. We also construct L as a subfield of MM
    EEp = [f for f in MM.subfields() if f[0].degree() == 4 and f[0].is_isomorphic(E)]
    LLp = [f for f in MM.subfields() if f[0].degree() == 4 and not f[0].is_isomorphic(E)]
    E, sigmaE, _ = EEp[0]
    L, sigmaL, _ = LLp[1]

    # take the unit of norm 1 in E
    sigma = E.automorphisms()[1]
    found = False
    for uu in E.units():
        if uu.absolute_minpoly().degree() == 4 and uu*sigma(uu) == 1:
            u = uu
            found = True
    if found == False:
        raise RuntimeError('did not find a unit of relative norm one')
    # from u construct gamma and its galois conjugate gamma_p
    gL = L.primitive_element()
    Gal_M_L = [t for t in MM.automorphisms() if t(sigmaL(gL)) == sigmaL(gL)]
    gE = E.primitive_element()
    Gal_M_E = [t for t in MM.automorphisms() if t(sigmaE(gE)) == sigmaE(gE)]
    Nu = Gal_M_L[1](sigmaE(u)) * sigmaE(u) # compute N_{M/L}(u)
    L_K = L.relativize(K.embeddings(L)[0], names = 'gL_K')
    gL_K = L_K.gen()
    u_L_K = Nu.minpoly().roots(L_K)[0][0] # this is Nu viewed as an element of the relative extension L_K
    a, c = u_L_K.vector()
    b, d = (u_L_K * gL_K).vector()
    gamma = Matrix(2,2,[a,b,c,d])
    # this matrix is constructed using the K-basis <1, gL_K>, but gL_K is not a generator of the ring of integers of L/K, so we need to change basis
    _, gOLK = module_generators(L_K) # gOLK is a generator of O_L_K over O_K
    O = L_K.maximal_order()
    X = Matrix([(1,0), gOLK.vector()]).transpose() # change of basis matrix
    # we check that gOLK is a generator of O_L_K over O_K, by checking that any generator of O_L_K when written in terms of <1, gOLK> has integral coefficients
    for gen in O.gens():
        if not all([o.is_integral() for o in (X.inverse()* vector(gen.vector())).list() ]):
            raise RuntimeError('it seems that we do not have an embedding of the maximal order')
    # now we conjugate gamma to express it in terms of the basis <1, gOLK>
    gamma = X**-1 * gamma * X
    Ms = [Matrix(2,2,[1,0,0,1]), Matrix(2,2,[0,1,-1,0]), Matrix(2,2,[2,1,1,1]), Matrix(2,2,[1+i,1,i,1]), Matrix(2,2,[1+i,1,1,1]), Matrix(2,2,[i,1,0,1]), Matrix(2,2,[i,1,2*i,1])]
    found = False
    for mm in Ms:
        aa, bb, cc, dd = (mm * gamma * ~mm).list()
        if cc.real() % 2 == 0 and cc.imag() % 2== 0:
            found = True
            gamma = mm * gamma * ~mm
            break
    if not found:
        raise RuntimeError('maybe there is no embedding into Gamma0(2)?')
    gamma_p = gamma.apply_morphism(a.parent().automorphisms()[1])
    assert u_L_K.minpoly() == gamma.minpoly()

    # Now we compute the matrix M and the scalar l such that gamma * M * gamma_p^{-1} = l * M
    A = PolynomialRing(K, names='xx,yy,bb,cc,ll')
    xx, yy, bb, cc, ll = A.gens()
    alpha = xx + i*yy
    alpha_p = xx-i*yy
    M = Matrix(A,2,2,[alpha, -bb,cc,-alpha_p])
    B = PolynomialRing(E, names='x,y,b,c,l')
    x,y,b,c,l = B.gens()
    fAB = A.hom([x,y,b,c,l], B, check = False)
    s = K.automorphisms()[1]
    Meq = gamma*M - ll*M*gamma_p
    eqns =[fAB(M.det())] +  [fAB((z + z.map_coefficients(s))/2) for z in Meq.list()] + [fAB((z - z.map_coefficients(s))/(2*i)) for z in Meq.list()]  + [c-1]
    I = B.ideal(eqns)
    pt = I.variety()[0]
    alpha = sigmaE(pt[x]) + K.embeddings(MM)[0](i) * sigmaE(pt[y])
    alpha_p = sigmaE(pt[x]) - K.embeddings(MM)[0](i) * sigmaE(pt[y])
    M = Matrix(MM, 2, 2, [alpha, -sigmaE(pt[b]), sigmaE(pt[c]), -alpha_p])
    Mp = Matrix(MM, 2, 2, [alpha_p, -sigmaE(pt[b]), sigmaE(pt[c]), -alpha])
    iKM = K.embeddings(MM)[0]
    # assert that M is fiexed under the action of gamma
    assert gamma.apply_map(iKM) * M * gamma_p.apply_map(iKM)**-1 == sigmaE(pt[l]) * M
    tau1 = M[0][0]
    tau2 = -M[1][1]
    # Now we compute the imge of tau1 and tau2 in Qp**2
    Kp = phi.codomain()
    L = Qq(p**2, Kp.precision_cap(),names='b')
    roots = [a[0] for a in MM.defining_polynomial().roots(L)]
    if len(roots) == 0:
        raise RuntimeError('tau does not live in Qp**2')
    # we find an embedding of MM into Qp**2 that extends phi
    found_emb = False
    for r in roots:
        iMML = MM.hom([r],L)
        if iMML(iKM(K.gen())) == L(phi(K.gen())):
            found_emb = True
            break
    if found_emb == False:
        raise RuntimeError('No embedding of MM into Qp**2 found that extends phi')
    return gamma, iMML(tau1), iMML(tau2)


def initial_seed(F, v, p):
    if isinstance(v, str):
        v = functional_from_label(v)
    v = tuple(v)
    L0 = []
    L1 = []
    for i, vi in enumerate(v, start=1):
        if vi == 0:
            continue
        Li = vectors_in_lattice(F, p, i)
        Lpi = vectors_in_lattice(F, p, p * i)
        L0.extend([(o, vi * sgn) for o, sgn in Li])
        L1.extend([(o, vi * sgn) for o, sgn in Lpi])
    assert sum(e for o, e in L0) == 0
    assert sum(e for o, e in L1) == 0
    return L0, L1


def find_value_one(maxD, cycle_type='smallCM'):
    stats = {'0':[], '1':[], 'oo':[], '?':[], 'err':[]}
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
                stats['1'].append(D)
            elif x1 == 0:
                print(f'{D = } yields zero')
                stats['0'].append(D)
            elif x1.valuation() < 40:
                print(f'{D = } yields infinity')
                stats['oo'].append(D)
            else:
                print(f'{D = } yields {x1}')
                stats['?'].append(D)
        except (TypeError, ValueError, RuntimeError) as e:
            stats['err'].append((D,str(e)))
        except PrecisionError:
            stats['oo'].append(D)
    return stats
