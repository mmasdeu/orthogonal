#/usr/bin/python
# -*- coding: utf-8 -*-
# File: dglpoint.py

from fire import Fire
from stopit import ThreadingTimeout
from sage.all import SageObject, srange, QuadraticField, load, Integer, parallel, GF, walltime, vector, save, fundamental_discriminant, Qq,  ModularForms, IntegralLattice, MatrixSpace
from sage.rings.padics.precision_error import PrecisionError
from multiprocessing import cpu_count
from concurrent import futures
from sage.misc.timing import cputime
from collections import defaultdict
from darmonpoints.divisors import Divisors
from darmonpoints.util import *
from dglutil import *

from time import sleep
import stopit

from xml.parsers.expat import ExpatError
from sage.rings.padics.precision_error import PrecisionError
from cysignals.signals import SignalError
from builtins import ModuleNotFoundError

G = None

x = QQ['x'].gen()
valid_ATR_fields = dict([(41, x**4 + 4*x**3 + x**2 - 6*x - 8),
                         (113, x**4 - 28*x**2 - 256),
                         (137, x**4 + 4*x**3 + 105*x**2 + 202*x - 224),
                         (313, x**4 + 4*x**3 - 111*x**2 - 230*x - 3032),
                         (337, x**4 + 4*x**3 - 3*x**2 - 14*x - 72),
                         (409, x**4 - 12*x**2 - 1600),
                         (457, x**4 + 4*x**3 - 15*x**2 - 38*x - 24),
                         (521, x**4 - 44*x**2 - 1600),
                         (569, x**4 + 4*x**3 - 111*x**2 - 230*x - 8216),
                         (593, x**4 - 92*x**2 - 256),
                         (809, x**4 + 4*x**3 + x**2 - 6*x - 200),
                         (857, x**4 + 4*x**3 - 255*x**2 - 518*x - 584),
                         (881, x**4 + 4*x**3 - 219*x**2 - 446*x - 5408),
                         (1153, x**4 + 4*x**3 - 27*x**2 - 62*x - 48),
                         (1201, x**4 + 4*x**3 - 219*x**2 - 446*x - 11888),
                         (1217, x**4 - 124*x**2 - 1024),
                         (1249, x**4 - 60*x**2 - 4096),
                         (1321, x**4 + 4*x**3 - 39*x**2 - 86*x - 26288),
                         (1553, x**4 - 92*x**2 - 4096),
                         (1657, x**4 - 68*x**3 - 1276*x**2 + 135712*x + 1478656),
                         (1777, x**4 - 156*x**2 - 1024),
                         (1889, x**4 + 4*x**3 - 11*x**2 - 30*x - 416),
                         (1993, x**4 - 172*x**2 - 576),
                         (2113, x**4 + 4*x**3 - 291*x**2 - 590*x - 21032), (2129, x**4 - 92*x**2 - 6400), (2273, x**4 + 4*x**3 + 53*x**2 + 98*x + 32), (2377, x**4 + 4*x**3 - 183*x**2 - 374*x - 39392), (2521, x**4 - 140*x**2 - 5184), (2593, x**4 + 4*x**3 - 147*x**2 - 302*x - 46808), (2617, x**4 + 4*x**3 + 465*x**2 + 922*x + 136), (2657, x**4 + 4*x**3 - 435*x**2 - 878*x - 5624), (2689, x**4 + 4*x**3 - 27*x**2 - 62*x - 432), (2729, x**4 + 4*x**3 - 39*x**2 - 86*x - 54800), (2833, x**4 - 92*x**2 - 9216), (2953, x**4 + 4*x**3 - 47*x**2 - 102*x - 88), (3001, x**4 - 204*x**2 - 1600), (3089, x**4 + 4*x**3 + 501*x**2 + 994*x - 800), (3217, x**4 + 4*x**3 - 3*x**2 - 14*x - 792), (3361, x**4 - 60*x**2 - 12544), (3433, x**4 - 100*x**3 - 2684*x**2 + 369056*x + 6718464), (3593, x**4 + 4*x**3 - 47*x**2 - 102*x - 248), (3761, x**4 + 4*x**3 - 219*x**2 - 446*x - 63728), (3881, x**4 - 236*x**2 - 1600), (3929, x**4 - 140*x**2 - 10816), (4049, x**4 - 220*x**2 - 4096), (4073, x**4 + 4*x**3 - 327*x**2 - 662*x - 55088), (4177, x**4 + 4*x**3 - 75*x**2 - 158*x - 83024), (4273, x**4 + 4*x**3 - 507*x**2 - 1022*x - 21248), (4289, x**4 + 4*x**3 - 59*x**2 - 126*x - 80), (4513, x**4 - 188*x**2 - 9216), (4657, x**4 - 148*x**3 - 476*x**2 + 589472*x + 8856576), (4721, x**4 + 4*x**3 - 219*x**2 - 446*x - 83168), (4793, x**4 + 4*x**3 - 111*x**2 - 230*x - 93752), (4801, x**4 + 4*x**3 - 59*x**2 - 126*x - 208), (4817, x**4 + 4*x**3 - 363*x**2 - 734*x - 63872), (4969, x**4 + 4*x**3 - 327*x**2 - 662*x - 73232), (4993, x**4 - 244*x**3 + 13348*x**2 + 347168*x + 589824), (5233, x**4 - 28*x**2 - 20736), (5393, x**4 + 4*x**3 - 651*x**2 - 1310*x - 1952), (5641, x**4 + 4*x**3 + 81*x**2 + 154*x + 72), (5657, x**4 + 4*x**3 - 543*x**2 - 1094*x - 39752), (5801, x**4 + 4*x**3 + x**2 - 6*x - 1448), (5897, x**4 - 44*x**2 - 23104), (6073, x**4 + 4*x**3 - 71*x**2 - 150*x - 112), (6217, x**4 + 4*x**3 - 183*x**2 - 374*x - 117152), (6329, x**4 + 4*x**3 - 71*x**2 - 150*x - 176), (6353, x**4 + 4*x**3 - 651*x**2 - 1310*x - 21392), (6449, x**4 - 28*x**2 - 25600), (6473, x**4 - 164*x**3 - 2172*x**2 + 936608*x + 19784704), (6529, x**4 + 4*x**3 - 579*x**2 - 1166*x - 47240), (6689, x**4 + 4*x**3 - 11*x**2 - 30*x - 1616), (7001, x**4 - 140*x**2 - 23104), (7121, x**4 - 220*x**2 - 16384), (7193, x**4 - 260*x**3 + 12036*x**2 + 862496*x + 5914624), (7321, x**4 + 4*x**3 - 543*x**2 - 1094*x - 73448), (7369, x**4 + 4*x**3 - 759*x**2 - 1526*x - 3680), (7393, x**4 - 188*x**2 - 20736), (7417, x**4 - 76*x**2 - 28224), (7489, x**4 + 4*x**3 - 291*x**2 - 590*x - 129896), (7793, x**4 - 28*x**2 - 30976), (7841, x**4 - 316*x**2 - 6400), (8009, x**4 + 4*x**3 - 79*x**2 - 166*x - 280), (8089, x**4 - 268*x**2 - 14400), (8161, x**4 + 4*x**3 - 75*x**2 - 158*x - 480), (8209, x**4 - 220*x**2 - 20736), (8273, x**4 - 84*x**3 - 13532*x**2 + 907168*x + 58491904), (8297, x**4 + 4*x**3 + 825*x**2 + 1642*x + 496), (8329, x**4 - 300*x**2 - 10816), (8369, x**4 + 4*x**3 - 19*x**2 - 46*x - 1960), (8377, x**4 - 204*x**2 - 23104), (8521, x**4 + 4*x**3 - 759*x**2 - 1526*x - 27008), (8609, x**4 - 188*x**2 - 25600), (8681, x**4 - 364*x**2 - 1600), (9137, x**4 - 284*x**2 - 16384), (9257, x**4 - 228*x**3 + 1924*x**2 + 1558432*x + 30647296), (9377, x**4 - 316*x**2 - 12544), (9433, x**4 + 4*x**3 - 831*x**2 - 1670*x - 16712), (9473, x**4 + 4*x**3 - 867*x**2 - 1742*x - 2168), (9497, x**4 + 4*x**3 - 55*x**2 - 118*x - 1504), (9601, x**4 - 380*x**2 - 2304), (9689, x**4 - 132*x**3 - 12284*x**2 + 1408288*x + 69222400), (9697, x**4 + 4*x**3 - 75*x**2 - 158*x - 864), (9817, x**4 + 4*x**3 + 105*x**2 + 202*x + 96), (9929, x**4 + 4*x**3 - 759*x**2 - 1526*x - 55520), (10009, x**4 - 12*x**2 - 40000), (10169, x**4 + 4*x**3 - 111*x**2 - 230*x - 202616), (10177, x**4 - 124*x**2 - 36864), (10337, x**4 - 316*x**2 - 16384), (10369, x**4 - 252*x**2 - 25600), (10433, x**4 + 4*x**3 - 91*x**2 - 190*x - 352), (10513, x**4 + 4*x**3 - 651*x**2 - 1310*x - 105632)])

def get_p_prec(f):
    fsp = f.split('_')
    p = ZZ(fsp[1])
    M = ZZ(fsp[-1].split('.')[0])
    return p, M

def get_label(f):
    fsp = f.split('_')
    return fsp[2] + '_' + fsp[3]

def label_from_functional(func, sep="·"):
    pos = ""
    neg = ""
    for i, o in enumerate(func):
        if o > 0:
            pos = pos + o * (sep + str(i + 1))
        elif o < 0:
            neg = neg + (-o) * (sep + str(i + 1))
    return pos[1:] + "_" + neg[1:]


def functional_from_label(label, sep="·"):
    func = {}
    pos, neg = label.split("_")
    for i in pos.split(sep):
        try:
            func[int(i) - 1] += 1
        except KeyError:
            func[int(i) - 1] = 1
    for i in neg.split(sep):
        try:
            func[int(i) - 1] -= 1
        except KeyError:
            func[int(i) - 1] = -1
    lfunc = max(func.keys()) + 1
    ans = [0 for _ in range(lfunc)]
    for ky, val in func.items():
        ans[ky] = val
    return tuple(ans)

def random_point_on_A0(K):
    p = K.prime()
    L = Qq(p**2, K.precision_cap(),names='b')    
    a = L(1)
    while any((a - i).valuation() > 0 for i in range(p)):
        try:
            a = L.random_element().unit_part()
        except ValueError:
            continue
    b = a.trace() - a
    return a,b

def get_alpha(D):
    try:
        f = valid_ATR_fields[D]
    except KeyError:
        raise ValueError(f'No bigRM cycle computed for discriminant {D}')

    K = NumberField(x*x - D, names ='t')
    E = NumberField(f, names = 'gE')
    EK = E.relativize(K.embeddings(E)[0], 'gEF')
    r0 = f.roots(EK)[0][0]
    alpha = r0.minpoly().discriminant()
    assert alpha.relative_norm().squarefree_part() == -1
    return alpha


class DGLGroup(SageObject):
    def __init__(self, F, P, level=2, parallelize=True):
        self.F = F
        (self.prime, ) = F.ideal(P).gens_reduced()
        self.p = ZZ(self.prime.norm())
        assert self.p.is_prime()
        self.level = F.ideal(level)
        self.parallelize = parallelize
        self.reduced_reps = {}
        coset_list = self.get_coset_list()
        self.reduced_rep_list = [coset_list[0]]
        # coset_list contains a dictionary mapping
        # g0 |-> [ (hi, ni, gi) ]_i such that
        # {g0 0 -> g0 oo} = sum_i n_i {higi 0 -> higi oo}
        # with hi in Gamma0(level) and gi reduced reps.
        for g0 in coset_list:
            try:
                ans = self.good_matrices(g0, coset_list=self.reduced_rep_list)
                self.reduced_reps[g0] = ans
            except RuntimeError:
                # g0 becomes a new reduced_rep
                self.reduced_rep_list.append(g0)
                one = coset_list[0]
                assert one == 1
                self.reduced_reps[g0] = [(one, 1, g0)]

    @cached_method
    def get_coset_list(self):
        r'''
        Find all cosets of the level.
        '''
        F = self.F
        index = self.level.norm()
        for P, _ in self.level.factor():
            index *= (1 + 1/P.norm())
        id = Matrix(F,2,2,1)
        id.set_immutable()
        ans = [id]
        bound = 1       
        while len(ans) < index:
            bound *= 2
            for m in self.enumerate_matrices(bound):
                if self.coset_rep(m, ans) is None:
                    m.set_immutable()
                    ans.append(m)
                if len(ans) == index:
                    break
        return ans

    def coset_rep(self, m, coset_list=None):
        assert m.determinant() == 1
        if coset_list is None:
            coset_list = self.get_coset_list()
        minv = m.adjugate()
        for c in coset_list:
            if (c * minv)[1,0] in self.level:
                return c
        return None

    def good_matrices(self, m, coset_list=None):
        assert m.determinant() == 1
        if coset_list is None:
            coset_list = self.get_coset_list()
        a, b, c, d = m.list()
        F = m.parent().base_ring()
        one = coset_list[0]
        assert one == 1
        r = F.gen()
        I = self.level
        if F == QQ or F == ZZ:
            raise NotImplementedError
        if c in I:
            assert m[1,0] in I
            return [ (m, 1, one) ]
        elif d in I:
            m1 = matrix([[b, -a], [d, -c]])
            assert m1[1,0] in I
            return [ (m1, -1, one) ]
        
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
                assert m1[1,0] in I
                assert m2[1,0] in I
                return [(m1, -1, one),(m2, 1, one)]
        g = self.coset_rep(m, coset_list=coset_list)
        if g is not None:
            gamma = m * g.adjugate()
            assert gamma[1,0] in I
            return [ (gamma, 1, g) ]            
        m = matrix([[b, -a], [d, -c]])
        g = self.coset_rep(m, coset_list=coset_list)
        if g is not None:
            gamma = m * g.adjugate()
            assert gamma[1,0] in I
            return [ (gamma, -1, g) ]
        raise RuntimeError(f'No good matrix found for m = {m.list()}!')

    def enumerate_matrices(self, bound=10):
        for a in self.F.elements_of_bounded_height(bound=bound):
            if not a.is_integral():
                continue
            for c in self.F.elements_of_bounded_height(bound=bound):
                if not c.is_integral():
                    continue
                gg, d, mb = xgcd_quadratic(a, c)
                if gg.norm() != 1:
                    continue
                M = Matrix(2,2,[a,-mb / gg,c,d / gg])
                assert M.determinant() == 1
                M.set_immutable()
                yield M


class Cocycle(SageObject):
    def __init__(self, functional, prec, initialize=True):
        self.prec = prec
        self.p = G.p
        self.G = G
        p = self.p
        self.F = G.F
        Fgen = self.F.gen()
        self._keylist = tuple((i,j) for i in range(self.p+1) for j in range(self.p+1))
        if isinstance(functional, str):
            functional = functional_from_label(functional)
        self.functional = tuple(functional)
        self.Fp = Qp(p,2*prec, type='floating-point')
        self.Rp = Zmod(p**prec)
        phi = self.F.hom([Fgen.minpoly().roots(self.Fp)[0][0]],check=False)
        if phi(G.prime).valuation() == 0:
            phi = self.F.hom([Fgen.minpoly().roots(self.Fp)[1][0]],check=False)
        self.phi = phi
        self.Ruv = PolynomialRing(PolynomialRing(self.Rp,'u'),'v')
        self.map_poly = self.Ruv.flattening_morphism().inverse()
        self.L0 = None
        self.datalist = [None]
        if initialize:
            self._initialize_data()

    def _initialize_data(self, check=True):
        Div = Divisors(MatrixSpace(self.F,2,2))
        v = self.functional
        self.G.get_coset_list()
        reps = self.G.reduced_rep_list
        L0 = {g : Div([]) for g in reps}
        L1 = {g : Div([]) for g in reps}
        for i, vi in enumerate(v, start=1):
            if vi == 0:
                continue
            for g in reps:
                Li = vectors_in_lattice(self.F, self.p, i, g)
                Lpi = vectors_in_lattice(self.F, self.p, self.p * i, g)
                L0[g] += Div([(vi * sgn, o) for o, sgn in Li])
                L1[g] += Div([(vi * sgn, o) for o, sgn in Lpi])
        degrees = [L0[g].degree() for g in reps]
        assert all(d == 0 for d in degrees), degrees
        degrees = [L1[g].degree() for g in reps]
        assert all(d == 0 for d in degrees), degrees
        if check:
            for g in reps:
                r, s = act_flt(g, Infinity), act_flt(g, self.F(0))
                for M, sgn in L0[g]:
                    alphaconj, mb, c, _ = M.list()
                    alpha = alphaconj.conjugate()
                    b = - mb
                    sqradius = (alpha.norm() - b*c) / (c*c)
                    center = alphaconj / c
                    qf = lambda x : 1 if x == Infinity else (x-center).norm() - sqradius
                    assert qf(r) * qf(s) < 0
        self.L0 = L0
        self.L1 = L1
        newdata = {}        
        for g in reps:
            lev1 = self.level1(L1[g])
            for ky, val in lev1.items():
                newdata[(g, *ky)] = val
        self.datalist.append(newdata)
        return

    # {x -> y}
    def manin_trick(self, x, y, check=True):
        F = self.F
        G = self.G
        try:
            y = F(y)
            ynum = F(y.numerator())
            yden = F(y.denominator())
        except TypeError:
            ynum = F(1)
            yden = F(0)
            y = Infinity
        try:
            x = F(x)
            xnum = F(x.numerator())
            xden = F(x.denominator())
        except TypeError:
            xnum = F(1)
            xden = F(0)
            x = Infinity
        ylist = matrices_for_unimodular_path(ynum, yden) # {oo -> y}
        xlist = matrices_for_unimodular_path(xnum, xden) # {oo -> x}
        ans = []
        for g in ylist:
            gi = G.coset_rep(g)
            m = g * gi.adjugate()
            assert m * gi == g
            for h, sgn, gj in G.reduced_reps[gi]:
                assert h[1,0] in F.ideal(2) # DEBUG
                ans.append((m * h, -sgn, gj))
        for g in xlist:
            gi = G.coset_rep(g)
            m = g * gi.adjugate()
            assert m * gi == g
            for h, sgn, gj in G.reduced_reps[gi]:
                assert h[1,0] in F.ideal(2) # DEBUG
                ans.append((m * h, sgn, gj))
        if check:
            # Check that the {x -> y} (y - x) is the same as the returned sum applied to {oo -> 0}
            res = defaultdict(ZZ)
            res[x] = 1
            res[y] = -1
            # We will check that the returned value equals {x -> y}. So we have initialized res to (x) - (y)
            for h, sgn, gi in ans:
                m = h * gi
                a, b, c, d = m.list()
                # Store in res the contribution {a/c -> b/d}
                if c == 0:
                    res[Infinity] -= sgn
                else:
                    res[a / c] -= sgn
                if d == 0:
                    res[Infinity] += sgn
                else:
                    res[b / d] += sgn
            assert all(o == 0 for o in res.values()), f'Error in Manin trick: {x = }, {y = } {res = }\n{ans = }'
        return ans

    def get_item(self, r, s, level=None, lift=False):
        # Write r->s as sum g0 -> goo using Manin trick
        # phigamma = lambda M : gamma * M * ~gamma.apply_map(lambda x : F(x).conjugate())
        F = self.F
        Div = Divisors(MatrixSpace(F,2,2))
        D = Div([])
        ans = {ky : 1 for ky in self._keylist}
        for h, exponent, gi in self.manin_trick(r,s):
            res = self._eval_translate_inf_zero(h, gi, level=level)
            phi_h = lambda M : h * M * ~h.apply_map(lambda x : F(x).conjugate())
            D += ZZ(exponent) * self.L0[gi].apply_map(phi_h)
            for ky, val in res.items():
                if exponent > 0:
                    ans[ky] *= val**exponent
                else:
                    ans[ky] = mul_inv(ans[ky], val**(-exponent))
        if lift:
            R0 = self.Ruv.base_ring()
            psi = R0.hom([ZZ['u'].gen()],base_map=lambda x:x.lift(),check=False)
            ans = {ky : f.change_ring(psi) for ky, f in ans.items()}        
        return ans, D
    
    def evaluate(self, r, s, tau, prec):
        # Write r->s as sum g0 -> goo using Manin trick
        # phigamma = lambda M : gamma * M * ~gamma.apply_map(lambda x : F(x).conjugate())
        F = self.F
        # Div = Divisors(MatrixSpace(F,2,2))
        # D = Div([])
        # ans = {ky : 1 for ky in self._keylist}
        ans = 1
        L = tau[0].parent()
        for h, exponent, gi in self.manin_trick(r,s):
            hinv = h.adjugate()
            hinv_conj = hinv.apply_map(lambda x : x.conjugate())
            hinv = hinv.apply_morphism(self.phi).change_ring(L)
            hinv_conj = hinv_conj.apply_morphism(self.phi).change_ring(L)
            hinv_tau = [act_flt(hinv, tau[0]), act_flt(hinv_conj, tau[1])]
            ans *= (self.eval_L0(gi, hinv_tau) * self.eval_PS(gi, hinv_tau, prec))**ZZ(exponent)
        return ans
    
    # Evaluate the power series part of the cocycle at h gi oo -> h gi 0
    def _eval_translate_inf_zero(self, h, gi, level=None, lift=False):
        S = self.Ruv
        R = S.base_ring()
        z1 = R.gen()
        z2 = S.gen()

        h.set_immutable()
        # Evaluate m * data[gi]
        mconj = h.apply_map(lambda x : x.trace() - x)
        A = h.apply_morphism(self.phi)
        Aconj = mconj.apply_morphism(self.phi)
        tmp_dict = {}
        for i in range(self.p+1):
            tmp_dict[i, 0] = ApplySingle(A, i, z1, self.prec, check=False)
            tmp_dict[i, 1] = ApplySingle(Aconj, i, z2, self.prec, check=False)

        transformation_list = {}
        for inky0 in range(self.p+1):
            outky0, s1 = tmp_dict[inky0, 0]
            if s1 is None:
                continue
            for inky1 in range(self.p+1):
                inky = (gi, inky0, inky1)
                outky1, s2 = tmp_dict[inky1, 1]
                if s2 is None:
                    continue
                outky = (outky0, outky1)
                try:
                    transformation_list[outky].append((inky, s1, s2))
                except KeyError:
                    transformation_list[outky] = [(inky, s1, s2)]
        ans = {}
        for outky in transformation_list:
            res = self.Ruv(1)
            for inky, s1, s2 in transformation_list[outky]:
                if level is None:
                    f = self.data[inky]
                else:
                    f = self.datalist[level][inky]
                newres = sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] \
                    for fi, i in zip(f.coefficients(), f.exponents()))
                res *= newres
            ans[outky] = res
        if lift:
            R0 = self.Ruv.base_ring()
            psi = R0.hom([ZZ['u'].gen()],base_map=lambda x:x.lift(),check=False)
            ans = {ky : f.change_ring(psi) for ky, f in ans.items()}
        return ans

    def check_cocycle_property_pseries(self, level=2, gamma=None):
        F = self.F        
        prec = self.Fp.precision_cap()
        if gamma is None:
            gamma = matrix(F,2,2,[-1,-1,2,1])
        A1 = 1
        A2 = 1
        B1 = 1
        B2 = 1
        for gi in self.G.reduced_rep_list:
            for level0 in range(2,level+2, 2):                
                LHS = self._eval_translate_inf_zero(gamma, gi, level=level0, lift=True)
                r = act_flt(gamma * gi, Infinity)
                s = act_flt(gamma * gi, 0)
                RHS = self.get_item(r, s, level0, lift=True)
                tau1 = random_point_on_A0(self.Fp)
                tau2 = random_point_on_A0(self.Fp)
                assert tau1 != tau2
                J = LHS
                A1 *= self.eval_PS(tau1, prec)
                A2 *= self.eval_PS(tau2, prec)
                J = RHS
                B1 *= self.eval_PS(tau1, prec)
                B2 *= self.eval_PS(tau2, prec)
            print(f'Quotient: {(A1 / B1) / (A2 / B2)}')


    def check_cocycle_property(self, gamma=None):
        F = self.F        
        if gamma is None:
            gamma = matrix(F,2,2,[-1,-1,2,1])
        ans = {}
        for gi in self.G.reduced_rep_list:
            r = act_flt(gamma * gi, Infinity)
            s = act_flt(gamma * gi, 0)
            Div = Divisors(MatrixSpace(self.F,2,2))
            J = self.L0
            phigamma = lambda M : gamma * M * ~gamma.apply_map(lambda x : F(x).conjugate())
            LHS = J[gi].apply_map(phigamma)
            mlist = self.manin_trick(r, s)
            RHS = Div([])
            for h, exponent, gj in mlist:
                phih = lambda M : h * M * ~h.apply_map(lambda x : F(x).conjugate())
                RHS += ZZ(exponent) * J[gj].apply_map(phih)
            if not (RHS - LHS).is_zero():
                print(f'Cocycle condition failed for {gi}!')
                print('LHS:', LHS)
                print('RHS:', RHS)
            ans[gi] = (RHS,LHS)
        return ans

    @cached_method
    def calculate_Tp_matrices(self, g=ZZ(1), Tplist=None):
        P = self.G.prime
        p = ZZ(P.norm())
        assert p.is_prime(), (P, p)
        if Tplist is None:
            Tplist = [matrix(2,2,[P, a, 0, 1]) for a in range(P.norm()) ] + [matrix(2,2,[1,0,0,P])]
        MS = []
        ginv = ~g
        for m0 in Tplist:
            m = ginv * m0
            m.set_immutable()
            # x = m0^{-1} * g · oo
            x = -m[1,1] / m[1,0] if m[1,0] != 0 else Infinity
            # y = m0^{-1} * g · 0
            y = -m[0,1] / m[0,0] if m[0,0] != 0 else Infinity
            mlist = self.manin_trick(x, y)
            if any((m * h).apply_map(lambda x : x.trace() - x).apply_morphism(self.phi).determinant().valuation() != 0 for h,_,_ in mlist):
                raise RuntimeError('Problem with embedding, try the other one.')
            for h, _, t in mlist:
                assert (m0 * h)[1,0] in self.G.level, m0 * h
            MS.extend((m0 * h, sgn, t) for h, sgn, t in mlist) # t is a restricted coset representative, m0 * h is what needs to act
        for m, _, t in MS:
            m.set_immutable()
            t.set_immutable()
        return MS

    def initialize_input_list(self):
        M = self.prec
        S = self.Ruv
        R = S.base_ring()
        z1 = R.gen()
        z2 = S.gen()
        input_list = {(h,i,j) : list() for h in self.G.reduced_rep_list for i in range(self.p+1) for j in range(self.p+1)}
        cnt = 0
        MS = [(gi, self.calculate_Tp_matrices(gi)) for gi in self.G.reduced_rep_list]
        num_matrices = sum(len(ms) for g, ms in MS)
        apply_single_dict = {}        
        for gi, ms in MS:
            for m, sgn, t in ms:
                update_progress(float(cnt+1)/num_matrices, msg="InitInput")        
                cnt += 1                
                mconj = m.apply_map(lambda x : x.trace() - x)
                A = m.apply_morphism(self.phi)
                Aconj = mconj.apply_morphism(self.phi)
                for i in range(self.p+1):
                    apply_single_dict[m,i,0] = ApplySingle(A, i, z1, M, check=False)
                    apply_single_dict[m,i,1] = ApplySingle(Aconj, i, z2, M, check=False)
                for inky0 in range(self.p+1):
                    outky0, s1 = apply_single_dict[m, inky0, 0]
                    if s1 is None:
                        continue
                    for inky1 in range(self.p+1):
                        inky = (t, inky0, inky1)
                        outky1, s2 = apply_single_dict[m, inky1, 1]
                        if s2 is None:
                            continue
                        outky = (gi, outky0, outky1)
                        input_list[outky].append((inky, s1, s2, sgn))
        self.G.input_list = input_list
        return input_list

    def eval_L0(self, gi, tau):
        vv = vector([-1, tau[0]])
        ww = vector([tau[1], 1])
        KK = vv.parent().base_ring()
        return prod((vv * P.apply_morphism(self.phi).apply_map(lambda o : KK(o)) * ww)**ZZ(n) for P, n in self.L0[gi])


    def eval_PS(self, gi, tau, prec):
        t0, t1 = tau
        p = t0.parent().prime()
        ans = 1
        x0s = {i : [1] for i in range(p+1)}
        x1s = {i : [1] for i in range(p+1)}
        for i in range(p+1):
            x0 = t0 if i == p else 1/(t0 - i)
            x1 = t1 if i == p else 1/(t1 - i)
            assert x0 == 0 or x0.valuation() >= 0
            assert x1 == 0 or x1.valuation() >= 0
            for _ in range(prec+1):
                x0s[i].append(x0 * x0s[i][-1])
                x1s[i].append(x1 * x1s[i][-1])
        for i0 in range(p+1):
            for j0 in range(p+1):
                ky = (gi, i0, j0)
                f = self.data[ky]
                x00 = x0s[i0]
                x11 = x1s[j0]
                ans *= sum(sum(o * x00[j] for o, j in zip(a.coefficients(), a.exponents()) if j < prec) * x11[i] \
                    for a, i in zip(f.coefficients(), f.exponents()) if i < prec)
        return ans

    def level1(self, V):
        p = self.p
        try:
            ncpus = ncpus
        except NameError:
            ncpus = cpu_count()
        res0 = {(i,j) : [1, 1] for i in range(p+1) for j in range(p+1)}
        res = {(i,j) : 1 for i in range(p+1) for j in range(p+1)}
        dg =  {(i,j) : 0 for i in range(p+1) for j in range(p+1)}
        input_vec = []
        for A, exponent in V:
            Ap = A.apply_map(lambda x : self.Rp(self.phi(x).lift()))
            input_vec.append((A, Ap, exponent))
        for A, Ap, exponent in input_vec:
            i, j, FF = self.compute_level1_contribution(self.matrix_to_poly(A))
            if exponent > 0:
                res0[i,j][0] *= FF**ZZ(exponent)
            else:
                res0[i,j][1] *= FF**ZZ(-exponent)
            dg[i,j] += exponent
        with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
            # Calculate inverses
            print('Calculating inverses...')
            future = {executor.submit(mul_inv, v0, v1) : ky for ky, (v0, v1) in res0.items()}
            for fut in futures.as_completed(future):
                ky = future[fut]
                val = fut.result()
                res[ky] = val
        assert all(val == 0 for val in dg.values()), list(dg.values())
        return res

    def matrix_to_poly(self, A):
        Rpol = PolynomialRing(self.F,2,names='u,v')
        t1, t2 = Rpol.gens()
        return vector(Rpol, [-1,t1]) * A * vector(Rpol, [t2,1])
    
    @parallel(ncpus=cpu_count())
    def compute_level1_contribution(self, pol):
        p = self.p
        t1, t2 = pol.parent().gens()
        polp = pol.change_ring(self.phi).change_ring(GF(p))
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
        ans = (pol.parent()(pol(tau1,tau2) * t1**d1 * t2**d2)).change_ring(self.phi)
        while j1 < 0:
            j1 += p
        while j2 < 0:
            j2 += p
        ans = self.map_poly(ans.change_ring(self.Rp))        
        return j1, j2, ans

    def next(self, **kwargs):
        global gF
        gF = self.datalist[-1]
        timing = kwargs.get('timing', False)
        p = kwargs.get('p', None)
        if p is None:
            p = self.p
        parallelize = kwargs.get('parallelize', True)
        ncpus = kwargs.get('ncpus', cpu_count())
        res = {}
        times = vector([0, 0, 0])
        if parallelize:
            with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
                # Iteration
                future_dict = {executor.submit(Transform, outky) : outky for outky in self.G.input_list}
                for fut in futures.as_completed(future_dict):
                    outky = future_dict[fut]
                    ans, (t1p, t2p, t3p) = fut.result()
                    times += vector([t1p, t2p, t3p])
                    res[outky] = ans
        else:
            assert 0
            for outky in self.G.input_list:
                res[outky], (t1p, t2p, t3p) = Transform(outky)
                times += vector([t1p, t2p, t3p])
        self.datalist.append(res)
        if timing:
            return self, tuple(times)
        else:
            return self

    def RMC(self, **kwargs):
        coset_list = kwargs.get('coset_list', None)
        parallelize = kwargs.get('parallelize', True)
        if coset_list is None:
            coset_list = G.reduced_rep_list
        p = kwargs.get('p', None)        
        if p is None:
            p = self.p
        try:
            ncpus = ncpus
        except NameError:
            ncpus = cpu_count()
        d = [-1]
        j = 0
        while sum(d) != 0:
            j += 1
            print(f'Iteration {j}')
            t = walltime()
            try:
                FF0, t0 = self.next(timing=True, **kwargs)
                FF, t1 = FF0.next(timing=True, **kwargs)
            except AttributeError:
                FF0, t0 = Cocycle.next(self, timing=True, **kwargs)
                FF, t1 = Cocycle.next(FF0, timing=True, **kwargs)
            t1, t2, t3 = vector(t0) + vector(t1)
            d = degrees(FF.datalist[-1])
            print(f'..done in {walltime(t)} seconds. {t1 = :.2f}, {t2 = :.2f}, {t3 = :.2f}.') # Degrees: {d}')
            print(f'{d}')
        print(f'Now computing product...')
        res = self.datalist[2::2]
        t = walltime()
        ans = multiply_dicts(res, parallelize=parallelize, ncpus=ncpus)
        print(f'..done in {walltime(t)} seconds.')
        R0 = self.Ruv.base_ring()
        psi = R0.hom([ZZ['u'].gen()],base_map=lambda x:x.lift(),check=False)
        self.data = {ky : f.change_ring(psi) for ky, f in ans.items()}
        return self

    def smallCMcycle(self, D):
        A = compute_gamma_tau(D).change_ring(self.F)
        t0, t1 = fixed_point(A, self.phi)
        hE = QuadraticField(-D,'w').class_number()
        return A, [t0, t1], hE

    def smallRMcycle(self, D):
        A = compute_gamma_tau(D).change_ring(self.F)
        t0, t1 = fixed_point(A, self.phi)
        hE = QuadraticField(D,'w').class_number()
        return A, [t0, t0], hE # not a typo!

    def bigRMcycle(self, D, n=1, alpha=None, version='old'):
        if version not in ['old', 'new']:
            raise ValueError('version must be either "old" or "new"')
        if version == 'old':
            if alpha is None:
                alpha = ATR_alpha(D, n)
            else:
                assert n == 1, 'n would not be used'
        else:
            assert alpha is None and n == 1
            alpha = get_alpha(D)
        try:
            A, tau0, tau1 = compute_gamma_tau_ATR(self.phi, alpha)
            E = tau_ATR_field(alpha)
        except IndexError:
            raise ValueError(f'No bigRM cycle computed for discriminant {D} (IndexError...)')
        return A, [tau0, tau1], E.class_number()

def multiply_one_dict(x,y):
    return {ky : x[ky] * y[ky] for ky in x}


def multiply_dicts(dict_list, parallelize=True, ncpus=4):
    r'''
    Multiply a list of dictionaries with the same keys.
    '''
    res = dict_list
    if parallelize:
        with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
            while len(res) > 1:
                print(len(res), 'dictionaries to multiply...')
                future_dict = {executor.submit(multiply_one_dict, res[i], res[i+1]) : i for i in range(0,len(res)-1, 2) }
                res = [] if len(res) % 2 == 0 else [res[-1]]
                for fut in futures.as_completed(future_dict):
                    res.append(fut.result())
        ans = res[0]
    else:
        ans = {ky : prod(r[ky] for r in res) for ky in res[0] }
    return ans

def get_predicted_field_and_prime_list(F, D, n, typ, char, names='z', prime_bound=600, version='old'):
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
        if version == 'new':
            alpha = get_alpha(D)
        else:
            alpha = ATR_alpha(D, n)
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
    if typ == 'bigRM':
        prime_list = [p for p in prime_range(prime_bound)]
    else:
        prime_list = [p for p in prime_range(prime_bound) if len(L.ideal(p).factor()) < L.degree()]
    if typ != 'bigRM':
        assert all(len(M.ideal(p).factor()) < M.degree() for p in prime_list)
    return H, prime_list

# Related to Manin Trick
def quo_rem_quadratic(a,b): # DEBUG
    if b == 0:
        raise ZeroDivisionError
    K = a.parent()
    (w,) = K.maximal_order().ring_generators()
    q0 = ((a/b)[0]).floor() + ((a/b)[1]).floor()*K.gen()
    for e0, e1 in product([0,1], repeat=2):
        q = q0 + e0 + e1*w
        r = a - q*b
        if r.norm().abs() < b.norm().abs():
            assert q.is_integral() and r.is_integral()
            return q, r
    raise RuntimeError

def quo_rem_quadratic_old(a,b):
    if b == 0:
        raise ZeroDivisionError
    K = a.parent()
    q = ((a/b)[0]).round('away') + ((a/b)[1]).round('away')*K.gen()
    r = a - q*b
    assert r.norm() < b.norm()
    return (q, r)

def xgcd_quadratic(a, b):
    F = a.parent()
    v1 = vector([F(1),F(0)])
    v2 = vector([F(0),F(1)])
    if b == 0:
        return (a, v1[0], v1[1])
    if a == 0:
        return (b, v2[0], v2[1])
    a0 = a
    b0 = b
    q, _ = quo_rem_quadratic(a0,b0)
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
        q, _ = quo_rem_quadratic(a0,b0)
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
    return p, q

# given the convergents [p_i] and [q_j] of a continued fraction for alpha, computes the matrices that Mi such that {Infinity, alpha} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def compute_Ms(p, q):
    Ms = []
    m = Matrix(2,2,[-p[0], 1, -q[0], 0]) # j = 0
    m.set_immutable()
    Ms.append(m)
    for j in range(1, len(p)):
        m = Matrix(2,2,[(-1)**(j-1)*p[j],p[j-1],(-1)**(j-1)*q[j],q[j-1]])
        m.set_immutable()
        Ms.append(m)
    assert all([A.det() == 1 for A in Ms])
    return Ms


# given a, c in Z[i] returns the matrices that Mi such that {Infinity, a/c} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def matrices_for_unimodular_path(a,c):
    if c == 0:
        return []
    return compute_Ms(*compute_convergents(continued_fraction(a,c)))


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
# the path g·∞ -> g·0
def vectors_in_lattice(F, p, n, g=1):
    g = MatrixSpace(F,2,2)(g)
    res = []
    r, s = act_flt(g, Infinity), act_flt(g, F(0))
    for mac in range(1,n+1):
        RHS = n - mac
        betas = all_elements_of_norm(F, RHS)
        for mc in divisors(mac):
            a = mac // mc
            for beta in betas:
                tmp = Matrix([[beta, mc], [a, -beta.conjugate()]])
                gm = g * tmp * g.adjugate().apply_map(lambda x:x.conjugate())
                if (gm[0,0] / p).is_integral() and ZZ(gm[1,0]) % p == 0 and ZZ(gm[0,1]) % p == 0:
                    continue
                aa = tmp[1,0]
                assert aa.trace() == 2 * aa
                aa = ZZ(aa.trace() / 2)
                assert aa > 0
                if aa % 4 == 1:
                    new_mat = gm
                    sgn = ZZ(1) if ZZ(new_mat[1,0].trace() / 2) % 4 == 1  else -ZZ(1)
                elif aa % 4 == 3:
                    new_mat = -gm
                    sgn = -ZZ(1) if ZZ(new_mat[1,0].trace() / 2) % 4 == 1 else ZZ(1)
                else:
                    continue
                alphaconj, mb, c, _ = gm.list()
                alpha = alphaconj.conjugate()
                sqradius = (alpha.norm() + mb*c) / (c*c)
                center = alpha / c
                center = center.conjugate()
                qf = lambda x : 1 if x == Infinity else (x-center).norm() - sqradius                    
                assert qf(r) * qf(s) < 0
                res.extend([(new_mat, sgn)])
    return res

def mul_inv(num, FF):
    try:
        Rp = FF.parent().base().base()
        a0inv = ~Rp(FF(0)(0))
    except AttributeError:
        assert FF == 1
        a0inv = 1
    y1 = 0
    y = a0inv
    while y != y1:
        y1 = y
        y = y1 * (2 - y1 * FF)
    return num * y

@cached_function
def ApplySingle(A, i, z, M, check=True):
    Rp = z.parent().base_ring()
    if not Rp.is_finite():
        Rp = Rp.base_ring()
    ((p,_),) = ZZ(Rp.cardinality()).factor()
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
    return ii, list_powers(substx, M)


def Transform(outky):
    global gF
    res = 1
    resinv = 1
    t1 = 0
    t2 = 0
    for inky, s1, s2, sgn in G.input_list[outky]:
        assert sgn == 1 or sgn == -1
        f = gF[inky]
        t = cputime()
        newres = sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] for fi, i in zip(f.coefficients(), f.exponents()))
        t1 += cputime(t)
        t = cputime()
        if sgn > 0:
            res *= newres
        else:
            resinv *= newres
        t2 += cputime(t)
    t = cputime()
    ans = mul_inv(res, resinv)
    t3 = cputime(t)
    return ans, (t1,t2,t3)

def degrees(res):
    deg = lambda x : (x.degree() if hasattr(x, 'degree') else 0)
    dumax = [deg(ff) for ky, ff in res.items()]
    dvmax = [max([deg(o) for o in ff.coefficients()]) for ky, ff in res.items()]
    return dumax + dvmax

def fixed_point(g, phi):
    a, b, c, d = g.list()
    K = g.parent()
    Kp = phi.codomain()
    x = Kp['x'].gen()
    f = phi(c)*x**2 + phi(d-a)*x - phi(b)
    p = Kp.prime()
    L = Qq(p**2, Kp.precision_cap(),names='b')
    return solve_quadratic(f.change_ring(phi).change_ring(L), L, return_all=True)




def RMCEval(J, D, cycle_type, prec, alpha=None, n=1, return_class_number=False):
    try:
        ncpus = ncpus
    except NameError:
        ncpus = cpu_count()
    if cycle_type == 'smallCM':
        A, tau0, hE = J.smallCMcycle(D)
    elif cycle_type == 'smallRM':
        A, tau0, hE = J.smallRMcycle(D)
    else:
        if cycle_type != 'bigRM':
            raise NotImplementedError('Cycle type should be either "smallCM" or "smallRM" or "bigRM"')
        A, tau0, hE = J.bigRMcycle(D, alpha=alpha, n=n, version='old')

    if any(t.trace() == 2 * t for t in tau0):
        raise ValueError(f'Tuple {(D, alpha) = } is not admissible for {cycle_type} cycle: the resulting tau is not in Hp.')

    ans = J.evaluate(Infinity, A[0,0]/A[1,0], tau0, prec)
    # if parallelize:
    #     res1 = 1
    #     with futures.ProcessPoolExecutor(max_workers=ncpus) as executor:
    #         future_dict = {executor.submit(Eval, act_flt(m.apply_morphism(phi).adjugate(), tau0), prec) : sgn for m, sgn in mlist}
    #         for fut in futures.as_completed(future_dict):
    #             sgn = future_dict[fut]
    #             res1 *= fut.result()**sgn
    # else:
    #     res1 = prod(Eval(act_flt(m.apply_morphism(phi).adjugate(), tau0), prec)**sgn for m, sgn in mlist)
    ans = ans.add_bigoh(prec + ans.valuation())
    if return_class_number:
        return ans, hE
    else:
        return ans

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
    F = phi.domain()
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
        raise RuntimeError('No embedding of MM into Qp^2 found that extends phi')
    return gamma, iMML(tau1), iMML(tau2)


def cocycle(q : int, label : str, M : int, fname=None, parallelize=True):
    global G
    F = QuadraticField(-1, names='r')
    p = ZZ(q)
    print(f'{p = }')
    print(f'{label = }')
    w = F.elements_of_norm(p)[0]
    G = DGLGroup(F, w, 2, parallelize=parallelize)
    if fname is None:
        fname = f'L0Jtuple_{p}_{label}_{M}.sobj'
    print('Entering Level1')
    F1 = Cocycle(label, M)
    F1.initialize_input_list()
    ### Next
    print('Next...')
    t = walltime()
    F2 = F1.next(parallelize=parallelize)
    print(f'..done next (in {walltime(t)} seconds)')
    ### RMC
    t = walltime()
    F2 = F2.RMC(parallelize=parallelize)
    del F2.datalist
    del F2.G.input_list
    print(f'Saving to {fname}...')
    save(F2, fname)
    print(f'Finished in {walltime(t)} seconds')
    return # F2, F1, G

def recognize(q : int, label : str, D, cycle_type : str, M : int, fname=None, timeout=20, outfile=None, logfile=None, max_degree=16):
    Jtau, hE = evaluate(q, label, D, cycle_type, M, fname)
    if isinstance(D,tuple):
        D, n = D
    else:
        n = 1
    success = False
    for deg in [2**i for i in range(1, 10) if 2**i <= max_degree]:
        with ThreadingTimeout(timeout) as to_ctx_mgr:
            ans = recognize_DGL_algdep(Jtau, deg, outfile=None)
        if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
            fwrite('Someone got impatient...', logfile)
            ans = (None, None, None)
            break
        if ans[0] is not None:
            success = True
            fwrite(f'Computed Jtau for {D = } {n = } {hE = }...', outfile)
            fwrite('SUCCESS: ' + str(ans[0].add_bigoh(10)) + str(ans[1:]), outfile)
            break
    if not success:
        fwrite(f'Computed Jtau for {D = } {n = } {hE = }... (not recognized)', outfile)
    if __name__ != '__main__':
        return Jtau, ans

def evaluate(q : int, label : str, D, cycle_type : str, M : int, fname=None, J=None):
    p = q

    if fname is None:
        if not isinstance(label, str):
            label = label_from_functional(label)
        fname = f'L0Jtuple_{p}_{label}_{M}.sobj'
    if isinstance(D,tuple):
        D, n = D
    else:
        n = 1
    if J is None:
        J = load(fname)
    Jtau, hE = RMCEval(J, D, cycle_type, M, alpha=None, n=n, return_class_number=True)
    if __name__ == '__main__':
        print(Jtau)
    return Jtau, hE

def functional(p, N = 20, MFs = None):
    # We only consider d with are not a norm from Z[i]
    sum_of_squares = lambda n : \
        all(e % 2 == 0 for p, e in ZZ(n).factor() if p % 4 == 3)
    valid_ds = [i for i in range(1, N) if not sum_of_squares(i)]

    if MFs is None:
        MFs = ModularForms(4*p).gens()
    A = Matrix([[QQ(g[o]) for o in valid_ds] for g in MFs])
    for i in range(A.nrows()):
        l = lcm([QQ(o).denominator() for o in A.row(i)])
        A.rescale_row(i, l)
    A = A.change_ring(ZZ)
    # print(A)
    # vectors in the kernel correspond to functionals that vanish on g and on the Eisenstein series
    # the first position of the vector in the kernel corresponds to a_1, and so on
    L = IntegralLattice(ZZ(len(valid_ds))).sublattice(A.right_kernel().basis())
    length_cap = 8
    all_vectors = sum(L.short_vectors(length_cap),[])
    ans = []
    ans_expanded = []
    W = L.submodule(ans)
    i = 0
    while W.rank() < L.rank():
        verbose(W.rank(),L.rank())
        try:
            while L.submodule(ans + [all_vectors[i]]).rank() == W.rank():
                i += 1
        except IndexError:
            length_cap *= QQ(3)/2
            length_cap = ZZ(length_cap.ceil())
            all_vectors = sum(L.short_vectors(length_cap),[])
        try:
            v = all_vectors[i]
        except IndexError:
            continue
        verbose(f'Adding {v = }')
        ans.append(v)
        vexp = [0 for _ in range(1,N)]
        for val, idx in zip(v, valid_ds):
            vexp[idx-1] = val
        ans_expanded.append(vexp)
        W = L.submodule(ans)
    return [label_from_functional(o) for o in ans_expanded]

# Act on power series
def act_pseries(cocycle, h, Phi, level=None, lift=False):
    S = cocycle.Ruv
    R = S.base_ring()
    z1 = R.gen()
    z2 = S.gen()
    p = cocycle.p
    prec = cocycle.prec
    Ruv = cocycle.Ruv
    h.set_immutable()
    mconj = h.apply_map(lambda x : x.trace() - x)
    A = h.apply_morphism(cocycle.phi)
    Aconj = mconj.apply_morphism(cocycle.phi)
    tmp_dict = {}
    for i in range(p+1):
        tmp_dict[i, 0] = ApplySingle(A, i, z1, prec, check=False)
        tmp_dict[i, 1] = ApplySingle(Aconj, i, z2, prec, check=False)

    transformation_list = {}
    for inky0 in range(p+1):
        outky0, s1 = tmp_dict[inky0, 0]
        if s1 is None:
            continue
        for inky1 in range(p+1):
            inky = (inky0, inky1)
            outky1, s2 = tmp_dict[inky1, 1]
            if s2 is None:
                continue
            outky = (outky0, outky1)
            try:
                transformation_list[outky].append((inky, s1, s2))
            except KeyError:
                transformation_list[outky] = [(inky, s1, s2)]
    ans = {}
    for outky in transformation_list:
        res = Ruv(1)
        for inky, s1, s2 in transformation_list[outky]:
            f = Phi[inky]
            newres = sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] \
                for fi, i in zip(f.coefficients(), f.exponents()))
            res *= newres
        ans[outky] = res
    if lift:
        R0 = Ruv.base_ring()
        psi = R0.hom([ZZ['u'].gen()],base_map=lambda x:x.lift(),check=False)
        ans = {ky : f.change_ring(psi) for ky, f in ans.items()}
    return ans

def get_p_prec(f):
    fsp = f.split('_')
    p = ZZ(fsp[1])
    M = ZZ(fsp[-1].split('.')[0])
    return p, M

def get_label(f):
    fsp = f.split('_')
    return fsp[2] + '_' + fsp[3]

# This command can run all of it
# for file in L0Jtuple_*.sobj;do tmux new-session -d -s `basename $file | sed 's/·//g' | sed 's/\.//g'` "conda run -n sage sage find_all_types.sage $file";done;

def eval_and_recognize(fname, typ = None, Dmin=1, Dmax=1000, outdir='outfiles', log='output.log', prime_bound=600, J=None):
    logfile = outdir + '/' + log
    if typ is None or typ == 'all':
        cycle_types = ['smallCM', 'smallRM', 'bigRM']
    elif typ == 'small':
        cycle_types = ['smallCM', 'smallRM']
    else:
        cycle_types = [typ]
    F = QuadraticField(-1, names='r')

    Dvalues = [D for D in srange(ZZ(Dmin),ZZ(Dmax)) if D.is_fundamental_discriminant()]
    nvalues = [1]

    fwrite(f'Doing cocycle {fname}...', logfile)
    p, M = get_p_prec(fname)
    if J is None:
        J = load(fname)
    label = get_label(fname)
    for D in Dvalues:
        for cycle_type in cycle_types:
            outfile = outdir + '/' + f'points_{cycle_type}_{p}_{label}_{M}.txt'
            nvalues = [1]
            for n in nvalues:
                try:
                    try:
                        fwrite(f'Computing {fname} {D},{n},{cycle_type}...', logfile)
                        Jtau0, hE = RMCEval(J, D, cycle_type, 10, alpha=None, n=n, return_class_number=True)
                        if Jtau0 == 1:
                            fwrite(f'Computed Jtau = 1 for {D = } {n = } {hE = }...', outfile)
                            continue
                        Jtau0, hE = RMCEval(J, D, cycle_type, M, alpha=None, n=n, return_class_number=True)
                        Cp = Jtau0.parent()
                    except (ValueError, NotImplementedError, TypeError, RuntimeError, PrecisionError, SignalError, KeyboardInterrupt, ModuleNotFoundError) as e:
                        if 'no elements of norm -1' not in str(e):
                            fwrite(f'Skipping {fname} {D},{n},{cycle_type}...({str(e)})', logfile)
                        continue
                    fwrite(f'Computed Jtau for {fname} {D} {n} {cycle_type}.', logfile)
                    if Jtau0 == 1:
                        fwrite(f'Computed Jtau = 1 for {D = } {n = } {hE = }...', outfile)
                        continue
                    for Jtau, char in [(Jtau0.norm(), 'triv'), (Jtau0**2 / Jtau0.norm(), 'conj')]:
                        success = False
                        try:
                            with stopit.ThreadingTimeout(120) as to_ctx_mgr:
                                H, prime_list = get_predicted_field_and_prime_list(F, D, n, cycle_type, char, names='z', prime_bound=prime_bound)
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                raise ExpatError
                        except ExpatError:
                            fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }... ({char} not recognized since H could not be computed)', outfile)
                            H = None
                            prime_list = None
                        for deg in [2,4,8]:
                            with stopit.ThreadingTimeout(30) as to_ctx_mgr:
                                ans = recognize_DGL_algdep(Jtau, deg, outfile=None, roots_of_unity=[1])
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                fwrite('Someone got impatient...', logfile)
                                ans = None
                                break
                            if ans[0] is not None:
                                success = True
                                fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }...', outfile)
                                support = [q for q, _ in ans[-1] if q in ZZ]
                                dta = list(ans[1:])
                                dta.append('?' if H is None else H.defining_polynomial())
                                msg = f'SUCCESS {char} with algdep: ' + str(ans[0].unit_part().add_bigoh(10)) + str(tuple(dta))
                                if prime_list is None:
                                    msg += ' warning: primes not checked!'
                                else:
                                    if any(q not in prime_list for q in support):
                                        msg += ' warning: ' + str([q for q in support if q not in prime_list])
                                fwrite(msg, outfile)
                                break
                        if not success and prime_list is not None and ('big' not in cycle_type or H.degree() <= 16):
                            with stopit.ThreadingTimeout(3600) as to_ctx_mgr:
                                ans = recognize_DGL_lindep(Jtau, H, prime_list = prime_list, outfile=None, recurse_subfields=True, degree_bound=8, algorithm='pari')
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                fwrite(f'Someone got impatient...  {cycle_type = } {p = } {label = } {D = }, {n = }', logfile)
                                ans = None
                            if ans is not None:
                                success = True
                                ffpoly, rec, d, clist_ans, field = ans
                                fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }...', outfile)
                                fwrite(f'SUCCESS {char} with lindep: ' + str(Jtau.unit_part().add_bigoh(10)) + f'({ffpoly}, {field}, {rec}, {d}, 0, 1, {clist_ans}, {H.defining_polynomial()})', outfile)

                        if not success:
                            fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }... ({char} not recognized)', outfile)
                    fwrite('...', outfile)
                except KeyboardInterrupt:
                    fwrite(f'WARNING! Keyboard interrupt so skipping {cycle_type = } {p = } {label = } {D = }, {n = }', logfile)
                    sleep(1)
                except Exception as e:
                    fwrite(f'WARNING! Unhandled exception so skipping {cycle_type = } {p = } {label = } {D = }, {n = } : {str(e)}', logfile)

if __name__ == '__main__':
  Fire({'cocycle': cocycle,'evaluate': evaluate,'functional':functional, 'recognize':recognize,
         'eval_and_recognize': eval_and_recognize})

