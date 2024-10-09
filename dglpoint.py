from util import *
from fire import Fire
import stopit

from sage.all import QuadraticField, load, Integer, RealNumber, parallel, GF, walltime, vector, save, fundamental_discriminant, Qq, srange, ModularForms, IntegralLattice
# Calculate DGL periods

parallelize = True
F = QuadraticField(-1, names='r')

def label_from_functional(func, sep='·'):
    pos = ''
    neg = ''
    for i, o in enumerate(func):
        if o > 0:
            pos = pos + o * (sep + str(i+1))
        elif o < 0:
            neg = neg + (-o) * (sep + str(i+1))
    return pos[1:] + '_' + neg[1:]

def functional_from_label(label, sep='·'):
    func = {}
    pos, neg = label.split('_')
    for i in pos.split(sep):
        try:
            func[int(i)-1] += 1
        except KeyError:
            func[int(i)-1] = 1
    for i in neg.split(sep):
        try:
            func[int(i)-1] -= 1
        except KeyError:
            func[int(i)-1] = -1
    lfunc = max(func.keys()) + 1
    ans = [0 for _ in range(lfunc)]
    for ky, val in func.items():
        ans[ky] = val
    return tuple(ans)

def cocycle(q : int, label : str, M : int, fname=None):
    global L0, p, J, F, parallelize, phi, Fp, Ruv, Rp, map_poly, inv_map_poly
    p = ZZ(q)
    print(f'{p = }')
    print(f'{label = }')
    w = F.elements_of_norm(p)[0]
    Fp = Qp(p,2*M, type='floating-point')
    Rp = Zmod(p**M)
    phi = F.hom([F.gen().minpoly().roots(Fp)[0][0]],check=False)
    if phi(w).valuation() == 0:
        phi = F.hom([F.gen().minpoly().roots(Fp)[1][0]],check=False)
    Ruv = PolynomialRing(PolynomialRing(Rp,'u'),'v')
    inv_map_poly = Ruv.flattening_morphism()
    map_poly = inv_map_poly.inverse()
    load('orthogonal.sage')

    if fname is None:
        fname = f'L0Jtuple_{p}_{label}_{M}.sobj'
    if isinstance(label, str):
        label = functional_from_label(label)
    L0, V = initial_seed(label, p)
    calculate_Tp_matrices(w, M)
    print('Entering Level1')
    t = walltime()
    F1, dg = Level1(V, M)
    print(f'..done Level1 (in {walltime(t)} seconds)')
    ### Next
    print('Next...')
    t = walltime()
    F2 = Next(F1)
    print(f'..done next (in {walltime(t)} seconds)')
    ### RMC
    t = walltime()
    J = RMC(F2, M)
    print(f'Saving to {fname}...')
    save((L0, J), fname)
    print(f'Finished in {walltime(t)} seconds')

def recognize(q : int, label : str, D, cycle_type : str, M : int, fname=None, timeout=20, outfile=None, logfile=None, max_degree=16):
    global L0, J, p, F, parallelize, phi, Fp, Ruv, Rp, map_poly, inv_map_poly
    p = q
    if fname is None:
        if not isinstance(label, str):
            label = label_from_functional(label)
        fname = f'L0Jtuple_{p}_{label}_{M}.sobj'
    if isinstance(D,tuple):
        D, n = D
    else:
        n = None
    w = F.elements_of_norm(p)[0]
    Fp = Qp(p,2*M, type='floating-point')
    Rp = Zmod(p**M)
    phi = F.hom([F.gen().minpoly().roots(Fp)[0][0]],check=False)
    if phi(w).valuation() == 0:
        phi = F.hom([F.gen().minpoly().roots(Fp)[1][0]],check=False)
    Ruv = PolynomialRing(PolynomialRing(Rp,'u'),'v')
    inv_map_poly = Ruv.flattening_morphism()
    map_poly = inv_map_poly.inverse()
    load('orthogonal.sage')
    L0, J = load(fname)
    Jtau, hE = RMCEval(D, cycle_type, M, n, return_class_number=True)
    success = False
    for deg in [2**i for i in range(1, 10) if 2**i <= max_degree]:
        with stopit.ThreadingTimeout(timeout) as to_ctx_mgr:
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

def evaluate(q : int, label : str, D, cycle_type : str, M : int, fname=None):
    global L0, J, p, F, parallelize, phi, Fp, Ruv, Rp, map_poly, inv_map_poly
    p = q
    if fname is None:
        if not isinstance(label, str):
            label = label_from_functional(label)
        fname = f'L0Jtuple_{p}_{label}_{M}.sobj'
    if isinstance(D,tuple):
        D, n = D
    else:
        n = 1
    w = F.elements_of_norm(p)[0]
    Fp = Qp(p,2*M, type='floating-point')
    Rp = Zmod(p**M)
    phi = F.hom([F.gen().minpoly().roots(Fp)[0][0]],check=False)
    if phi(w).valuation() == 0:
        phi = F.hom([F.gen().minpoly().roots(Fp)[1][0]],check=False)
    Ruv = PolynomialRing(PolynomialRing(Rp,'u'),'v')
    inv_map_poly = Ruv.flattening_morphism()
    map_poly = inv_map_poly.inverse()
    load('orthogonal.sage')
    L0, J = load(fname)
    x = RMCEval(D, cycle_type, M, n)
    if __name__ == '__main__':
        print(x)
    else:
        return x

def sum_of_squares(n):
    return all(e % 2 == 0 for p, e in ZZ(n).factor() \
               if p % 4 == 3)

def functional(p, N = 20):
    # We only consider d with are not a norm from Z[i]
    valid_ds = [i for i in range(1, N) if not sum_of_squares(i)]

    MFs = ModularForms(4*p).gens()
    A = Matrix([[QQ(g[o]) for o in valid_ds] for g in MFs])
    for i in range(A.nrows()):
        l = lcm([QQ(o).denominator() for o in A.row(i)])
        A.rescale_row(i, l)
    A = A.change_ring(ZZ)

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


def check_cocycle(L0, J, M, version=2):
    if version == 2:
        mlist = [Matrix(QQ,2,2,1), Matrix(QQ,2,2,[0,-1,1,0])]
    else:
        U = Matrix(QQ,2,2,[0,1,-1,1])
        mlist = [Matrix(QQ,2,2,1), U, U**2]
    apply_single_dict = {}
    S = Ruv
    R = S.base_ring()
    z1 = R.gen()
    z2 = S.gen()
    for m in mlist:
        m.set_immutable()
        A = m.apply_morphism(phi)
        for i in range(p+1):
            ii, subst1 = ApplySingle(A, i, z1, M, check=False)
            jj, subst2 = ApplySingle(A, i, z2, M, check=True)
            apply_single_dict[m,i,0] = (ii, list_powers(subst1,M))
            apply_single_dict[m,i,1] = (jj, list_powers(subst2,M))

    input_list = {(i,j) : list() for i in range(p+1) for j in range(p+1)}
    for m in mlist:
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
                input_list[outky].append((inky, s1, s2))
    res = Ruv(1)
    for inky, s1, s2 in input_list[outky]:
        f = gF[inky]
        res *= sum(sum(aij * s1[j] for aij, j in zip(fi.coefficients(), fi.exponents())) * s2[i] for fi, i in zip(f.coefficients(), f.exponents()))
    return res

if __name__ == '__main__':
  Fire({'cocycle': cocycle,'evaluate': evaluate,'functional':functional, 'recognize':recognize})

