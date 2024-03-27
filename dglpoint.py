from util import *
from fire import Fire
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
    w = 1 + 2 * F.gen()
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
    print('..saving...')
    label = label_from_functional(label)
    save((L0, J), fname)
    print(f'Finished in {walltime(t)} seconds')

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
    w = 1 + 2 * F.gen()
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

    MM = []
    MFs = ModularForms(4*p).gens()
    for g0 in MFs:
        g = list(g0.q_expansion(N+1))
        l = lcm([QQ(o).denominator() for o in g])
        MM.append([l * o for o in g])
    A = Matrix([[QQ(g[o]) for o in valid_ds] for g in MFs])
    l = lcm([QQ(o).denominator() for o in A.list()])
    A = (l * A).change_ring(ZZ)
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

if __name__ == '__main__':
  Fire({'cocycle': cocycle,'evaluate': evaluate,'functional':functional})
