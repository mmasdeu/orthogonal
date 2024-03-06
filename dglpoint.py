# from sage.all_cmdline import *
from util import *
from fire import Fire
from sage.all import QuadraticField, load, Integer, RealNumber, parallel, GF, walltime, vector, save, fundamental_discriminant, Qq, srange
# Calculate DGL periods

ncpus = 128
parallelize = True
p = 5
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

def cocycle(label : str, M : int, fname=None):
    global L0, J, p, F, ncpus, parallelize, phi, Fp, Ruv, Rp, map_poly, inv_map_poly
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
        fname = f'L0Jtuple_{label}_{M}.sobj'
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
    save((L0, J), f'L0Jtuple_{label}_{M}.sobj')
    print(f'Finished in {walltime(t)} seconds')

def evaluate(label : str, D, cycle_type : str, M : int, fname=None):
    global L0, J, p, F, ncpus, parallelize, phi, Fp, Ruv, Rp, map_poly, inv_map_poly
    if fname is None:
        if not isinstance(label, str):
            label = label_from_functional(label)
        fname = f'L0Jtuple_{label}_{M}.sobj'
    L0, J = load(fname)
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
    x = RMCEval(D, cycle_type, M, n)
    if __name__ == '__main__':
        print(x)
    else:
        return x

if __name__ == '__main__':
  Fire({'cocycle': cocycle,'evaluate': evaluate})
