from util import *
# Calculate DGL periods
ncpus = 64
parallelize = True
p = 5
D = -1
x = QQ['x'].gen()
F.<r> = QuadraticField(D)
Fp = Qp(p,2*M, type='floating-point')
Rp = Zmod(p**M)
eps = F.unit_group().gens()[0]
units = [F(eps**i) for i in range(eps.multiplicative_order())]
phi = F.hom([F.gen().minpoly().roots(Fp)[0][0]],check=False)
if phi(1+2*r).valuation() != 0:
    phi = F.hom([F.gen().minpoly().roots(Fp)[1][0]],check=False)
Ruv = PolynomialRing(PolynomialRing(Rp,'u'),'v')
inv_map_poly = Ruv.flattening_morphism()
map_poly = inv_map_poly.inverse()
load('orthogonal.sage')

def cocycle(label : str, p : int):
    global J
    if instance(label, str):
        label = functional_from_label(label)
    L0, V = initial_seed(label, p)
    calculate_Tp_matrices(1+2*r)
    print('Entering Level1')
    t = walltime()
    F1, dg = Level1(V)
    print(f'..done Level1 (in {walltime(t)} seconds)')
    ### Next
    print('Next...')
    t = walltime()
    F2 = Next(F1)
    print(f'..done next (in {walltime(t)} seconds)')
    ### RMC
    t = walltime()
    J = RMC(F2)
    print(f'Finished in {walltime(t)} seconds')
    label = label_from_functional(label)
    save((L0, J), f'L0Jtuple_{label}_{M}.sobj')
    return (L0, J)

def eval(label : str, D : int, cycle_type : str, prec : int):
    global L0, J, M
    M = prec
    L0, J = load(f'L0Jtuple_{label}_{M}.sobj')
    x = RMCEval(D, cycle_type, prec)
    return x

