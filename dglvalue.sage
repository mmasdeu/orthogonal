from util import *
import sys
# Calculate DGL periods

if __name__ == "__main__":
    global L0, J, M
    label = str(sys.argv[1])
    D = int(sys.argv[2])
    cycle_type = str(sys.argv[3])
    M = int(sys.argv[4])
    print(f'label = {label}, D = {D}, cycle_type = {cycle_type}, precision = {M}')
    ncpus = 64
    parallelize = True
    p = 5
    x = QQ['x'].gen()
    F.<r> = QuadraticField(-1)
    Fp = Qp(p,2*M, type='floating-point')
    Rp = Zmod(p**M)
    eps = F.unit_group().gens()[0]
    units = [F(eps**i) for i in range(eps.multiplicative_order())]
    phi = F.hom([F.gen().minpoly().roots(Fp)[0][0]],check=False)
    if phi(1+2*r).valuation() == 0:
        phi = F.hom([F.gen().minpoly().roots(Fp)[1][0]],check=False)
    Ruv = PolynomialRing(PolynomialRing(Rp,'u'),'v')
    inv_map_poly = Ruv.flattening_morphism()
    map_poly = inv_map_poly.inverse()
    load('orthogonal.sage')
    L0, J = load(f'L0Jtuple_{label}_{M}.sobj')
    x = RMCEval(D, cycle_type, M)
    print(x)

