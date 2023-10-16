from util import *
# Calculate DGL periods

## Speed
# M = 20, nonparallel: about 30s total
# M = 100, parallel: Level1 170s, Next 248s, first iteration 247s, second iter 201s, third 180s,...
# M = 100, parallel with nonparallel Apply: Level1 174s, Next 260s, first iteration
# M = 300, parallel: Level1 3910s, Next 16332s, first iteration 17594, second 14447,

### Initialize
p = 5
M = 20
parallelize = True
ncpus = 2
embedding_idx = 1

D = -1
x = QQ['x'].gen()
F.<r> = QuadraticField(D)
Fp = Qp(p,2*M, type='floating-point')
Rp = Zmod(p**M)
eps = F.unit_group().gens()[0]
units = [F(eps**i) for i in range(eps.multiplicative_order())]
phi = F.hom([F.gen().minpoly().roots(Fp)[embedding_idx][0]],check=False)
print(phi)
Ruv = PolynomialRing(PolynomialRing(Rp,'u'),'v')
inv_map_poly = Ruv.flattening_morphism()
map_poly = inv_map_poly.inverse()
load('orthogonal.sage')
MS, input_list = calculate_Tp_matrices(1+2*r)

### Find initial seed
L0, V = initial_seed((0, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0), p)
### Level1
def run_all():
    print('Entering Level1')
    t = walltime()
    F1, dg = Level1(V)
    print(f'..done Level1 (in {walltime(t)} seconds)')
    ### Next
    print('Next...')
    t = walltime()
    F1 = Next(F1)
    print(f'..done next (in {walltime(t)} seconds)')
    ### RMC
    t = walltime()
    J = RMC(F1)
    print(f'Finished in {walltime(t)} seconds')
    A = compute_gamma_tau(17*9)
    x1, x10, x11 = RMCEval(L0, J, A)

    I = our_sqrt(Qp(p,160, type='floating-point')(-1))

    xnum =  19^2 * (1+I)^8 * (1-4*I)^5 * (5+2*I)^3 * (5+4*I) * (5+6*I)^2
    xden = (1+2*I)^5 * (1-2*I)^6 * (1-6*I)^2 * (9-4*I)^2 * (7-8*I)^2 *  (2-13*I)^2

    print('This should be large: ', (x1.parent()(xnum/xden) - x1).valuation()) # should be large

    return J

### Save and load
# To save:
# save((L0, J), 'L0Jtuple.sobj')
# L0, J = load('L0Jtuple.sobj')

