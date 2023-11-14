from util import *
# Calculate DGL periods
# Example discriminants up to 300 [17, 24, 44, 96, 104, 116, 153, 156, 164, 200, 204, 224, 257, 273, 297]
# h = 1: [17, 24, 44, 96, 116, 153, 164, 297]

## Speed
# M = 20, nonparallel: about 30s total
# M = 100, parallel: Level1 170s, Next 248s, first iteration 247s, second iter 201s, third 180s,...
# M = 100, parallel with nonparallel Apply: Level1 174s, Next 260s, first iteration
# M = 300, parallel: Level1 3910s, Next 16332s, first iteration 17594, second 14447,


### M = 100 - antz5
# Finished in 658.7516839504242 seconds
# CPU times: user 7min 28s, sys: 1min 57s, total: 9min 25s
# Wall time: 15min 7s

### M = 300 - antz5
# Entering Level1
# Calculating inverses...
# Looking at as_completed
# ..done Level1 (in 513.5984840393066 seconds)
# Next...
# ..done next (in 3244.2414684295654 seconds)
# Iteration 1
# ..done in 4430.110747098923 seconds.
# Iteration 2
# ..done in 3239.9868392944336 seconds.
# Iteration 3
# ..done in 2680.5017473697662 seconds.
# Iteration 4
# ..done in 2409.334589958191 seconds.
# Iteration 5
# ..done in 2226.1840937137604 seconds.
# Iteration 6
# ..done in 2097.554993391037 seconds.
# Iteration 7
# ..done in 2003.137940645218 seconds.
# Iteration 8
# ..done in 1912.285739660263 seconds.
# Iteration 9
# ..done in 1854.8608467578888 seconds.
# Iteration 10
# ..done in 1777.2119376659393 seconds.
# Iteration 11
# ..done in 1687.341745853424 seconds.
# Iteration 12
# ..done in 1619.5119152069092 seconds.
# Iteration 13
# ..done in 1556.705640554428 seconds.
# Iteration 14
# ..done in 1508.1171824932098 seconds.
# Iteration 15
# ..done in 1474.0365979671478 seconds.
# Iteration 16
# ..done in 1411.3740196228027 seconds.
# Iteration 17
# ..done in 1382.1933391094208 seconds.
# Iteration 18
# ..done in 1351.6515498161316 seconds.
# Iteration 19
# ..done in 1295.6327142715454 seconds.
# Iteration 20
# ..done in 1224.034026145935 seconds.
# Iteration 21
# ..done in 1166.423781633377 seconds.
# Iteration 22
# ..done in 1118.7835175991058 seconds.
# Iteration 23
# ..done in 1113.985188484192 seconds.
# Iteration 24
# ..done in 1129.1612219810486 seconds.
# Iteration 25
# ..done in 1147.2020468711853 seconds.

### Initialize
# M = 50
# parallelize = True
ncpus = 64
# embedding_idx = 1

p = 5
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
MS = calculate_Tp_matrices(1+2*r)

### Find initial seed
L0, V = initial_seed((0, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0), p)
### Level1
def run_all():
    global J
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
    save((L0, J), 'L0Jtuple.sobj')
    x1 = RMCEval(L0, 17*9)

    I = our_sqrt(Qp(p,300, type='floating-point')(-1))

    xnum =  19^2 * (1+I)^8 * (1-4*I)^5 * (5+2*I)^3 * (5+4*I) * (5+6*I)^2
    xden = (1+2*I)^5 * (1-2*I)^6 * (1-6*I)^2 * (9-4*I)^2 * (7-8*I)^2 *  (2-13*I)^2

    print('This should be large: ', (x1.parent()(xnum/xden) - x1).valuation()) # should be large

    return J

### Save and load
# To save:
# save((L0, J), 'L0Jtuple.sobj')
# L0, J = load('L0Jtuple.sobj')

