from util import *
# Calculate DGL periods

### Initialize
p = 5
M = 20
parallelize = False
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

MS = calculate_Tp_matrices(1+2*r)
functional_list = vanishing_functional(10)


###
zero_discs = [32, 88, 92, 152, 153, 168, 172, 188, 232, 248, 252, 268, 272, 308, 312, 328, 332, 348, 368]
found_discs = {D : [] for D in zero_discs}
for i, v in enumerate(functional_list):
    L0, V = initial_seed(v, p)
    F1, dg = Level1(V)
    F1 = Next(F1, MS)
    J = RMC(F1)
    for D in zero_discs:
        A = compute_gamma_tau(D)
        try:
            x1, x10, x11 = RMCEval(L0, J, A)
        except PrecisionError:
            print(f'Discriminant {D} cannot be evaluated for cocycle {i}')
            continue
        if (x1 - 1).valuation() < 2:
            print(f'Discriminant {D} gives nonzero for cocycle {i}')
            found_discs[D].append(v)


### Find initial seed
L0, V = initial_seed((0, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0), p)

### Level1
F1, dg = Level1(V)
### Next
F1 = Next(F1, MS)
### RMC
%time J = RMC(F1)

A = compute_gamma_tau(17*9)
x1, x10, x11 = RMCEval(L0, J, A)

xg = (19^2 * (1+r)^8 * (1-4*r)^5 * (5+2*r)^3 * (5+4*r) * (5+6*r)^2) / \
    ((1+2*r)^5 * (1-2*r)^6 * (1-6*r)^2 * (9-4*r)^2 * (7-8*r)^2 *  (2-13*r)^2)

xgp = (19^2 * (1+I)^8 * (1-4*I)^5 * (5+2*I)^3 * (5+4*I) * (5+6*I)^2) / \
    ((1+2*I)^5 * (1-2*I)^6 * (1-6*I)^2 * (9-4*I)^2 * (7-8*I)^2 *  (2-13*I)^2)



prime_list = [ell for ell in prime_range(200) if len(QuadraticField(-17,'y').ideal(ell).factor()) == 1]
# [ell for ell, e in xg.norm().factor()]

I = our_sqrt(Qp(p,160, type='floating-point')(-1))

xnum =  19^2 * (1+I)^8 * (1-4*I)^5 * (5+2*I)^3 * (5+4*I) * (5+6*I)^2
xden = (1+2*I)^5 * (1-2*I)^6 * (1-6*I)^2 * (9-4*I)^2 * (7-8*I)^2 *  (2-13*I)^2

print('This should be large: ', (x1.parent()(xnum/xden) - x1).valuation()) # should be large

Fp2 = Qq(p^2,80,names='y')
x1 = Fp2(xnum/xden)
recognize_DV_lindep(x1, F, prime_list, Cp = Fp2)




# ### Check
# _ = gp.eval('\\r cocycle.gp')

# res1 = {ky : inv_map_poly(val) for ky, val in F1.items()}

# x1, x2 = res1[0,0].parent().gens()
# y1, y2 = PolynomialRing(Zmod(p**M),names='y1,y2').gens()
# psi = y1.parent().hom([x1,x2],codomain=x1.parent(),check=False)


# ###################

# for i, j in J.keys():
#     fij = sage_eval(gp.eval(f'J[2][{i+1},{j+1}]'), locals={'x1':y1, 'x2':y2}).change_ring(psi)
#     gij = J[i,j]
#     if fij != gij:
#         print(f'!!! {i}, {j}')



