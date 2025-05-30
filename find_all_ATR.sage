global L0, J, M
label = '37_6'
cycle_type = 'bigRM'
M = 300
fname = f'L0Jtuple_{label}_{M}.sobj'
ncpus = 64
parallelize = True
outfile = 'Jlist_37_6.sage'
###

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
L0, J = load(fname)

Dvalues = [D for D in srange(1000) if D.is_fundamental_discriminant()]
nvalues = [n for n in srange(-7,7) if n.is_squarefree()]
fwrite(f'L = Qq({p}**2, {M}, names="b"', outfile)
fwrite('Jlist = [\\', outfile)
for D, n in [(D, n) for D in Dvalues for n in nvalues]:
    try:
        Jtau, hE = RMCEval(D, cycle_type, M, n, return_class_number=True)
    except (ValueError, NotImplementedError, TypeError, RuntimeError) as e:
        print(f'Skipping {D},{n}...({str(e)})')
        continue
    if hE > 4:
        continue
    fwrite(str((D, n, Jtau, hE)) + ',\\', outfile)
fwrite(']', outfile)
