global L0, J, M
label = '37_6'
cycle_type = 'smallCM'
M = 300
fname = f'L0Jtuple_{label}_{M}.sobj'
ncpus = 64
parallelize = True
outfile = 'algebraicsmallCM2.txt'
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

Dvalues = [D for D in srange(10000) if D.is_fundamental_discriminant()]
n = 1
for D in Dvalues:
    prime_list = prime_range(2*D)
    try:
        Jtau, hE = RMCEval(D, cycle_type, M, return_class_number=True)
    except (ValueError, NotImplementedError, TypeError, RuntimeError, PrecisionError) as e:
        print(f'Skipping {D},{n}...({str(e)})')
        continue
    fwrite(f'Computed Jtau for {D = } {hE = }...', outfile)
    fwrite('Jtau = %s'%Jtau.add_bigoh(10), outfile)
    E.<zrel> = F.extension(x^2 - D)
    Eabs.<z> = E.absolute_field()
    try:
        _ = Eabs._pari_init_()
    except:
        pass
    H = Eabs.extension(Eabs['x'](Eabs._pari_nf.bnfinit().bnrinit(1).bnrclassfield()[0].sage(locals={'x':x,'y':z})),names='c0').absolute_field(names='z')
    plist = [p for p in prime_list if len(QuadraticField(-D,names='tt').ideal(p).factor()) == 1]
    ans = recognize_DGL_lindep(Jtau, H, prime_list=plist, outfile=outfile, algorithm='pari', recurse_subfields=True)
    if ans is not None:
        fwrite(str(ans[:-2]),outfile)
    fwrite('...', outfile)

