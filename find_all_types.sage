from glob import glob
from sage.rings.padics.precision_error import PrecisionError
load('orthogonal.sage')
cocycle_fnames = sorted(glob('L0Jtuple*.sobj'))

def get_p_prec(f):
    fsp = f.split('_')
    p = ZZ(fsp[1])
    M = ZZ(fsp[-1].split('.')[0])
    return p, M

def get_label(f):
    fsp = f.split('_')
    return fsp[2] + '_' + fsp[3]


global L0, J, M
cycle_types = ['smallCM', 'smallRM', 'bigRM']
# fname = f'L0Jtuple_{label}_{M}.sobj'
ncpus = 64
parallelize = True
F = QuadraticField(-1, names='r')
logfile='outfiles/output.log'
# outfile = {ct : f'algebraic_points_{ct}.txt' for ct in cycle_types}
###

Dvalues = [D for D in srange(1000) if D.is_fundamental_discriminant()]
nvalues = [1]

for fname in cocycle_fnames:
    fwrite(f'Doing cocycle {fname}...', logfile)
    p, M = get_p_prec(fname)
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
    label = get_label(fname)
    for D in Dvalues:
        for cycle_type in cycle_types:
            if cycle_type == 'bigRM':
                nvalues = [n for n in range(-7,7) if ZZ(n).is_squarefree()]
            else:
                nvalues = [1]
            for n in nvalues:
                outfile = f'outfiles/points_{cycle_type}_{p}_{label}_{M}.txt'
                try:
                    Jtau, hE = RMCEval(D, cycle_type, M, n, return_class_number=True)
                except (ValueError, NotImplementedError, TypeError, RuntimeError, PrecisionError) as e:
                    fwrite(f'Skipping {D},{n},{cycle_type}...({str(e)})', logfile)
                    continue
                if Jtau == 1:
                    fwrite(f'Computed Jtau = 1 for {D = } {n = } {hE = }...', outfile)
                    continue
                success = False
                for deg in [2,4,8,16]:
                    ans = recognize_DGL_algdep(Jtau,deg, outfile=None)
                    if ans[0] is not None:
                        success = True
                        fwrite(f'Computed Jtau for {D = } {n = } {hE = }...', outfile)
                        fwrite('SUCCESS: ' + str(ans[0].add_bigoh(10)) + str(ans[1:]), outfile)
                        break
                if not success:
                    fwrite(f'Computed Jtau for {D = } {n = } {hE = }... (not recognized)', outfile)
                fwrite('...', outfile)

