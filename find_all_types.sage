from glob import glob
from time import sleep
from fire import Fire
import stopit

from sage.rings.padics.precision_error import PrecisionError
from cysignals.signals import SignalError
from builtins import ModuleNotFoundError

load('orthogonal.sage')
# cocycle_fnames = sorted(glob('L0Jtuple*.sobj'))

def get_p_prec(f):
    fsp = f.split('_')
    p = ZZ(fsp[1])
    M = ZZ(fsp[-1].split('.')[0])
    return p, M

def get_label(f):
    fsp = f.split('_')
    return fsp[2] + '_' + fsp[3]

# Ths command can run all of it
# for file in L0Jtuple_*.sobj;do tmux new-session -d -s `basename $file | sed 's/·//g' | sed 's/\.//g'` "conda run -n sage sage find_all_types.sage $file";done;

F.<i> = QuadraticField(-1)
ncpus = 4
parallelize = True

def evaluate_cocycle(fname, typ = None, Dmin=1, Dmax=1000):
    global L0, J, M, F, p, Rp, phi
    if typ is None or typ == 'all':
        cycle_types = ['smallCM', 'smallRM', 'bigRM']
    elif typ == 'small':
        cycle_types = ['smallCM', 'smallRM']
    else:
        cycle_types = [typ]
    F = QuadraticField(-1, names='r')
    logfile='outfiles/output.log'

    Dvalues = [D for D in srange(ZZ(Dmin),ZZ(Dmax)) if D.is_fundamental_discriminant()]
    nvalues = [1]

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
            outfile = f'outfiles/points_{cycle_type}_{p}_{label}_{M}.txt'
            if cycle_type == 'bigRM':
                nvalues = [-7, -6, -5, -3, -3, -1, 1, 2, 3, 5, 6, 7]
            else:
                nvalues = [1]
            for n in nvalues:
                try:
                    try:
                        Jtau0, hE = RMCEval(D, cycle_type, M, n, return_class_number=True)
                        Cp = Jtau0.parent()
                    except (ValueError, NotImplementedError, TypeError, RuntimeError, PrecisionError, SignalError, KeyboardInterrupt, ModuleNotFoundError) as e:
                        fwrite(f'Skipping {D},{n},{cycle_type}...({str(e)})', logfile)
                        continue
                    if Jtau0 == 1:
                        fwrite(f'Computed Jtau = 1 for {D = } {n = } {hE = }...', outfile)
                        continue
                    for Jtau, char in [(Cp(Jtau0.norm()), 'triv'), (Cp(Jtau0**2 / Jtau0.norm()), 'conj')]:
                        success = False
                        for deg in [2,4,8]:
                            with stopit.ThreadingTimeout(20) as to_ctx_mgr:
                                ans = recognize_DGL_algdep(Jtau, deg, outfile=None, roots_of_unity=[1])
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                fwrite('Someone got impatient...', logfile)
                                ans = None
                                break
                            if ans[0] is not None:
                                success = True
                                fwrite(f'Computed Jtau for {D = } {n = } {hE = }...', outfile)
                                fwrite(f'SUCCESS {char} with algdep: ' + str(ans[0].unit_part().add_bigoh(10)) + str(ans[1:]), outfile)
                                break
                        if not success and cycle_type in ['smallRM', 'smallCM']:
                            H, prime_list = get_predicted_field_and_prime_list(F, D, n, cycle_type, char, names='z')
                            with stopit.ThreadingTimeout(1800) as to_ctx_mgr:
                                ans = recognize_DGL_lindep(Jtau, H, prime_list = prime_list, outfile=None, recurse_subfields=True)
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                fwrite(f'Someone got impatient...  {cycle_type = } {p = } {label = } {D = }, {n = }', logfile)
                                ans = None
                            if ans is not None:
                                success = True
                                fwrite(f'Computed Jtau for {D = } {n = } {hE = }...', outfile)
                                ffpoly, rec, clist_ans, field = ans
                                fwrite(f'SUCCESS {char} with lindep: ' + str(Jtau.unit_part().add_bigoh(10)) + f'({ffpoly}, {field}, {rec}, 0, 1, {clist_ans})', outfile)

                        if not success:
                            fwrite(f'Computed Jtau for {D = } {n = } {hE = }... ({char} not recognized)', outfile)
                    fwrite('...', outfile)
                except KeyboardInterrupt:
                    fwrite(f'WARNING! Keyboard interrupt so skipping {cycle_type = } {p = } {label = } {D = }, {n = }', logfile)
                    sleep(1)
                except Exception as e:
                    fwrite(f'WARNING! Unhandled exception so skipping {cycle_type = } {p = } {label = } {D = }, {n = } : {str(e)}', logfile)



if __name__ == '__main__':
  Fire(evaluate_cocycle)

