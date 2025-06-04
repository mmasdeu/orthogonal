from glob import glob
from time import sleep
from fire import Fire
from xml.parsers.expat import ExpatError
import stopit

from sage.rings.padics.precision_error import PrecisionError
from cysignals.signals import SignalError
from builtins import ModuleNotFoundError

global J

J = None
ncpus = 4
parallelize = True

from dglpoint import *
# cocycle_fnames = sorted(glob('L0Jtuple*.sobj'))

def get_p_prec(f):
    fsp = f.split('_')
    p = ZZ(fsp[1])
    M = ZZ(fsp[-1].split('.')[0])
    return p, M

def get_label(f):
    fsp = f.split('_')
    return fsp[2] + '_' + fsp[3]

# This command can run all of it
# for file in L0Jtuple_*.sobj;do tmux new-session -d -s `basename $file | sed 's/·//g' | sed 's/\.//g'` "conda run -n sage sage find_all_types.sage $file";done;

def evaluate_cocycle(fname, typ = None, Dmin=1, Dmax=1000, outdir='outfiles', log='output.log', prime_bound=600):
    global L0, J, M, F, p, Rp, phi
    logfile = outdir + '/' + log
    if typ is None or typ == 'all':
        cycle_types = ['smallCM', 'smallRM', 'bigRM']
    elif typ == 'small':
        cycle_types = ['smallCM', 'smallRM']
    else:
        cycle_types = [typ]
    F = QuadraticField(-1, names='r')

    Dvalues = [D for D in srange(ZZ(Dmin),ZZ(Dmax)) if D.is_fundamental_discriminant()]
    nvalues = [1]

    fwrite(f'Doing cocycle {fname}...', logfile)
    p, M = get_p_prec(fname)

    # w = F.elements_of_norm(p)[0]
    # Fp = Qp(p,2*M, type='floating-point')
    # Rp = Zmod(p**M)
    #phi = F.hom([F.gen().minpoly().roots(Fp)[0][0]],check=False)
    # if phi(w).valuation() == 0:
    #    phi = F.hom([F.gen().minpoly().roots(Fp)[1][0]],check=False)
    # Ruv = PolynomialRing(PolynomialRing(Rp,'u'),'v')
    # inv_map_poly = Ruv.flattening_morphism()
    # map_poly = inv_map_poly.inverse()
    # load('orthogonal.sage')
    J = load(fname)
    label = get_label(fname)
    for D in Dvalues:
        for cycle_type in cycle_types:
            outfile = outdir + '/' + f'points_{cycle_type}_{p}_{label}_{M}.txt'
            nvalues = [1]
            for n in nvalues:
                try:
                    try:
                        fwrite(f'Computing {fname} {D},{n},{cycle_type}...', logfile)
                        Jtau0, hE = RMCEval(J, D, cycle_type, M, n=n, return_class_number=True)
                        Cp = Jtau0.parent()
                    except (ValueError, NotImplementedError, TypeError, RuntimeError, PrecisionError, SignalError, KeyboardInterrupt, ModuleNotFoundError) as e:
                        if 'no elements of norm -1' not in str(e):
                            fwrite(f'Skipping {fname} {D},{n},{cycle_type}...({str(e)})', logfile)
                        continue
                    fwrite(f'Computed Jtau for {fname} {D} {n} {cycle_type}.', logfile)
                    if Jtau0 == 1:
                        fwrite(f'Computed Jtau = 1 for {D = } {n = } {hE = }...', outfile)
                        continue
                    for Jtau, char in [(Jtau0.norm(), 'triv'), (Jtau0**2 / Jtau0.norm(), 'conj')]:
                        success = False
                        try:
                            with stopit.ThreadingTimeout(120) as to_ctx_mgr:
                                H, prime_list = get_predicted_field_and_prime_list(F, D, n, cycle_type, char, names='z', prime_bound=prime_bound)
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                raise ExpatError
                        except ExpatError:
                            fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }... ({char} not recognized since H could not be computed)', outfile)
                            H = None
                            prime_list = None
                        for deg in [2,4,8]:
                            with stopit.ThreadingTimeout(30) as to_ctx_mgr:
                                ans = recognize_DGL_algdep(Jtau, deg, outfile=None, roots_of_unity=[1])
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                fwrite('Someone got impatient...', logfile)
                                ans = None
                                break
                            if ans[0] is not None:
                                success = True
                                fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }...', outfile)
                                support = [q for q, _ in ans[-1] if q in ZZ]
                                dta = list(ans[1:])
                                dta.append('?' if H is None else H.defining_polynomial())
                                msg = f'SUCCESS {char} with algdep: ' + str(ans[0].unit_part().add_bigoh(10)) + str(tuple(dta))
                                if prime_list is None:
                                    msg += ' warning: primes not checked!'
                                else:
                                    if any(q not in prime_list for q in support):
                                        msg += ' warning: ' + str([q for q in support if q not in prime_list])
                                fwrite(msg, outfile)
                                break
                        if not success and prime_list is not None and ('big' not in cycle_type or H.degree() <= 16):
                            with stopit.ThreadingTimeout(3600) as to_ctx_mgr:
                                ans = recognize_DGL_lindep(Jtau, H, prime_list = prime_list, outfile=None, recurse_subfields=True, degree_bound=8, algorithm='pari')
                            if to_ctx_mgr.state in [to_ctx_mgr.TIMED_OUT, to_ctx_mgr.INTERRUPTED]:
                                fwrite(f'Someone got impatient...  {cycle_type = } {p = } {label = } {D = }, {n = }', logfile)
                                ans = None
                            if ans is not None:
                                success = True
                                ffpoly, rec, d, clist_ans, field = ans
                                fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }...', outfile)
                                fwrite(f'SUCCESS {char} with lindep: ' + str(Jtau.unit_part().add_bigoh(10)) + f'({ffpoly}, {field}, {rec}, {d}, 0, 1, {clist_ans}, {H.defining_polynomial()})', outfile)

                        if not success:
                            fwrite(f'Computed Jtau{char} = {(Jtau.trace()/2).lift()} for {D = } {n = } {hE = }... ({char} not recognized)', outfile)
                    fwrite('...', outfile)
                except KeyboardInterrupt:
                    fwrite(f'WARNING! Keyboard interrupt so skipping {cycle_type = } {p = } {label = } {D = }, {n = }', logfile)
                    sleep(1)
                # except Exception as e:
                #     fwrite(f'WARNING! Unhandled exception so skipping {cycle_type = } {p = } {label = } {D = }, {n = } : {str(e)}', logfile)



if __name__ == '__main__':
  Fire(evaluate_cocycle)

