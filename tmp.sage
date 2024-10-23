# set_random_seed(1)
from dglpoint import *
load('orthogonal.sage')
p = 5
M = 20
label = '3·7_6'
cocycle(p, label, M)
xprec = evaluate(p, label, 9*17, 'smallCM', M)
print(xprec)
I = our_sqrt(Qp(p,M, type='floating-point')(-1))
xnum =  19^2 * (1+I)^8 * (1-4*I)^5 * (5+2*I)^3 * (5+4*I) * (5+6*I)^2
xden = (1+2*I)^5 * (1-2*I)^6 * (1-6*I)^2 * (9-4*I)^2 * (7-8*I)^2 *  (2-13*I)^2
print('This should be large: ', ((xprec.parent()(xnum/xden) / xprec) -1).valuation()) # should be large

for D in range(1,1000):
    try:
        J = evaluate(p, label, D, 'smallCM', M)
    except (ValueError, NotImplementedError, TypeError, RuntimeError, PrecisionError, SignalError) as e:
            print(f'Skipping {D},{n}...({str(e)})')
            if "No good matrix" in str(e):
                break
            continue
    print(f'{J = }')
