from util import *
outfile='success_6_12_300.txt'
Nunits = 20
L.<b> = Qq(5**2, 300)
load('Jlist_6_12_300.sage')
for D, n, J, hE in Jlist:
    F.<r> = QuadraticField(D)
    alpha = F.elements_of_norm(-1)[0]
    x = QQ['x'].gen()
    E.<z> = F.extension(x^2-(alpha*n))
    Eabs.<zabs> = E.absolute_field()
    u1, u2 = Eabs.units()
    philist = [Eabs.hom([o]) for o, e in Eabs.defining_polynomial().change_ring(L).roots()]
    for phi in philist:
        for i in range(-Nunits, Nunits+1):
            for j in range(-Nunits, Nunits+1):
                J1 = phi(u1)**i * phi(u2)**j * J
                ans = recognize_DGL_algdep(J1, 2)
                if ans != (None,None,None):
                    fwrite(str(D, n, i, j, 2), outfile)
                    ans = recognize_DGL_algdep(J1, 4)
                if ans != (None,None,None):
                    fwrite(str(D, n, i, j, 4), outfile)
