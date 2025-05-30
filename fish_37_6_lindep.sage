from util import *
outfile='success_37_6_300_lindep.txt'
Nunits = 20
L.<b> = Qq(5**2, 300)
load('Jlist_37_6_300.sage')
prime_list = prime_range(1000)
for D, n, J, hE in Jlist:
    if hE > 1:
        continue
    F.<r> = QuadraticField(D)
    alpha = F.elements_of_norm(-1)[0]
    x = QQ['x'].gen()
    E.<z> = F.extension(x^2-(alpha*n))
    Eabs.<zabs> = E.absolute_field()
    plist = [p for p in prime_list if len(Eabs.ideal(p).factor()) == 1]
    print(f'{len(plist) = }')
    ans = recognize_DGL_lindep(J, Eabs, prime_list=prime_list,class_number=hE, outfile=outfile, algorithm='pari')
    print(D, n, ans)
