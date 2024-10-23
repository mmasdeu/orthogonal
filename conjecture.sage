from pysat.examples.hitman import Hitman

def quadform_support(Dlist, max_prime, D):
    ans = sorted(list(reduce(lambda x, y : x.union(y), (set(quadform_support_single(D, max_prime)) for D in Dlist), set([]))))
    if D is not None:
        K.<d> = QuadraticField(D)
        return [ell for ell in ans if len(K.ideal(ell).factor()) < 2]
    else:
        return ans

@cached_function
def quadform_support_single(D, max_prime=500):
    ells = []
    for Q in BinaryQF_reduced_representatives(D):
        Qpol = Q.polynomial()
        ells.extend([ ell for ell in prime_range(max_prime+1) if ((Qell := Qpol.change_ring(GF(ell))) != 0) and sum(e for _, e in Qell.factor()) == 2])
    return sorted(list(set(ells)))

def covering_discriminants(S, discriminants):
    ans = []
    for D in discriminants:
        try:
            if all(ell in quadform_support_single(D, max_prime=ell) for ell in S):
                ans.append(D)
        except ValueError:
            pass
    return ans


D = -17
p = 5
S = [2,5,17,19,29,37,41,73]

D = -377
p = 5
S = [2,5,11,17,19,23,29,31,47,61,101,107,109,113,149,157,197,337,349,389,397,401,433]

D = -553
p = 5
S = [2,11,23,29,37,89,97,101,181,241,257,313,433,593]

D = -41
p = 13
S = [2,13,17,23,29,31,41,53,101,109]

D = -177
p = 5
S = [13,17,29,37,41,53,61,197,257,293,449,461]

D = -33
p = 5
S = [5,13,53,61,113,137]

D = -105
p = 17
S = [3,5,7,17,23,61,173]

D = 12
p = 5
S = [19]

D = 137
p = 5
S = [3,23,43,47,67,71,79,127,163,179,191,199,223,271,307,331]

Dmax = D+1
S0 = sorted(list(set(S).difference(quadform_support([QuadraticField(D).discriminant()],max_prime=max(S), D=None))))
discriminants = sorted(list(set([delta for D0 in srange(-D, 0) if (delta := QuadraticField(D0).discriminant()).abs() < Dmax])))
# discriminants = [D0 for D0 in srange(Dmax) if D0.is_fundamental_discriminant()]

cover_sets = [covering_discriminants([p], discriminants=discriminants) for p in S0]


res = []
with Hitman(bootstrap_with=cover_sets, htype='sorted') as hitman:
    for hs in hitman.enumerate():
        supp = quadform_support(tuple(hs), max_prime=max(S0), D=D)
        score = RR(len(supp)) - len(S0)
        newres = (score, [D.factor() for D in hs])
        print(newres, set(supp).difference(set(S0)))
        res.append(newres)
        print(min(res, key=lambda x:x[0]))
        if len(res) > 10**3:
            break
res = sorted(res,key=lambda x:x[0])

