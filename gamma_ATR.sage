# Related to Manin Trick

def quo_rem_in_gaussian_integers(a,b):
    if b == 0:
        raise ZeroDivisionError
    K = a.parent()
    q = ((a/b)[0]).round() + ((a/b)[1]).round()*K.gen()
    r = a - q*b
    assert r.norm() < b.norm()
    return (q, r)

# computes the continued fraction of a/c, where a and c belong to Z[i]
def continued_fraction(a, c):
    a_over_c = a/c
    cf = []
    while c != 0:
        q, r = quo_rem_in_gaussian_integers(a,c)
        cf.append(q)
        a = c
        c = r
    assert a_over_c == eval_cf(cf)
    return cf

# given a continued fraction, computes its value as an element of Q(i); only for debugging
def eval_cf(cf):
    if len(cf) == 1:
        return cf[0]
    if len(cf) == 2:
        return cf[0] + 1/cf[1]
    new_cf = [cf[i+1] for i in range(len(cf)-1)]
    return eval_cf([cf[0],eval_cf(new_cf)])

# computes the convergents of a continued fraction, following Hardy-Wright theorem 149
def compute_convergents(cf):
    N = len(cf)
    if N == 1:
        return([[cf[0]], [1]])
    p = [cf[0], cf[1] * cf[0] +1]
    q = [1, cf[1]]
    for n in range(2, N):
        p.append(cf[n] * p[n-1] + p[n-2])
        q.append(cf[n] * q[n-1] + q[n-2])
    assert all([p[n]*q[n-1] - p[n-1]*q[n] == (-1)^(n-1) for n in range(1,N)])
    #assert p[-1]/q[-1] == eval_cf(cf)
    return(p, q)

# given the convergents [p_i] and [q_j] of a continued fraction for alpha, computes the matrices that Mi such that {Infinity, alpha} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def compute_Ms(p, q):
    Ms = []
    # Ms.append(Matrix(2,2,[1,0,0,1])) # j = -1, not necessary since we do (Infinity, alpha) instead of (0, alpha) 
    Ms.append(Matrix(2,2,[-p[0], 1, -q[0], 0])) # j = 0
    for j in range(1, len(p)):
        Ms.append(Matrix(2,2,[(-1)^(j-1)*p[j],p[j-1],(-1)^(j-1)*q[j],q[j-1]]))
    assert all([A.det() == 1 for A in Ms])
    return Ms

# given a, c in Z[i] returns the matrices that Mi such that {Infinity, a/c} = {M0(0), M0(Infinity)} + {M1(0), M1(Infinity)} + ... as in display (2.1.8) of Cremona's book
def matrices_for_unimodular_path(a,c):
    cf = continued_fraction(a, c)
    p, q = compute_convergents(cf)
    return compute_Ms(p, q)


def good_matrices(m):
    a, b, c, d = m.list()
    F = m.parent().base_ring()
    r = F.gen()
    I = F.ideal(2)
    if F == QQ or F == ZZ:
        raise NotImplementedError
    if c in I:
        return [ (m,-1) ]
    elif d in I:
        m = matrix([[b, a], [d, c]])
        return [ (m, 1) ]
    else:
        for d1, d2 in cartesian_product_iterator([[1,-1,r,-r] for _ in range(2)]):
            delta = -(m.determinant())
            A, B, C, D = m.list()
            y = matrix([[-C, d1], [D, d2]]).determinant()/delta
            if y in I:
                x = matrix([[d1, A], [d2, -B]]).determinant()/delta
                m1 = matrix([[A,x],[C,y]])
                m2 = matrix([[x,B],[y,D]])
                m1.rescale_col(1, -d1**-1)
                m2.rescale_col(0, d2**-1)
                m1.swap_columns(0,1)
                assert m1.determinant() == 1
                assert m2.determinant() == 1
                return [(m1,1),(m2,-1)]
    raise RuntimeError('No matrix found!')


def module_generators(K):
    x=var('x')
    y=var('y')
    F=K.base_field()
    f=F.polynomial()
    g=K.relative_polynomial()
    a=F.gen()
    b=K.gen()

    # Equivalent pari objects
    FP=F.pari_bnf().subst(x,y)
    fP=pari(f)
    KP=K.pari_rnf()
    gP=KP[0]
    BP=gp.rnfhnfbasis(FP,gP)
    E=[gp.matbasistoalg(FP,BP.vecextract(1)).lift(),gp.matbasistoalg(FP,BP.vecextract(2)).lift()]

    A=Matrix(F,2,2,0)
    for jj in range(2):
        for ii in [1,2]:
            tmp=E[jj][ii,1].Vec().sage()
            if(len(tmp)==2):
                A[ii-1,jj]=tmp[0]*a+tmp[1]
            else:
                A[ii-1,jj]=tmp[0]
    return (Matrix(K,1,2,[1,b])*A).list()



# it accepts a real quadratic field F and an element alpha in F of norm -1; then E = F(sqrt(alpha) is the ATR extension and K = Q(i) is contained in the galois closure of E
def compute_gamma_tau_ATR(F, alpha):
    x = QQ['x'].gen()
    R.<y> = PolynomialRing(F)
    E.<w> = NumberField(y^2-alpha)
    K.<i> = NumberField(x^2 + 1)
    MM.<gMM> = E.galois_closure()
    # Now we redefine E, because we want to view it as a subfield of MM. We also construct L as a subfield of MM
    EEp = [f for f in MM.subfields() if f[0].degree() == 4 and f[0].is_isomorphic(E)]
    LLp = [f for f in MM.subfields() if f[0].degree() == 4 and not f[0].is_isomorphic(E)]
    E, sigmaE, _ = EEp[0]
    L, sigmaL, _ = LLp[1]

    # take the unit of norm 1 in E
    sigma = E.automorphisms()[1]
    found = False
    for uu in E.units():
        if uu.absolute_minpoly().degree() == 4 and uu*sigma(uu) == 1:
            u = uu
            found = True
    if found == False:
        raise RuntimeError('did not find a unit of relative norm one')
    # from u construct gamma and its galois conjugate gamma_p 
    gL = L.primitive_element()
    Gal_M_L = [t for t in MM.automorphisms() if t(sigmaL(gL)) == sigmaL(gL)]
    gE = E.primitive_element()
    Gal_M_E = [t for t in MM.automorphisms() if t(sigmaE(gE)) == sigmaE(gE)]
    Nu = Gal_M_L[1](sigmaE(u)) * sigmaE(u) # compute N_{M/L}(u)
    L_K.<gL_K> = L.relativize(K.embeddings(L)[0])
    u_L_K = Nu.minpoly().roots(L_K)[0][0] # this is Nu viewed as an element of the relative extension L_K
    a, c = u_L_K.vector()
    b, d = (u_L_K * gL_K).vector()
    gamma = Matrix(2,2,[a,b,c,d])
    # this matrix is constructed using the K-basis <1, gL_K>, but gL_K is not a generator of the ring of integers of L/K, so we need to change basis
    _, gOLK = module_generators(L_K) # gOLK is a generator of O_L_K over O_K
    O = L_K.maximal_order()
    X = Matrix([(1,0), gOLK.vector()]).transpose() # change of basis matrix
    # we check that gOLK is a generator of O_L_K over O_K, by checking that any generator of O_L_K when written in terms of <1, gOLK> has integral coefficients
    for gen in O.gens():
        if not all([o.is_integral() for o in (X.inverse()* vector(gen.vector())).list() ]):
            raise RuntimeError('it seems that we do not have an embedding of the maximal order')
    # now we conjugate gamma to express it in terms of the basis <1, gOLK>    
    gamma = X^-1 * gamma * X
    Ms = [Matrix(2,2,[1,0,0,1]), Matrix(2,2,[0,1,-1,0]), Matrix(2,2,[2,1,1,1]), Matrix(2,2,[1+i,1,i,1]), Matrix(2,2,[1+i,1,1,1]), Matrix(2,2,[i,1,0,1]), Matrix(2,2,[i,1,2*i,1])]
    found = False
    for M in Ms:
        aa, bb, cc, dd = (M*gamma*M.inverse()).list()
        if cc.real() % 2 == 0 and cc.imag() % 2== 0:
            found = True
            gamma = M*gamma*M.inverse()
            break
    if not found:
        raise RuntimeError('maybe there is no embedding into Gamma0(2)?')
    gamma_p = gamma.apply_morphism(a.parent().automorphisms()[1])
    assert u_L_K.minpoly() == gamma.minpoly()

    # Now we compute the matrix M and the scalar l such that gamma * M * gamma_p^{-1} = l * M
    A.<xx,yy,bb,cc,ll> = PolynomialRing(K)
    alpha = xx + i*yy
    alpha_p = xx-i*yy
    M = Matrix(A,2,2,[alpha, -bb,cc,-alpha_p])
    B.<x,y,b,c,l> = PolynomialRing(E)
    fAB = A.hom([x,y,b,c,l], B, check = False)
    s = K.automorphisms()[1]
    Meq = gamma*M - ll*M*gamma_p
    eqns =[fAB(M.det())] +  [fAB((z + z.map_coefficients(s))/2) for z in Meq.list()] + [fAB((z - z.map_coefficients(s))/(2*i)) for z in Meq.list()]  + [c-1]
    I = B.ideal(eqns)
    pt = I.variety()[0]
    alpha = sigmaE(pt[x]) + K.embeddings(MM)[0](i) * sigmaE(pt[y])
    alpha_p = sigmaE(pt[x]) - K.embeddings(MM)[0](i) * sigmaE(pt[y])
    M = Matrix(MM, 2, 2, [alpha, -sigmaE(pt[b]), sigmaE(pt[c]), -alpha_p])
    Mp = Matrix(MM, 2, 2, [alpha_p, -sigmaE(pt[b]), sigmaE(pt[c]), -alpha])
    iKM = K.embeddings(MM)[0]
    # assert that M is fiexed under the action of gamma
    assert gamma.apply_map(iKM) * M * gamma_p.apply_map(iKM)**-1 == sigmaE(pt[l]) * M
    tau1 = M[0][0]
    tau2 = - M[1][1]
    return gamma, gamma_p, tau1, tau2
    
x = QQ['x'].gen()
F.<r> = NumberField(x^2 - 13)
alpha = F.elements_of_norm(-1)[0]
gamma, gamma_p, tau1, tau2 = compute_gamma_tau_ATR(F, 6*alpha)
print(gamma)
# print(gamma.minpoly())
# mlist0 = matrices_for_unimodular_path(gamma[0,0], gamma[1,0])
# mlist = sum((good_matrices(m) for m in mlist0),[])
# print(mlist)
