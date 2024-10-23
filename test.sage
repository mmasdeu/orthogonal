
def recognize_DGL_lindep(J, L, prime_list, Cp = None, units=None, outfile=None, **kwargs):
    r'''
    '''
    recurse_subfields = kwargs.pop('recurse_subfields', False)
    degree_bound = kwargs.pop('degree_bound', Infinity)
    print(f'Doing DGL-lindep with {L = }')
    if recurse_subfields:
        kwargs['phi'] = None
        ans = None
        for L0, _, _ in L.subfields():
            if L0.degree() <= degree_bound:
                try:
                    ans = recognize_DGL_lindep(J, L0, prime_list, Cp=Cp, units=units, outfile=outfile, **kwargs)
                except ValueError:
                    pass
                if ans is not None:
                    return ans
        return None
    K = J.parent()
    p = K.prime()
    prec = J.precision_relative()
    # embed L in a field containing K
    if Cp is None:
        Cp = K
        K_to_Cp = K.hom([K.gen()])
    else:
        K_to_Cp = K.hom([polynomial_roots(K._exact_modulus(), Cp)[0]])
    embeddings = polynomial_roots(L.polynomial(), Cp)
    if len(embeddings) == 0:
        raise ValueError(f'L (={L}) does not embed into Cp (={Cp})')
    phi_list = kwargs.get('phi_list', None)
    if phi_list is None:
        phi_list = [L.hom([rt]) for rt in embeddings]

    for phi in phi_list:
        V = [None]
        Vlogs = [K_to_Cp(J.log(0))]
        W = ['J']
        lp = 0
        hL = 1
        glist = []
        for ell in prime_list:
            for pp, _ in L.ideal(ell).factor():
                verbose(f'(Factoring {ell}')
                hL0 = 0
                gens = [None, None]
                while len(gens) > 1:
                    hL0 += hL
                    gens = (pp**hL0).gens_reduced(proof=False)
                hL = hL0
                glist.append((gens[0], hL, ell))
        glist = [(g**ZZ(hL / e), ell) for g, e, ell in glist if phi(g).valuation() == 0]
        norm_one = J.norm() == 1
        while len(glist) > 0:
            g0, ell = glist[0]
            phiv = phi(g0)
            if norm_one:
                g1 = next((o for o, ellp in glist[1:] if ellp == ell and (phiv / phi(o)).norm().log() == 0), 1)
            else:
                g1 = 1
            V.append(g0 / g1)
            Vlogs.append((phiv / phi(g1)).log(0))
            W.append(ell)
            glist = [(o,ell) for o, ell in glist if o != g0 and o != g1]

        # Add units
        if units is None:
            units = list(L.units(proof=False))
        else:
            units = [L(o) for o in units]
        for i, u in enumerate(units):
            V.append(u)
            Vlogs.append(phi(u).log(0))
            W.append(f'u{i}')

        # Truncate precision if prec is specified
        prec = kwargs.get('prec',prec)
        if prec is not None:
            Vlogs = [o.add_bigoh(prec) for o in Vlogs]

        # OK now cross fingers...
        clist = kwargs.get('clist', None)
        if clist is not None:
            assert len(clist) == len(Vlogs), f'Provided clist is of incorrect length (should be {len(Vlogs)})'
            return sum([c * o for c, o in zip(clist, Vlogs)])
        verbose(f'Running lindep with {len(Vlogs) = }')
        if kwargs.get('debug_Vlogs', None) is not None:
            return Vlogs
        clist = our_lindep(Vlogs, **kwargs)
        verbose(f'clist = {clist}')
        verbose(f'W = {W}')
        verbose(f'Should be zero : {sum([c * o for c, o in zip(clist, Vlogs)])}')
        if clist[0] < 0:
            clist = [-o for o in clist]
        if clist[0] == 0:
            verbose(f'Not recognized: clist[0] = {clist[0]}')
            return None
        ht = 2 * sum((1+RR(o).abs()).log(p) for o in clist)
        verbose(f'(confidence factor: { ht / prec})')
        if ht > prec:
            verbose(f'Not recognized (confidence factor: ht / prec = {ht / prec}): clist = {clist}')
            return None
        clist_ans = [(u,v) for u,v in zip(W, clist) if v != 0]
        fwrite("# SUCCESS!", outfile)
        fwrite(f'# {clist_ans}', outfile)
        algebraic = kwargs.get('algebraic', True)
        if not algebraic:
            return clist_ans
        else:
            verbose(str(clist_ans))
            if not clist[0] > 0:
                verbose(f'Redundant set of primes?')
                return None
            # fact = Factorization([(u,-a) for u, a in zip(V[1:],clist[1:])])
            assert len(V) == len(clist)
            # assert clist[0] % hL == 0
            hL = 1 # DEBUG
            J_alg = prod(u**-a for u, a in zip(V[1:], clist[1:]))
            # J_alg = fact.prod() # DEBUG # (L['x'].gen()**hL - fact.prod()).roots(L)[0][0]
            remainder = clist[0]# // hL
            try:
                check = (phi(J_alg) / K_to_Cp(J)**remainder).log(lp).valuation() >= prec
                if not check:
                    print('Did not pass check! Returning value anyway...')
            except ValueError as e:
                print('Did not pass check because it errored! (error = %s)'%str(e))
            ffsmall = J_alg.parent().polynomial()
            try:
                ffsmall = sage_eval(str(pari('polredabs(%s)'%ffsmall)), locals={'x': QQ['x'].gen()})
            except PariError:
                pass
            return J_alg, J_alg.minpoly(), remainder, clist_ans, ffsmall  # DEBUG - didn't check that they match...
