import pandas as pd
import time
import os
import glob
import fire
import re
import datetime
from yattag import Doc, indent
from copy import deepcopy
from dglpoint import cycle_support
doc, tag, text, line = Doc().ttl()


class str_or_int(SageObject):
    def __init__(self, x):
        self.val = x
    def _latex_(self):
        try:
            x = ZZ(self.val)
            return latex(x)
        except TypeError:
            x = self.val
            if x == 'J':
                return '\\text{J}'
            else:
                return f'{x[0]}_{{{x[1]}}}'

    def _repr_(self):
        if self.val not in ZZ:
            return str(self.val)
        p = ZZ(self.val)
        if p == 1:
            color = None
        elif p == 2:
            color = 'green'
        elif p % 4 == 3:
            color = 'blue'
        else:
            color = 'red'
        if color is None:
            return str(p)
        else:
            return(f'<span style="color:{color}">{p}</span>')


def print_factorization(fact):
    f1 = [(str_or_int(b),e) for b, e in fact[1:]]
    try:
        Jpow = fact[0][0]
        if isinstance(Jpow, str):
            Jpow = fact[0][1]
    except IndexError:
        Jpow = '1'
    Jpow = f'J^{Jpow}' if Jpow != 1 else 'J'
    ans_color = f'{Jpow} = ' + str(Factorization(f1, sort=False, simplify=False)).replace('(', '').replace(')','').replace('1 * ','')
    ans_raw = f'{Jpow} = ' + latex(Factorization(fact[1:])).replace('(', '').replace(')','').replace('1 * ','')
    return ans_color, ans_raw

def get_Dnh(ln):
    ans = list(re.search('D = (.[0-9]*) n = (.[0-9]*) hE = (.[^\.]*)...', ln).groups())
    ans[0] = int(ans[0])
    ans[1] = int(ans[1])
    return ans
    # return tuple(int(o) for o in re.search('D = (.[0-9]*) n = (.[0-9]*) hE = (.[^\.]*)...', ln).groups())

def get_poly_and_field(ln):
    x = ZZ['x'].gen()
    data_list = re.search('.*O\(.*\^10\)\(([^,]*), ([^,]*), ([^,]*), ([^,]*), ([^,]*), ([^,]*), (.*), ([^\)]*)\)(.*)', ln).groups()
    vals = []
    for o in data_list:
        try:
            newval = sage_eval(o, locals={'x':x})
        except (NameError, SyntaxError):
            newval = copy(o)
        vals.append(newval)
    poly, field, nrm, d, i, denom, pw, hpoly, msg = vals
    if type(pw) == list:
        pw, pw_raw = print_factorization(pw)
    if 'algdep' in ln:
        pw += ' (a)'
        pw_raw += ' (a)'
    else: # lindep
        pw += ' (l)'
        pw_raw += ' (l)'
    return poly, field, d, pw, pw_raw, hpoly

def get_header(fname):
    tp, p, label, prec = get_file_data(fname)
    return f'{tp} points. p = {str(p)}, label = {str(label)}, prec = {str(prec)}'

def get_file_data(fname):
    fsp = (fname.split('/')[-1]).strip('.txt').split('_')[1:]
    if len(fsp) == 5:
        tp, p, lbl1, lbl2, prec = fsp
        parity = 'even'
    else:
        tp, p, lbl1, lbl2, parity, prec = fsp
    label = lbl1 + '_' + lbl2 + '_' + parity
    return tp, p, label, prec

def make_table(fname):
    tp, p, label, prec = get_file_data(fname)
    val = {'p': p, 'label': label, 'type': tp}
    J = '?'
    maxD = 1
    with open(fname, 'r') as f:
        ds = []
        for ln in f:
            if ln[:3] == '...':
                continue
            if 'Jtriv' in ln and 'Jconj' in ln:
                Jtriv, Jconj = tuple(QQ(o) for o in re.search('Jtriv = (.[0-9/]*), Jconj = (.[0-9/]*)', ln).groups())
            if 'COMMENTS:' in ln:
                val['comments'] = ln.split('COMMENTS:')[1].strip()
            else:
                val['comments'] = ''
            if 'Computed' in ln:
                Dn, nn, hE = get_Dnh(ln)
                maxD = max(maxD, Dn)
                val['D'] = Dn
                val['n'] = nn         
                val['hE'] = hE
                    
                if 'Jtau = 1' in ln:
                    v1 = deepcopy(val)
                    v1.update({'field' : 'x', 'factor' : 'J = 1', 'poly' : 'x - 1', 'hpoly' : '-', 'J' : '1', 'trivial' : True, 'recognized' : True, 'o(ζ)' : 1})
                    v2 = deepcopy(v1)
                    v1['char'] = 'triv'
                    v2['char'] = 'conj'
                    ds.extend([v1,v2])
                else:
                    try:
                        J = QQ(re.search('J[a-z]* = (.[0-9/]*)', ln).groups()[0])
                    except AttributeError:
                        J = '?'
                if 'not recognized' in ln:
                    v1 = deepcopy(val)
                    char = 'triv' if 'triv' in ln else 'conj'
                    v1.update({'char' : char, 'field' : '?', 'factor' : '?', 'poly' : '?', 'hpoly' : '?', 'J' : J, 'trivial' : False, 'recognized' : False, 'o(ζ)' : '?'})
                    ds.append(v1)
            elif 'SUCCESS' in ln:
                poly, field, d, factor, factor_raw, hpoly = get_poly_and_field(ln)
                char = 'triv' if 'triv' in ln else 'conj'
                v1 = deepcopy(val)
                if v1['comments'] == '':
                    v1['comments'] = str(primes)
                trivial = True if str(poly) == "x - 1" else False
                v1.update({'char' : char, 'field' : field, 'factor' : factor, 'factor_raw': factor_raw, 'poly' : poly, 'hpoly' : hpoly, 'J' : J, 'trivial' : trivial, 'recognized' : True, 'o(ζ)' : d})
                ds.append(v1)
        ds = sorted(list(map(dict, frozenset(frozenset(i.items()) for i in ds))), key = lambda o : o['D'])
        return ds, (p, label, tp, maxD)

def make_latex(tbl, stats, last_update=None, p=None, typ=None):
    x = QQ['x'].gen()
    tbl = [o for o in tbl if o['recognized'] and not o['trivial'] and 'J = 1' not in o['factor_raw']]
    def add_dollars(o):
        ans = {}
        for entry in o:
            if entry in ['p', 'D', 'n']:
                ans[entry] = f'${latex(sage_eval(str(o[entry]), locals={"x": x}))}$'
            elif entry in ['field', 'poly']:
                ans[entry] = f'${latex(sage_eval(str(o[entry]), locals={"x": x}))}$'
            elif entry == 'factor_raw':
                ans['factor'] = f'${str(o[entry])[:-3]}$ {o[entry][-3:]}'
            elif entry == 'factor':
                continue
            else:
                ans[entry] = str(o[entry]).replace('_','\_')
        return ans

    tbl = [add_dollars(o) for o in tbl]
    df = pd.DataFrame(tbl, dtype=pd.StringDtype())
    if p is not None:
        df = df[['label', 'D', 'n', 'field', 'factor']]
        column_format = 'p{5em}p{3em}p{.5em}p{7cm}p{6cm}'
        caption = f'Recognized points of type {typ} for p = {p}'
        label = f'tbl:table-p{p}'
    else:
        df = df[['p', 'label', 'D', 'n', 'field', 'factor']]
        column_format = 'p{1em}p{5em}p{3em}p{.5em}p{7cm}p{6cm}'
        caption = f'Recognized points of type {typ} points for all primes'
        label = f'tbl:table-all'
    print(df.style.format_index(axis=1).hide(axis=0).to_latex(environment='longtable',\
        column_format=column_format, caption=caption, label=label, hrules=True))
    return

def make_tables(path='outfiles/', typ = 'all', p=None, format='html'):
    if typ == 'all':
        type_list = ['smallCM', 'smallRM', 'bigATR']
    else:
        type_list = [typ]
    flist0 = glob.glob(path + '/*.txt')
    if p is not None:
        type_list_p = [f'{t}_{p}' for t in type_list]
    else:
        type_list_p = type_list
    flist = [f for f in flist0 if any(t in f for t in type_list_p)]
    tbl = []
    stats = []
    for f in flist:
        hdr = get_header(f)
        typ = next(typ for typ in type_list if typ in hdr)
        newtbl, newstats = make_table(f)
        tbl.extend(newtbl)
        stats.append(newstats)
    tbl = sorted(tbl, key = lambda v : ZZ(v['D']))
    last_update = time.ctime(max(os.path.getmtime(os.path.join(root, f)) for root,_,files in os.walk(path) for f in files))
    if format == 'html':
        make_html(tbl, stats, last_update, p, typ)
    elif format == 'latex':
        make_latex(tbl, stats, last_update, p, typ)
    else:
        raise NotImplementedError(f'Format {format} not implemented')

def make_html(tbl, stats, last_update=None, p=None, typ=None):
    df = pd.DataFrame(tbl, dtype=pd.StringDtype())    
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            doc.asis('<meta charset="utf-8">')
            with tag('title'):
                doc.asis('Tables for DGL points')
            doc.asis('''
 <link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css'>
<link href="https://cdn.datatables.net/v/dt/jq-3.7.0/jszip-3.10.1/dt-2.0.6/b-3.0.2/b-colvis-3.0.2/b-html5-3.0.2/b-print-3.0.2/r-3.0.2/sb-1.7.1/datatables.min.css" rel="stylesheet">
<link rel='stylesheet' href='https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css'>
<link rel='stylesheet' href='https://cdn.datatables.net/buttons/1.2.2/css/buttons.bootstrap.min.css'>
<link rel='stylesheet' href='https://cdn.datatables.net/searchbuilder/1.7.1/css/searchBuilder.dataTables.css'>
<link rel='stylesheet' href='https://cdn.datatables.net/v/dt/dt-2.0.7/fh-4.0.1/datatables.min.css'>
<link rel='stylesheet' href='https://cdn.datatables.net/v/dt/dt-2.0.7/fh-4.0.1/datatables.min.css'>
            ''')

        doc.asis('<link rel="stylesheet" type="text/css" href="df_style.css"/>')
        with tag('body'):
            line('h1', 'Tables of DGL points')
            if last_update is not None:
                line('p', f'Last update: {last_update}')
            line('h2', 'Experimental observations')
            with tag('ul'):
                line('li', 'smallCM + conj is always 1: smallCM points are always defined over Qp')
                line('li', 'smallCM + triv is 1 when D = 1 (mod 8)')
                line('li', 'smallRM + triv is 1 when D = 0 (mod 4)')
                line('li', 'smallRM + triv is very often defined over Q(i), but we have some exceptions!')
            line('h2', 'Recovering data')
            with tag('ul'):
                line('li', 'For "triv": Cp.<t> = Qq(p**2, prec); J = Cp(2*J0)')
                line('li', 'For "conj": Cp.<t> = Qq(p**2, prec); a = Cp(J0); t0 = t - t.trace()/2; b = ((1-a**2)/t0.norm()).sqrt(); J = a + b * t0')
            line('h2', 'Summary')
            with tag('ul'):
                tot_trivial = len([o for o in tbl if o["factor"] == "J = (1)"])
                tot_unr = len([o for o in tbl if o["factor"] == "?"])
                tot_rows = len(tbl)
                recog = [o for o in tbl if not o['trivial'] and o['recognized']]
                total_recog = len(recog)
                recog_smallCM = len([o for o in recog if o['type'] == 'smallCM'])
                recog_smallRM = len([o for o in recog if o['type'] == 'smallRM'])
                recog_bigATR = len([o for o in recog if o['type'] == 'bigATR'])
                tot_recog = total_recog
                line('li', f'Total completely trivial: {tot_trivial} (simultaneously for both characters)')
                line('li', f'Total unrecognized: {tot_unr}')
                line('li', f'Total recognized (nontrivial): {tot_recog} (smallCM: {recog_smallCM}, smallRM: {recog_smallRM}, bigATR: {recog_bigATR})')
                line('li', f'Total rows: {tot_rows}')
                line('li', f'Largest computed computed discriminants:')
                with tag('ul'):
                    for p, label, typ, maxD in stats:
                        line('li', f'For p = {p}, label = {label}, type = {typ}: {maxD}')
            line('h2', 'Column description')
            with tag('ul'):
                line('li', 'p : prime')
                line('li', 'label : the cocycle being used. The numbers before the "_" have coefficient +1, while the numbers after have coefficient -1. Repeated numbers mean that the corresponding coefficient is larger.')
                line('li', 'D : absolute value of discriminant')
                line('li', 'type : either smallCM or smallRM')
                line('li', 'char : "triv" means J·J\' and "conj" means J/J\'')
                line('li', 'recognized : whether the point has been recognized as algebraic')
                line('li', 'field : polynomial defining the field HJ generated by Jtriv or Jconj')
                line('li', 'factor : a representation for the factorization of the ideal (Jtriv) or (Jconj). Only the prime below each prime in HJ is written. The list gives the primes in the support. Primes = 3 (mod 4) are represented in blue, while those = 1 (mod 4) are represented in red. The symbols (a) and (l) indicate which recognition algorithm was used.')
                line('li', 'poly : minimal polynomial of ζ·Jtriv or ζ·Jconj for a certain root of unity ζ')
                line('li', 'hpoly : polynomial defining the extension over which lindep was carried out')
                line('li', 'o(ζ) : order of ζ')
                line('li', 'J : an integer-encoded form for Jtriv or Jconj. See above for how to recover it')
            
            
            
            columns = ['p', 'label', 'D', 'n', 'type', 'char', 'trivial', 'recognized', 'field', 'factor', 'poly', 'hpoly', 'J', 'hE', 'o(ζ)', 'comments']
            doc.asis(df.to_html(classes=['table', 'table-striped', 'table-bordered' ],\
                        table_id='dgltable', sparsify=False,\
                        escape=False, index=False,\
                        columns=columns)[:-8]+"<tfoot>\n" + " ".join(["<th>"+ i +"</th>\n" for i in columns])+"</tr>\n  </tfoot></table>")

            doc.asis('''
            <!-- partial -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/pdfmake.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/vfs_fonts.js"></script>
<script src="https://cdn.datatables.net/v/dt/jq-3.7.0/jszip-3.10.1/dt-2.0.6/b-3.0.2/b-colvis-3.0.2/b-html5-3.0.2/b-print-3.0.2/r-3.0.2/sb-1.7.1/datatables.min.js"></script>
<script src='https://cdn.datatables.net/plug-ins/2.0.5/dataRender/ellipsis.js'></script>
<script src='https://cdn.datatables.net/fixedheader/4.0.1/js/dataTables.fixedHeader.min.js'></script>
            <script  src="./script.js"></script>
            ''')
    print((doc.getvalue()))

if __name__ == "__main__":
    fire.Fire(make_tables)
