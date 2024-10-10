import pandas as pd
import time
import os
import glob
import fire
import re
import datetime
from yattag import Doc, indent
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
        return str(self.val)


def print_factorization(fact):
    f1 = [(str_or_int(b),e) for b, e in fact[1:]]
    Jpow = fact[0][0]
    if isinstance(Jpow, str):
        Jpow = fact[0][1]
    plist0 = sorted(list(set([b for b, e in fact[1:] if not isinstance(b, str)])))
    if len(plist0) == 0:
        plist = ['1']
    else:
        plist = ['[']
        for p in plist0:
            if p == 2:
                color = 'green'
            elif p % 4 == 3:
                color = 'blue'
            else:
                color = 'red'
            plist.append(f'<span style="color:{color}">{p}</span>,')
        plist[-1] = plist[-1][:-1]
        plist.append(']')
    Jpow = f'J^{Jpow}' if Jpow != 1 else 'J'
    return f'{Jpow} = ' + str(Factorization(f1, sort=False, simplify=False))[6:] + \
        f' ({"".join(plist)})'

def get_Dn(ln):
    return tuple(int(o) for o in re.search('D = (.[0-9]*) n = (.[0-9]*)', ln).groups())

def get_poly_and_field(ln):
    x = ZZ['x'].gen()
    data_list = re.search('.*O\(.*\^10\)\(([^,]*), ([^,]*), ([^,]*), ([^,]*), ([^,]*), ([^,]*), (.*), (.*)\)(.*)', ln).groups()
    vals = []
    for o in data_list:
        try:
            newval = sage_eval(o, locals={'x':x})
        except (NameError, SyntaxError):
            newval = copy(o)
        vals.append(newval)
    poly, field, nrm, d, i, denom, pw, hpoly, msg = vals
    if type(pw) == list:
        pw = print_factorization(pw)
    if 'algdep' in ln:
        pw += ' (a)' + msg
    else: # lindep
        pw += ' (l)' + msg
    return poly, field, d, pw, hpoly

def get_header(fname):
    fsp = (fname.split('/')[-1]).strip('.txt').split('_')[1:]
    tp, p, lbl1, lbl2, prec = fsp
    label = lbl1 + '_' + lbl2
    return f'{tp} points. p = {str(p)}, label = {str(label)}, prec = {str(prec)}'

def get_file_data(fname):
    fsp = (fname.split('/')[-1]).strip('.txt').split('_')[1:]
    tp, p, lbl1, lbl2, prec = fsp
    label = lbl1 + '_' + lbl2
    return tp, p, label, prec

def make_table(fname):
    # print(f'Processing {fname}...')
    tp, p, label, prec = get_file_data(fname)
    val = {'p': p, 'label': label, 'type': tp}
    J = '?'
    with open(fname, 'r') as f:
        ds = []
        for ln in f:
            if ln[:3] == '...':
                continue
            if 'Jtriv' in ln and 'Jconj' in ln:
                Jtriv, Jconj = tuple(QQ(o) for o in re.search('Jtriv = (.[0-9/]*), Jconj = (.[0-9/]*)', ln).groups())
            if 'Computed' in ln:
                Dn, nn = get_Dn(ln)
                val['D'] = Dn
                if 'Jtau = 1' in ln:
                    v1 = copy(val)
                    v1.update({'field' : 'x', 'factor' : 'J = (1)', 'poly' : 'x - 1', 'hpoly' : '-', 'J' : '1', 'trivial' : True, 'recognized' : True, 'o(ζ)' : 1})
                    v2 = copy(v1)
                    v1['char'] = 'triv'
                    v2['char'] = 'conj'
                    ds.extend([v1,v2])
                else:
                    try:
                        J = QQ(re.search('J[a-z]* = (.[0-9/]*)', ln).groups()[0])
                    except AttributeError:
                        J = '?'
                if 'not recognized' in ln:
                    v1 = copy(val)
                    J, char = (J, 'triv') if 'triv' in ln else (J, 'conj')
                    v1.update({'char' : char, 'field' : '?', 'factor' : '?', 'poly' : '?', 'hpoly' : '?', 'J' : J, 'trivial' : False, 'recognized' : False, 'o(ζ)' : '?'})
                    ds.append(v1)
            elif 'SUCCESS' in ln:
                poly, field, d, factor, hpoly = get_poly_and_field(ln)
                J, char = (J, 'triv') if 'triv' in ln else (J, 'conj')
                v1 = copy(val)
                trivial = True if str(poly) == "x - 1" else False
                v1.update({'char' : char, 'field' : field, 'factor' : factor, 'poly' : poly, 'hpoly' : hpoly, 'J' : J, 'trivial' : trivial, 'recognized' : True, 'o(ζ)' : d})
                ds.append(v1)
        return ds

def main(path='outfiles/'):
    type_list = ['smallCM', 'smallRM', 'bigRM']
    flist0 = glob.glob(path + '/*.txt')
    flist = []
    for typ in type_list:
        flist.extend([f for f in flist0 if typ in f])

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
            tbl = []
            for f in flist:
                hdr = get_header(f)
                typ = next(typ for typ in type_list if typ in hdr)
                tbl.extend(make_table(f))
            tbl = sorted(tbl, key = lambda v : ZZ(v['D']))
            df = pd.DataFrame(tbl, dtype=pd.StringDtype())


            line('h1', 'Tables of DGL points')
            last_update = time.ctime(max(os.path.getmtime(root) for root,_,_ in os.walk(path)))
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
                tot_recog = len([o for o in tbl if not o['trivial'] and o['recognized']])
                line('li', f'Total completely trivial: {tot_trivial} (simultaneously for both characters)')
                line('li', f'Total unrecognized: {tot_unr}')
                line('li', f'Total recognized (nontrivial): {tot_recog}')
                line('li', f'Total rows: {tot_rows}')
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
                line('li', 'hpoly : polynomial defining the compositum of Q(i) with the Hilbert class field of real/imaginary field of absolute discriminant D.')
                line('li', 'o(ζ) : order of ζ')
                line('li', 'J : an integer-encoded form for Jtriv or Jconj. See above for how to recover it')


            columns = ['p', 'label', 'D', 'type', 'char', 'trivial', 'recognized', 'field', 'factor', 'poly', 'hpoly', 'o(ζ)', 'J']
            doc.asis(df.to_html(classes=['table', 'table-striped', 'table-bordered' ], table_id='dgltable', sparsify=False, escape=False, index=False, columns=columns)[:-8]\
+"<tfoot>\n" + " ".join(["<th>"+ i +"</th>\n" for i in columns])+"</tr>\n  </tfoot></table>")

            doc.asis('''
            <!-- partial -->
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/pdfmake.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/vfs_fonts.js"></script>
<script src="https://cdn.datatables.net/v/dt/jq-3.7.0/jszip-3.10.1/dt-2.0.6/b-3.0.2/b-colvis-3.0.2/b-html5-3.0.2/b-print-3.0.2/r-3.0.2/sb-1.7.1/datatables.min.js"></script>
<script src='https://cdn.datatables.net/plug-ins/2.0.5/dataRender/ellipsis.js'></script>
<script src='https://cdn.datatables.net/fixedheader/4.0.1/js/dataTables.fixedHeader.min.js'></script>
            <script  src="./script.js"></script>
            ''')
    print(indent(doc.getvalue()))

if __name__ == "__main__":
    fire.Fire(main)
