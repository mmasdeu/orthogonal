import pandas as pd
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
    data_list = re.search('.*O\(.*\^10\)\(([^,]*), ([^,]*), ([^,]*), ([^,]*), ([^,]*), (.*)\)', ln).groups()
    vals = []
    for o in data_list:
        try:
            newval = sage_eval(o, locals={'x':x})
        except (NameError, SyntaxError):
            newval = copy(o)
        vals.append(newval)
    poly, field, nrm, i, denom, pw = vals
    if type(pw) == list:
        try:
            pw = print_factorization(pw)
        except:
            pw = [(a, b) for b, a in pw] # DEBUG: for backward compatibility
            pw = print_factorization(pw)
    if 'algdep' in ln:
        pw += ' (a)'
    else:
        pw += ' (l)'
    return poly, field, pw

def get_header(fname):
    fsp = (fname.split('/')[-1]).strip('.txt').split('_')[1:]
    tp, p, lbl1, lbl2, prec = fsp
    label = lbl1 + '_' + lbl2
    return f'{tp} points. p = {str(p)}, label = {str(label)}, prec = {str(prec)}'

def make_table(fname):
    # print(f'Processing {fname}...')
    with open(fname, 'r') as f:
        ds = {}
        trivial_list = []
        unrecognized_list = []
        for ln in f:
            if ln[:3] == '...':
                continue
            if 'Computed' in ln:
                Dn, nn = get_Dn(ln)
                if 'Jtau = 1' in ln:
                    trivial_list.append((Dn, 'triv'))
                elif 'not recognized' in ln:
                    if 'triv' in ln:
                        unrecognized_list.append((Dn,'triv'))
                    else:
                        assert 'conj' in ln
                        unrecognized_list.append((Dn, 'conj'))
                continue
            elif 'SUCCESS' in ln:
                poly, field, pw = get_poly_and_field(ln)
                if 'triv' in ln:
                    ds[(Dn,'triv')] = {'poly' : poly, 'field' : field, 'factor' : pw} # only save D
                else:
                    assert 'conj' in ln
                    ds[(Dn,'conj')] = {'poly' : poly, 'field' : field, 'factor' : pw} # only save D

        return ds, trivial_list, unrecognized_list

def main(path='outfiles/'):
    type_list = ['smallCM', 'smallRM', 'bigRM']
    flist0 = glob.glob(path + '/*.txt')
    flist = []
    for typ in type_list:
        flist.extend([f for f in flist0 if typ in f])

    doc.asis('<!DOCTYPE html>')
    total_triv = 0
    total_unr = 0
    total_succ = 0
    with tag('html'):
        with tag('head'):
            doc.asis('<meta charset="utf-8">')
            with tag('title'):
                doc.asis('Tables for DGL points')
            # doc.asis('''<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script>''')
            # doc.asis('''<script src="https://cdn.datatables.net/2.0.5/js/dataTables.js"></script>''')
            doc.asis('''<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script><script type="text/javascript" id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>''')
            # doc.asis('''<script type="text/javascript">$(document).ready( function () {
            #  new DataTable('table.dgltable');
            #  } );</script>''')
            for typ in type_list:
                t1, t2 = [o for o in type_list if o != typ]
                doc.asis('''
                <script>
                function toggle_{typ}() {{
		var x = document.getElementsByClassName("{typ}");
                var y = document.getElementsByClassName("{t1}");
                var z = document.getElementsByClassName("{t2}");
		for(var i =0;i<x.length;i++){{
                    x[i].style.display = 'block'
                }}
		for(var i =0;i<y.length;i++){{
                    y[i].style.display = 'none'
                }}
		for(var i =0;i<z.length;i++){{
                    z[i].style.display = 'none'
                }}
	        }}
                </script>'''.format(typ=typ, t1=t1, t2=t2))
            doc.asis('''
            <script>
            function view_all() {
            var x = document.getElementsByClassName("smallCM");
            var y = document.getElementsByClassName("smallRM");
            var z = document.getElementsByClassName("bigRM");
            for(var i =0;i<x.length;i++){
                x[i].style.display = 'block'
            }
            for(var i =0;i<y.length;i++){
                y[i].style.display = 'block'
            }
            for(var i =0;i<z.length;i++){
                z[i].style.display = 'block'
            }
            }
            </script>''')

        doc.asis('<link rel="stylesheet" type="text/css" href="df_style.css"/>')
        with tag('body'):
            dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            line('p', f'Last update: {dt}')
            for typ in type_list:
                doc.asis(f'<button onclick="toggle_{typ}()">View {typ}</button>')
            doc.asis(f'<button onclick="view_all()">View All</button>')
            for f in flist:
                hdr = get_header(f)
                typ = next(typ for typ in type_list if typ in hdr)
                with tag('div', klass = typ):
                    line('h1', hdr)
                    tbl, triv, unr = make_table(f)
                    total_triv += len(triv)
                    total_unr += len(unr)
                    total_succ += len(tbl)
                    tbl = pd.DataFrame(tbl, index=['field', 'poly', 'factor']).transpose()
                    line('h2', 'Trival values')
                    text(str(triv))
                    line('h2', 'Unrecognized values')
                    if len(unr) == 0:
                        text('None. Hooray!')
                    else:
                        text(str(unr))
                    line('h2', 'Successes')
                    doc.asis(tbl.to_html(classes='dgltable', table_id='dgltable', escape=False))
            line('h1', "Summary")
            with tag('ul'):
                line('li', f'Total trivial values: {total_triv}')
                line('li', f'Total unrecognized values: {total_unr}')
                line('li', f'Total success values: {total_succ}')
    print(indent(doc.getvalue()))

if __name__ == "__main__":
    fire.Fire(main)
