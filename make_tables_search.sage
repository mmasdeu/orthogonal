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

def get_file_data(fname):
    fsp = (fname.split('/')[-1]).strip('.txt').split('_')[1:]
    tp, p, lbl1, lbl2, prec = fsp
    label = lbl1 + '_' + lbl2
    return tp, p, label, prec

def make_table(fname):
    # print(f'Processing {fname}...')
    tp, p, label, prec = get_file_data(fname)
    val = {'p': p, 'label': label, 'type': tp}
    with open(fname, 'r') as f:
        ds = []
        for ln in f:
            if ln[:3] == '...':
                continue
            if 'Computed' in ln:
                Dn, nn = get_Dn(ln)
                val['D'] = Dn
                if 'Jtau = 1' in ln:
                    v1 = copy(val)
                    v1.update({'field' : 'x', 'factor' : 'J = (1)', 'poly' : 'x-1'})
                    v2 = copy(v1)
                    v1['char'] = 'triv'
                    v2['char'] = 'conj'
                    ds.extend([v1,v2])
                elif 'not recognized' in ln:
                    v1 = copy(val)
                    v1.update({'char' : 'triv' if 'triv' in ln else 'conj', 'field' : '?', 'factor' : '?', 'poly' : '?'})
                    ds.append(v1)
                continue
            elif 'SUCCESS' in ln:
                poly, field, pw = get_poly_and_field(ln)
                char = 'triv' if 'triv' in ln else 'conj'
                v1 = copy(val)
                v1.update({'char' : char, 'field' : field, 'factor' : pw, 'poly' : poly})
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
<link rel='stylesheet' href='https://cdn.datatables.net/1.10.12/css/dataTables.bootstrap.min.css'>
<link rel='stylesheet' href='https://cdn.datatables.net/buttons/1.2.2/css/buttons.bootstrap.min.css'>''')
        doc.asis('<link rel="stylesheet" type="text/css" href="df_style.css"/>')
        with tag('body'):
            dt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            line('p', f'Last update: {dt}')
            tbl = []
            for f in flist:
                hdr = get_header(f)
                typ = next(typ for typ in type_list if typ in hdr)
                tbl.extend(make_table(f))
            df = pd.DataFrame(tbl, dtype=pd.StringDtype())# .transpose()
            columns = ['p', 'label', 'D', 'type', 'char', 'field', 'factor', 'poly']
            doc.asis(df.to_html(classes=['table', 'table-striped', 'table-bordered' ], table_id='dgltable', sparsify=False, escape=False, index=False, columns=columns)[:-8]\
+"<tfoot>\n" + " ".join(["<th>"+ i +"</th>\n" for i in columns])+"</tr>\n  </tfoot></table>")

            doc.asis('''
<!-- partial -->
  <script src='https://cdnjs.cloudflare.com/ajax/libs/jquery/3.1.0/jquery.min.js'></script>
<script src='https://cdn.datatables.net/1.10.12/js/jquery.dataTables.min.js'></script>
<script src='https://cdn.datatables.net/buttons/1.2.2/js/dataTables.buttons.min.js'></script>
<script src='https://cdn.datatables.net/buttons/1.2.2/js/buttons.colVis.min.js'></script>
<script src='https://cdn.datatables.net/buttons/1.2.2/js/buttons.html5.min.js'></script>
<script src='https://cdn.datatables.net/buttons/1.2.2/js/buttons.print.min.js'></script>
<script src='https://cdn.datatables.net/1.10.12/js/dataTables.bootstrap.min.js'></script>
<script src='https://cdn.datatables.net/buttons/1.2.2/js/buttons.bootstrap.min.js'></script>
<script src='https://cdnjs.cloudflare.com/ajax/libs/jszip/2.5.0/jszip.min.js'></script>
<script src='https://cdn.rawgit.com/bpampuch/pdfmake/0.1.18/build/vfs_fonts.js'></script>
<script src='https://cdn.rawgit.com/bpampuch/pdfmake/0.1.18/build/pdfmake.min.js'></script><script  src="./script.js"></script>
            ''')
    print(indent(doc.getvalue()))

if __name__ == "__main__":
    fire.Fire(main)
