#!/bin/bash
for typ in smallCM smallRM bigATR all; do
    if [ $1 == "outfiles" ]; then
        fname="html/dgltable_$typ.html"
    else
        fname="html/dgltable_$1_$typ.html"
    fi
    conda run -n darmonpoints-dev sage ~/orthogonal/make_tables_search.sage ~/orthogonal/$1 $typ --p=None --format=html > $fname
    scp -i /home/masdeu/.ssh/id_rsa $fname masdeu@dixie.mat.uab.cat:~/www/
done

