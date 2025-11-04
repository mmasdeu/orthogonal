#!/bin/bash
for typ in smallCM smallRM bigATR all; do
    conda run -n darmonpoints-dev sage ~/orthogonal/make_tables_search.sage ~/orthogonal/$1 $typ $prime html > dgltable_$1_$typ_$prime.html
    scp -i /home/masdeu/.ssh/id_rsa dgltable_$1_$typ.html masdeu@dixie.mat.uab.cat:~/www/
done

