#!/bin/bash
conda run -n darmonpoints-dev sage ~/orthogonal/make_tables_search.sage ~/orthogonal/$1 > dgltable_$1.html
scp -i /home/masdeu/.ssh/id_rsa dgltable_$1.html masdeu@dixie.mat.uab.cat:~/www/

