#!/bin/bash
conda run -n darmonpoints-dev sage ~/orthogonal/make_tables_search.sage ~/orthogonal/august > dgltable_august.html
scp dgltable_august.html masdeu@dixie.mat.uab.cat:www/

