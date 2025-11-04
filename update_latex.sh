#!/bin/bash
for prime in 5 13 17; do
    conda run -n darmonpoints-dev sage ~/orthogonal/make_tables_search.sage ~/orthogonal/$1 bigATR $prime latex > dgltable_$1_${prime}_bigATR.tex
done

