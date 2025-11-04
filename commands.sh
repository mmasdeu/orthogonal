#!/bin/bash
for file in L0Jtuple_*.sobj; do
   for typ in all; do
      tmux new-session -d -s `basename $file | sed 's/·//g' | sed 's/\.//g'`$typ \
      "conda run -n darmonpoints-dev sage dglpoint.py eval_and_recognize $file $typ 1 3000 --outdir=evenodd --max_n=3 --max_hE=20 --max_lindep_primes=100 --degree_bound=16";
   done
done
