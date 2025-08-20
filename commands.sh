#!/bin/bash
for file in L0Jtuple_*.sobj; do
    tmux new-session -d -s `basename $file | sed 's/·//g' | sed 's/\.//g'` \
    "conda run -n darmonpoints-dev sage dglpoint.py eval_and_recognize $file bigRM 1 5000 --max_size=10000 --outdir=darmonfest";
done
