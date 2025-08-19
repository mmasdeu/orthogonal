for file in L0Jtuple_*.sobj; do
    tmux new-session -d -s `basename $file | sed 's/·//g' | sed 's/\.//g'` \
    "conda run -n darmonpoints-dev sage dglpoint.py eval_and_recognize $file bigRM 1 10000 --max_size=9000 --outdir=august";
done
