for file in L0Jtuple_*.sobj; do
    tmux new-session -d -s `basename $file | sed 's/·//g' | sed 's/\.//g'` \
    "conda run -n sage sage find_all_types.sage $file bigRM 1 10000 newpoints";
done
