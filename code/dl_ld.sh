#!/bin/bash
# resumable parallel range download of the 1000G EUR LD reference
cd $HOME/rc_nafld_mr/work
U=http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz; T=1567805147; N=8; S=$(( (T+N-1)/N ))
for i in $(seq 0 $((N-1))); do
  a=$((i*S)); b=$(( (i+1)*S-1 )); [ $b -ge $T ] && b=$((T-1))
  have=$(stat -c%s part$i 2>/dev/null || echo 0)
  need=$(( b-a+1 ))
  [ $have -ge $need ] && continue
  timeout 160 curl -sS -r $((a+have))-$b http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz >> part$i 2>/dev/null &
done
wait
tot=0; for i in $(seq 0 $((N-1))); do tot=$((tot+$(stat -c%s part$i 2>/dev/null || echo 0))); done
echo "downloaded $tot of $T bytes"
