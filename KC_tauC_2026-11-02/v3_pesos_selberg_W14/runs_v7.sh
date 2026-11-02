#!/bin/bash
# tag scheme profile param
while read TAG SCH PROF PAR; do
  [ -z "$TAG" ] && continue
  if [ ! -f nf7_${TAG}_phi.bin ]; then python3 nufft_phi_v7.py zeross.bin nf7_$TAG $SCH $PROF $PAR 32 > nf7_${TAG}.log 2>&1 || { echo "phi fail $TAG"; continue; }; fi
  ./tauC_fromphi_v5 zeross.bin nf7_$TAG full7_$TAG > full7_${TAG}.out 2> full7_${TAG}.err
  python3 analyze_KC_v1.py full7_$TAG > analysis7_${TAG}.txt
  rm -f nf7_${TAG}_phi.bin nf7_${TAG}_dphi.bin
  echo "done $TAG"
done <<LIST
selb_X200 wfixed selberg 200
w14_X200 wfixed w14 200
selb_th0.5 wsize selberg 0.5
w14_th0.5 wsize w14 0.5
selb_th0.7 wsize selberg 0.7
w14_th0.7 wsize w14 0.7
LIST
echo ALLDONE
