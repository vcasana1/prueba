#!/bin/bash
for TH in 1.0 1.2 1.4; do
  python3 nufft_phi_v5.py zeross.bin nf_size$TH size $TH 32 > nf_size${TH}.log 2>&1 || { echo "phi fail $TH"; continue; }
  ./tauC_fromphi_v5 zeross.bin nf_size$TH full_size$TH > full_size${TH}.out 2> full_size${TH}.err
  python3 analyze_KC_v1.py full_size$TH > analysis_size${TH}_v5.txt
  rm -f nf_size${TH}_phi.bin nf_size${TH}_dphi.bin
  echo "done $TH"
done
echo ALLDONE
