#!/bin/bash
TH=1.4
python3 nufft_phi_v6.py zeross.bin nf_size$TH size $TH 32 > nf_size${TH}.log 2>&1 || { echo "phi fail"; exit 1; }
./tauC_fromphi_v5 zeross.bin nf_size$TH full_size$TH > full_size${TH}.out 2> full_size${TH}.err
python3 analyze_KC_v1.py full_size$TH > analysis_size${TH}_v5.txt
rm -f nf_size${TH}_phi.bin nf_size${TH}_dphi.bin nf_size1.2_phi.bin nf_size1.2_dphi.bin nf_size1.0_dphi.bin
echo ALLDONE
