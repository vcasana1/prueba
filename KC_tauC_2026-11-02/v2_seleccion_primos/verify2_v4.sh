#!/bin/bash
for P in 1.1 1.2 1.4; do SCHEME=size PARAM=$P WIN=100000,101000 ./tauC_v4 zeross.bin w100000_size$P 32 2>/dev/null; done
for P in 1.1; do SCHEME=size PARAM=$P WIN=1000000,1002000 ./tauC_v4 zeross.bin w1000000_size$P 32 2>/dev/null; done
echo ALLDONE
