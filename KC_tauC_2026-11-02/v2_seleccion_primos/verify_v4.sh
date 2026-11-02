#!/bin/bash
# tau^C exacto en ventanas para varias reglas de seleccion
for W in "100000,101000" "1000000,1002000"; do
  for S in "paper 1" "count 4" "count 16" "size 0.8" "size 0.9" "size 1.0"; do
    set -- $S
    tag="w${W%%,*}_${1}${2}"
    SCHEME=$1 PARAM=$2 WIN=$W ./tauC_v4 zeross.bin $tag 32 2>/dev/null
  done
done
echo ALLDONE
