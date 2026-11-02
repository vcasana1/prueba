#!/bin/bash
# usage: run_alto_v8.sh T0 W tag
T0=$1; W=$2; TAG=$3
python3 zwin_v8.py $T0 $W 256 z_$TAG > z_$TAG.log 2>&1 || exit 1
python3 phiwin_v8.py $T0 $W 32 pw_$TAG > pw_$TAG.log 2>&1 || exit 1
python3 - <<PY
import numpy as np, subprocess
d=np.load('z_$TAG'+'_zeros.npz'); K0=int(d['K0']); fr0=float(d['fr0'])
subprocess.run(['./subnivel_rel_v8','pw_$TAG','$T0','$W','32',str(K0),repr(fr0),'tc_$TAG.txt'])
PY
rm -f pw_${TAG}_phi.bin pw_${TAG}_dphi.bin
echo DONE
