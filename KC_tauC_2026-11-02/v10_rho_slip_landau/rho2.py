import numpy as np, sys
from rho import run
T=float(sys.argv[1]); h=float(sys.argv[2])
L=[("0",0.0),("0.5",0.5),("1",1.0),("2",2.0),("log2",np.log(2)),("log3",np.log(3)),("log4",np.log(4)),("log5",np.log(5)),("log6",np.log(6)),("log7",np.log(7))]
r,sd=run(T,[v for _,v in L],h)
print("T=%g"%T, " ".join("%s:%.5f"%(n,abs(r[v])) for n,v in L)); sys.stdout.flush()
