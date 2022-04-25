import numpy as np, pandas as pd, matplotlib.pyplot as plt, glob
from mpl_toolkits import mplot3d

files=glob.glob("*.txt")
dfs=pd.read_csv(files[2], sep='\s+',header=None)
dfx = pd.read_csv(files[3], sep='\s+',header=None) 
dfy = pd.read_csv(files[4], sep='\s+',header=None)

x=np.asarray(dfx)
y=np.asarray(dfy)
S=np.asarray(dfs)

X,Y = np.meshgrid(y,x)

plt.contourf(Y,X,S,levels=100, cmap='jet')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("solution.png", dpi=1200, bbox_inches='tight')
plt.show()

dfi=pd.read_csv(files[0], sep='\s+',header=None)
dfr=pd.read_csv(files[1], sep='\s+',header=None)

plt.semilogy(dfi,dfr)
plt.xlabel('Iterations')
plt.ylabel('Residual')
plt.savefig("residualSemilog.png", dpi=1200, bbox_inches='tight')
plt.show()

plt.loglog(dfi,dfr)
plt.xlabel('Iterations')
plt.ylabel('Residual')
plt.savefig("residualLogLog.png", dpi=1200, bbox_inches='tight')
plt.show()
