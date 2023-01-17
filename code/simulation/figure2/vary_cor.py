import numpy as np
import knockpy
import warnings
import random
from knockpy.knockoff_filter import KnockoffFilter
import scipy
import os
import statsmodels.api as sm

### Set working directory to 'figure2' folder, e.g., 
# os.chdir("~/code/simulation/figure2")

n = 500 # number of data points
p = 60  # number of features
### Change the replicate index
#replicate = int(os.getenv('SLURM_ARRAY_TASK_ID'))
#random.seed(replicate)
### Change the correlation
#rho = float(os.getenv("att"))
delta = 6.5
p1 = 30
q = 0.1


Sigma = np.eye(p)
for i in range(0,p):
    for j in range(0,p):
        Sigma[i,j] = np.power(rho, abs(i-j))


nonzero = random.sample(range(0, p), p1)
zero = np.setdiff1d(range(0, p), nonzero)
beta = np.zeros(p)
for nzr in nonzero:
    beta[nzr] = random.sample([-1,1],1)[0]*delta


# Sample X
X = np.random.multivariate_normal(mean=np.zeros(p), cov= Sigma, size=(n,))*1/np.sqrt(n)
mu = np.dot(X, beta)
def f(x):
    return np.exp(x)/(1+np.exp(x))
y = np.random.binomial(1, f(mu))


start = time.clock()
kfilter = KnockoffFilter(ksampler='gaussian', knockoff_kwargs={'method':'mvr'})
_ = kfilter.forward(X=X, y=y, fdr=q)
Xk = kfilter.Xk
# Fit MLE with X and Xk
binomial_model = sm.GLM(y, np.hstack((X, Xk)), family=sm.families.Binomial())
binomial_results = binomial_model.fit()
Z = binomial_results.params
tau = 1/np.sqrt(np.diag(np.linalg.inv(np.matmul(np.hstack((X, Xk)).T, np.hstack((X, Xk))))))
# Get feature importance
W = abs(Z[:p]*tau[:p]) - abs(Z[p:]*tau[p:])
kfilter.W = W
# Make selections
rejections = kfilter.make_selections(kfilter.W, fdr = 0.1)
#Calculate FDP and power
if rejections.sum()!=0:
	power = np.dot(rejections, beta != 0) / (beta != 0).sum()
	fdp = np.dot(rejections, beta == 0) / rejections.sum()
else:
	power = 0
	fdp = 0

end = time.clock()
knockoff_time = end-start


data_save = [fdp, power, knockoff_time]






