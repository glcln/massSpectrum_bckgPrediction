import numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson
 
# données observées
n_obs = np.array([5, 2])
lam0 = np.array([3.0, 1.5])
a = np.array([0.2, 0.2])  # impact systématique par canal
 
# NLL (Poisson + nuisances indépendantes)
def nll(theta, n):
    lam = lam0*np.exp(a*theta)        # theta est maintenant un vecteur
    pois = np.sum(lam - n*np.log(lam))
    gauss = 0.5*np.sum(theta**2)      # contraintes indépendantes
    return pois + gauss
 
# fit profilé
def fit(n):
    res = minimize(lambda t: nll(t, n), np.zeros_like(n))
    return res.fun
 
# statistique observée
nll_obs = fit(n_obs)
nll_nom = nll(np.zeros_like(n_obs), n_obs)
q_obs = 2*(nll_nom - nll_obs)
 
print("q_obs =", q_obs)
 
# Monte Carlo
Nmc = 50000
q_mc = []
 
for _ in range(Nmc):
    sim = poisson.rvs(lam0)
    nll_sim = fit(sim)
    nll_nom_sim = nll(np.zeros_like(sim), sim)
    q = 2*(nll_nom_sim - nll_sim)
    q_mc.append(q)
 
q_mc = np.array(q_mc)
p_value = np.mean(q_mc >= q_obs)
 
print("p-value =", p_value)
