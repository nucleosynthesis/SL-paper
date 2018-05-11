# -*- python -*-

from __future__ import division
import numpy as np

"""
Calculation of next-to-leading order simplified likelihood coefficients

SL backgrounds b_i for signal region i are parametrised as
  b_i = A_i + B_i*th_i + C_i/2 * th_i^2,
where the th_i are distributed as a multivariate normal.

TODO:
 * Check marginalisation functions
 * Skew argument -> m3=None?
"""

def getCoeffCi_fast(m2ii, m3iii):
    "Compute c_i from the minimal necessary moment data"
    if not m3iii:
        return 0.0
    rtn = -np.sign(m3iii) * np.sqrt(2*m2ii)
    arg = 8 * m2ii**3 / m3iii**2 - 1
    if arg < 0:
        raise Exception("Negative sqrt term in c_III coeff extraction!")
    rtn *= np.cos(4*np.pi/3 + np.arctan(np.sqrt(arg))/3)
    return rtn

def getCoeffBi_fast(m2ii, ci):
    "Compute b_i from the second-order moment data and a precomputed c_i"
    return np.sqrt(m2ii - 2*ci**2)

def getCoeffAi_fast(m1i, ci):
    "Compute a_i from the first-order moment data and a precomputed c_i"
    return m1i - ci


def getCoeffCi(m1i, m2ii, m3iii):
    "Get c_i from a uniform interface of all moment data"
    return getCoeffCi_fast(m2ii, m3iii)

def getCoeffBi(m1i, m2ii, m3iii):
    "Get b_i from a uniform interface of all moment data (internally recomputes c_i)"
    return getCoeffBi_fast(m2ii, getCoeffCi_fast(m2ii, m3iii))

def getCoeffAi(m1i, m2ii, m3iii):
    "Get a_i from a uniform interface of all moment data (internally recomputes c_i)"
    return getCoeffAi_fast(m1i, getCoeffCi_fast(m2ii, m3iii))


def getCoeffsABCi(m1i, m2ii, m3iii, skew=True):
    "Efficiently compute all three SL coefficients for bin i"
    ci = getCoeffCi_fast(m2ii, m3iii) if skew else 0.0
    return getCoeffAi_fast(m1i, ci), getCoeffBi_fast(m2ii, ci), ci


def getCoeffsABC(m1, m2, m3=None, skew=True):
    "Efficiently compute all three SL coefficients for all bins"
    As, Bs, Cs = [], [], []
    if m3 is not None:
        #m3 = np.zeros(len(m1))
        assert len(m1) == len(m3)
    assert m2.shape == tuple([len(m1), len(m1)])
    for i in range(len(m1)):
        A, B, C = getCoeffsABCi(m1[i], m2[i][i], m3[i] if m3 is not None else None, skew)
        As.append(A)
        Bs.append(B)
        Cs.append(C)
    return As, Bs, Cs


def getRhoIJ(m1i, m2ii, m3iii, m1j, m2jj, m3jjj, m2ij, skew=True):
    "Compute the SL correlation matrix element rho_ij"
    if skew and (m3iii or m3jjj): #< TODO: not correct if only one m3 coeff is zero: FIX!
        epsilon = 1e-10
        ci = getCoeffCi_fast(m2ii, m3iii)
        cj = getCoeffCi_fast(m2jj, m3jjj)
        ci += epsilon if ci >= 0 else -epsilon
        cj += epsilon if cj >= 0 else -epsilon
        cicj = ci*cj
        # print ci, cj, cicj
        #
        bi = getCoeffBi_fast(m2ii, ci)
        bj = getCoeffBi_fast(m2jj, cj)
        bibj = bi*bj
        #
        discr1 = bibj**2
        discr2 = 8*cicj*m2ij
        discr = discr1 + discr2
        #
        rtn = (np.sqrt(abs(discr)) - bibj) / 4 / cicj
    else:
        rtn = m2ij / np.sqrt(m2ii*m2jj)
    return rtn

def getRho(m1, m2, m3, skew=True):
    "Compute the SL correlation matrix rho"
    # if m3 is None or m3 == 0:
    #     bg_diagsig = np.sqrt(bg_m2.diagonal())
    #     bg_sig = np.diag(bg_diagsig)
    #     bg_invsig = np.diag(np.reciprocal(bg_diagsig))
    #     bg_corr = reduce(np.matmul, [bg_invsig, BG_M2, bg_invsig])
    #
    # TODO: there must be a better, vectorised way to calculate this!
    rho = np.ones(m2.shape)
    for i in range(m2.shape[0]):
        for j in range(m2.shape[1]):
            rhoij = getRhoIJ(m1[i], m2[i,i], m3[i] if m3 is not None else 0.0,
                             m1[j], m2[j,j], m3[j] if m3 is not None else 0.0,
                             m2[i,j], skew)
            rho[i,j] = rhoij
            if i != j:
                rho[j,i] = rhoij
    return rho


class SLParams(object):
    """Unified container for SL background moments, observed data, signal
    prediction, computed coefficients, and statistical functions"""

    def __init__(self, bg_m1, bg_m2, bg_m3=None, obs=None, sig=None, skew=True):
        import scipy.stats as st
        self.bg_m1 = np.array(bg_m1)
        self.bg_m2 = np.array(bg_m2).reshape([self.size,self.size])
        self.bg_m3 = np.array(bg_m3) if bg_m3 is not None else None
        self.obs = np.array(obs) if obs is not None else None
        self.sig = np.array(sig) if sig is not None else None
        #
        self.aparams, self.bparams, self.cparams = getCoeffsABC(self.bg_m1, self.bg_m2, self.bg_m3, skew=skew)
        self.rhoparams = getRho(self.bg_m1, self.bg_m2, self.bg_m3, skew=skew)
        self.thetadbn = st.multivariate_normal(np.zeros(self.size), self.rhoparams)
        #
        self.llmax = None #< for caching of unconditional max LL

    @property
    def size(self):
        "Number of bins / signal regions"
        return len(self.bg_m1)

    def expbi(self, i, theta_i):
        "Expected background rate in bin i, for SL nuisance component theta_i"
        return self.aparams[i] + self.bparams[i]*theta_i + self.cparams[i]*theta_i**2

    def expbs(self, thetas):
        "Expected background rates in all bins, for SL nuisance vector thetas"
        bs = np.array([self.expbi(i, th) for (i, th) in enumerate(thetas)])
        bs[bs <= 0] = 1e-3 #< force any negatives to ~0
        return bs

    def expsi(self, i, mu):
        "Expected signal rate in bin i, for signal strength mu"
        return mu*self.sig[i]

    def expss(self, mu):
        "Expected signal rates in all bins, for signal strength mu"
        return mu*self.sig

    def expsbi(self, i, mu, theta_i):
        "Expected s+b rate in bin i, for signal strength mu and SL nuisance theta_i"
        return self.expsi(i, mu) + self.expbi(i, theta_i)

    def expsbs(self, mu, thetas):
        "Expected s+b rates in all bins, for signal strength mu and SL nuisance vector thetas"
        return self.expss(mu) + self.expbs(thetas)

    def sampleth(self, nsamp):
        "Draw nsamp theta vectors from the SL distribution"
        return self.thetadbn.rvs(nsamp)

    def sampleb(self, nsamp):
        "Draw nsamp background arrays from the SL distribution and parametrisation"
        ths = self.thetadbn.rvs(nsamp)
        bs = np.zeros_like(ths)
        for n in range(nsamp):
            bs_n = self.expbgs(ths[n,:])
            bs[n,:] = bs_n
        return bs

    def samplesb(self, mu, nsamp):
        "Draw nsamp s+b arrays from the SL distribution and parametrisation, with signal strength mu"
        bs = self.sampleb(nsamp)
        sbs = mu*self.sig + bs #< check broadcasting
        return sbs

    def like(self, mu, thetas):
        """Calculate the simplified likelihood for the given (mu, {theta}) params

        TODO:
         * Broadcast across multidimensional numpy arrays of thetas and mus
        """
        import scipy.stats as st
        assert self.obs is not None
        assert mu == 0 or self.sig is not None
        p_th = self.thetadbn.pdf(thetas)
        counts = self.expsbs(mu, thetas)
        ps_poiss = st.poisson.pmf(self.obs, counts)
        p_poiss = np.prod(ps_poiss)
        p_tot = p_th * p_poiss
        return p_tot

    def loglike(self, mu, thetas):
        """Calculate the simplified log-likelihood for the given (mu, {theta}) params

        TODO:
         * Broadcast across multidimensional numpy arrays of thetas and mus
        """
        import scipy.stats as st
        assert self.obs is not None
        assert mu == 0 or self.sig is not None
        ll_th = self.thetadbn.logpdf(thetas)
        counts = self.expsbs(mu, thetas)
        lls_poiss = st.poisson.logpmf(self.obs, counts)
        ll_poiss = np.sum(lls_poiss)
        ll_tot = ll_th + ll_poiss
        return ll_tot

    def maxloglike(self, mu=None, rtnparams=False):
        """Calculate the conditional or unconditional maximum log likelihood

        If mu=None, the unconditional (i.e. absolute) maximum LL will be computed
        over (mu, {theta}). Otherwise, the likelihood conditional on mu wil be
        computed over {theta}.

        TODO:
         * Speed up minimisation with analytic gradient function
        """
        from scipy.optimize import minimize
        # def calc_gradll(mu_thetas):
        #     pass
        if mu is None:
            if self.llmax is not None: #< return cached value if available
                return self.llmax
            def calc_optll_unconditional(mu_thetas):
                return -self.loglike(mu_thetas[0], mu_thetas[1:])
            minres = minimize(calc_optll_unconditional, [0 for _ in range(self.size+1)])
            llopt = -calc_optll_unconditional(minres.x)
            self.llmax = llopt
        else:
            def calc_optll_conditional(thetas, mu): #< capture mu?
                return -self.loglike(mu, thetas)
            minres = minimize(calc_optll_conditional, np.zeros(self.size), args=(mu,))
            llopt = -calc_optll_conditional(minres.x, mu)
        # print llopt
        if rtnparams:
            return llopt, minres.x
        else:
            return llopt

    def tmu(self, mu):
        """Calculate profile likelihood chi2-distributed test statistic, t_mu = -2 ln lambda_mu"""
        return -2 * (self.maxloglike(mu) - self.maxloglike()) #< NB. use of cached unconditional maxLL

    def lambdamu(self, mu):
        """Calculate profile likelihood ratio, lambda_mu = L(mu, theta_hathat)/L(mu_hat, theta_hat)"""
        return np.exp(self.tmu(mu))

    def avgloglike(self, mu, nsamp=10000):
        """Calculate marginalised LL = ln <p_SL>, using nsamp samples

        TODO:
         * Sampling convergence is *very* slow...
         * Automatic convergence checking (use doubling scheme)
        """
        ths = self.sampleth(nsamp)
        p_avg = 0
        for n in range(nsamp):
            p_avg += self.like(mu, ths[n,:])/nsamp
        ll_avg = np.log(p_avg)
        return ll_avg
