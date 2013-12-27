"""
:NAME:
dist

:DESC:
statistical functions
-Poisson
-Binomial

binomial recursive formula:
           P[x] * p(N-x)
  P[x+1] = -------------
              q(x+1)

written by matthieu defrance
received by defrance@bigre.ulb.ac.be
"""

from math import *

def pbinom_right_left(n, N, p):
    """Compute P(x>=n) if n >= N*p or P(x<=n) else (binomial)
    n -- number of success
    N -- number of trials
    p -- prob of success
    """
    if n >= N*p:
        return pbinom(n, N, p)
    else:
        return pbinom_left(n, N, p)


def pbinom(n, N, p):
    """Compute P(x>=n) (binomial)
    n -- number of success
    N -- number of trials
    p -- prob of success

    P[x] =  P[x-1] * p (N-x+1) / (q x)
    """
    assert(p > 0 and p <= 1.0)
    assert(n <= N)

    if n == 0 or p == 1:
        return 1.0

    logp = log(p)
    logq = log(1-p)
    logbin = N * logq

    S = 0.0
    for x in range(1, N+1):
        #try:
        #    logbin += LOG[N - x + 1] - LOG[x] + logp - logq
        #except KeyError:
        logbin += log(N - x + 1) - log(x) + logp - logq
        if x >= n:
            S += exp(logbin)
        if x > N*p and exp(logbin) <= 0.0:
            break

    return max(min(S, 1.0), 0.0)


def pbinom_left(n, N, p):
    """Compute P(x<=n) (binomial)
    n -- number of success
    N -- number of trials
    p -- prob of success
    """
    assert(p > 0 and p <= 1.0)
    assert(n <= N)

    if n == 0 or p == 1:
        return 1.0

    logp = log(p)
    logq = log(1-p)
    logbin = N * logq

    S = exp(logbin)
    for x in range(1, n+1):
        logbin += log(N - x + 1) - log(x) + logp - logq
        S += exp(logbin)

    return max(min(S, 1.0), 0.0)


def ppois(n, lamb):
    """Compute P(x>=n) (Poisson)
    """
    assert(n >= 1)
    assert (lamb >= 0)
    
    loglamb = log(lamb)
    logpois = -lamb + loglamb

    S = 0.0
    if n == 1:
        S += exp(logpois)

    x = 2
    while not (x > lamb and exp(logpois) <= 0):
        logpois += loglamb - log(x)
        if x >= n:
            S += exp(logpois)
        x += 1
    #print x, lamb, exp(logpois)
    return max(min(S, 1.0), 0.0)


def sum_of_binomials(proba, trials, a, b):
    logproba = log(proba)
    q = max(0, 1 - proba)
    logq = log(q)
    sum_of_bin = 0
    logbin = 0
    prev_value = -1
    expected = trials * proba
    logbin = trials * log(q)
    if a == 0:
	    sum_of_bin = exp(logbin)
	    loop_start = 1
    else:
	    loop_start = a
	    for x in range(1,a):
	        logbin += log(trials - x + 1) - log(x) + logproba - logq

    for x in range(loop_start, b+1):
	    logbin += log(trials - x + 1) - log(x) + logproba - logq
	    sum_of_bin += exp(logbin)
        
	    if (x > expected) and (sum_of_bin <= prev_value):
	        break
	    prev_value = sum_of_bin

    return max(min(sum_of_bin, 1.0), 0.0)


if __name__ == '__main__':
    pbinom(1, 3, 0.2)
    pbinom_right_left(1, 3, 0.2)
    pbinom_left(1, 3, 0.2)        
    ppois(4, 12.4)



