# csm
An R package for Causal State Models

# The Proposal
An R package “csm” with a framework for using Causal State Models.

#Info
Causal State Models are equivalent of Hidden Markov Models with a property, that state transition probability and observation probability are united. They are the unique minimal sufficient statistics for a sequence of random variables Xt drawn from a discrete alphabet A.
CSSR (Causal-State Splitting Reconstruction) algorithm is the way to estimate Causal State Model from the data using statistical tests.
Here is a paper that describes algorithm and the theory in details

http://arxiv.org/abs/cs/0210025

There are also an improvement, that gives algorithm ability to be even more automatic ( find the best values of some of the parameters for CSSR algorithm, described in the following paper.

http://www.sciencedirect.com/science/article/pii/S0167865509001597

Here is an extension that allows to determine the statistical significance level of the inferred model.

http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=6193099&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D6193099

	By implementing all of these algorithms and techniques the framework will be able to find all the statistically significant patterns from discrete data time series with a minimal input and bias from a user.

