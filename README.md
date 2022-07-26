# eiViaPBApprox
Code for estimating individual-level voter preferences using aggregate election data, making use of a normal approximation to the Poisson Binomial (see Rosenman &amp; Viswanathan, 2018). 

These scripts give a simple method for fitting ecological inference models using the method described in Rosenman & Viswanathan, 2018 and Rosenman, 2019. The training data is obtained at the voting precinct level. With each precinct, we associate a model matrix that contains the relevant covariates for every individual who participated in the election, along with the total number of votes received by a major party candidate in that precinct. The given method is for learning binary classifiers, so we assume we are dealing with a two-candidate election. 

Our goal is to learn a vector, $\beta$, which relates the individual-level covariates 