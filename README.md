# eiViaPBApprox
Code for estimating individual-level voter preferences using aggregate election data, making use of a normal approximation to the Poisson Binomial. 

These scripts give a simple method for fitting ecological inference models using the method described in Rosenman & Viswanathan (2018) and Rosenman (2019). The training data is obtained at the voting precinct level. With each precinct, we associate a model matrix that contains the relevant covariates for every individual who participated in the election, along with the total number of votes received by a major party candidate in that precinct. The given method is for learning binary classifiers, so we assume we are dealing with a two-candidate election. 

Our goal is to learn a vector, $\beta$, which relates the individual-level covariates to the probability of voting for the given candidate via the standard logistic regression formulation. For individual $i$ with covariates $x_i$, denote as $D_i$ the indicator that individual $i$ voted for the Democratic candidate in the election. Our model is thus

$$ P(D_i = 1) = \frac{1}{1 + exp(- \beta^T x_i)}.$$

Dummy data is provided in the repo in the R data file "trainData.rData". This demonstrates the format needed for training the model. In "Fit EI Model.R", we provide a very simple script for training the model via Newton's method. Note that the log likelihood need not be concave, so this procedure can fail if there are too many covariates included in the model, or too few precincts from which to learn. We have found that in reasonably large settings (~1,000 precincts, ~20-50 covariates), the models tend to converge quickly. 

For details about how we have approximated the true Poisson Binomial likelihood with a Gaussian, see Rosenman & Viswanathan (2018) and Rosenman (2019). The file "Helper Functions.R" contains the relevant functions for computing the approximate gradient and Hessian of the likelihood at each value of $\beta$. 