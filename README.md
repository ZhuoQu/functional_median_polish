# functional_median_polish
This documents include the one-way and two-way median polish methods for univariate and multivariate functional data

### please run test_Australia.R directly. 
The required function is loaded from the start.
one-way median polish has been tested by the mortality rate of Australian females (univariate functional data)
and the mortality rate of Australian females and males jointly (bivariate functional data).


## If we use different years as a category effect, then each subgroup only contains one observation.
In the one-way model,
according to Y_{i,j}=m+a_{i} + epsilon_{i,j}, epsilon_{i,j}=0. 
