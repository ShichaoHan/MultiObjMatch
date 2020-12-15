# MultiObjMatch: Multi-objective Matching for R


The `MultiObjMatch` offers a user-friendly R package that implements matching of two groups of subjects to treated and control groups in observational studies. This package allows users to form matches that achieve a specified balance among the three objectives: the number of treated units matched, the total variation imbalance on the marginal distribution of key categorical variables, and sum of within-pair distance,  Researchers are allowed to form matches that meet user-specified design goals for matching problems in observational studies. More detailed discussion can be found in [Pimentel and Kelz (2020)](https://www.tandfonline.com/doi/pdf/10.1080/01621459.2020.1720693?casa_token=ubLCKouur94AAAAA:LiiihGbqOwfWHhb2UpxZYXqpKsCQWPB5u8OgyjETfIq9ucrM5OIgLq_OZWYz2DgEc2wxIWrWAoMq).    
  
Besides the main matching algorithm, the package also contains useful functions for generate numeric and graphical diagnostics. 
   
   
## Using MultiObjMatch  
  
This section provides a guide to the package usage. Before running the example, users can install the package from github.  


```r
library(devtools)
install_github("ShichaoHan/MultiObjMatch", ref="main")
library(MultiObjMatch)
```   
  
In the demo below, the dataset "lalonde" is loaded from the package MatchIt.   
  
```r
library(cobalt)
data("lalonde", package="cobalt")
```

  

### Matching 

After data pre-processing, users can use the main matching function __multiObjMatch__ to generate a set of possible matches. Users need to specify the input data frame, the name of treatment variable, the name of the outcome variable, the list of columns for measuring pairwise distance, the column for exact matching, the column for measuring the balance for, and the list of  the coefficients for objective functions.    
```r
psCols <- c("age", "educ", "married", "nodegree")
treatVal <- "treat"
responseVal <- "re78"  
pairDistVal <- c("age", "married", "educ", "nodegree")
exactVal <- c("educ") 
myBalVal <- c("race")
r1s <- c(0.01, 0.02, 0.03, 0.04, 0.1, 0.3, 0.5, 0.7, 0.9)
r2s <- c(1)
matchResult <- multiObjMatch(lalonde, treatVal, responseVal, pairDistVal, 
    exactVal, myBalVal, rho1=r1s, rho2=r2s, propensityCols = psCols, 
    pScores = NULL, idCol = NULL, maxUnMatched = 0.1, caliperOption=0.25, 
    toleranceOption=1e-1, maxIter=3, rho.max.f = 10)
```
   
### Numeric Diagnostics  
  
Users can use the main function __compare_matching__ on specified covariates to compare the covariate balance across different matching. 
```r

compare_matching(matchResult, covList=c("age", "educ", "race", "married", "nodegree"))
```  
  
The number of matched units and percentage of matched units can be automatically generated using the helper function __getUnmatched__:  

```r
getUnmatched(matchResult)
```
  
### Graphical Diagnostics  
  
There are three helper functions that generate the graphical diagnostics. 

```r
generate_tv_graph(matchResult)
generate_pairdistance_graph(matchResult)
generate_pairdistance_balance_graph(matchResult)

```
  
### Get Matched Data   
  
Users can use the helper function __matched_data__ to obtain the dataframe containing only the matched treated and control pairings by passing in the result from the main matching function and the index of the match.     

```r
matched_data(matchResult, 13)
```  
  
Then, users can use outcome analysis of their choice upon the matched data. 
