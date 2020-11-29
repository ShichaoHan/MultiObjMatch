# MultiObjMatch: Multi-objective Matching for R


The `MultiObjMatch` offers a user-friendly R package that implements matching of two groups of subjects to treated and control groups in observational studies. With user-specified emphasis on three design goals: maximizing closeness in pair-wise measured by Mahalanobis distances, maximizing the number of matched treated units, and minimizing the distance of balance variable's distributions in both groups as measured by total variation distance.    
  
Besides the main matching algorithm, the package also contains useful functions for generate numeric and graphical diagnostics. 
   
   
## Using MultiObjMatch  
  
This section provides a guide to the package usage. Before running the example, users can install the package from github.  


```r
library(devtools)
install_github("ShichaoHan/MultiObjMatch")
library(MultiObjMatch)
```   
  
In the demo below, the dataset "lalonde" is loaded from the package MatchIt.   
  
```r
library(MatchIt)
data("lalonde")
```

  

### Matching 

After data pre-processing, users can use the main matching function __multiObjMatch__ to generate a set of possible matches. Users need to specify the input data frame, the name of treatment variable, the name of the outcome variable, the list of columns for measuring pairwise distance, the column for exact matching, the column for measuring the balance for, and the list for $$\rho_1$$ and $$\rho_2$$ -  the coefficients for $$f_1$$ and $$f_2$$.    
```r
psCols <- c("age", "educ", "married", "nodegree")
treatVal <- "treat"
responseVal <- "re78"  
pairDistVal1 <- c("age","married","educ", "nodegree")
exactVal <- c("educ") 
myBalVal <- c("race")
r1s <- c(0.01,0.02,0.03,0.1,0.5,1,1.5,2,2.5,4,5,6,7,8,9,10)
r2s <- c(1)
matchResult <- multiObjMatch(lalonde, treatVal, responseVal, pairDistVal, 
    exactVal, myBalVal, rho1=r1s, rho2=r2s, propensityCols = psCols, 
    pScores = NULL, idCol = NULL, maxUnMatched = 0.1, caliperOption=0.25, 
    toleranceOption=1e-1, maxIter=3, rho.max.f = 10)
```
   
### Numeric Diagnostics  
  
Users can use the main function __compare_matching__ on specified covariates to compare the covariate balance across different matching. 
```r

compare_matching(mathResult, covList=c("age", "educ", "race", "married", "nodegree"))
```  
  
The number of matched units and percentage of matched units can be automatically generated using the helper function __getUnmatched__:  

```r
getUnmatched(matchResult)
```


