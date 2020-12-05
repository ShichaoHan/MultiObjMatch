
library(dplyr)
library(optmatch)
library(rcbalance)
library(matchMulti)
library(cobalt)
# library(vioplot)
library(cobalt)

rho_proposition <- function(paircosts.list, rho.max.factor=10, rho1old, rho2old){
  rho.min <- 1e-2
  max.dist <- max(paircosts.list)
  min.dist <- min(paircosts.list)

  rho1 <- c(0, rho.min,0.5,1,2,
            4,5, 8,
            max(min.dist, quantile(paircosts.list)[2]), quantile(paircosts.list)[4], rho.max.factor*max.dist)
  rho2 <- seq(0, max(rho1), max(rho1) / 5 )

  if(is.null(rho1old) & is.null(rho2old)){
    result <- generate_rhos(rho1, rho2)
  } else if(is.null(rho1old)){
    result <- generate_rhos(rho1, rho2old)
  } else if(is.null(rho2old)){
    result <- generate_rhos(rho1old, rho2)
  } else{
    result <- generate_rhos(rho1old, rho2old)
  }
  return(result)

}

generate_rhos <- function(rho1.list, rho2.list){
  result <- list()
  ind = 1
  for(i in 1:length(rho1.list)){
    for(j in 1:length(rho2.list)){
      result[[ind]] <- c(rho1.list[i], rho2.list[j])
      ind = ind + 1
    }
  }
  return(result)

}


obj.to.match <- function(out.elem, already.done = NULL, prev.obj = NULL){
  tcarcs <- length(unlist(out.elem$net$edgeStructure))
  edge.info <- extractEdges(out.elem$net)
  one.sol <- function(sol){
    x <- sol[1:tcarcs]
    match.df <- data.frame(treat = as.factor(edge.info$startn[1:tcarcs]), x = x, control = edge.info$endn[1:tcarcs])
    matched.or.not <- daply(match.df, .(match.df$treat), function(treat.edges) c(as.numeric(as.character(treat.edges$treat[1])),
                                                                                 sum(treat.edges$x)), .drop_o = FALSE)
    if (any(matched.or.not[, 2] == 0)) {
      match.df <- match.df[-which(match.df$treat %in% matched.or.not[which(matched.or.not[,
                                                                                          2] == 0), 1]), ]
    }
    match.df$treat <- as.factor(as.character(match.df$treat))
    matches <- as.matrix(daply(match.df, .(match.df$treat), function(treat.edges) treat.edges$control[treat.edges$x ==
                                                                                                        1], .drop_o = FALSE))
    matches - length(out.elem$net$treatedNodes)
  }
  if(is.null(already.done)) return(llply(out.elem$solutions, one.sol))
  new.ones <- setdiff(1:length(out.elem$solutions), already.done)
  out.list <- list()
  out.list[already.done] <- prev.obj
  out.list[new.ones] <- llply(out.elem$solutions[new.ones], one.sol)
  return(out.list)
}



#' Fit propensity scores using the given covariates
#'
#' @param df dataframe that contains a column named "treat", the treatment vector, and columns of covariates specified.
#' @param covs factor of column names of covariates used for fitting a propensity score model.
#'
#' @return factor of fitted propensity scores.
#' @export
getPropensityScore <- function(df, covs){
  pscoreModel <- glm(treat ~.,
                     data = df[c('treat', covs)], family = binomial("logit"))
  return(predict(pscoreModel, type="response"))
}


#' Get a factor of string representing the exact variables.
#'
#' @param dat dataframe containing all the variables in exactList
#' @param exactList factor of names of the variables on which we want exact or close matching.
#'
#' @return factor of concatenated variable values for the variables we want exact or close matching on.
#' @export
getExactOn <- function(dat, exactList){
  if(length(exactList) == 1){
    exactOn = dat[,exactList[1]]
  } else {
    exactOn = paste(dat[,exactList[1]],dat[,exactList[2]])
    if(length(exactList) >=3){
      for(i in 3:length(exactList)){
        exactOn = paste(exactOn, dat[,exactList[i]])
      }
    }
  }
  return(exactOn)
}

#' The main matching function that can return a matching object
#'
#' @param df data frame that contain columns indicating treatment, outcome and covariates
#' @param treatCol character of name of the column indicating treatment assignment
#' @param responseCol character of name of the column indicating the outcome variable
#' @param distList factor of the names of the variables used for calculating within-pair distance
#' @param exactlist factor of the names of the variables that we want exact matching on
#' @param myBalCol charactor of column name of the variable that we want to evaluate balance on
#' @param rho1 factor of values of rho1. Default value is c(1).
#' @param rho2 factor of values of rho2. Default value is c(1,2,3)
#' @param propensityCols factor of names of columns that users want to fit a propensity score model
#' @param pScores character of the name of the column that indicate the propensity score; default is NULL
#' @param idCol character of the name of the column that indicate the id for each unit; default is NULL
#' @param maxUnMatched double of the maximum proportion of unmatched unit that can be accepted; default is 0.25
#' @param caliperOption double of the caliper value; default is 0.25
#' @param toleranceOption double of tolerance of close match distance; default is 1e-2
#' @param maxIter interger of the maximum number of iterations to search for (rho1, rho2) pair that improve the matching
#' @param rho.max.f double of the scaling factor used in proposal for rhos; default is 10
#'
#' @return a list whose elements are (1) "rhoList": list of rhos for each match (2) "matchList": list of matches indexed by number (3) "treatmentCol": character of treatment variable (4) "covs": factor of names of the variables used for calculating within-pair distance (5) "exactCovs": factor of names of variables that we want exact or close match on
#' (6) "idMapping": factor of row index for each observation in the sorted data frame for internal use (6) "stats": data frame of important statistics (total variation distance) for variable on which marginal balance is measured (7)"b.var": character of the variable on which marginal balance is measured (8) "dataTable": data frame sorted by treatment value
#' (10) "df": data frame of input data (11) "pair_cost": list of pair-wise distance sum for each match
multiObjMatch <- function(df, treatCol, responseCol, distList, exactlist, myBalCol, rho1=c(1), rho2=c(1,2,3),propensityCols = NULL, pScores = NULL, idCol = NULL, maxUnMatched = 0.25, caliperOption=0.25, toleranceOption=1e-2, maxIter=0, rho.max.f = 10){
  ## 0. Data preprocessing
  dat = df
  dat$originalID = 1:nrow(df)
  dat$tempSorting = df[,treatCol]
  dat = dat[order(-dat$tempSorting),]
  rownames(dat) <- as.character(1:nrow(dat))
  columnNames <- colnames(df)
  xNames <- columnNames[!columnNames %in% c(treatCol, responseCol, "originalID", "tempSorting")]

  dat['response'] = df[responseCol]
  dat['treat'] = df[treatCol]

  result <- structure(list(), class = 'multiObjMatch')

  ## 1. Fit a propensity score model if the propensity score is not passed in
  if(is.null(pScores)){
    if(is.null(propensityCols)){
      pscore = getPropensityScore(dat,xNames)
    } else {
      pscore = getPropensityScore(dat, propensityCols)
    }
  } else{

    pscore = as.numeric(list(dat[pScores]))
  }


  ## 2. Construct nets
  base.net <- netFlowMatch(dat$treat)

  ## 3. Define distance and cost


  dist.i <- build.dist.struct(
    z = dat$treat, X = dat[distList],
    exact = getExactOn(dat, exactlist), calip.option = "user",
    calip.cov = pscore, caliper = caliperOption)




  net.i <- addExclusion(makeSparse(base.net, dist.i))
  my.bal <- apply(dat[,match(myBalCol, colnames(dat)), drop = FALSE], 1, function(x) paste(x, collapse = "."))
  net.i <- addBalance(net.i, treatedVals =
                        my.bal[dat$treat==1], controlVals = my.bal[dat$treat==0])
  paircost.v <- unlist(dist.i)
  rho.max.factor = rho.max.f
  ### Create edge costs for each of the objective functions we will trade off
  pair.objective.edge.costs <- pairCosts(dist.i, net.i)
  excl.objective.edge.costs <-  excludeCosts(net.i,1)
  bal.objective.edge.costs <- balanceCosts(net.i, max(paircost.v)*rho.max.factor)

  f1list = list(pair.objective.edge.costs)
  f2list = list(excl.objective.edge.costs)
  f3list = list(bal.objective.edge.costs)


  ## Propose possible rho values for multi-objective optimization problem
  rho_list <- list()
  rho_list <- rho_proposition(paircost.v, rho.max.f, rho1, rho2)

  # #### EXPENSIVE PART ####
  #
  solutions <- list()
  solution.nets <- list()
  rho_counts = 1
  for(rho in rho_list){

    temp.sol <- solveP1(net.i,
                         f1list, f2list, f3list, rho1 = rho[1], rho2 = rho[2], tol = toleranceOption)

    solutions[[as.character(rho_counts)]] <- temp.sol$x
    solution.nets[[as.character(rho_counts)]] <- temp.sol$net
    print(paste('Matches finished for rho1 = ', rho[1], " and rho2 = ", rho[2]))
    rho_counts = rho_counts + 1
  }

  match.list <- obj.to.match(list('net' = net.i, 'costs' = pair.objective.edge.costs, 'balance' = bal.objective.edge.costs, 'rho.v' = 1:length(rho_list), 'solutions' = solutions))
  match.list <- match.list[order(as.numeric(names(match.list)))]

  ## remove the matching results that result in zero matches
  temp.length <- laply(match.list, nrow)
  temp.names <- names(match.list)
  nonzero.names <- temp.names[temp.length != 0]
  zero.names <- temp.names[temp.length == 0]

  if(length(zero.names) == length(temp.names)){
    stop('Error: Initial rho values result in non-matching. Please try a new set of initial rho values.')
  }
  match.list <- match.list[names(match.list) %in% nonzero.names == TRUE]


  ## Iterative Search for Best Rhos (if necessary)
  numIter = 1


  tempRhos <- as.numeric(names(match.list))
  percentageUnmatched <- 1 - as.vector(laply(match.list, nrow) / sum(dat[,treatCol]))
  minInd <- which.min(percentageUnmatched)
  bestRhoInd<- tempRhos[minInd]
  bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
  bestRho1 <- rho_list[[bestRhoInd]][1]
  bestRho2 <- rho_list[[bestRhoInd]][2]
  ind = length(rho_list) + 1
  while((bestPercentageSoFar > maxUnMatched) & (numIter <= maxIter)){

    if(is.finite(log10(bestRho1)) & is.finite(log10(bestRho2))){
    proposedNewRho1s <- c(0.5* (-c((abs(rnorm(1))* 0.01 + bestRho1/10):floor(log10(bestRho1)))), 0.5*(-c(abs(rnorm(1))* 0.01 +bestRho1:log10(bestRho1 + 0.01* abs(rnorm(1)))*10)))

    proposedNewRho2s <- c(0.5* (-c((abs(rnorm(1))* 0.01 + bestRho2/10):floor(log10(bestRho2)))), 0.5*(-c(abs(rnorm(1))* 0.01 +bestRho2:log10(bestRho2+ 0.01* abs(rnorm(1)))*10)))
    } else {
      proposedNewRho1s <- c(rnorm(1) * 1e-1, rnorm(1) * 1e2, rnorm(1) * 1e3)
      proposedNewRho2s <- c(rnorm(1) * 1e-1, rnorm(1) * 1e2, rnorm(1) * 1e3)
     }


    for(r1 in proposedNewRho1s){
      rho_list[[ind]] = c(r1, bestRho2)
      temp.sol <- solveP1(net.i,
                          f1list, f2list, f3list, rho1 = r1, rho2 = bestRho2, tol = toleranceOption)

      solutions[[as.character(ind)]] <- temp.sol$x
      solution.nets[[as.character(ind)]] <- temp.sol$net
      print(paste('Matches finished for rho1 = ', r1, " and rho2 = ", bestRho2))
      ind = ind + 1
    }

    for(r2 in proposedNewRho2s){
      rho_list[[ind]] = c(bestRho1, r2)
      temp.sol <- solveP1(net.i,
                          f1list, f2list, f3list, rho1 = bestRho1, rho2 = r2, tol = toleranceOption)

      solutions[[as.character(ind)]] <- temp.sol$x
      solution.nets[[as.character(ind)]] <- temp.sol$net
      print(paste('Matches finished for rho1 = ', bestRho1, " and rho2 = ", r2))
      ind = ind + 1
    }

    match.list <- obj.to.match(list('net' = net.i, 'costs' = pair.objective.edge.costs, 'balance' = bal.objective.edge.costs, 'rho.v' = 1:length(solutions), 'solutions' = solutions))
    match.list <- match.list[order(as.numeric(names(match.list)))]

    ## remove the matching results that result in zero matches
    temp.length <- laply(match.list, nrow)
    temp.names <- names(match.list)
    nonzero.names <- temp.names[temp.length != 0]
    match.list <- match.list[names(match.list) %in% nonzero.names == TRUE]


    tempRhos <- as.numeric(names(match.list))
    percentageUnmatched <- 1 - as.vector(laply(match.list, nrow) / sum(dat[,treatCol]))
    oldMinInd <- minInd
    minInd <- which.min(percentageUnmatched[-oldMinInd])
    if(minInd >= oldMinInd){
      minInd = minInd + 1
    }
    bestRhoInd<- tempRhos[minInd]
    bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
    if(sum(percentageUnmatched==bestPercentageSoFar) > 1){
      bestRhoInd <- sample(which(percentageUnmatched == bestPercentageSoFar), size = 1)
    }
    bestRho1 <- rho_list[[bestRhoInd]][1]
    bestRho2 <- rho_list[[bestRhoInd]][2]

    numIter = numIter + 1
  }

  my.stats1 <- t(laply(match.list, descr.stats_general, df=dat, treatCol = treatCol, b.vars = myBalCol, pair.vars = distList, extra = TRUE))

  solutions.old <- solutions
  solution.nets.old <- solution.nets

  pair_cost_sum <- c()
  for(ind in names(match.list)){
    x <- solutions[[ind]]
    pair_cost_sum[[ind]] <- sum(paircost.v*x[1:length(paircost.v)])
  }

  result$rhoList <- rho_list
  result$matchList <-  match.list
  ## Store some information about treatment column, balance covariates and exact match column
  result$treatmentCol = treatCol
  result$covs = distList
  result$exactCovs = exactlist
  result$idMapping = dat$originalID
  result$stats = my.stats1
  result$b.var = myBalCol
  result$dataTable = dat
  result$t = treatCol
  result$df = df
  result$pair_cost = pair_cost_sum
  return(result)
}


generate_balance_table <- function(matchingResult, covList=NULL, display.All=FALSE, statList=c("mean.diffs")){
  originalDF <- matchingResult$dataTable

  numTreated = sum(originalDF$treat)
  percentageUnmatched <-  1- as.vector(laply(matchingResult$matchList, nrow) / numTreated)
  sortedMatchList <- matchingResult$matchList[order(percentageUnmatched)]

  matchIndex <- names(sortedMatchList)
  if(display.All==FALSE & length(matchIndex) >= 5){
    matchIndex <- matchIndex[as.integer(quantile(1:length(matchIndex)))]
  }
  balanceTable <- list()

  for(ind in matchIndex){

    matchTab <- matchingResult$matchList[[ind]]
    treatedUnits <- as.numeric(rownames(matchTab))
    controlUnits <- as.vector(matchTab[,1]) + numTreated
    if(!is.null(covList)){
      covariates <- originalDF[c(treatedUnits, controlUnits), covList]
    } else{
      covariates <- originalDF[c(treatedUnits, controlUnits), matchingResult$covs]
    }
    balanceTable[[ind]] <- bal.tab(covariates, s.d.denom= "pooled", treat = as.vector(originalDF[c(treatedUnits, controlUnits), 'treat']), stats = statList)
  }
  return(balanceTable)
}


compare_tables <- function(balanceTable){
  inds <- names(balanceTable)
  rnames <- rownames(balanceTable[[inds[1]]]$Balance)
  result <- data.frame(balanceTable[[inds[1]]]$Balance[, 1])
  colnames(result) <- c("type")
  count = 1
  for(ind in inds){
    result[ind] <- balanceTable[[ind]]$Balance[rnames,2]
    count = count + 1
  }
  rownames(result) <- rnames
  return(result)
}

#' Generate a table that compare the covariate balance across different possible matches
#'
#' @param matchingResult an object returned by the main matching function multiObjMatch
#' @param covList factor of names of covariates that we want to evaluate covariate balance on; default is NULL. When set to Null, the program will compare the covariates that have been used to construct a propensity model
#' @param display_option boolean value of whether to display all the matches; default is FALSE, where matches at each quantile is displayed
#' @param stat character of the name of the statistic used for measuring covariate balance; default is "mean.diff". This argument is the same as used in "cobalt" package
#'
#' @return a dataframe that shows covariate balance in different matches
#' @export
compare_matching <- function(matchingResult, covList=NULL, display_option=FALSE, stat="mean.diff"){
  return(compare_tables(generate_balance_table(matchingResult, covList, display.All=display_option, c(stat))))
}


get_five_index <- function(matchingResult){
  df <- compare_matching(matchingResult)
  matches <- colnames(df)[2:length(colnames(df))]
  return(matches)

}




#' Plotting function that generate the total variation imbalance v.s. number of unmatched treated units
#'
#' @param matchingResult an object returned by the main matching function multiObjMatch
#'
#' @return NULL
#' @export
generate_tv_graph <- function(matchingResult){
  inds <- get_five_index(matchingResult)
  graph_labels <- names(matchingResult$matchList)
  for(i in 1:length(graph_labels)){
    if(sum(as.numeric(inds) == as.numeric(graph_labels[i])) == 0){
      graph_labels[i] = " "
    }
  }
  treatedSize = sum(matchingResult$dataTable$treat)
  samp.size <- laply(matchingResult$matchList, nrow)
  f1 <- matchingResult$stats[1,] * samp.size
  f2 <- treatedSize  - samp.size
  plot(f2,f1, pch = 20, xlab = 'Treated Units Unmatched', ylab = 'Total Variation Imbalance', cex = 1.2, xlim = c(0, treatedSize), ylim = c(0,max(f1)), cex.lab = 1.2, cex.axis= 1.2, col = rgb(0,0,0,0.3))
  abline(h=0, lty = 'dotted')
  points(f2, f1, col = rgb(0,0,1,1),pch = 20)
  text(f2, f1, labels = graph_labels, pos=2)
}

#' Plotting function that generate sum of pair-wise distance v.s. number of unmatched treated units
#'
#' @param matchingResult an object returned by the main matching function multiObjMatch
#'
#' @return NULL
#' @export
generate_pairdistance_graph <- function(matchingResult){
  inds <- get_five_index(matchingResult)
  graph_labels <- names(matchingResult$matchList)
  for(i in 1:length(graph_labels)){
    if(sum(as.numeric(inds) == as.numeric(graph_labels[i])) == 0){
      graph_labels[i] = " "
    }
  }
  treatedSize = sum(matchingResult$dataTable$treat)
  samp.size <- laply(matchingResult$matchList, nrow)
  f1 <- matchingResult$pair_cost
  f2 <- treatedSize  - samp.size
  plot(f2,f1, pch = 20, xlab = 'Treated Units Unmatched', ylab = 'Pair-wise Distance Sum', cex = 1.2, xlim = c(0, treatedSize), ylim = c(0,max(f1)), cex.lab = 1.2, cex.axis= 1.2, col = rgb(0,0,0,0.3))
  abline(h=0, lty = 'dotted')
  points(f2, f1, col = rgb(0,0,1,1),pch = 20)
  text(f2, f1, labels = graph_labels, pos=2)
}

#' Plotting function that generate sum of pair-wise distance v.s. total variation imbalance on balance variable
#'
#' @param matchingResult an object returned by the main matching function multiObjMatch
#'
#' @return NULL
#' @export
generate_pairdistance_balance_graph <- function(matchingResult){
  inds <- get_five_index(matchingResult)
  graph_labels <- names(matchingResult$matchList)
  for(i in 1:length(graph_labels)){
    if(sum(as.numeric(inds) == as.numeric(graph_labels[i])) == 0){
      graph_labels[i] = " "
    }
  }
  treatedSize = sum(matchingResult$dataTable$treat)
  samp.size <- laply(matchingResult$matchList, nrow)
  f1 <- matchingResult$pair_cost
  f2 <- matchingResult$stats[1,] * samp.size
  plot(f2,f1, pch = 20, xlab = 'Total Variation Imbalance on Balance Variable', ylab = 'Pair-wise Distance Sum', cex = 1.2, xlim = c(0, max(f2)), ylim = c(0,max(f1)), cex.lab = 1.2, cex.axis= 1.2, col = rgb(0,0,0,0.3))
  abline(h=0, lty = 'dotted')
  points(f2, f1, col = rgb(0,0,1,1),pch = 20)
  text(f2, f1, labels = graph_labels, pos=2)
}

convert_index <- function(matchingResult){
  idMap <- matchingResult$idMapping
  numTreated = sum(matchingResult$dataTable$treat)
  matchIndex <- names(matchingResult$matchList)
  result <- list()

  for(ind in matchIndex){

    matchTab <- matchingResult$matchList[[ind]]
    treatedUnits <- as.numeric(rownames(matchTab))
    controlUnits <- as.vector(matchTab[,1]) + numTreated
    treatedOriginalID <- idMap[treatedUnits]
    controlOriginalID <- idMap[controlUnits]
    tempMatch <- list()
    tempMatch$treated = treatedOriginalID
    tempMatch$control = controlOriginalID
    result[[ind]] <- tempMatch
  }
  return(result)
}

matched_index <- function(matchingResult){
  return(convert_index(matchingResult))
}


#' A function that returns the dataframe that contains only matched pairs from the original data frame with specified match index
#'
#' @param matchingResult an object returned by the main matching function multiObjMatch
#' @param match_num Integer index of match that the user want to extract paired observations from
#'
#' @return dataframe that contains only matched pair data
#' @export
matched_data <- function(matchingResult, match_num){
  matched_idx <- matched_index(matchingResult)
  treated_idx <- matched_idx[[as.character(match_num)]]$treated
  control_idx <- matched_idx[[as.character(match_num)]]$control
  return(matchingResult$df[c(treated_idx, control_idx),])
}



#' A function that generate the percentage of unmatched units for each match.
#'
#' @param matchingResult matchingResult object that contains information for all matches
#'
#' @return data frame that contain number of matched units as well as the percentage of matched units
#' @export
#'
#' @examples
getUnmatched <- function(matchingResult){
  matchingIndex <- names(matchingResult$matchList)
  matchedUnits <- laply(matchingResult$matchList, nrow)
  matchedPercentage <- matchedUnits / sum(matchingResult$dataTable$treat)
  result <- data.frame(matchingIndex, matchedUnits, matchedPercentage)
  colnames(result) <- c("Matching Index", "Number of Matched Units", "Percentage of Matched Units")
  return(result)

}

ate_calculate <- function(df, responseColName){
  calculate_ate <- function(matches){
    treatedInd <- matches$treated
    controlInd <- matches$control
    return(mean(df[treatedInd, responseColName]) - mean(df[controlInd, responseColName]))
  }

}

getATE <- function(df, responseCol, matchingResult){
  converted_ind = convert_index(matchingResult)
  ateList = laply(converted_ind, ate_calculate(df, responseCol))
  matchingIndex <- names(matchingResult$matchList)
  result <- data.frame(matchingIndex, ateList)
  colnames(result) <- c("Matching Index", "ATE")
  return(result)
}














































