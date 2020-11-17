require(plyr)
require(gtools)

#assumes table is built control x treated with zero outcomes in upper left and events in bottom right.
test.ae <- function(tab2x2, a = 0, gam = 1, null.less = TRUE, exact = TRUE){
	#if attributable effect is negative, transform problem so it is positive
	if(a < 0){
		tab2x2 <- t(tab2x2)
		null.less <- !null.less
		a <- -a
	}
	#if a is too large retain null immediately
	if(null.less){
		if(a > sum(tab2x2[,2])) return(1)

		#shift table to account for a if necessary
		shifttab <- tab2x2
		if(a > 0){
			shift1.1 <- min(a,tab2x2[2,2])
			shifttab[2,2] <- shifttab[2,2] - shift1.1
			shifttab[2,1] <- shifttab[2,1] + shift1.1
			if(a > tab2x2[2,2]){
				shifttab[1,2] <- shifttab[1,2] - (a - shift1.1)
				shifttab[1,1] <- shifttab[1,1] + (a - shift1.1)
			}
		}
		#set values for 2x2 table
		lambda_tab <- rbind(c(0,gam/(1 + gam)),c(gam/(1 + gam),1))

		#calculate test statistic
		if(exact) 	return(1 - pbinom(shifttab[1,2],shifttab[2,1]+shifttab[1,2],lambda_tab[1,2]))
		test.stat <- (sum(shifttab[,2]) - sum(lambda_tab*shifttab))/sqrt(sum(shifttab*lambda_tab*(1 - lambda_tab)))
		#test.stat <- (sum(shifttab[,2]) - a - sum(lambda_tab*shifttab))/sqrt(sum(shifttab*lambda_tab*(1 - lambda_tab)))
		return(1 - pnorm(test.stat))
	}
	if(!null.less){
		if(a > sum(tab2x2[,2])) return(0)
		
		#shift table to account for a if necessary
		shifttab <- tab2x2
		if(a > 0){
			shift1.1 <- min(a,tab2x2[2,2])
			shifttab[2,2] <- shifttab[2,2] - shift1.1
			shifttab[2,1] <- shifttab[2,1] + shift1.1
			if(a > tab2x2[2,2]){
				shifttab[1,2] <- shifttab[1,2] - (a - shift1.1)
				shifttab[1,1] <- shifttab[1,1] + (a - shift1.1)
			}
		}
		#set values for 2x2 table
		lambda_tab <- rbind(c(0,1/(1 + gam)),c(1/(1 + gam),1))

		#calculate test statistic
		if(exact) return(pbinom(shifttab[1,2],shifttab[2,1]+shifttab[1,2],lambda_tab[1,2]))
		test.stat <- (sum(tab2x2[,2]) - a - sum(lambda_tab*shifttab))/sqrt(sum(shifttab*lambda_tab*(1 - lambda_tab)))
		return(pnorm(test.stat))
	}
}


#what do I want to produce?
#2x2 table, confidence interval for attributable effect.
conf.ae <- function(tab2x2, gam = 1, searchlims = c(-10000,10000), coverage = 0.95){
	target.qt <- (1 - coverage)/2
	ll <- binsearch(function(a){test.ae(tab2x2, a = a, gam=gam,null.less=TRUE)}, range = searchlims, target = target.qt)$where
	if(length(ll) > 1) ll <- ll[2]
	ul <- binsearch(function(a){test.ae(tab2x2, a = a, gam=gam,null.less=FALSE)}, range = c(ll,searchlims[2]), target = target.qt)$where
	if(length(ul) > 1) ul<- ul <- ul[1] 
	return(c(ll,ul))
}

#get a Hodges-Lehmann point estimate
#for now require gamma = 1; this becomes an interval with sensitivity analysis but I shouldn't need that for the current paper.
pt.ae <- function(tab2x2, gam = 1, searchlims = c(-10000,10000)){
	#need to  find a for which test statistic equals null expectation; equivalent to a one-sided p-value of 0.5
	target.qt <- 0.5
	hl <- binsearch(function(a){test.ae(tab2x2, a = a, gam=gam,null.less=TRUE)}, range = searchlims, target = target.qt)$where
	if(length(hl) > 1) hl <- sum(hl)/2
	if(gam > 1){
		hl2 <- binsearch(function(a){test.ae(tab2x2, a = a, gam=gam,null.less=FALSE)}, range = searchlims, target = target.qt)$where
		if(length(hl2) > 1) hl2 <- sum(hl2)/2
		hl <- (hl + hl2)/2
	}
	return(hl)
}


eq.test <- function(tab2x2, iota, gam = 1){
		pval1 <- test.ae(tab2x2, gam = gam, a = iota, null.less = FALSE)
		pval2 <- test.ae(tab2x2, gam = gam, a = -iota, null.less = TRUE)		
		max(pval1, pval2)
}







