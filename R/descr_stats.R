
dummy <- function (x, data = NULL, sep = "", drop = TRUE, fun = as.integer,
                   verbose = FALSE, name = NULL)
{
  if (is.null(data)) {
    if(is.null(name)) name <- as.character(sys.call(1))[2]
    name <- sub("^(.*\\$)", "", name)
    name <- sub("\\[.*\\]$", "", name)
  }
  else {
    if (length(x) > 1)
      stop("More than one variable provided to produce dummy variable.")
    name <- x
    x <- data[, name]
  }
  if (drop == FALSE && class(x) == "factor") {
    x <- factor(x, levels = levels(x), exclude = NULL)
  }
  else {
    x <- factor(x, exclude = NULL)
  }
  if (length(levels(x)) < 2) {
    if (verbose)
      warning(name, " has only 1 level. Producing dummy variable anyway.")
    return(matrix(rep(1, length(x)), ncol = 1, dimnames = list(rownames(x),
                                                               c(paste(name, sep, x[[1]], sep = "")))))
  }
  mm <- model.matrix(~x - 1, model.frame(~x - 1), contrasts = FALSE)
  colnames.mm <- colnames(mm)
  if (verbose)
    cat(" ", name, ":", ncol(mm), "dummy varibles created\n")
  mm <- matrix(fun(mm), nrow = nrow(mm), ncol = ncol(mm), dimnames = list(NULL,
                                                                          colnames.mm))
  colnames(mm) <- sub("^x", paste(name, sep, sep = ""), colnames(mm))
  if (!is.null(row.names(data)))
    rownames(mm) <- rownames(data)
  return(mm)
}


descr.stats <- function(matches, b.vars, pair.vars, extra = FALSE){
    treatment.status <- c(rep(1, nrow(matches)), rep(0,
        nrow(matches)))

	#total variation imbalance
    matched.info <- rbind(my.slice[my.slice$img == 1,][as.numeric(rownames(matches)),
        , drop = FALSE], my.slice[my.slice$img == 0,][as.vector(matches), , drop = FALSE])
    interact.factors.matched = as.factor(apply(matched.info[,
        match(b.vars, colnames(matched.info)), drop = FALSE],
        1, function(x) paste(x, collapse = ".")))
	fb.tab <- table('balance.variable' = interact.factors.matched, treatment.status)
	tv.sum <- sum(abs(fb.tab[,1] - fb.tab[,2]))
  tv.prop <-  sum(abs(fb.tab[,1] - fb.tab[,2]))/sum(fb.tab)
	names(tv.sum) <- paste(c('TV',b.vars), collapse = '.')
	names(tv.prop) <- paste(c('TV',b.vars), collapse = '.')

	#chi-squared test
	chisq.p <-  chisq.test(fb.tab)$p.value
	names(chisq.p) <- paste(c('Chisq',b.vars), collapse = '.')

	interact.pv = as.character(apply(matched.info[,
			match(pair.vars, colnames(matched.info)), drop = FALSE],
			1, function(x) paste(x, collapse = ".")))
	interact.pv[is.na(interact.pv)] <- 'NA.explicit'
	exp.treated <- interact.pv[matched.info$img==1]
	exp.control <- interact.pv[matched.info$img==0]
	pair.interact.exact <- mean(as.character(exp.treated) == as.character(exp.control))
	names(pair.interact.exact) <- paste(c('Exact',pair.vars), collapse = '.')


	#do standardized diffs
	#need to generate a vector of names of non-factor variables.
	sdiff.vars <- unique(c(b.vars, pair.vars))
	factor.vars <- sdiff.vars[aaply(sdiff.vars, 1, function(x) is.factor(my.slice[[x]]) || is.character(my.slice[[x]]))]

	#augment my.slice with dummified versions of factors
	temp.my.slice <- my.slice
	add.vars <- c()
	for(varn in factor.vars){
		varn.dummy <- dummy(my.slice[[varn]], name = varn)
		for(coln in colnames(varn.dummy)){
			temp.my.slice[[coln]] <- varn.dummy[,which(colnames(varn.dummy) == coln)]
		}
		add.vars <- c(add.vars, colnames(varn.dummy))
	}
	sdiff.vars <- setdiff(sdiff.vars, factor.vars)
	sdiff.vars <- unique(c(sdiff.vars, add.vars, 'exp', 'comorb','age'))

	z <- temp.my.slice$img
	marg.bal <- aaply(sdiff.vars, 1, function(x) {
			sdiff(x, 'img', orig.data = temp.my.slice, match.data = temp.my.slice[c(which(z==1)[as.numeric(rownames(matches))], which(z==0)[matches]),])[6]
		})
	names(marg.bal) <- paste('SD',sdiff.vars,sep ='.')


	if(!extra) return(c(tv.prop,  chisq.p,  pair.interact.exact, marg.bal))

	#row-wise tests
	pval.vector <- rep(NA, nrow(fb.tab))
	for(i in 1:nrow(fb.tab)){
		new.tab <- rbind(fb.tab[i,], colSums(fb.tab[-i,]))
		pval.vector[i] <- chisq.test(new.tab)$p.value
	}
	names(pval.vector) <- paste('Pval',rownames(fb.tab),sep ='.')

	#calculate rate of exact matching on pair.vars
	exact.prop <- rep(NA, length(pair.vars))
	for(i in 1:length(exact.prop)){
		col.idx <- match(pair.vars[i], colnames(matched.info))
		exp.treated <- as.character(my.slice[,col.idx][my.slice$img==1][as.numeric(rownames(matches))])
		exp.control <- as.character(my.slice[,col.idx][my.slice$img==0][as.vector(matches)])
		exp.treated[is.na(exp.treated)] <- 'NA.explicit'
		exp.control[is.na(exp.control)] <- 'NA.explicit'
		exact.prop[i] <- mean(as.character(exp.treated) == as.character(exp.control))
	}
	names(exact.prop) = paste('Exact', pair.vars, sep = '.')


	if(length(pair.vars) == 1) return(c( tv.prop,  chisq.p,  exact.prop, median(pval.vector)), marg.bal)


	return(c('TV' = tv.prop, 'Chi-squared pval' = chisq.p, pair.interact.exact, exact.prop, median(pval.vector), marg.bal))
}



descr.stats_general <- function(matches, df, treatCol, b.vars, pair.vars, extra = FALSE){
  treatment.status <- c(rep(1, nrow(matches)), rep(0,
                                                   nrow(matches)))

  #total variation imbalance
  matched.info <- rbind(df[df[,treatCol] == 1,][as.numeric(rownames(matches)),
                                                     , drop = FALSE], df[df[,treatCol] == 0,][as.vector(matches), , drop = FALSE])
  interact.factors.matched = as.factor(apply(matched.info[,
                                                          match(b.vars, colnames(matched.info)), drop = FALSE],
                                             1, function(x) paste(x, collapse = ".")))
  fb.tab <- table('balance.variable' = interact.factors.matched, treatment.status)
  tv.sum <- sum(abs(fb.tab[,1] - fb.tab[,2]))
  tv.prop <-  sum(abs(fb.tab[,1] - fb.tab[,2]))/sum(fb.tab)
  names(tv.sum) <- paste(c('TV',b.vars), collapse = '.')
  names(tv.prop) <- paste(c('TV',b.vars), collapse = '.')

  #chi-squared test
  chisq.p <-  chisq.test(fb.tab)$p.value
  names(chisq.p) <- paste(c('Chisq',b.vars), collapse = '.')

  interact.pv = as.character(apply(matched.info[,
                                                match(pair.vars, colnames(matched.info)), drop = FALSE],
                                   1, function(x) paste(x, collapse = ".")))
  interact.pv[is.na(interact.pv)] <- 'NA.explicit'
  exp.treated <- interact.pv[matched.info[,treatCol]==1]
  exp.control <- interact.pv[matched.info[,treatCol]==0]
  pair.interact.exact <- mean(as.character(exp.treated) == as.character(exp.control))
  names(pair.interact.exact) <- paste(c('Exact',pair.vars), collapse = '.')


  #do standardized diffs
  #need to generate a vector of names of non-factor variables.
  sdiff.vars <- unique(c(b.vars, pair.vars))
  factor.vars <- sdiff.vars[aaply(sdiff.vars, 1, function(x) is.factor(df[[x]]) || is.character(df[[x]]))]

  #augment my.slice with dummified versions of factors
  temp.df <- df
  add.vars <- c()
  for(varn in factor.vars){
    varn.dummy <- dummy(df[[varn]], name = varn)
    for(coln in colnames(varn.dummy)){
      temp.df[[coln]] <- varn.dummy[,which(colnames(varn.dummy) == coln)]
    }
    add.vars <- c(add.vars, colnames(varn.dummy))
  }
  sdiff.vars <- setdiff(sdiff.vars, factor.vars)
  ## SMALL Change in line 134
  ## sdiff.vars <- unique(c(sdiff.vars, add.vars, 'exp', 'comorb','age'))

  z <- temp.df[,treatCol]
  marg.bal <- aaply(sdiff.vars, 1, function(x) {
    sdiff(x, treatCol, orig.data = temp.df, match.data = temp.df[c(which(z==1)[as.numeric(rownames(matches))], which(z==0)[matches]),])[6]
  })
  names(marg.bal) <- paste('SD',sdiff.vars,sep ='.')


  if(!extra) return(c(tv.prop,  chisq.p,  pair.interact.exact, marg.bal))

  #row-wise tests
  pval.vector <- rep(NA, nrow(fb.tab))
  for(i in 1:nrow(fb.tab)){
    new.tab <- rbind(fb.tab[i,], colSums(fb.tab[-i,]))
    pval.vector[i] <- chisq.test(new.tab)$p.value
  }
  names(pval.vector) <- paste('Pval',rownames(fb.tab),sep ='.')

  #calculate rate of exact matching on pair.vars
  exact.prop <- rep(NA, length(pair.vars))
  for(i in 1:length(exact.prop)){
    col.idx <- match(pair.vars[i], colnames(matched.info))
    exp.treated <- as.character(df[,col.idx][df[,treatCol]==1][as.numeric(rownames(matches))])
    exp.control <- as.character(df[,col.idx][df[,treatCol]==0][as.vector(matches)])
    exp.treated[is.na(exp.treated)] <- 'NA.explicit'
    exp.control[is.na(exp.control)] <- 'NA.explicit'
    exact.prop[i] <- mean(as.character(exp.treated) == as.character(exp.control))
  }
  names(exact.prop) = paste('Exact', pair.vars, sep = '.')


  if(length(pair.vars) == 1) return(c( tv.prop,  chisq.p,  exact.prop, median(pval.vector)), marg.bal)


  return(c('TV' = tv.sum, 'Chi-squared pval' = chisq.p, pair.interact.exact, exact.prop, median(pval.vector), marg.bal))
}
