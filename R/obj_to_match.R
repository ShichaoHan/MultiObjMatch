obj.to.match <- function(out.elem, already.done = NULL, prev.obj = NULL){
	tcarcs <- length(unlist(out.elem$net$edgeStructure))
	edge.info <- extractEdges(out.elem$net)
	one_sol <- function(sol){
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
	if(is.null(already.done)) return(llply(out.elem$solutions, one_sol))
	new.ones <- setdiff(1:length(out.elem$solutions), already.done)
	out.list <- list()
	out.list[already.done] <- prev.obj
	out.list[new.ones] <- llply(out.elem$solutions[new.ones], one_sol)
	return(out.list)
}
