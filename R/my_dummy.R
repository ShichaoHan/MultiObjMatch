#this is a modified version of the function "dummy" from the R package dummies
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
