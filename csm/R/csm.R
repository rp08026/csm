#  File csm.R
#  Part of the csm package, http://CRAN.R-project.org/package=csm
#
#  Copyright (C) 2014 Rudolfs Petrovs
#  Copyright (C) 1995-2014 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# ------------------------------------------------------------------------------
#' @export
#' @aliases csm
#' @title Create a CSM object
#' @description Creates a Causal State Model object of class "csm" or "csmlike".
#' @usage csm(trans.df)
#' @param trans.df a dataframe with transitions.
#' @details 'trans.df' row specifies one transition. It is c("From","To","Symbol", "Probability") - first two are state names,
#' then emitted symbol and then this transition probability.
#' @examples
#' coin <- csm(data.frame(rep("COIN", 2), rep("COIN", 2), c("H", "T"), rep(0.5, 2)))
csm <- function (trans.df)
{
    # Check the argument
    if(!is.data.frame(trans.df))
        stop("trans.df is not a dataframe object.")
    if(!is.numeric(trans.df[[4]]))
        stop("trans.df[[4]], the probabilities column, is not numeric.")
    # Make sure the type of arguments
    from <- as.factor(trans.df[[1]])
    to <- as.factor(trans.df[[2]])
    symb <- as.factor(trans.df[[3]])
    prob <- as.numeric(trans.df[[4]])
    # Prepare indexing
    states <- unique(c(levels(from), levels(to)))
    symbols <- levels(symb)
    stl <- length(states)
    syl <- length(symbols)
    statesInd <- 1:stl
    names(statesInd) <- states
    symbolsInd <- 1:syl
    names(symbolsInd) <- symbols
    transInd <- cbind(statesInd[from],statesInd[to], symbolsInd[symb], prob)
    # Create trans.array
    trans.array <- array(0L, c(stl, stl, syl))
    trans.array[transInd[,1:3]] <- transInd[,4]
    dimnames(trans.array) <- list(from = states, to = states, symbol = symbols)
    # Check trans.array
    sumStateOutput <- rowSums(trans.array)
    names(sumStateOutput) <- NULL
    if(!isTRUE(all.equal(sumStateOutput, rep(1, stl))))
        stop("Sum of all output transitions of all states does not add up to 1.")
    ans <- list(states = states, symbols = symbols, trans.array = trans.array)
    class(ans) <- c("csm")
    if(!all(sapply(statesInd, function(a) {isTRUE(length(transInd[which(transInd[,1]==a),3]) == length(unique(transInd[which(transInd[,1]==a),3])))
    }))) {
        class(ans) <- c("csmlike")
        warning("In true Causal State Model given the current state and the next symbol, there is a unique successor state.\nThis model does not have this property so its class will be 'csmlike'.")
    }
    ans
}
# ------------------------------------------------------------------------------
csm2 <- function (states, symbols, trans.array)
{
    # Check the arguments.
#     if (missing(states) || missing(symbols) || missing(trans.array))
#         stop("essential arguments are missing")
    if (!is.vector(states)) {
        warning("'states' argument is not a vector, coercing to a vector")
        states <- as.factor(states)
    }
    if (!is.vector(symbols)) {
        warning("'symbols' argument is not a vector, coercing to a vector")
        symbols <- as.factor(symbols)
    }
    stl <- length(states)
    syl <- length(symbols)
    if (!is.array(trans.array)) {
        warning("'trans.array' argument is not an array, coercing to an array")
        trans.array <- array(trans.array, c(stl, stl, syl))
    }
    if (length(levels(as.factor(states))) != length(states)) {
        stop("'states' argument has duplicate states")
    }
    if (length(levels(as.factor(symbols))) != length(symbols)) {
        stop("'symbols' argument has duplicate symbols")
    }
    # Check if there is only one state or only one symbol.
    if(stl != 1 && syl != 1 ) {
        if (length(dim(trans.array)) != 3)
            stop("'trans.array' argument does not have 3 dimensions")
        if (dim(trans.array)[1] != dim(trans.array)[2])
            stop("'trans.array' first two dimensions do not form a square transition matrix")
        if (dim(trans.array)[1] != stl)
            stop("'trans.array' first two dimensions are not equal to the number of states")
        if (syl != dim(trans.array)[3])
            stop("'trans.array' third dimension is not equal to the number of symbols")
        sumStateOutput <- rowSums(trans.array)
        names(sumStateOutput) <- NULL
        if(!isTRUE(all.equal(sumStateOutput, rep(1, stl))))
            stop("Sum of all output transitions of all states does not add up to 1.")
    }
    
    dimnames(trans.array) <- list(from = states, to = states, symbol = symbols)
    ans <- list(states = as.character(states), symbols = as.character(symbols), trans.array = trans.array)
    class(ans) <- c("csm")
    ans
}
# ------------------------------------------------------------------------------
#' @importFrom stats simulate
#' @export
#' @aliases simulate.csm simulate.csmlike
#' @name simulate
#' @title Simulate a CSM and CSMlike objects
#' @description Simulate an observable sequence from the Causal State Model
#' @method simulate csm
#' @method simulate csmlike
#' @param object an object of class "csm" or "csmlike".
#' @param nsim a number of observations to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number generator should be initialized ('seeded').
#' @param \dots additional optional arguments.
#' @details
#' The result of simulation is a sequence of hidden states under $states 
#' and a sequence of observable emitted symbols under $symbols. Additional information
#' about used indexing is available under $statesInd, $symbolsInd and $transInd.
#' @examples 
#' even <- csm(data.frame(c("A", "A", "B"),c("A", "B", "A"), c(0, 1, 1), c(0.5, 0.5, 1)))
#' simulate(even, 1000)
#' 
## Simulates observations of csm (Causal State Model)
simulate.csm <- function(object, nsim = 1, seed = NULL, ...)
{
    # Random Number Generator initialization
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)                     # initialize the RNG if necessary
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    # Vectors for simulation data
    statesInd <- rep(NA, nsim) 
    transInd <- rep(NA, nsim)
    symbolsInd <- rep(NA, nsim)
    
    states <- rep(NA, nsim)
    symbols <- rep(NA, nsim)
    # Sampling initial state
    statesTotal <- length(object$states)
    symbolsTotal <- length(object$symbols)
    statesInd[1] <- sample(1:statesTotal, 1)
    states[1] <- object$states[statesInd[1]]
    # Complete all the data vectors except last index
    for (i in 1:nsim) {
        # Sampling the transition between all possible states and emission symbols 
        # according to probabilities in trans.array
        transInd[i] <- sample(statesTotal * symbolsTotal, 1, prob = object$trans.array[statesInd[i],,])
        # Writing the emiited symbol
        symbolsInd[i] <- ceiling(transInd[i] / statesTotal)
        symbols[i] <- object$symbols[symbolsInd[i]]
        # Writing the next state
        statesInd[i+1] <- (transInd[i] - ((symbolsInd[i]-1) * statesTotal))
        states[i+1] <- object$states[statesInd[i+1]]
    }
    # Make the return list
    val <- list(states = states[-nsim-1], symbols = symbols, statesInd = statesInd[-nsim-1], 
                symbolsInd = symbolsInd,transInd = transInd )
    attr(val, "seed") <- RNGstate
    val
}
# ------------------------------------------------------------------------------
#' @export
simulate.csmlike <- function(object, nsim = 1, seed = NULL, ...) 
    simulate.csm(object, nsim, seed, ...)
# ------------------------------------------------------------------------------
#' @export
print.csm <- function(x, ...) 
{
    cat("Causal State Model specification:\n")
    cat("States:\n")
    print(x$states)
    cat("Symbols:\n")
    print(x$symbols)
    cat ("Transition probability matrix:\n")
    print(x$trans.array)
    return(invisible(x))
}
# ------------------------------------------------------------------------------
#' @export
print.csmlike <- function(x, ...)
{
    cat("Causal State-like Model specification:\n")
    cat("States:\n")
    print(x$states)
    cat("Symbols:\n")
    print(x$symbols)
    cat ("Transition probability matrix:\n")
    print(x$trans.array)
    return(invisible(x))
}
# ------------------------------------------------------------------------------
#' @export
summary.csm <- function (object, ...) 
{
    return(invisible(print.csm(object)))
}
# ------------------------------------------------------------------------------
#' @export
summary.csmlike <- function (object, ...) 
{
    return(invisible(print.csmlike(object)))
}


# Common models
#' @export
coin <- csm(data.frame(rep("COIN", 2), rep("COIN", 2), c("H", "T"), rep(0.5, 2)))
dice <- csm(data.frame(rep("DICE", 6), rep("DICE", 6), c(1:6), rep(1/6, 6)))
even <- csm(data.frame(c("A", "A", "B"),c("A", "B", "A"), c(0, 1, 1), c(0.5, 0.5, 1)))
FRETlike <- csm(data.frame(rep(c("A", "B"), each= 2), rep(c("A", "B"), times = 2), c(0, 1, 0, 1), c(0.8, 0.2, 0.1, 0.9))) 
# FRETlike <- csm2(c("A", "B"), c(0, 1), array(c(0.8, 0.1, 0, 0, 0, 0, 0.2, 0.9), c(2, 2, 2)))
# even <- csm2(c("A", "B"), c(0, 1), array(c(0.5, 0, 0, 0, 0, 1, 0.5, 0), c(2, 2, 2)))
# FRETlike <- csm2(c("A", "B"), c(0, 1), array(c(0.8, 0.1, 0, 0, 0, 0, 0.2, 0.9), c(2, 2, 2)))
# inferredFRET <- CSSR(simulate(FRETlike, 10000)$symbols, 3, 0.05)