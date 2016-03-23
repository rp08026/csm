#  File CSSR.R
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
# Check if data for CSSR is appropriate
checkCSSRData <- function (data) UseMethod("checkCSSRData")

checkCSSRData.integer <- function (data) data

checkCSSRData.numeric <- function (data) 
{
    warning("CSSR function currently supports only numeric integer data, data coerced to integer.")
    as.integer(data)
}
checkCSSRData.logical <- function (data) 
{
    warning("CSSR function uses logical data as numeric, data coerced to integer.")
    as.integer(data)
}
checkCSSRData.character <- function (data) return(data)

checkCSSRData.list <- function (data)
{
    data <- lapply(data, checkCSSRData)
    dataClass <- sapply(data, class)
    if(length(levels(factor(dataClass))) != 1)
        stop("CSSR function only supports lists with same classes for all objects (which are distinct time-series from the same source).")
    data
}
checkCSSRData.matrix <- function (data)
{
    dataVector <- as.vector(data)
    dataClass <- class(dataVector)
    dataVector <- checkCSSRData(dataVector)
    matrix(dataVector, nrow(data), ncol(data))
}
checkCSSRData.default <- function (data) stop("CSSR function currently does not support this class/type of input data.")

# ------------------------------------------------------------------------------
#     Specification:
#     It has to:
#     1. Split data, that contains "NA" values to a list without any "NA" values.
#     2. Sequence direction must be preserved.
#     3. Supported data types: discrete vectors, matrices, lists
#     4. Should NOT support tables, functions, arrays.
#     5. Does not currently support dataframes, factors.

splitbyNA <- function (x, ...) UseMethod("splitbyNA")

splitbyNA.vector <- function(x, ...) {
    if(!anyNA(x))
        return(x)
    data.na <- is.na(x)
    data <- !data.na
    data.id <- cumsum(data.na) + 1
    result <- split(x[data], data.id[data], ...)
    names(result) <- NULL
    result
}

splitbyNA.integer <- function(x, ...) splitbyNA.vector(x, ...)
splitbyNA.numeric <- function(x, ...) splitbyNA.vector(x, ...)
splitbyNA.logical <- function(x, ...) splitbyNA.vector(x, ...)
splitbyNA.character <- function(x, ...) splitbyNA.vector(x, ...)

splitbyNA.list <- function(x, ...) {
    id.na <- sapply(x, anyNA)
    if(!any(id.na))
        return(x)
    newsplit <- unlist(lapply(x[id.na], splitbyNA.vector, ...), recursive = FALSE)
    c(x[!id.na], newsplit)
}

splitbyNA.matrix <- function(x, MARGIN, ...) {
    if(!anyNA(x))
        return(x)
    dims <- dim(x)
    if(dims[1] > dims[2] && MARGIN == 1)
        warning("Length of rows is smaller than length of columns, but rows are used as distinct time series. If this is incorrect, use argument MARGIN = 2.")
    if(dims[2] > dims[1] && MARGIN == 2)
        warning("Length of columns is smaller than length of rows, but columns are used as distinct time series. If this is incorrect, use argument MARGIN = 1.")
    if(dims[1] ==  dims[2] && MARGIN == 1)
        warning("Length of rows is equal to length of columns, CHECK if it is correct that rows are used as distinct time series. If this is incorrect, change argument MARGIN.")
    if(dims[1] ==  dims[2] && MARGIN == 2)
        warning("Length of rows is equal to length of columns, CHECK if it is correct that columns are used as distinct time series. If this is incorrect, change argument MARGIN.")
    split <- apply(x, MARGIN, splitbyNA.vector, ...)
    id.list <- sapply(split, is.list)
    c(split[!id.list], unlist(split[id.list], recursive = FALSE))
}
# ------------------------------------------------------------------------------
parser <- function(data, L, MARGIN, cut.end = TRUE)
    UseMethod("parser")

parser.numeric <- function(data, L, MARGIN, cut.end = TRUE) parser.vector(as.integer(data), L, cut.end)
parser.character <- function(data, L, MARGIN, cut.end = TRUE) parser.vector(data, L, cut.end)

parser.vector <- function(vector, L, cut.end, alphabet = levels(factor(data)))
{
    # Check input values
    dl <- length(vector)
    if(dl<L)
        stop("Vector's length is smaller than L")
    # Find alphabet values
    data <- as.character(vector)
    al <- length(alphabet)
    # Create all possible permutation of words.
    pn <- al*(al^L-1)/(al-1) # permutation number
    word.names <- character(pn)
    for(j in 1:L)
    {
        a <- matrix("", nrow = al^j, ncol = j)
        for(i in 1:j)
        {
            a[,i] <-rep(alphabet, each=al^(j-i), times = al^(i-1))
        }
        word.names[((al^j-1)/(al-1)):(al*(al^j-1)/(al-1))] <- apply(a, 1, paste0, collapse ="")
    }
    # Check if there are same words made by different symbol sequences
    if(length(unique(word.names)) != pn)
        stop("Symbols contain same characters, parser is unable to calculate frequencies correctly")
    # Create matrix of words
    word.matrix <- matrix("", nrow = dl, ncol = L)
    for(k in 1:L)
    {
        word.matrix[1:(dl+1-k),k] <- data[k:dl]
    }
    if(cut.end) {
        dl <- dl-L+1
        word.matrix <- word.matrix[1:dl,]
    }
    # Get the vector with all words
    if( length(word.matrix) > L )
        word.vector <- apply(word.matrix, 1, paste0, collapse ="")
    else
        word.vector <- paste0(word.matrix, collapse ="")
    # Calculate word frequencies
    final.vector <- vector("integer", pn)
    names(final.vector) <- word.names
    word.table <- table(word.vector)
    names.expr <- apply(cbind(rep("^", pn), word.names), 1, paste0, collapse ="")
    sum.list <- lapply(names.expr, grep,  x=word.names, value=TRUE)
    for(l in 1:pn)
    {
        final.vector[l] <- sum(word.table[sum.list[[l]]], na.rm =TRUE)
    }
    attr(final.vector, "L") <- L
    attr(final.vector, "alphabet") <- alphabet
    attr(final.vector, "alph.len") <- al
    attr(final.vector, "length") <- dl
    attr(final.vector, "cut.end") <- cut.end
    attr(final.vector, "class") <- "csm.parsed.vector"
    final.vector
}
parser.list <- function(data, L, MARGIN, cut.end = TRUE)
{   
    list.lengths <- lapply(data, length)
    data2 <- data[which(list.lengths >= L)]
    factors <- lapply(data, factor)
    alphabet <- levels(unlist(factors))
    x <- lapply(data2, parser.vector, L, alphabet, cut.end)
    dl <- sapply(x,attr, "length")
    x.split <- split(unlist(x), names(x[[1]]))
    final.vector <- sapply(x.split, sum)
    final.vector <- final.vector[names(x[[1]])]
    attr(final.vector, "L") <- L
    attr(final.vector, "alphabet") <- alphabet
    attr(final.vector, "alph.len") <- length(alphabet)
    attr(final.vector, "length") <- sum(dl)
    attr(final.vector, "cut.end") <- cut.end
    attr(final.vector, "class") <- "csm.parsed.vector"
    final.vector
}


parser.matrix <- function(data, L, MARGIN, cut.end = TRUE)
{
    dims <- dim(data)
    if(dims[1] > dims[2] && MARGIN == 1)
        warning("Length of rows is smaller than length of columns, but rows are used as distinct time series. If this is incorrect, use argument MARGIN = 2.")
    if(dims[2] > dims[1] && MARGIN == 2)
        warning("Length of columns is smaller than length of rows, but columns are used as distinct time series. If this is incorrect, use argument MARGIN = 1.")
    if(dims[1] ==  dims[2] && MARGIN == 1)
        warning("Length of rows is equal to length of columns, CHECK if it is correct that rows are used as distinct time series. If this is incorrect, change argument MARGIN.")
    if(dims[1] ==  dims[2] && MARGIN == 2)
        warning("Length of rows is equal to length of columns, CHECK if it is correct that columns are used as distinct time series. If this is incorrect, change argument MARGIN.")
    if(any(apply(data, MARGIN, length) < L))
        stop("Length of obesrvable sequence is smaller than L. Make sure argument MARGIN is correct.")
    split <- apply(data, MARGIN, parser.vector, L, cut.end)
    #list.lengths <- lapply(data, length)
    #data2 <- data[which(list.lengths >= L)]
    factors <- lapply(data, factor)
    alphabet <- levels(unlist(factors))
    #x <- lapply(data2, parser.vector, L, cut.end, alphabet)
    #dl <- sapply(x,attr, "length")
    # x.split <- split(unlist(x), names(x[[1]]))
    final.vector <- apply(split, 1, sum)
    attr(final.vector, "L") <- L
    attr(final.vector, "alphabet") <- alphabet
    attr(final.vector, "alph.len") <- length(alphabet)
    attr(final.vector, "length") <- sum(final.vector[1:length(alphabet)])
    attr(final.vector, "cut.end") <- cut.end
    attr(final.vector, "class") <- "csm.parsed.vector"
    final.vector
}
# ------------------------------------------------------------------------------
#' @export
#' @title Use CSSR algorithm to infer Causal State Model
#' @description Infers Causal State Model from data
#' @usage CSSR(data, L, alpha = 0.05, MARGIN = 1, verbose = FALSE)
#' @param data a vector, matrix or list with data.
#' @param L length of longest history examined.
#' @param alpha a significance level used to decide whether two distributions differ.
#' @param MARGIN an argument only relevant if "data" is matrix, needed for apply function. Defaults to 1(rows are time-series sequences).
#' @param verbose flag to print out intermediate information. Defaults to FALSE.
#' @details Accepted data at the moment is limited to integer, character or coered to integer numeric vectors,
#' matrixes and lists with same-class vectors.
#' Below is the abstract from help from authors of algorithm. 
#' See http://bactra.org/CSSR/ for the source.
#' 
#' CSSR tries to infer the minimal Markovian model capable of generating a
#' time-series, or set of time-series from the same source.  The program
#' implements the algorithm proposed in the paper "Blind Construction of Optimal
#' Nonlinear Recursive Predictors for Discrete Sequences", hereafter BC. [...]
#' We won't describe the algorithm in any detail here (see BC for that), but the
#' next two paragraphs say a little about what it produces and how it does it.
#' 
#' The output of the algorithm is a set of states which form a Markov chain.  Each
#' state has a certain probability of emitting any of the symbols in the original
#' time series.  The current state and the symbol it emits fix the next state.
#' (The states are "future-resolving", if you're from nonlinear dynamics, or
#' "deterministic", if you're from automata theory.)  Each state, moreover,
#' corresponds to a distinct set of strings, in the following sense.  If the state
#' A contains a string w, and at time t the time-series ends with w, then at time
#' t the Markov chain is in state A.  The set of states, their transition
#' probabilities and connections, is called the state machine.
#' 
#' The algorithm uses a recursive inference procedure to find the simplest set of
#' states with the above properties that can reproduce the statistical properties
#' of the data.  If we could give the algorithm an infinitely long time series,
#' and let it consider infinitely long sub-strings, it would produce the causal
#' states of the process, which are its ideal predictors (see BC for a formal
#' definition).  Since we have only finite data, there is always some probability
#' that the inferred or estimated states are not the true causal states.
#' Nonetheless, for the rest of this file, when we say "causal states", we mean
#' the estimated causal states.
#' 
#' Some Suggestions About Parameters
#' 
#' It is always good to use as much data as you can.  While it is generally good
#' practice to hold back some data for testing or cross-validation, we recommend
#' that this be minimized.  High-entropy processes are especially data-hungry.
#' (See BC.)  For reference, let us call the number of data-points N.
#' 
#' The two key parameters of the program are the maximum history length, L, and
#' the significance level used in the test, s. For any given process, there is a
#' minimum history length M, such that the true states cannot be found if L < M.
#' The number of states returned may be less than the correct number _or higher.
#' If L >= M, and there is enough data, there will generally be a "plateau" of
#' values of L where the correct number of states is returned.  For fixed N, if we
#' keep increasing L, then past a certain point there are not enough examples of
#' each string in the data.  This tends to erroneously create new states, which
#' spawn others through determinization.  Thus there is generally a "blow-up" when
#' L is too large (relative to N and s).  A rough guide-line is to limit L to no
#' more than log(N)/log(k), where k is the alphabet size (see BC for
#'                                                        details).
#' 
#' In general, one should use as small an L as possible, since under-sampling,
#' even before the blow-up, will reduce the accuracy of many probability
#' estimates.  Blow-up can be delayed by reducing s --- that is, reducing the
#' probability of mistakenly splitting a state --- but this carries the risk of
#' failing to create valid new states.  We suggest exploring the data at low L and
#' high s initially, and then increasing L and lowering s.  If a stable
#' architecture is found, it should be recorded at the lowest possible L.
#' @examples
#' even <- csm(data.frame(c("A", "A", "B"),c("A", "B", "A"), c(0, 1, 1), c(0.5, 0.5, 1)))
#' data <- simulate(even, 1000)$symbols
#' acsm <- CSSR(data, 3, 0.001)
#' 

CSSR <- function(data, L, alpha, MARGIN = 1, cut.end = TRUE, verbose = FALSE)
{
    null.state = "*"
    # Check for essential arguments
    if (missing(data) || missing(L) || missing(alpha))
        stop("Essential arguments are missing.")
    L <- as.integer(L)
    if (L<1)
        stop("L must be a natural number (L >= 1).")
    # Check data
    significance <- 1L - alpha
    L <- as.integer(L)
    data <- checkCSSRData(data)
    # Parse the data
    word.count <- parser(splitbyNA(data, MARGIN), L+1, cut.end)
    alphabet <- attr(word.count, "alphabet")
    al <- attr(word.count, "alph.len")
    dl <- attr(word.count, "length")
    # Initialize
    if(any(names(word.count) == null.state))
        stop("Data contains a word for the null state : ", null.state)
    state.list <- list(A = null.state)
    state.morph <- list()
    state.morph[[1]] <- word.count[alphabet]/dl
    # 1. Homogenize
    # Rearrangement of original algorithm for vectorisation purposes.
    # 1.1. Find all suffix morphs
    sl <- al*(al^L-1)/(al-1)
    index <- 1:sl
    suffix.count <- matrix(word.count[(al+1):length(word.count)], nrow = sl, byrow=TRUE)
    rownames(suffix.count) <- names(word.count)[index]
    colnames(suffix.count) <- alphabet
    suffix.morph <- suffix.count/word.count[index]
    
    suffix.count.all <- matrix(c(word.count[alphabet], t(suffix.count)), nrow = sl+1, byrow=TRUE)
    rownames(suffix.count.all) <- c(null.state, names(word.count)[index])
    colnames(suffix.count.all) <- alphabet
    
    suffix.morph.left <- na.exclude(suffix.morph)
    suffix.count.left <- suffix.count[names(suffix.morph.left),]
    suffix.state <- rep(NA, length(rownames(suffix.count.left)))
    names(suffix.state) <- rownames(suffix.count.left)
    # 1.2. Assign suffix morphs to states
    last.round.states <- 1
    last.round.suffixes <- ""
    for (l in 0:(L-1)) {
        names.start <- al*(al^l-1)/(al-1)+1
        names.end <- al*(al^(l+1)-1)/(al-1)
        child.suffixes <- names(word.count)[names.start:names.end]
        child.suffixes.left <- child.suffixes[which(child.suffixes %in% rownames(suffix.morph.left))]
        state.morph <- lapply(state.list, function(x) 
            colSums(matrix(suffix.count.all[x,], ncol = al)/sum(colSums(matrix(suffix.count.all[x,], ncol = al)))))
        pvaluematrix <- matrix(NA, nrow = length(state.list), ncol = al^(l+1))
        pvaluematrix <- sapply(state.morph, function(x) 
            pchisq(colSums(na.omit((t(suffix.morph[child.suffixes.left,])-x)^2/x*rowSums(suffix.count[child.suffixes.left,]))), al-1))
        
        this.round.states <- rep(NA, al^(l+1))
        this.round.suffixes <- rep(NA, al^(l+1))
        for (m in 1:length(child.suffixes.left) ) {
            parent.suffix <- gsub("^.", "", child.suffixes.left[m])
            parent.state <- last.round.states[which(last.round.suffixes == parent.suffix)]
            if(pvaluematrix[child.suffixes.left[m],parent.state] < significance) {
                state.list[[parent.state]] <- c(state.list[[parent.state]], child.suffixes.left[m])
                suffix.state[child.suffixes.left[m]] <- parent.state
                this.round.states[m] <- parent.state
                this.round.suffixes[m] <- child.suffixes.left[m]
                state.morph <- lapply(state.list, function(x) 
                    colSums(matrix(suffix.count.all[x,], ncol = al)/sum(colSums(matrix(suffix.count.all[x,], ncol = al)))))
                pvaluematrix <- sapply(state.morph, function(x) 
                    pchisq(colSums(na.omit((t(suffix.morph[child.suffixes.left,])-x)^2/x*rowSums(suffix.count[child.suffixes.left,]))), al-1))
                next
            } 
            if(min(pvaluematrix[child.suffixes.left[m],]) < significance) {
                state.new <- which.min(pvaluematrix[child.suffixes.left[m],])
                state.list[[state.new]] <- c(state.list[[state.new]], child.suffixes.left[m])
                suffix.state[child.suffixes.left[m]] <- state.new
                this.round.states[m] <- state.new
                this.round.suffixes[m] <- child.suffixes.left[m]
                state.morph <- lapply(state.list, function(x) 
                    colSums(matrix(suffix.count.all[x,], ncol = al)/sum(colSums(matrix(suffix.count.all[x,], ncol = al)))))
                pvaluematrix <- sapply(state.morph, function(x) 
                    pchisq(colSums(na.omit((t(suffix.morph[child.suffixes.left,])-x)^2/x*rowSums(suffix.count[child.suffixes.left,]))), al-1))
                next                    
            } else {
                state.new <- length(state.list) + 1
                state.list[state.new] <- child.suffixes.left[m]
                suffix.state[child.suffixes.left[m]] <- state.new
                this.round.states[m] <- state.new
                this.round.suffixes[m] <- child.suffixes.left[m]
                state.morph <- lapply(state.list, function(x) 
                    colSums(matrix(suffix.count.all[x,], ncol = al)/sum(colSums(matrix(suffix.count.all[x,], ncol = al)))))
                pvaluematrix <- sapply(state.morph, function(x) 
                    pchisq(colSums(na.omit((t(suffix.morph[child.suffixes.left,])-x)^2/x*rowSums(suffix.count[child.suffixes.left,]))), al-1))
            }
        }
        last.round.states <- this.round.states
        last.round.suffixes <- this.round.suffixes
    } 
    names(state.list) <- LETTERS[1:length(state.list)]
    names(state.morph) <- LETTERS[1:length(state.morph)]
    # Comments in verbose mode
    if(verbose) {
        cat("Results after homogenisation:\n")
        cat("List of states and assigned morhps:\n")
        print(state.list)
        cat("List of morphs transition probabilities:\n")
        print(state.morph)
    }
    # 2. Determinize
    # 2.1. Get the transition structure
    transition.structure <- findTransitionMatrix(state.list, state.morph, suffix.state, alphabet, al, L)
    structure.names <- dimnames(transition.structure)
    # Delete states that do not have enough output transitions.
    transitionsSum <- apply(transition.structure, 1, sum)
    transitionsIndex <-  which(transitionsSum >= 1)
    while(length(which(transitionsSum < 0.99)) != 0) {
        transitionsIndex <- which(transitionsSum >= 1)
        transition.structure <- transition.structure[transitionsIndex,transitionsIndex,]
        if(length(transitionsIndex) > 1)
        transitionsSum <- apply(transition.structure, 1, sum)
        else {
            transition.structure <- array(transition.structure, c(1,1,length(alphabet)))
            dimnames(transition.structure) <-c(from = structure.names[[1]][transitionsIndex], to = structure.names[[2]][transitionsIndex], structure.names[3])
            transitionsSum <- sum(transition.structure)
        }
    }
    states.true <- names(transitionsIndex)
    states.fake <- names(state.list)[!(names(state.list) %in% states.true)]
    state.list <- state.list[names(transitionsIndex)]
    
    # Comments in verbose mode
    if(verbose) {
        if(length(states.fake) > 0)
            cat(paste("States", paste(states.fake, collapse = ", "),"were deleted, as were fake states - had no output transitions or their sum was less than 1.\n"))
        cat("Transition probability matrix:\n")
        print(transition.structure)
    }
    transitionsSum <- apply(transition.structure, 1, sum)
    # 2.2. Eliminate transient states
    transient <- eliminateTransient(transition.structure, state.list)
    transition.structure.rec <- transient[[1]]
    state.list.rec <- transient[[2]]
    # Comments in verbose mode
    if(verbose) {
        cat("Results after elimination of transient states:\n")
        cat("List of states:\n")
        print(state.list.rec)
        cat ("Transition probability matrix:\n")
        print(transition.structure.rec)
    }
    # 2.3.1 Successor states
    state.list.det <- state.list.rec
    transition.structure.det <- transition.structure.rec
    if(L > 1) {        
        for(s in 1:length(state.list.rec)) {
            histories <- state.list.rec[[s]][which(nchar(state.list.rec[[s]]) < L)]
            history.number <- length(histories)
            successor.suffixes <- matrix(paste0(rep(histories, each = al), alphabet), history.number, al, byrow=TRUE)
            successor.states <- matrix(mapply(successor.suffixes, FUN = function(y) which(lapply(state.list.rec, FUN = function(x) any(grep(paste0("^",y,"$"), x))) == TRUE)), history.number, al)
            for(t in 1:al) {
                suc.total <- length(unique(unlist(successor.states[,t])))
                if(suc.total > 1) {
                    new.state.histories.index <- lapply(levels(factor(as.vector(unlist(successor.states[,t])))), function(x) which(x == successor.states[,t]))
                    new.state.histories <- lapply(new.state.histories.index, function(x) histories[x])
                    state.list.det[[s]] <- c(new.state.histories[[1]], state.list.rec[[s]][which(nchar(state.list.rec[[s]]) == L)])
                    state.list.det[(length(state.list.det)+1):(length(state.list.det)+suc.total-1)] <- new.state.histories[2:suc.total]
                }
            }
        }
    }
    transition.structure.det <- findTransitionMatrix(state.list.det, state.morph, suffix.state, alphabet, al, L)
    # 2.3.2 Check transition structure
    
    # 2.4. Remove new transient states
    determinisedTransient <- FALSE   
    if(length(state.list.det)!=length(state.list.rec)) {
        determinisedTransient <- TRUE
        transient2 <- eliminateTransient(transition.structure.det, state.list.det)
        transition.structure.rec2 <- transient2[[1]]
        state.list.rec2 <- transient2[[2]]
    } else {
        state.list.rec2 <- state.list.det
        transition.structure.rec2 <- transition.structure.det
    }
        
    
    # Comments in verbose mode
    if(verbose) {
        cat("Results after determinisation without elimination of new transient states:\n")
        cat("List of states:\n")
        print(state.list.det)
        cat("Transition structure:\n")
        print(transition.structure.det)
    }
    if(determinisedTransient && verbose) {
        cat("Results after elimination of transient states from determinised model:\n")
        cat("List of states:\n")
        print(state.list.rec2)
        cat("Transition structure:\n")
        print(transition.structure.rec2)
    } else cat("No transient states after determinisation.\n")
    csm2(names(state.list.rec2), alphabet, transition.structure.rec2)
}
# ------------------------------------------------------------------------------

eliminateTransient <- function (transition.structure, state.list, sn = length(state.list))
{
    transition.total <- matrix(rowSums(matrix(transition.structure, nrow = sn^2)), sn, sn)
    transition.total.list <- apply(transition.total, 2, function(x) which(x!=0))
    transient.states <- numeric()
    recurrent.states <- numeric()
    for(q in 1:sn) {
        needed.states <- (1:sn)[-q]
        visited.states <- c(transition.total.list[[q]], q)
        if(all(needed.states %in% visited.states)) {
            recurrent.states <- append(recurrent.states, q)
            next
        }
        if(any(visited.states %in% recurrent.states)) {
            recurrent.states <- append(recurrent.states, q)
            next
        }
        chains <- as.list(transition.total.list[[q]])
        while(length(chains) > 0) {
            cn <- length(chains)
            for(r in 1:cn) {
                chain.continuation <- transition.total.list[[chains[[r]]]]
                if(any(visited.states %in% chain.continuation)) {
                    visited.states <- append(visited.states, chain.continuation[!any(chain.continuation %in% visited.states)])
                    chains <- c(chains, as.list(chain.continuation[!any(chain.continuation %in% visited.states)]))
                }
                if(all(needed.states %in% visited.states)) {
                    recurrent.states <- append(recurrent.states, q)
                    break
                    next
                }
                if(any(visited.states %in% recurrent.states)) {
                    recurrent.states <- append(recurrent.states, q)
                    break
                    next
                }
            }
            chains[1:cn] <- NULL
        }
        
    }
    transition.structure.rec <- transition.structure[recurrent.states, recurrent.states, ]
    state.list.rec <- state.list[recurrent.states]
    list(transition.structure.rec, state.list.rec)
}
# ------------------------------------------------------------------------------
findTransitionMatrix <- function (state.list, state.morph, suffix.state, alphabet, al, L)
{
    sn <- length(state.list)
    trnum <- sn*al # Transition number (if no transtitions are forbidden)
    transition.structure <- array(0, dim = c(sn, sn, al))
    transition.counter <- lapply(state.list, FUN = function(x) alphabet)
    
    states <- 1
    state.suffixes <- ""
    
    checked.states <- numeric()
    next.states <- numeric()
    next.suffixes <- character()
    for(p in 0:L) {
        
        for(o in 1:length(states)) {
            for(n in 1:al) {
                suf <- paste0(state.suffixes[o], alphabet[n])
                transition.structure[states[o],suffix.state[suf],n] <- state.morph[[states[o]]][n]
                transition.counter[[states[o]]][n] <- NA
                next.states <- c(next.states, suffix.state[suf])
                next.suffixes <- c(next.suffixes, suf)
            }
            checked.states <- c(checked.states, states)
        }
        states <- next.states[!is.na(next.states)]
        state.suffixes <- next.suffixes[!is.na(next.states)]
        next.states <- numeric()
        next.suffixes <- character()
        if(sum(is.na(unlist(transition.counter))) == trnum) break
    }
    dimnames(transition.structure) <- list(from = names(state.list), to = names(state.list), symbol = alphabet)
    transition.structure
}