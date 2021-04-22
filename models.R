# Federico Marotta (federico.marotta@edu.unito.it)
# Apr 2021

library(mvtnorm)

# The posterior distribution and its gradient
ppost <- function(beta, y, geno1, geno2, sigmasq = 1, b) {
    # browser()
    w <- log(exp(geno1 %*% beta) + exp(geno2 %*% beta))
    dmvnorm(y, w, sigmasq*diag(length(y)), log = T) - b/(2*sigmasq) * sum(beta^2)
}
gradppost <- function(beta, y, geno1, geno2, sigmasq = 1, b) {
    e1 <- exp(geno1 %*% beta)
    e2 <- exp(geno2 %*% beta)
    w <- c(log(e1 + e2))
    wprime <- geno1 * c(e1 / (e1 + e2)) + geno2 * c(e2 / (e1 + e2))
    colSums((y - w) * wprime) / sigmasq - b/sigmasq * beta
}

NestedCV <- function(d.train, params.grid,
                     outer.foldid = NULL, inner.foldid = NULL,
                     n.outer.folds = 5, n.inner.folds = 10,
                     verbose = 1) {

    # Define the folds
    if (is.null(outer.foldid)) {
        outer.foldid <- cut(1:nrow(d.train), breaks = n.outer.folds, labels = FALSE)
        outer.foldid <- sample(outer.foldid)
        # for (i in 1:n.outer.folds) {
        #     cat("Fold n.", i, ": ")
        #     cat(d.train$IID[foldid == i], "\n")
        # }
    }

    # Find the dimension of the parameter space
    n_params <- length(params.grid)

    nested.cv <- lapply(1:n.outer.folds, function(f) {
        if (verbose > 0)
            message("Outer fold", f, "--- Performing the CV on 4/5ths of the data")

        inner <- CV(d.train[outer.foldid != f], params.grid,
                    foldid = inner.foldid, n.folds = n.inner.folds,
                    fit.best = TRUE, verbose = verbose - 1)
        if (is.null(inner$fit))
            return(NULL)

        if (verbose > 0)
            message("Outer fold", f, "--- Testing on 1/5th of the data")

        outer <- Test(inner$fit, d.train[outer.foldid == f], verbose = verbose - 2)
        inner$fit <- NULL
        outer[, seq_len(n_params)] <- NULL

        list(inner = inner, outer = outer)
    })

    if (any(sapply(nested.cv, is.null))) {
        message(warnings())
        message(str(nested.cv))
        return(NULL)
    }

    # nested.cv is a list of five lists, each with two elements: inner and outer.
    # inner is itself, for each fold, a list of three elements: full, summ and fit.
    inner.cv <- lapply(nested.cv, "[[", "inner")
    outer.cv <- lapply(nested.cv, "[[", "outer")

    inner.cv.full <- lapply(seq_along(inner.cv), function(i) {
        inner.cv[[i]]$full$outer_fold <- i
        inner.cv[[i]]$full
    })
    inner.cv.full <- rbindlist(inner.cv.full)
    inner.cv.summ <- lapply(seq_along(inner.cv), function(i) {
        inner.cv[[i]]$summ$outer_fold <- i
        inner.cv[[i]]$summ
    })
    inner.cv.summ <- rbindlist(inner.cv.summ)

    outer.cv.full <- MergeFolds(outer.cv)
    outer.cv.summ <- AggregateFolds(outer.cv.full)

    # Return a list of two elements. each element is a list of two data.tables
    list(inner.cv = list(full = inner.cv.full, summ = inner.cv.summ),
         outer.cv = list(full = outer.cv.full, summ = outer.cv.summ))
}


CV <- function(d.train, params.grid,
               foldid = NULL, n.folds = 5,
               fit.best = FALSE, verbose = 1) {

    # Define the folds
    if (is.null(foldid)) {
        foldid <- cut(1:nrow(d.train), breaks = n.folds, labels = FALSE)
        foldid <- sample(foldid)
    }

    n_params <- length(params.grid)

    # For each fold, loop through all the parameters, then merge
    # everything at the end. This is more efficient than having the loop
    # over the parameters as the outer one. The loop over the parameters
    # must be done inside the user-specified function.
    kfold.cv <- lapply(1:n.folds, function(f) {
        if (verbose > 0)
            message("Fold", f, "--- Doing CV")

        Train(d.train[foldid != f],
              d.train[foldid == f],
              params.grid,
              verbose = verbose - 1)
    })

    if (any(sapply(kfold.cv, is.null))) {
        message(warnings())
        message(str(kfold.cv))
        return(NULL)
    }

    kfold.cv.full <- MergeFolds(kfold.cv)
    l <- list(full = kfold.cv.full, summ = AggregateFolds(kfold.cv.full))

    if (fit.best) {
        if (verbose > 0)
            message("CV --- Fitting the best model")

        # Train the model with the best parameters for this gene
        kfold.cv.params <- as.list(as.numeric(l$summ[which.min(pred_perf_pval), 1:n_params])) # TODO: use user-specified loss() to find best params
        kfold.cv.fit <- Train(d.train,
                              NULL,
                              kfold.cv.params,
                              return.fit = TRUE,
                              verbose = verbose - 1)
        l <- append(l, list(fit = kfold.cv.fit))
    }

    l
}


# TODO: check that d.test is null when return.fit

# TODO: check also that the number of params is one when return.fit

Train <- function(d.train, d.test = NULL, params.grid,
                  return.fit = FALSE, verbose = 0) {
    # browser()
    y <- d.train$TPM
    g1 <- grep("@0", names(d.train), value = T)
    geno1 <- as.matrix(d.train[, ..g1])
    g2 <- grep("@1", names(d.train), value = T)
    geno2 <- as.matrix(d.train[, ..g2])

    fit <- matrix(NA, nrow = length(params.grid[[1]]), ncol = ncol(geno1))
    tryCatch({
    for (i in seq_along(params.grid[[1]])) {
        map <- optim(rep(0, ncol(geno1)), ppost,
                          method = "BFGS", gr = gradppost,
                          control = list(fnscale = -1, maxit = 500),
                          y = y,
                          geno1 = geno1, geno2 = geno2,
                          b = params.grid[[1]][i])
        if (map$convergence != 0)
            warning("optim did not converge\n")
        fit[i, ] <- map$par
    }
    names(fit) <- g1
    fit <- c(list(fit), params.grid)

    if (!is.null(d.test))
        return(Test(fit, d.test, verbose = verbose))

    if (return.fit)
        return(fit)
    }, error = function(e) {
        message(e$message)
        return(NULL)
    })
}


Test <- function(fit, d.test, verbose = 0, return.pred = F) {
    # browser()
    y <- d.test$TPM
    g1 <- grep("@0", names(d.test), value = T)
    geno1 <- as.matrix(d.test[, ..g1])
    g2 <- grep("@1", names(d.test), value = T)
    geno2 <- as.matrix(d.test[, ..g2])

    if (nrow(fit[[1]]) == 1) {
        beta <- matrix(fit[[1]], nrow = 1)
    } else {
        beta <- fit[[1]]
    }
    names(beta) <- g1

    pearsons <- lapply(seq_len(nrow(fit[[1]])), function(p) {
        pred_y <- c(log(exp(geno1 %*% beta[p, ]) + exp(geno2 %*% beta[p, ])))

        pearson <- compute_pearson(y, pred_y, by = d.test$GENE)
        cbind(
            b = as.character(fit[[2]][p]), # These need to be characters otherwise the merge doesn't work
            data.name = pearson$data.name,
            mse = cvmse(y, pred_y),
            rsq = cvrsq(y, pred_y),
            pearson[, c("parameter.df", "estimate.cor", "p.value")]
        )
    })
    rbindlist(pearsons)
}


Coef <- function(fit) {
    matrix(fit[[1]], nrow = 1)[1, ]
}


Predict <- function(fit, newx) {
    g1 <- grep("@0", names(newx), value = T)
    geno1 <- as.matrix(newx[, ..g1])
    g2 <- grep("@1", names(newx), value = T)
    geno2 <- as.matrix(newx[, ..g2])
    beta <- Coef(fit)
    c(log(exp(geno1 %*% beta) + exp(geno2 %*% beta)))
}
