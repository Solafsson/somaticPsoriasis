library(foreach)
library(doMC)

registerDoMC(cores=10)

library(lpSolve)

label.switching <- function (method, zpivot, z, K, prapivot, p, complete, mcmc, 
    sjwinit, data, constraint, groundTruth, thrECR, thrSTE, thrSJW, 
    maxECR, maxSTE, maxSJW, userPerm) 
{
    runStart <- proc.time()
    cat("\n")
    L <- length(method)
    for (l in 1:L) {
        if ((method[l] %in% c("ECR", "ECR-ITERATIVE-1", "ECR-ITERATIVE-2", 
            "STEPHENS", "SJW", "PRA", "DATA-BASED", "AIC", "USER-PERM")) == 
            FALSE) {
            stop(cat(paste("method:", method[l], "is not recognised"), 
                "\n"))
        }
    }
    if (missing(z) == FALSE) {
        m <- dim(z)[1]
        nCheck <- dim(z)[2]
        kMinCheck <- max(z)
        kCheck <- kMinCheck
    }
    else {
        if (missing(mcmc) == FALSE) {
            m <- dim(mcmc)[1]
            kCheck <- dim(mcmc)[2]
            kMinCheck <- kCheck
        }
        else {
            if (missing(p) == FALSE) {
                m <- dim(p)[1]
                kCheck <- dim(p)[3]
                nCheck <- dim(p)[2]
                kMinCheck <- kCheck
            }
            else {
                stop(cat(paste("    [ERROR]: At least one of z, mcmc, or p should be provided"), 
                  "\n"))
            }
        }
    }
    if (missing(p) == FALSE) {
        nCheck <- dim(p)[2]
    }
    if (missing(p) == FALSE) {
        kCheck <- dim(p)[3]
    }
    if (missing(mcmc) == FALSE) {
        kCheck <- dim(mcmc)[2]
    }
    if (missing(data) == FALSE) {
        if (is.null(dim(data)) == TRUE) {
            nCheck = length(data)
        }
        else {
            nCheck <- dim(data)[1]
        }
    }
    if (missing(K) == TRUE) {
        K = kCheck
        cat(paste("    [WARNING]: K is not provided. According to input it is assumed that K = ", 
            K, ".", sep = ""), "\n")
    }
    if (max(z) < K) {
        cat(paste("    [WARNING]: max sampled latent allocation = ", 
            max(z), " < ", "K = ", K, ". This indicates that the MCMC sampler has not converged and/or the presence of redundant components.", 
            sep = ""), "\n")
    }
    if (missing(z) == FALSE) {
        if (m != dim(z)[1]) {
            stop(cat(paste("    [ERROR]: MCMC iterations are not equal to dim(z)[1]."), 
                "\n"))
        }
    }
    if (missing(mcmc) == FALSE) {
        if (m != dim(mcmc)[1]) {
            stop(cat(paste("    [ERROR]: MCMC iterations are not equal to dim(mcmc)[1]."), 
                "\n"))
        }
    }
    if (missing(p) == FALSE) {
        if (m != dim(p)[1]) {
            stop(cat(paste("    [ERROR]: MCMC iterations are not equal to dim(p)[1]."), 
                "\n"))
        }
    }
    if (missing(mcmc) == FALSE) {
        if (missing(kCheck)) {
            kCheck = dim(mcmc)[2]
        }
        if (kCheck != dim(mcmc)[2]) {
            stop(cat(paste("    [ERROR]: K is not equal to dim(mcmc)[2]."), 
                "\n"))
        }
    }
    if (missing(p) == FALSE) {
        if (missing(kCheck)) {
            dim(p)[3]
        }
        if (kCheck != dim(p)[3]) {
            stop(cat(paste("    [ERROR]: K is not equal to dim(p)[3]."), 
                "\n"))
        }
    }
    if (missing(K) == FALSE) {
        if (kCheck != K) {
            stop(cat(paste("    [ERROR]: Number of components is not consistent with the input."), 
                "\n"))
        }
        if (kMinCheck > K) {
            stop(cat(paste("    [ERROR]: Number of components should be at least equal to ", 
                kMinCheck, ",", sep = ""), "\n"))
        }
    }
    if (missing(z) == FALSE) {
        if (nCheck != dim(z)[2]) {
            stop(cat(paste("    [ERROR]: Number of observations is not equal to dim(z)[2]."), 
                "\n"))
        }
    }
    if (missing(p) == FALSE) {
        if (nCheck != dim(p)[2]) {
            stop(cat(paste("    [ERROR]: Number of observations is not equal to dim(p)[2]."), 
                "\n"))
        }
    }
    if (missing(zpivot) == FALSE) {
        if (is.null(dim(zpivot)) == TRUE) {
            if (nCheck != length(zpivot)) {
                stop(cat(paste("    [ERROR]: Number of observations is not equal to length(zpivot)."), 
                  "\n"))
            }
        }
        else {
            if (nCheck != dim(zpivot)[2]) {
                stop(cat(paste("    [ERROR]: Number of observations is not equal to dim(zpivot)[2]."), 
                  "\n"))
            }
        }
    }
    else {
        zpivot = array(1)
    }
    if (missing(prapivot) == FALSE) {
        if (kCheck != dim(prapivot)[1]) {
            stop(cat(paste("    [ERROR]: K is not equal to dim(prapivot)[1]."), 
                "\n"))
        }
    }
    if (missing(prapivot) == FALSE) {
        if (dim(mcmc)[3] != dim(prapivot)[2]) {
            stop(cat(paste("    [ERROR]: J is not equal to dim(prapivot)[2]."), 
                "\n"))
        }
    }
    if (missing(data) == FALSE) {
        if (is.null(dim(data)) == TRUE) {
            nX <- length(data)
        }
        else {
            nX <- dim(data)[1]
        }
        if (nCheck != nX) {
            stop(cat(paste("    [ERROR]: data length is not compatible with the input."), 
                "\n"))
        }
    }
    if (missing(groundTruth) == FALSE) {
        if (length(groundTruth) != nCheck) {
            stop(cat(paste("    [ERROR]: length(groundTruth) is not equal to number of observations."), 
                "\n"))
        }
        if (all(groundTruth == floor(groundTruth)) == FALSE) {
            stop(cat(paste("    [ERROR]: non-integer groundTruth entries are not allowed."), 
                "\n"))
        }
    }
    if (("ECR" %in% method) == TRUE) {
        if (missing(z) == TRUE) 
            stop(cat(paste("    [ERROR]: z is required for ECR."), 
                "\n"))
        if (missing(zpivot) == TRUE) 
            stop(cat(paste("    [ERROR]: zpivot is required for ECR."), 
                "\n"))
    }
    if (("ECR-ITERATIVE-1" %in% method) == TRUE) {
        if (missing(z) == TRUE) 
            stop(cat(paste("    [ERROR]: z is required for ECR-ITERATIVE-1."), 
                "\n"))
    }
    if (("ECR-ITERATIVE-2" %in% method) == TRUE) {
        if (missing(z) == TRUE) 
            stop(cat(paste("    [ERROR]: z is required for ECR-ITERATIVE-2."), 
                "\n"))
    }
    if (("ECR-ITERATIVE-2" %in% method) == TRUE) {
        if (missing(p) == TRUE) 
            stop(cat(paste("    [ERROR]: p is required for ECR-ITERATIVE-2."), 
                "\n"))
    }
    if (("STEPHENS" %in% method) == TRUE) {
        if (missing(p) == TRUE) 
            stop(cat(paste("    [ERROR]: p is required for STEPHENS."), 
                "\n"))
    }
    if (("SJW" %in% method) == TRUE) {
        if (missing(data) == TRUE) 
            stop(cat(paste("    [ERROR]: data is required for SJW."), 
                "\n"))
    }
    if (("SJW" %in% method) == TRUE) {
        if (missing(complete) == TRUE) 
            stop(cat(paste("    [ERROR]: complete is required for SJW."), 
                "\n"))
    }
    if (("SJW" %in% method) == TRUE) {
        if (missing(z) == TRUE) 
            stop(cat(paste("    [ERROR]: z is required for SJW."), 
                "\n"))
    }
    if (("SJW" %in% method) == TRUE) {
        if (missing(mcmc) == TRUE) 
            stop(cat(paste("    [ERROR]: mcmc is required for SJW."), 
                "\n"))
    }
    if (("AIC" %in% method) == TRUE) {
        if (missing(mcmc) == TRUE) 
            stop(cat(paste("    [ERROR]: mcmc is required for AIC."), 
                "\n"))
    }
    if (("DATA-BASED" %in% method) == TRUE) {
        if (missing(z) == TRUE) 
            stop(cat(paste("    [ERROR]: z is required for DATA-BASED."), 
                "\n"))
    }
    if (("DATA-BASED" %in% method) == TRUE) {
        if (missing(data) == TRUE) 
            stop(cat(paste("    [ERROR]: data is required for DATA-BASED."), 
                "\n"))
    }
    if (kCheck < 2) {
        stop(cat(paste("    [ERROR]: K should be at least equal to 2."), 
            "\n"))
    }
    if (nCheck < 2) {
        stop(cat(paste("    [ERROR]: n should be at least equal to 2."), 
            "\n"))
    }
    if (("PRA" %in% method) && (kCheck > 8)) {
        cat(paste("    [WARNING]: PRA is not suggested for ", 
            kCheck, "components"), "\n")
    }
    if (("SJW" %in% method) && (kCheck > 8)) {
        cat(paste("    [WARNING]: SJW is not suggested for ", 
            kCheck, "components"), "\n")
    }
    if (missing(thrECR)) {
        thrECR <- 10^(-6)
    }
    if (missing(thrSTE)) {
        thrSTE <- 10^(-6)
    }
    if (missing(thrSJW)) {
        thrSJW <- 10^(-6)
    }
    if (missing(maxECR)) {
        maxECR <- 100
    }
    if (missing(maxSTE)) {
        maxSTE <- 100
    }
    if (missing(maxSJW)) {
        maxSJW <- 100
    }
    minThreshold <- 1e-12
    thrECR <- max(minThreshold, thrECR)
    thrSTE <- max(minThreshold, thrSTE)
    thrSJW <- max(minThreshold, thrSJW)
    if (maxECR < 1) {
        maxECR <- 100
    }
    if (maxSTE < 1) {
        maxSTE <- 100
    }
    if (maxSJW < 1) {
        maxSJW <- 100
    }
    if (("USER-PERM" %in% method) == TRUE) {
        if (missing(K) == TRUE) {
            stop(cat(paste("   [ERROR]: K is not supplied."), 
                "\n"))
        }
        if (is.list(userPerm) == FALSE) {
            temp <- vector("list", length = 1)
            temp[[1]] <- userPerm
            userPerm <- temp
            temp <- 0
            userLength <- 1
        }
        else {
            userLength <- length(userPerm)
            if (is.null(names(userPerm)) == TRUE) {
                names(userPerm) <- paste("user", (1:userLength), 
                  sep = "-")
            }
        }
        names(userPerm) <- paste("USER", (1:userLength), sep = "-")
        cat(paste("    Checking user-supplied permutations for consistency..."))
        userStatus <- rep(1, userLength)
        for (j in 1:userLength) {
            if (dim(z)[1] != dim(userPerm[[j]])[1]) {
                userStatus[j] = 0
                cat("\n")
                cat(paste("    [ERROR]: number of MCMC samples should be equal to number of permutations."))
            }
            sss <- 1:K
            if (K != dim(userPerm[[j]])[2]) {
                userStatus[j] = 0
                cat("\n")
                cat(paste("    [ERROR]: Number of components should be equal to permutation size."))
            }
            for (i in 1:dim(userPerm[[j]])[1]) {
                myCheck <- table(match(unique(userPerm[[j]][i, 
                  ]), sss))
                if (length(myCheck) != K) {
                  userStatus[j] = 0
                  cat(paste("\n"))
                  cat(paste("    problem in line ", i, " of supplied permutation set ", 
                    j, ":", sep = ""), "\n")
                  cat(paste("   "), paste(userPerm[[j]][i, ]), 
                    paste("is not a permutation of {1,...,", 
                      K, "}"), "\n")
                  stop(cat(paste("    [ERROR]: user-defined input is not valid"), 
                    "\n"))
                }
            }
        }
        cat(paste("done."), "\n")
    }
    fr <- 1
    dimname <- L
    if ((is.array(zpivot) == TRUE) && (("ECR" %in% method) == 
        TRUE)) {
        fr <- dim(zpivot)[1]
        dimname <- L + fr - 1
    }
    if ("ECR" %in% method) {
        ind <- which(method == "ECR")
        if (length(ind) > 1) {
            stop(paste("ECR appearing more than 1 times"))
        }
        nam <- numeric(dimname)
        if (ind == 1) {
        }
        else {
            for (i in 1:(ind - 1)) {
                nam[i] <- method[i]
            }
        }
        for (i in 1:fr) {
            nam[ind + i - 1] <- paste("ECR", i, sep = "-")
        }
        if ((ind + fr) <= dimname) {
            for (i in (ind + fr):(dimname)) {
                nam[i] <- method[i - fr + 1]
            }
        }
        if (fr == 1) {
            nam[ind] <- "ECR"
        }
    }
    else {
        nam <- method
    }
    nams1 <- nam
    if ((missing(constraint) == FALSE) && (length(constraint) == 
        1) && (constraint == "ALL")) {
        J <- dim(mcmc)[3]
        constraint = 1:J
    }
    fr <- 1
    if ((missing(constraint) == FALSE) && (length(constraint) > 
        1) && (("AIC" %in% method) == TRUE)) {
        fr <- length(constraint)
        dimname <- dimname + fr - 1
    }
    if ("AIC" %in% method) {
        ind <- which(nams1 == "AIC")
        if (length(ind) > 1) {
            stop(paste("AIC appearing more than 1 times"))
        }
        if (ind == 1) {
        }
        else {
            for (i in 1:(ind - 1)) {
                nam[i] <- nams1[i]
            }
        }
        for (i in 1:fr) {
            nam[ind + i - 1] <- paste("AIC", i, sep = "-")
        }
        if ((ind + fr) <= dimname) {
            for (i in (ind + fr):(dimname)) {
                nam[i] <- nams1[i - fr + 1]
            }
        }
        if (fr == 1) {
            nam[ind] <- "AIC"
        }
    }
    nams2 <- nam
    if ("USER-PERM" %in% method) {
        fr <- userLength
        dimname <- dimname + fr - 1
        ind <- which(nams2 == "USER-PERM")
        if (length(ind) > 1) {
            stop(paste("USER-PERM appearing more than 1 times"))
        }
        if (ind == 1) {
        }
        else {
            for (i in 1:(ind - 1)) {
                nam[i] <- nams2[i]
            }
        }
        for (i in 1:fr) {
            nam[ind + i - 1] <- names(userPerm)[i]
        }
        if ((ind + fr) <= dimname) {
            for (i in (ind + fr):(dimname)) {
                nam[i] <- nams2[i - fr + 1]
            }
        }
        if (fr == 1) {
            nam[ind] <- "USER"
        }
    }
    myPrettyPrint <- function(gap0, gap1, gap2, gap3, word1, 
        word2, word3) {
        nch1 = nchar(as.character(word1))
        nch2 = nchar(as.character(word2))
        nch3 = nchar(as.character(word3))
        cat(rep("", gap0), ".", word1, paste(rep("", max(gap1 - 
            nch1, 1))), word2, rep("", max(1, gap2 - nch2)), 
            word3, rep("", max(1, gap3 - nch3)), ".", "\n")
    }
    permutations <- vector("list", length = dimname)
    names(permutations) <- nam
    timings <- numeric(dimname)
    names(timings) <- nam
    f <- 0
    fic <- 0
    fuser <- 0
    t <- 1
    userIterator <- 0
    gap0 <- 4
    gap1 <- 30
    gap2 <- 20
    gap3 <- 30
    cat(paste("    ......................................................................................\n"))
    myPrettyPrint(gap0, gap1, gap2, gap3, "Method", "Time (sec)", 
        "Status    ")
    cat(paste("    ......................................................................................\n"))
    for (l in 1:L) {
        if (method[l] == "ECR") {
            if (is.array(zpivot) == TRUE) {
                while (f < dim(zpivot)[1]) {
                  f <- f + 1
                  if (dim(zpivot)[2] != dim(z)[2]) {
                    stop(paste("length(zpivot) and number of columns of z are not equal"))
                  }
                  if (K < max(z)) {
                    stop(paste("K should be at least equal to", 
                      max(z)))
                  }
                  tpm <- proc.time()
                  permutations[[t]] <- ecr(zpivot[f, ], z, K)$permutations
                  time <- proc.time() - tpm
                  time <- round(as.numeric(time[3]), 3)
                  voutsas <- paste(method[l], " (pivot ", f, 
                    " of ", dim(zpivot)[1], ")", sep = "")
                  myPrettyPrint(gap0, gap1, gap2, gap3, voutsas, 
                    time, "OK")
                  timings[t] <- time
                  t <- t + 1
                }
    		
            }
            else {
                if (missing(zpivot)) {
                  stop(paste("zpivot is missing"))
                }
                else {
                  if (length(zpivot) != dim(z)[2]) {
                    stop(paste("length(zpivot) and number of columns of z are not equal"))
                  }
                  if (K < max(z)) {
                    stop(paste("K should be at least equal to", 
                      max(z)))
                  }
                  tpm <- proc.time()
                  permutations[[t]] <- ecr(zpivot, z, K)$permutations
                  time <- proc.time() - tpm
                  time <- round(as.numeric(time[3]), 3)
                  myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                    time, "OK")
                  timings[t] <- time
                  t <- t + 1
                }
            }
        }
        if (method[l] == "ECR-ITERATIVE-1") {
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- ecr.iterative.1(z, K, threshold = thrECR, 
                maxiter = maxECR)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]), 3)
            myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                time, hold$status)
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "ECR-ITERATIVE-2") {
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(p)) {
                stop(paste("p is missing"))
            }
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- ecr.iterative.2(z, K, p, thrECR, maxECR)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]), 3)
            myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                time, hold$status)
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "PRA") {
            tpm <- proc.time()
            permutations[[t]] <- pra(mcmc, prapivot)$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]), 3)
            myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                time, "OK")
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "STEPHENS") {
            if (missing(p)) {
                stop(paste("p is missing"))
            }
            tpm <- proc.time()
            hold <- stephens(p, thrSTE, maxSTE)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]), 3)
            myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                time, hold$status)
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "SJW") {
            if (missing(mcmc)) {
                stop(paste("mcmc is missing"))
            }
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(complete)) {
                stop(paste("Complete log-likelihood function is missing"))
            }
            if (missing(data)) {
                stop(paste("Data is missing"))
            }
            if (missing(sjwinit)) {
                sjwinit = 0
            }
            tpm <- proc.time()
            hold <- sjw(mcmc, z, complete, x = data, sjwinit, 
                thrSJW, maxSJW)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]), 3)
            myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                time, hold$status)
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "DATA-BASED") {
            if (missing(z)) {
                stop(paste("z is missing"))
            }
            if (missing(data)) {
                stop(paste("Data is missing"))
            }
            if (K < max(z)) {
                stop(paste("K should be at least equal to", max(z)))
            }
            tpm <- proc.time()
            hold <- dataBased(x = data, K, z)
            permutations[[t]] <- hold$permutations
            time <- proc.time() - tpm
            time <- round(as.numeric(time[3]), 3)
            myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                time, "OK")
            timings[t] <- time
            t <- t + 1
        }
        if (method[l] == "AIC") {
            if ((missing(constraint) == FALSE) && (length(constraint) > 
                1)) {
                if (missing(mcmc)) {
                  stop(paste("MCMC is missing"))
                }
                while (fic < length(constraint)) {
                  fic <- fic + 1
                  tpm <- proc.time()
                  permutations[[t]] <- aic(mcmc, constraint[fic])$permutations
                  time <- proc.time() - tpm
                  time <- round(as.numeric(time[3]), 3)
                  voutsas <- paste(method[l], " (constraint ", 
                    fic, " of ", length(constraint), ")", sep = "")
                  myPrettyPrint(gap0, gap1, gap2, gap3, voutsas, 
                    time, "OK")
                  timings[t] <- time
                  t <- t + 1
                }
            }
            else {
                if (missing(mcmc)) {
                  stop(paste("MCMC is missing"))
                }
                if (missing(constraint)) {
                  constraint = 1
                }
                tpm <- proc.time()
                hold <- aic(mcmc, constraint)
                permutations[[t]] <- hold$permutations
                time <- proc.time() - tpm
                time <- round(as.numeric(time[3]), 3)
                myPrettyPrint(gap0, gap1, gap2, gap3, method[l], 
                  time, "OK")
                timings[t] <- time
                t <- t + 1
            }
        }
        if (method[l] == "USER-PERM") {
            while (userIterator < userLength) {
                tpm <- proc.time()
                permutations[[t]] <- userPerm[[userIterator + 
                  1]]
                userPerm[[userIterator + 1]] <- 0
                time <- proc.time() - tpm
                time <- round(as.numeric(time[3]), 3)
                if (userIterator > 1) {
                  voutsas <- paste(method[l], " (input ", userIterator + 
                    1, " of ", userLength, ")", sep = "")
                }
                else {
                  voutsas <- paste(method[l])
                }
                if (userStatus[userIterator + 1] == 1) {
                  myPrettyPrint(gap0, gap1, gap2, gap3, voutsas, 
                    "NA", "OK")
                }
                else {
                  myPrettyPrint(gap0, gap1, gap2, gap3, voutsas, 
                    "NA", "FAIL")
                }
                timings[t] <- time
                t <- t + 1
                userIterator <- userIterator + 1
            }
        }
    }
    cat(paste("    ......................................................................................\n"))
    cat(paste("\n"))
    if ((missing(z) == FALSE)) {
        if (missing(groundTruth) == TRUE) {
            cat(paste("    Relabelling all methods according to method", 
                nam[1], "..."))
            temp <- z
            for (i in 1:m) {
                temp[i, ] <- order(permutations[[1]][i, ])[z[i, 
                  ]]
            }
            zpivot <- apply(temp, 2, function(y) {
                uy <- unique(y)
                uy[which.max(tabulate(match(y, uy)))]
            })
            myComparison <- compare.clust(zpivot, permutations, 
                z, K)
            best.clusterings <- myComparison$clusters
            similarity.matrix <- myComparison$similarity[1:dimname, 
                1:dimname]
            rownames(best.clusterings) <- nam
            permutations <- myComparison$permutations
        }
        else {
            if (length(groundTruth) != dim(z)[2]) {
                stop(paste("groundTruth size not compatible"))
            }
            cat(paste("    Relabelling all methods according to ground truth ..."))
            temp <- z
            myComparison <- compare.clust(groundTruth, permutations, 
                z, K)
            best.clusterings <- myComparison$clusters
            similarity.matrix <- myComparison$similarity
            rownames(best.clusterings) <- nam
            permutations <- myComparison$permutations
        }
        results <- list(permutations, best.clusterings, timings, 
            similarity.matrix)
        names(results) <- c("permutations", "clusters", "timings", 
            "similarity")
    }
    cat(paste(" done!\n"))
    cat(paste("    Retrieve the", dimname, "permutation arrays by typing:\n"))
    for (i in 1:dimname) {
        cat(paste("        [...]$permutations$\"", nam[i], "\"", 
            sep = "", "\n"))
    }
    cat(paste("    Retrieve the", dimname, "best clusterings: [...]$clusters\n"))
    cat(paste("    Retrieve the", dimname, "CPU times: [...]$timings\n"))
    moustakas <- 0
    if (missing(groundTruth) == FALSE) {
        moustakas <- 1
    }
    cat(paste("    Retrieve the", dimname + moustakas, "X", dimname + 
        moustakas, "similarity matrix: [...]$similarity\n"))
    runFinish <- proc.time() - runStart
    runFinish <- round(as.numeric(runFinish[3]), 1)
    cat(paste("    Label switching finished. Total time: ", runFinish, 
        " seconds.", sep = ""), "\n")
    return(results)
}

ecr <- function (zpivot, z, K)
{
    if (K < max(z)) {
        stop("K should be at least equal to max(z)")
    }
    if (dim(z)[2] != length(zpivot)) {
        if (K < max(z)) {
            stop("zpivot has not the same length with simulated z's")
        }
    }
    n <- length(zpivot)
    m <- dim(z)[1]
    k <- K
    st <- 1:k
    perm <- array(data = NA, dim = c(m, K))
    cost.matrix <- matrix(numeric(k * k), nrow = k, ncol = k)
    s <- 1:n
    out = foreach(iter = 1:m, .combine=rbind) %dopar% {
        alloc <- z[iter, ]
        for (i in 1:k) {
            ind <- which(alloc == i)
            so <- zpivot[ind]
            l <- length(ind)
            for (j in 1:k) cost.matrix[i, j] <- l - length(so[so ==
                j])
        }
        matr <- lp.assign(cost.matrix)$solution
        for (i in 1:k) {
            perm[iter, i] <- st[matr[, i] > 0]
        }

	return(perm[iter,])
    }
    results <- list(out)
    names(results) <- c("permutations")
    return(results)
}

