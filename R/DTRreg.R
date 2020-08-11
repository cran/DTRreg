# Doubly-robust dynamic treatment regimen estimation via dynamic weighted ordinary least squares (dWOLS), assuming linear blip functions and logistic regression for treatment model
#' @export
DTRreg <- function(outcome, blip.mod, treat.mod, tf.mod, data=NULL, method = "gest", weight = "default", var.estim="none", B = 200, M = 0, truncate = 0, verbose = FALSE, interrupt = FALSE, treat.range = NULL, missing = "default", interactive = FALSE, treat.mod.man = NULL, type = "DTR") {
  # if interactive mode chosen, build standardized input from user input
  if (interactive==TRUE & var.estim != "bs.quiet") {
    cat("DTR estimation interactive mode!\n")
    outcome <- get(readline("Enter outcome variable: "))
    K <- readline("Enter number of stages: ")
    for (j in 1:K) {
      var <- c()
      new.var <- ""
      cat("Enter each stage",j,"blip covariate (type STOP when done)\n")
      while (new.var != "STOP") {
        new.var <- readline(": ")
        if (new.var != "STOP") {var <- append(var,new.var)}
      }
      # convert input to formula form
      if (length(var) == 0) {var <- "1"}
      if (length(var) == 1) {blip.mod[[j]] <- paste("~",var)}
      else {blip.mod[[j]] <- paste("~",paste(var,collapse="+"))}
      cat("Stage",j,"blip model:",blip.mod[[j]],"\n")
      blip.mod[[j]] <- as.formula(blip.mod[[j]])
    }
    for (j in 1:K) {
      var <- c()
      new.var <- ""
      cat("Enter each stage",j,"treatment-free covariate (type STOP when done)\n")
      while (new.var != "STOP") {
        new.var <- readline(": ")
        if (new.var != "STOP") {var <- append(var,new.var)}
      }
      if (length(var) == 0) {var <- "1"}
      if (length(var) == 1) {tf.mod[[j]] <- paste("~",var)}
      else {tf.mod[[j]] <- paste("~",paste(var,collapse="+"))}
      cat("Stage",j,"treatment-free model:",tf.mod[[j]],"\n")
      tf.mod[[j]] <- as.formula(tf.mod[[j]])

    }
    for (j in 1:K) {
      var <- c()
      new.var <- ""
      cat("Enter stage",j,"treatment variable\n")
      out <- readline(": ")
      cat("Enter each stage",j,"treatment model covariate (type STOP when done)\n")
      while (new.var != "STOP") {
        new.var <- readline(": ")
        if (new.var != "STOP") {var <- append(var,new.var)}
      }
      if (length(var) == 0) {var <- "1"}
      if (length(var) == 1) {treat.mod[[j]] <- paste(out,"~",var)}
      else {treat.mod[[j]] <- paste(out,"~",paste(var,collapse="+"))}
      cat("Stage",j,"treatment model:",treat.mod[[j]],"\n")
      treat.mod[[j]] <- as.formula(treat.mod[[j]])
    }
    method <- readline("Dynamic WOLS (dwols), G-estimation (gest), or Q-learning (qlearn)?\n> ")
    var.estim <- readline("Variance estimation: none, bootstrap or sandwich?\n> ")
  }
  obj <- list()
  # checking input validity
  try(match.arg(method,c("dwols","gest","qlearn")))
  try(match.arg(var.estim,c("none","bootstrap","sandwich","bs.quiet")))
  # if q-learning selected, then set method to dwols and define qlearn to 1 for later
  qlearn <- 0
  if (method == "qlearn") {
    method <- "dwols"
    qlearn <- 1
  }
  if (var.estim == "sandwich" & method == "dwols") {stop("Sandwich variance estimation only available with G-estimation.\n")}
  if (class(treat.mod) == "formula") {treat.mod <- list(treat.mod); blip.mod <- list(blip.mod); tf.mod <- list(tf.mod)}
  if (length(treat.mod) != length(blip.mod) | length(treat.mod) != length(tf.mod)) {stop("-treat.mod-, -blip.mod- and -tf.mod- must all have same length.")}
  if (truncate < 0 | truncate > 0.5) {stop("-truncate- must lie between 0 and 0.5.")}
  # infer stages from treat.mod
  K <- length(treat.mod)
  obj$K <- K
  # if no data set, create one from the models
  if (is.null(data)) {
    # form list of variables
    model.list <- c(blip.mod, tf.mod, treat.mod)
    # extract variables: note we drop any intercept-only terms
    data <- cbind(outcome, do.call('cbind', lapply(model.list[which(model.list != '~1')], get_all_vars)))
    # remove duplicates
    data <- data[!duplicated(lapply(data,summary))]
    obj$outcome <- outcome
  } else {
        if (var.estim == 'bs.quiet') {
    obj$outcome <- outcome
	} else {
    obj$outcome <- eval(substitute(outcome),envir=data)
	}
  }
  # make sure data is a data frame (can't be a matrix)
  data <- data.frame(data)
  rownames(data) <- seq(1:nrow(data))
  obj$data <- data

  obj$blip.mod <- blip.mod
  obj$treat.mod <- treat.mod
  obj$tf.mod <- tf.mod
  Y <- obj$outcome
  if (any(is.na(Y))) {
    stop("Missing values in the outcome are not allowed.")
  }
  obj$obs.Y <- Y
  obj$blip.list <- obj$psi <- obj$covmat <- list()
  drop <- c()
  N <- length(Y)
  keep <- c(1:N)
  remove <- rep(0,N)
  # work in stages starting from final stage (K)
  for (j in K:1) {
    # store current outcome for use if missingness occurs
    Y.orig <- Y
    # get data: Hpsi = blip, Hbeta = treatment-free, Halpha = treatment
    A <- model.response(model.frame(treat.mod[[j]],data, na.action='na.pass'))
    Hpsi <- model.matrix(blip.mod[[j]],model.frame(blip.mod[[j]],data, na.action='na.pass'))
    Hbeta <- model.matrix(tf.mod[[j]],model.frame(tf.mod[[j]],data, na.action='na.pass'))
    Halpha <- model.matrix(treat.mod[[j]],model.frame(treat.mod[[j]],data, na.action='na.pass'))
    # list of treatment variables
    treat.var <- all.vars(treat.mod[[j]])[1]
    # store blip variable names for output
    obj$blip.list[[j]] <- colnames(Hpsi)
    # missing data: identify if missing in Hpsi, Hbeta, Halpha or A; drop and store sample size
    if (missing == "default" | missing == "ipcw") {
      # pi = vector weights for IPCW
      pi <- rep(1,N)
      # drop adds to current drop so future missingness is accommodated
      new.drop <- which(apply(is.na(cbind(A,Hpsi,Hbeta,Halpha)),1,any))
      # remove keeps track of which patients were dropped at what stage
      if(j==K) remove[new.drop] <- j
      if(j!=K) {
        dropk <- new.drop[!(new.drop %in% drop)]
        remove[dropk] <- j
      }
      drop <- unique(new.drop,drop)
      # if missing, later than stage 1, and requested, calculate IPCW
      if (length(new.drop) > 0 & missing == "ipcw" & j > 1) {
        # this works on 'keep' from previous stage: full population who *could* be missing now
        Miss <- rep(1,N)
        Miss[new.drop] <- 0
        Miss <- Miss[keep]
        # build full history from all previous stages
        H.list <- unique(c(unlist(sapply(blip.mod[(j-1):1],all.vars)),unlist(sapply(tf.mod[(j-1):1],all.vars)),unlist(sapply(treat.mod[(j-1):1],all.vars))))
        H <- data.matrix(data[,which(colnames(data) %in% H.list)])
        H <- H[keep,]
        pi[which(!apply(is.na(H),1,any))] <- 1 - fitted(glm(Miss~.,data=as.data.frame(H), family = 'binomial'))
        obj$ipcw[[j]] <- 1/pi
      }

	if (missing == 'ipcw') {
	keep <- c(1:N)
      keep <- keep[!(keep %in% new.drop)]				
	} else {
      keep <- keep[!(keep %in% drop)]		
	}

      # only work with those we want to kep at this stage
      pi <- pi[keep]
      A <- A[keep]
      Hpsi <- Hpsi[keep,]
      Hbeta <- Hbeta[keep,]
      Halpha <- Halpha[keep,]
      Y <- Y[keep]
    }
    # force matrices in case intercept-only models
    Hpsi <- as.matrix(Hpsi)
    Hbeta <- as.matrix(Hbeta)
    Halpha <- as.matrix(Halpha)
    obj$n[[j]] <- length(keep)
    n <- length(Y)
    # check if treatments are binary, if not, and dWOLS selected, abort
    A.bin <- as.numeric(length(unique(A)) <= 2)
    if (A.bin == 0 & method %in% c("dwols","dWOLS")) {
      stop("Non-binary treatment only suitable for G-estimation analysis.\n")
    }
    # if binary and not 0/1, recode (and print warning)
    if (A.bin & !all(unique(A) %in% c(0,1))) {
      A <- ifelse(A==min(A),0,1)
      # suppress output if bootstrapping
      if (var.estim!="bs.quiet") {
        cat("Warning: stage",j,"treatment recoded to 0/1.\n")
      }
    }
    # if no treat.range but outcome non-binary, set treat range to min/max of A
    if (A.bin == 0 & is.null(treat.range)) {
      treat.range <- c(min(A),max(A))
      if (var.estim!="bs.quiet") {
        cat("Warning: no treatment range specified, [min,max] of observed treatment used.\n")
      }
    }
    obj$treat.range[[j]] <- treat.range
    # treatment model: if binary logistic, if non-binary linear
    if (A.bin == 1) {
      obj$cts[[j]] <- "bin"
      # note that if user-specified we take specified P(A = 1) instead
      if (class(treat.mod.man) == "list") {
        if (class(treat.mod.man[[j]]) != "formula") {
          Ahat <- treat.mod.man[[j]]
        }
      } else {
        alpha <- glm(A~Halpha,binomial)
        obj$treat.mod.fitted[[j]] <- alpha
        Ahat <- fitted(alpha)
      }
    } else {
      alpha <- lm(treat.mod[[j]],data)
      obj$treat.mod.fitted[[j]] <- alpha
      Ahat <- fitted(alpha)
      if (treat.var %in% colnames(Hpsi)){
        # continuous with quadratic blip
        obj$cts[[j]] <- "cts.q"
      } else {
        # continuous with linear blip
        obj$cts[[j]] <- "cts.l"
      }
    }
    # if dWOLS
    if (method == "dwols") {
      # weights
      if (is.function(weight) == FALSE) {
        w <- abs(A - Ahat)
      } else {
        # weight function specified for pi = 1, therefore
        w <- weight(Ahat)
        w[A == 0] <- w[A==0]*Ahat[A==0]/(1-Ahat[A==0])
      }
      # if q-learning selected, set all weights to 1
      if (qlearn == 1) {w <- rep(1,n)}
      # update with IPCW
      w <- w*pi
      # blip parameter estimates extracted via dimensions of Hbeta and Hpsi
      est <- solve(t(cbind(Hbeta,A*Hpsi))%*%cbind(w*Hbeta,w*A*Hpsi)) %*% t(cbind(Hbeta,A*Hpsi)) %*% (w*Y)
      psi <- est[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])]
      beta <- est[1:dim(Hbeta)[2]]
    }
    # g-estimation
    if (method == "gest") {
      Hd <- cbind(Hbeta,A*Hpsi)
      if (A.bin == 1) {
        # binary treatment
        Hw <- cbind(Hbeta,Hpsi*(A - Ahat)*pi)
      } else {
        # continuous treatment, requires some fiddly cleanup
        # ***INPUT CLEANUP STARTS
        # separate Hpsi into A/not-A components (if applicable)
        p1.list <- p2.list <- c()
        for (p.var in colnames(Hpsi)) {
          # note there may be interaction terms, so have to separate on :
          if (treat.var %in% strsplit(p.var,split=":")[[1]]) {
            p2.list <- append(p2.list,which(colnames(Hpsi) == p.var))} 						else {
              p1.list <- append(p1.list,which(colnames(Hpsi) == p.var))
            }
        }
        Hpsi1 <- as.matrix(Hpsi[,p1.list])
        # if we have A terms need to extract what they're multiplied with
        # 3 options: A by itself, or an A interaction as :A or A:
        psi2.list <- c()
        for (p.var in colnames(Hpsi)[p2.list]) {
          # if p.var is treat.var, we do nothing (assumed intercept term)
          if (p.var != treat.var) {
            # if variable starts with A: then replace with nothing
            if (strsplit(p.var,split=":")[[1]][1] == treat.var) {psi2.list <- append(psi2.list,paste(strsplit(p.var,split=":")[[1]][-1],collapse=":"))}
            # if variable ends with :A then replace with nothing
            if (tail(strsplit(p.var,split=":")[[1]],n=1) == treat.var) {psi2.list <- append(psi2.list,paste(head(strsplit(p.var,split=":")[[1]],n=-1),collapse=":"))}
            # if variable has :A: in the middle, replace this with :
            if (strsplit(p.var,split=paste(":",treat.var,":",sep=""))[[1]][1] != p.var) {psi2.list <- append(psi2.list,paste(strsplit(p.var,split=paste(":",treat.var,":",sep=""))[[1]],collapse=":"))}
          }
        }
        # turn psi2.list into a model so we can get a new data frame
        # intercept always implied
        if (length(psi2.list) > 0) {psi2.mod <- paste("~1+",paste(psi2.list,collapse="+")); obj$b.list[[j]] <- psi2.list} else {psi2.mod <- "~1"}
        if (length(p1.list) > 1) {psi1.mod <- paste("~1+",paste(colnames(Hpsi)[p1.list[-1]],collapse="+"))} else {psi1.mod <- "~1"}
        # store both components of the blip model for prediction function
        obj$blip.list.cts[[j]] <- c(psi1.mod,psi2.mod)
        Hpsi2 <- model.matrix(as.formula(psi2.mod),data)
        # ***INPUT CLEANUP ENDS
        A2 <- A^2
        A2hat <- (summary(alpha)$sigma)^2 + Ahat^2
        # if no A^2 terms, just use Hpsi
        ifelse(length(p2.list) > 0,Hw <- cbind(Hbeta,Hpsi1*(A-Ahat),Hpsi2*(A2-A2hat)),Hw <- cbind(Hbeta,Hpsi*(A-Ahat)))
      }
      # get estimates
      est <- solve(t(Hw)%*%Hd) %*% t(Hw) %*% (Y)
      # blip parameters
      psi <- est[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])]
      # non-blip parameters
      beta <- est[1:dim(Hbeta)[2]]
    }
    # use estimates to identify optimal treatments
    # if treatment is binary:
    if (A.bin == 1) {opt <- as.numeric(Hpsi %*% psi > 0)}
    # if treatment is continuous
    if (A.bin == 0) {
      if (obj$cts[[j]] == "cts.l") {
        # if no treatment term within blip, then max/min based on sign
        opt <- as.numeric(Hpsi %*% psi > 0)*max(treat.range) + as.numeric(Hpsi %*% psi <= 0)*min(treat.range)
      } else {
        # otherwise optimal is found by maximizing wrt A, need to check second derivative
        # if sign(Hpsi2 %*% psi2) < 0, then blip is maxed by setting first derivative to zero
        # otherwise evaluate at max/min or treatment range and pick the larger
        # dif is blip(minA) - blip(maxA); so if negative take max, otherwise min
        sec.deriv <- Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]
        dif <- sign(min(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + min(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]) - max(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + max(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]))
        opt <- as.numeric(sec.deriv < 0)*(-(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])])/(2*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)])) + as.numeric(sec.deriv >= 0)*(as.numeric(dif >= 0)*min(treat.range) + as.numeric(dif < 0)*max(treat.range))
      }
      # make sure opt is between limits
      opt <- pmin(pmax(opt,min(treat.range)),max(treat.range))
    }
    # ensure opt is in vector form
    opt <- as.vector(opt)
    # store estimates, optimal treatments and blip parameters
    obj$beta[[j]] <- beta
    obj$psi[[j]] <- psi
    obj$opt.treat[[j]] <- opt
    # Sandwich variance estimator done stage-by-stage
    if (var.estim == "sandwich") {
      D <- A - Ahat
      # HD: what precedes the (Y-gamma-beta) term in the est eq.
      if (method == "gest") {HD <- Hpsi*(A-Ahat)}
      if (method == "dwols") {HD <- A*Hpsi*abs(A-Ahat)}
      # residuals are the (Y-gamma-beta) residuals
      residuals <- as.vector(Y - ((A*Hpsi)%*%psi) - Hbeta%*%beta)
      # estimating equation = (Y - gamma - E[Y - gamma|H])*(A - Ahat)*Hpsi
      U <- residuals * HD
      # dU/dpsi:
      M <- (1/n) * -crossprod(HD, A * Hpsi)
      # if binary
      if (A.bin == 1) {
        score.alpha <- D * Halpha
        score.varsig <- residuals * Hbeta
        U.adj <- U
        # dT.dalpha = Ahat(1-Ahat) Halpha Halpha
        dT.dalpha <- (1/n)*(-crossprod(Ahat*(1-Ahat)*Halpha, Halpha))
        dU.dalpha <- (1/n)*(-crossprod(Ahat*(1-Ahat)*residuals*Hpsi,Halpha))
        U.adj <- U.adj - t(dU.dalpha%*%solve(dT.dalpha) %*% t(score.alpha))
        # now with respect to beta
        dT.dvarsig <- (1/n) * (-crossprod(Hbeta))
        dU.dvarsig <- (1/n) * (-crossprod(HD, Hbeta))
        U.adj <- U.adj-t(dU.dvarsig%*%solve(dT.dvarsig)%*%t(score.varsig))
        if (j == K) {U.store <- M.store <- deriv.regret.store <- list()}
        U.store[[j]] <- matrix(nrow=N,ncol=dim(U)[2])
        deriv.regret.store[[j]] <- matrix(nrow=N,ncol=dim(Hpsi)[2])
        U.store[[j]][keep,] <- U
        M.store[[j]] <- M
        # derivative of regrets wrt psi
        deriv.regret.store[[j]][keep,] <- (opt - A) * Hpsi
        if (j < K) {
          for (l in (j+1):K) {
            dU2.dpsi2 <- -M.store[[l]]
            dU1.dpsi2<-(1/n)*(crossprod(HD,deriv.regret.store[[l]][keep,]))
            U.adj <- U.adj - t(dU1.dpsi2 %*% solve(dU2.dpsi2) %*% t(U.store[[l]][keep,]))
          }
        }
      }
      # if linear, then score = (1/s^2) * (A - Ahat)*Halpha
      if (A.bin == 0) {
        score.alpha <- D * Halpha
        score.varsig <- residuals * Hbeta
        U.adj <- U
        # dT.dalpha = -1/s^2 t(Halpha) * Halpha
        dT.dalpha <- (1/n) * (-crossprod(Halpha,Halpha))
        # dU/dalpha = (y - gamma - G)*Hpsi*[...], if no A^2 term, [...] is -Halpha
        # A^2 term then also -2*Halpha*alpha-hat - 2/(n-2)*sum(Halpha*(Ahat-A))
        # N.B. this slight ugliness is because of the var(A-hat|...) term
        if (obj$cts[[j]] == "cts.l") {
          dU.dalpha <- (1/n) * (-crossprod(residuals * Hpsi, Halpha))
        } else {
          Hpsi.cts <- cbind(Hpsi1,Hpsi2*2*(Ahat + (1/(dim(Hbeta)[1]-2))*sum(Ahat-A)))
          dU.dalpha <- (1/n) * (-crossprod(residuals * Hpsi.cts,Halpha))
        }
        U.adj <- U.adj - t(dU.dalpha %*% solve(dT.dalpha) %*% t(score.alpha))#
        dT.dvarsig <- (1/n) * (-crossprod(Hbeta))
        dU.dvarsig <- (1/n) * (-crossprod(HD, Hbeta))
        U.adj <- U.adj - t(dU.dvarsig %*% solve(dT.dvarsig) %*% t(score.varsig))
        # storage for earlier intervals
        if (j == K) {U.store <- M.store <- deriv.regret.store <- list()}
        U.store[[j]] <- U
        M.store[[j]] <- M
        # derivative of regrets wrt psi, depends on whether A term in reg
        if (obj$cts[[j]] == "cts.l") {deriv.regret.store[[j]] <- (opt - A) * Hpsi}
        if (obj$cts[[j]] == "cts.q") {deriv.regret.store[[j]] <- opt*(cbind(Hpsi1,opt*Hpsi2))-A*Hpsi}
        if (j < K) {
          for (l in (j+1):K) {
            dU2.dpsi2 <- -M.store[[l]]
            dU1.dpsi2 <- (1/n) * (crossprod(HD, deriv.regret.store[[l]]))
            U.adj <- U.adj - t(dU1.dpsi2 %*% solve(dU2.dpsi2) %*% t(U.store[[l]]))
          }
        }
      }
      IF <- t(solve(M, t(U.adj)))
      Sigma <- (1/n) * var(IF)
      obj$covmat[[j]] <- Sigma
    }
    # update pseudo-Y by adding regret, note that if missing we take Y = Y.orig
    if (obj$cts[[j]] %in% c("cts.l","bin")) {
      # for linear blip
      Y.orig[keep] <- Y.orig[keep] + (as.numeric(type == "DTR")*(opt) - A)*(Hpsi%*%psi)
      Y.orig[drop] <- Y.orig[drop]
      obj$regret[[j]] <- (opt - A)*(Hpsi%*%psi)
    } else {
      # for quadratic blip
      Y.orig[keep] <- Y.orig[keep] + as.numeric(type == "DTR")*(opt)*(Hpsi1 %*% psi[1:length(p1.list)] + as.numeric(type == "DTR")*(opt)*Hpsi2 %*% psi[(length(p1.list)+1):(length(psi))]) - A*(Hpsi%*%psi)
      Y.orig[drop] <- Y.orig[drop]
      obj$regret[[j]] <- opt*(Hpsi1 %*% psi[1:length(p1.list)] + opt*Hpsi2 %*% psi[(length(p1.list)+1):(length(psi))]) - A*(Hpsi%*%psi)
    }
    Y <- Y.orig
    # fitted values at each stage
    Y.fit <- Hbeta %*% beta + A*(Hpsi %*% psi)
    k <- j
    while (k < K) {
      k <- k + 1
      Y.fit[which(rownames(obj$regret[[k]])%in% keep)] <- Y.fit[which(rownames(obj$regret[[k]])%in% keep)] - obj$regret[[k]][which(rownames(obj$regret[[k]])%in% keep)] # correctly use keep indicator
    }
    obj$fitted[[j]] <- Y.fit
  }
  # keep track of which patients are removed at which stage
  obj$remove <- remove
  # optimal outcome according to blip estimates is current pseudo-Y
  obj$opt.Y <- Y
  obj$type <- type
  # standard errors
  if (var.estim=="bootstrap") {
    if (M == 0) {M <- nrow(data)}
    # note time for progress reports
    ptm <- proc.time()
    psi.boot <- list()
    # continue indicator
    cont <- "u"
    for (i in 1:B) {
      b.sample <- sample(1:M, replace=TRUE)
      data.boot <- data[b.sample,]
      # set outcome vector manually to avoid dependency on input format
      outcome.boot <- Y[b.sample]
      psi.boot[[i]] <- DTRreg(outcome.boot, blip.mod, treat.mod, tf.mod, data.boot, method, var.estim="bs.quiet", verbose=verbose,treat.mod.man=treat.mod.man)$psi
      # if verbose, display ETA and give option to abort
      if (verbose == T & i >= 10) {
        # only display if projected > 30s, then only display every 30s
        eta <- (B-i)*((proc.time() - ptm)[3]/i)
        # if very long (> 10 mins) give option to abort but if ignored continue
        if (i > 50 & eta > 600 & cont != "y" & interrupt == T) {
          cont <- readline(paste("Estimated run time",round(eta/60),"minutes. Continue? y/n: "))
          if (cont == "n" | cont == "no" | cont == "NO") {stop("Aborted.\n")} else {cont <- "y"}
        }
        if (i == 10) {last <- eta+31}
        if (eta > 30 & eta < (last-30)) {
          cat("Approximately",round(eta),"seconds remaining.\n")
          last <- eta
        }
      }
    }
    psi.boot <- do.call(function(...) mapply(rbind, ..., SIMPLIFY=FALSE), psi.boot)
    # if requested, truncate estimates based on specified percentile
    if (truncate > 0) {
      for (j in 1:K) {
        psi.boot[[j]] <- apply(psi.boot[[j]],2,DTRreg.trunc,truncate)
      }
    }
    covmat <- lapply(psi.boot, var)
    obj$covmat <- covmat
    obj$psi.boot <- psi.boot
  }
  # if variance estimated (and we're not currently bootstrapping, estimate non-regularity
  if (!(var.estim %in% c("none","bs.quiet")) & A.bin == 1) {
    for (j in K:1) {
      Hpsi <- model.matrix(blip.mod[[j]],data)
      # lower/upper limits on psi
      psi.l <- obj$psi[[j]] - 1.96*sqrt(diag(obj$covmat[[j]]))
      psi.u <- obj$psi[[j]] + 1.96*sqrt(diag(obj$covmat[[j]]))
      # look at max/min value of blip based on sign of covariates and max/min of parameter CIs
      Hpsi.p <- apply(Hpsi,c(1,2),FUN=function(x) max(x,0))
      Hpsi.n <- apply(Hpsi,c(1,2),FUN=function(x) min(x,0))
      # lower blip is positive values * lower CI + negative values * upper CI
      blip.l <- Hpsi.p %*% psi.l + Hpsi.n %*% psi.u
      blip.u <- Hpsi.p %*% psi.u + Hpsi.n %*% psi.l
      obj$nonreg[j] <- sum(blip.l < 0 & blip.u > 0)/dim(Hpsi)[1]
    }
  }
  # return
  class(obj) <- "DTRreg"
  obj
}

# printing
#' @export
print.DTRreg <- function(x,...) {
	K <- x$K
	# longest blip variable name lengths, if applicable, get length of blip excluding treatment terms
	max.length <- 0
	p.length <- c()
	for (j in 1:K) {
		max.length <- max(max.length,nchar(x$blip.list[[j]]))
		if (j <= length(x$blip.list.cts)) {if (is.null(x$blip.list.cts[[j]]) == F) {p.length[j] <- length(strsplit(x$blip.list.cts[[j]][1],split="\\+")[[1]])} else {p.length[j] <- 0}}
	}
	if (x$type == "DTR") {cat("DTR estimation over",K,"stages:\n\n")}
	cat("Blip parameter estimates")
	if (is.null(x$covmat[1])) {cat(" (standard errors not requested)")}
	cat("\n")
	# first blank space should be the length of longest variable + 1
	blank.paste <- function(n) {if (n > 0) {for (i in 1:n) {cat(" ")}}}
	blank.paste(max.length)
	# if no standard errors, just print estimates
	if (is.null(x$covmat[1][[1]])) {
		cat(" Estimate\n")
		for (j in 1:K) {
			cat("Stage ",j," (n = ",x$n[[j]],")\n",sep="")
			for (p in 1:length(x$psi[[j]])) {
				blank.paste(max.length-nchar(x$blip.list[[j]][p]))
				cat(x$blip.list[[j]][p])
				blank.paste(3 - as.numeric(x$psi[[j]][p] < 0))
				cat(paste(sprintf("%.4f",x$psi[[j]][p]),"\n",sep=""))
			}
			cat("\n")
		}
	}
	# otherwise, print standard errors and CIs
	if (is.null(x$covmat[1][[1]]) == F) {
		cat(" Estimate Std. Error    95% Conf. Int\n")
		for (j in 1:K) {
			cat("Stage ",j," (n = ",x$n[[j]],")\n",sep="")
			for (p in 1:length(x$psi[[j]])) {
				blank.paste(max.length-nchar(x$blip.list[[j]][p]))
				cat(x$blip.list[[j]][p])
				blank.paste(3 - as.numeric(x$psi[[j]][p] < 0))
				cat(paste(sprintf("%.4f",x$psi[[j]][p]),"     ",sprintf("%.4f",sqrt(x$covmat[[j]][p,p])),sep=""))
				blank.paste(2 - as.numeric(x$psi[[j]][p] < 0))
				cat(paste("[",sprintf("%.4f",x$psi[[j]][p]-1.96*sqrt(x$covmat[[j]][p,p])),",",sprintf("%.4f",x$psi[[j]][p]+1.96*sqrt(x$covmat[[j]][p,p])),"]\n",sep=""))
			}
			cat("\n")
		}
	}
	# print warnings if non-regularity concerns
	if (x$type == "DTR") {
	if (is.null(x$covmat[1][[1]]) == F) {
		for (j in 1:K) {
			if (x$cts[[j]] == "bin") {
			if (x$nonreg[j] > 0.05) {cat("Warning: possible non-regularity at stage ",j," (prop = ",round(x$nonreg[j],3),")\n",sep="")}
		}
		}
		cat("\n")
	}
	# print DTRregs explicitly
	if (!("cts.q" %in% x$cts)) {cat("Recommended dynamic treatment regimen:\n")}
	# if binary just treat/don't treat
	for (j in 1:K) {
		if (x$cts[[j]] == "bin") {
			cat("Stage ",j,": treat if ",sprintf("%.4f",x$psi[[j]][1]),sep="")
			if (length(x$psi[[j]]) > 1) {
				for (p in 2:length(x$psi[[j]])) {
					if(x$psi[[j]][p] >= 0) {cat(" + ")} else {cat(" - ")}
					cat(sprintf("%.4f",abs(x$psi[[j]][p])),x$blip.list[[j]][p])
				}
			}
		cat(" > 0\n")
		} else {
			# set to max/min of range depending on sign of blip
			if (x$cts[[j]] == "cts.l") {
				cat("Stage ",j,": maximum treatment if ",sprintf("%.4f",x$psi[[j]][1]),sep="")
				if (length(x$psi[[j]]) > 1) {
					for (p in 2:length(x$psi[[j]])) {
						if(x$psi[[j]][p] >= 0) {cat(" + ")} else {cat(" - ")}
						cat(sprintf("%.4f",abs(x$psi[[j]][p])),x$blip.list[[j]][p])
					}
				}
			cat(" > 0, otherwise minimum treatment.\n")
			}
		}
	}
	}
}
#' @export
summary.DTRreg <- function(object,...) {
	print.DTRreg(object)
}

# optimal outcome prediction for new data ASSUMING CORRECT MODELS
#' @export
predict.DTRreg <- function(object, newdata, treat.range=NULL, ...) {
	x <- object
	Hpsi <- model.matrix(x$blip.mod[[1]],newdata)
	Hbeta <- model.matrix(x$tf.mod[[1]],newdata)
	if (is.null(treat.range) & x$cts[[1]] != "bin") {treat.range <- x$treat.range[[1]]}
	# options: binary, cts with range, or cts without range
	# binary:
	if (x$cts[[1]] == "bin") {
		opt <- as.numeric(Hpsi%*%x$psi[[1]] > 0)
		obj <- Hbeta%*%x$beta[[1]] + opt*(Hpsi%*%x$psi[[1]])
	}
	# cts with no A^2 term in blip
	if (x$cts[[1]] == "cts.l") {
		opt <- as.numeric(Hpsi%*%x$psi[[1]] > 0)*max(treat.range) + as.numeric(Hpsi%*%x$psi[[1]] <= 0)*min(treat.range)
		obj <- Hbeta%*%x$beta[[1]] + opt*(Hpsi%*%x$psi[[1]])
	}
	# cts with A^2 term in blip
	if (x$cts[[1]] == "cts.q") {
		Hpsi1 <- model.matrix(as.formula(x$blip.list.cts[[1]][1]), newdata)
		Hpsi2 <- model.matrix(as.formula(x$blip.list.cts[[1]][2]), newdata)
		psi <- x$psi[[1]]
		sec.deriv <- Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]

		dif <- sign(min(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + min(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]) - max(treat.range)*(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])] + max(treat.range)*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)]))
		opt <- as.numeric(sec.deriv >= 0)*(-(Hpsi1 %*% psi[1:(dim(Hpsi1)[2])])/(2*Hpsi2 %*% psi[(dim(Hpsi1)[2]+1):length(psi)])) + as.numeric(sec.deriv < 0)*(as.numeric(dif >= 0)*min(treat.range) + as.numeric(dif > 0)*max(treat.range))

		opt <- pmin(pmax(opt,min(treat.range)),max(treat.range))
		obj <- Hbeta%*%x$beta[[1]] + opt*(Hpsi1%*%psi[1:dim(Hpsi1)[2]] + Hpsi2%*%psi[(dim(Hpsi1)[2]+1):length(psi)])
	}
	return(obj)
}
#' @export
coef.DTRreg <- function(object,...) {
	x <- object
	K <- x$K
	obj <- x$psi
	for (j in 1:K) {
		names(obj)[j] <- paste("stage",j,sep="")
		names(obj[[j]]) <- x$blip.list[[j]]
	}
	return(obj)
}
# confidence interval
#' @export
confint.DTRreg <- function(object, parm = NULL, level = 0.95, type = "se", ...) {
  x <- object
  K <- x$K
  obj <- x$psi
  if (is.null(x$covmat[1][[1]]) == F) {
    if(type == "se"){
      for (j in 1:K) {
        names(obj)[j] <- paste("stage",j,sep="")
        obj[[j]] <- matrix(nrow=length(x$psi[[j]]),ncol=2)
        rownames(obj[[j]]) <- x$blip.list[[j]]
        colnames(obj[[j]]) <- c(paste(100*((1 - level)/2),"%"), paste(100*((1 + level)/2),"%"))
        for (p in 1:length(x$psi[[j]])) {
          obj[[j]][p,] <- c(x$psi[[j]][p] + qnorm((1-level)/2) * sqrt(x$covmat[[j]][p,p]), x$psi[[j]][p] + qnorm((1+level)/2)*sqrt(x$covmat[[j]][p,p]))
        }
      }
    }else if(type == "percentile"){
      if(length(x$psi.boot) == 0){
        stop("Percentile confidence intervals only available with bootstrap.")
      }else{
        for (j in 1:K) {
          names(obj)[j] <- paste("stage",j,sep="")
          obj[[j]] <- matrix(nrow=length(x$psi[[j]]),ncol=2)
          rownames(obj[[j]]) <- x$blip.list[[j]]
          colnames(obj[[j]]) <- c(paste(100*((1 - level)/2),"%"), paste(100*((1 + level)/2),"%"))
          for (p in 1:length(x$psi[[j]])) {
            obj[[j]][p,] <- c(quantile(x$psi.boot[[j]][,p], probs = (1-level)/2), quantile(x$psi.boot[[j]][,p], probs = (1+level)/2))
          }
        }
      }
    }
    return(obj)
  } else {
    cat("confint() only available when DTRreg called with a variance estimation option.\n")
  }
}

# diagnostics
#' @export
plot.DTRreg <- function(x, method = "DTRreg", ...) {
  # original Y is in x$obs.Y
  # fitted Y (one for each stage) is in x$fitted
  if(method == "DTRreg"){
    cat("Y residuals (for treatment-free and blip models):\n")
    for (j in 1:x$K) {
      advance <- readline("Hit <Return> to see next plot:")
      obsK <- x$obs.Y[x$remove<j]
      fitK <- x$fitted[[j]]
      plot(fitK,obsK-fitK,xlab="Fitted values",ylab="Residuals",main=paste("Stage",j,"Residuals vs Fitted"))
      abline(h=0,lty=2,col=4)
      points(lowess(fitK,obsK - fitK),type="l",col=2,lwd=2)
      if (length(x$blip.list[[j]]) > 1) {
        for (var in x$blip.list[[j]][2:length(x$blip.list[[j]])]) {
          advance <- readline("Hit <Return> to see next plot:")
          plot(x$data[x$remove<j,which(colnames(x$data) == var)],obsK-fitK,xlab=paste(var),ylab="Residuals",main=paste("Stage",j,"Residuals vs",var))
          abline(h=0,lty=2,col=4)
          points(lowess(x$data[x$remove<j,which(colnames(x$data) == var)],obsK-fitK),type="l",col=2,lwd=2)
        }
      }
    }
    cat("Treatment models:\n")
    for (j in 1:x$K) {
      cat("Stage",j,"\n")
      plot(x$treat.mod.fitted[[j]],sub.caption=paste("Stage",j,"treatment"))
    }
  }else if(method == "DWSurv"){
    cat("Y residuals (for treatment-free and blip models):\n")
    for (j in 1:x$K) {
      advance <- readline("Hit <Return> to see next plot:")
      obsK <- x$log.obs[[j]]
      fitK <- x$log.fitted[[j]]
      plot(fitK,obsK-fitK,xlab="Fitted values",ylab="Residuals",main=paste("Stage",j,"Residuals vs Fitted"))
      abline(h=0,lty=2,col=4)
      points(lowess(fitK,obsK - fitK),type="l",col=2,lwd=2)
    }
    cat("Treatment models:\n")
    for (j in 1:x$K) {
      advance <- readline("Hit <Return> to see next plot:")
      cat("Stage",j,"\n")
      plot(x$treat.mod.fitted[[j]],sub.caption=paste("Stage",j,"treatment"))
    }
    cat("Censoring models:\n")
    for (j in 1:x$K) {
      advance <- readline("Hit <Return> to see next plot:")
      cat("Stage",j,"\n")
      plot(x$cens.mod.fitted[[j]],sub.caption=paste("Stage",j,"censoring"))
    }
  }
}

# truncation function for truncate option
DTRreg.trunc <- function(x,p) {
	l <- quantile(x,p); u <- quantile(x,1-p)
	x[x < l] <- l; x[x > u] <- u
	x
}

# get_all_vars modified to accommodate intercept-only components
get_all_vars_DTR <- function(x,n) {
	if (dim(get_all_vars(x))[1] == 0) {
		return(rep(1,n))
	} else {
		return(get_all_vars(x))
	}
}

#' @export
chooseM <- function(outcome, blip.mod, treat.mod, tf.mod, data=NULL, method = "gest", weight = "default", missing = "default", treat.mod.man = NULL, B1 = 500, B2 = 500)
{
	# fit initial DTRreg model to the data to get original estimates + phat
	mod1 <- DTRreg(outcome,blip.mod, treat.mod, tf.mod, data, method, weight, missing, treat.mod.man, var.estim="bootstrap", M = 0, truncate = 0, B = 50)
	blip.psi.1 <- mod1$psi[[1]]
	p <- mod1$nonreg[2]

	# make sure the data are in a data frame
	if (is.null(data)) {
    data <- cbind(outcome,do.call("cbind",lapply(blip.mod,get_all_vars)),do.call("cbind",lapply(treat.mod,get_all_vars)),do.call("cbind",lapply(tf.mod,get_all_vars)))
    data <- data[!duplicated(lapply(data,summary))]
	}
  	data <- data.frame(data)

	# deal with missing data here

	# initialize
	alpha <- 0.025
	n <- nrow(data)

	while(alpha <= 0.5) # loop until alpha reaches 0.5, with break when the coverage reaches the nominal 95% rate
	{
		# reset coverage
		coverage <- 0

		# reset matrix to save estimates
    	est <- matrix(NA, ncol = B2 + 3, nrow = 1)
    	for(j in 1:B1) # loop over B1 first stage bootstrap samples
    	{
     		# draw a n-out-of-n bootstrap sample
      		index <- sample(1:n, n, replace = TRUE)
      		boot1 <- as.data.frame(data[index,])

      		# if prespecified probability of treatment with treat.mod.man, adjust
      		treat.mod.manN <- NULL
      		if(class(treat.mod.man)=="list")
      		{
      			treat.mod.manN <- treat.mod.man
      			numT <- length(treat.mod.man)
      			for(t in 1:numT) treat.mod.manN[[t]] <- treat.mod.man[[t]][index]
      		}

      		# fit the model to b1-th bootstrap sample
      		res1 <- try(DTRreg(outcome, blip.mod, treat.mod, tf.mod, treat.mod.man = treat.mod.manN, method, var.estim = "bootstrap", data = boot1))
      		esb1 <- res1$psi[[1]][1] # only consider main effect of treatment
      		est[1] <- esb1

      		# estimate m for each b1 bootstrap sample
      		phat <- res1$nonreg[2]
      		est[2] <- phat

      		# resampling size
      		m <- n^((1 + alpha*(1-phat))/(1 + alpha))
      		est[3] <- m

           	for(k in 1:B2) # loop over B2 second stage bootstrap samples
      		{
        		# resample with replacement
        		index <- sample(1:n, floor(m), replace = TRUE)
        		boot2 <- boot1[index,]

        		# if prespecified probability of treatment with treat.mod.man
        		treat.mod.manM <- NULL
      			if(class(treat.mod.man)=="list")
      			{
      				treat.mod.manM <- treat.mod.manN
      				numT <- length(treat.mod.man)
      				for(t in 1:numT) treat.mod.manM[[t]] <- treat.mod.manN[[t]][index]
      			}

        		# fit the model to bootstrap sample
        		res2 <- try(DTRreg(outcome, blip.mod, treat.mod, tf.mod, treat.mod.man = treat.mod.manM, method, data = as.data.frame(boot2)))
        		esb2 <- res2$psi[[1]][1]

        		# save the (b1,b2) bootstrap estimates in the (k+3) column
        		est[k + 3] <- esb2
      		}
      		quan <- quantile(sqrt(m)*(est[4:(B2+3)]-esb1), probs = c(0.025,0.975))
      		coverage <- coverage + ((esb1 - quan[2]/sqrt(m) <= blip.psi.1[1]) & (esb1 - quan[1]/sqrt(m) >= blip.psi.1[1]))
    	}
		coverage <- coverage/B1
		# print coverage and alpha
		cat("alpha=",alpha,", coverage = ", coverage, "\n")
		if(coverage >= 0.95)
    	{
      		cat("Final value of alpha = ", alpha,"\n")
      		m <- n^((1 + alpha*(1-p))/(1 + alpha))
      		cat("Selected subsample size m = ", floor(m),"\n")
      		break
    	} else {
       		alpha <- alpha + 0.025
    	}
	}
	return(list(m=as.numeric(floor(m))))
}


# Doubly-robust dynamic treatment regimen estimation via dynamic weighted survival modeling (DWSurv), assuming linear blip and treatment-free functions and logistic regression for treatment and censoring models
#' @export
DWSurv <- function(time, blip.mod, treat.mod, tf.mod, cens.mod, data = NULL, weight = "default", var.estim = "none", asymp.opt = "adjusted", boot.opt = "standard", B = 500, optimization = "max", quiet = FALSE) {
  obj <- list()
  # for compatibility with other functions
  obj$type <- "DTR"
  # save formula
  obj$time.mod <- time
  obj$blip.mod <- blip.mod
  obj$treat.mod <- treat.mod
  obj$tf.mod <- tf.mod
  obj$cens.mod <- cens.mod
  # checking input validity
  try(match.arg(var.estim, c("none", "asymptotic", "bootstrap")))
  if(var.estim == "asymptotic") try(match.arg(asymp.opt, c("adjusted", "naive")))
  if(var.estim == "bootstrap") try(match.arg(boot.opt, c("standard", "empirical", "normal")))
  try(match.arg(optimization, c("min", "max")))
  if (length(treat.mod) != length(blip.mod) | length(treat.mod) != length(tf.mod) | length(treat.mod) != length(cens.mod) | length(treat.mod) != length(time)) {stop("-treat.mod-, -blip.mod-, -tf.mod-, -cens.mod- and -time- must all have same length.")}
  # infer stages from blip.mod
  K <- length(blip.mod)
  obj$K <- K
  obj$optimization <- optimization
  # checking input validity (continued)
  for(j in 1:K){
    if(class(treat.mod[[j]]) != "formula") stop(paste("The stage ", j, " treatment model must be supplied as a formula with the treatment variable on the left hand side.", sep = ""))
    if(class(cens.mod[[j]]) != "formula") stop("The stage ", j, " censoring model must be supplied as a formula with the censoring indicator on the left hand side.")
    if(class(time[[j]]) != "formula" | length(all.vars(time)) > 1) stop("The stage ", j, " survival time must be supplied as a formula with nothing on the left hand side and the survival time variable on the right hand side.")
  }
  #time <- list(time); treat.mod <- list(treat.mod); blip.mod <- list(blip.mod); tf.mod <- list(tf.mod); cens.mod <- list(cens.mod)
  
  # verify if cens.mod supplied as no censoring
  nocens <- FALSE
  if(length(cens.mod[[1]])==2) nocens <- TRUE
  
  # if no data set, create one from the models
  if (is.null(data)) {
    data <- cbind(do.call("cbind",lapply(time,get_all_vars)), do.call("cbind",lapply(blip.mod,get_all_vars)),do.call("cbind",lapply(treat.mod,get_all_vars)),do.call("cbind",lapply(tf.mod,get_all_vars)), do.call("cbind",lapply(cens.mod,get_all_vars)))
    data <- data[!duplicated(lapply(data,summary))]
    obj$time <- lapply(time,get_all_vars)
    obj$status <- do.call("cbind",lapply(cens.mod,get_all_vars))[,1]
  } else {
    obj$time <- list()
    for(j in 1:K){
      obj$time[[j]] <- data[,which(colnames(data) == all.vars(time[[j]]))]
    }
    obj$status <- data[,which(colnames(data) == all.vars(cens.mod[[1]])[1])] 
  }
  # make sure data is a data frame (can't be a matrix)
  data <- data.frame(data)
  data$id <- seq(1:nrow(data))
  obj$data <- data
  Y <- obj$time
  delta <- obj$status
  
  if (any(is.na(delta))) { # missing survival times will be normal if patients don't reach stage
    stop("Missing values in status are not allowed.")
  }
  obj$blip.list <- obj$psi <- obj$covmat <- obj$covmat.all <- list() # to output list of covariates in blip functions, list of estimates of psi, var-covar matrix at each stage
  N <- nrow(data)
  # calculate eta j i.e. whether an individual enters stage j and remove individuals with missing data
  eta <- list()
  drop <- c()
  for(j in 1:K) {
    if(j > 1){
      eta[[j]] <- ifelse(is.na(Y[[j]]) | Y[[j]] == 0 | eta[[j-1]] == 0, 0, 1)
    }else{
      eta[[j]] <- ifelse(is.na(Y[[j]]) | Y[[j]] == 0, 0, 1)
    }
    dataetaj <- data[eta[[j]] == 1,]
    A <- model.response(model.frame(treat.mod[[j]], dataetaj, na.action = 'na.pass'))
    Hpsi <- model.matrix(blip.mod[[j]], model.frame(blip.mod[[j]], dataetaj, na.action = 'na.pass'))
    Hbeta <- model.matrix(tf.mod[[j]], model.frame(tf.mod[[j]], dataetaj, na.action = 'na.pass'))
    Halpha <- model.matrix(treat.mod[[j]], model.frame(treat.mod[[j]], dataetaj, na.action = 'na.pass'))
    if(nocens == FALSE) {
      Hlambda <- model.matrix(cens.mod[[j]], model.frame(cens.mod[[j]], dataetaj, na.action = 'na.pass'))
      drop <- append(drop,dataetaj$id[which(apply(is.na(cbind(A, Hpsi, Hbeta, Halpha, Hlambda)), 1, any))])
    }else{
      drop <- append(drop,dataetaj$id[which(apply(is.na(cbind(A, Hpsi, Hbeta, Halpha)), 1, any))])
    }
    drop <- unique(drop)
    if(length(drop) > 0) {
      eta[[j]][which(data$id %in% drop)] <- 0
      Y[[j]][which(data$id %in% drop)] <- 0
    }
  }
  
  # keep track of id removed
  obj$drop <- drop
  # keep track of who enter stages
  obj$eta <- eta
  # store quantities for variance estimation
  A.store <- Hpsi.store <- U.adj.store <- HD.store <- Htheta.store <- Halpha.store <- ds.dalpha.store <- score.alpha.store <- list()
  delta.store <- w.store <- Ahat.store <- Hlambda.store <- ds.dlambda.store <- score.lambda.store <- Dhat.store <- list()
  # save bootstrap estimates across stages if necessage
  obj$psi.boot <- list()
  # work in stages starting from final stage (K)
  Ycumul <- rep(0, nrow(data))
  for (j in K:1) {
    Ycumul[eta[[j]] == 1] <- Ycumul[eta[[j]] == 1] + Y[[j]][eta[[j]] == 1]
    Yj <- Ycumul[eta[[j]] == 1]
    dataetaj <- data[eta[[j]] == 1, ]
    deltaj <- delta.store[[j]] <- delta[eta[[j]] == 1]
    
    # get data: Hpsi = blip, Hbeta = treatment-free, Halpha = treatment, Hlambda = censoring
    A <- A.store[[j]] <- as.numeric(model.response(model.frame(treat.mod[[j]], dataetaj, na.action = 'na.pass')))
    Hpsi <- model.matrix(blip.mod[[j]], model.frame(blip.mod[[j]], dataetaj, na.action = 'na.pass'))
    Hbeta <- model.matrix(tf.mod[[j]], model.frame(tf.mod[[j]], dataetaj, na.action = 'na.pass'))
    Halpha <- model.matrix(treat.mod[[j]], model.frame(treat.mod[[j]], dataetaj, na.action = 'na.pass'))
    if(nocens == FALSE) Hlambda <- model.matrix(cens.mod[[j]], model.frame(cens.mod[[j]], dataetaj, na.action = 'na.pass'))
    
    # list of censoring and treatment variables
    treat.var <- all.vars(treat.mod[[j]])[1]
    
    # store blip variable names for output
    obj$blip.list[[j]] <- colnames(Hpsi)
    
    # force matrices in case intercept-only models
    Hpsi <- Hpsi.store[[j]] <-  as.matrix(Hpsi)
    Hbeta <- as.matrix(Hbeta)
    Halpha <- Halpha.store[[j]] <- as.matrix(Halpha)
    if(nocens == FALSE) Hlambda <- Hlambda.store[[j]] <- as.matrix(Hlambda)
    obj$n[[j]] <- nrow(dataetaj)
    
    # check if treatments are binary, if not, and dWOLS selected, abort
    A.bin <- as.numeric(length(unique(A)) <= 2)
    if (A.bin == 0) {
      stop("Non-binary treatment not suitable for DWSurv analysis.\n")
    }
    # if binary and not 0/1, recode (and print warning)
    if (A.bin & !all(unique(A) %in% c(0, 1))) {
      A <- A.store[[j]] <- ifelse(A == min(A),0,1)
      # suppress output if bootstrapping
      if (quiet == FALSE) {
        cat("Warning: stage", j, "treatment recoded to 0/1.\n")
      }
    }
    # for compatibility with other functions
    obj$cts[[j]] <- "bin"
    
    # treatment model: binary logistic
    alpha <- glm(A ~ -1 + Halpha, family = binomial(link = "logit"))
    obj$treat.mod.fitted[[j]] <- alpha
    Ahat <- Ahat.store[[j]] <- fitted(alpha)
    obj$obs.treat[[j]] <- A
    
    # censoring model: binary logistic
    if(sum(deltaj) != length(deltaj) & nocens == FALSE){
      lambda <- glm(deltaj ~ -1 + Hlambda, family = binomial(link = "logit"))
      obj$cens.mod.fitted[[j]] <- lambda
      Dhat <- Dhat.store[[j]] <- fitted(lambda)
      obj$obs.cens[[j]] <- deltaj
    }else if(sum(deltaj) == length(deltaj) & nocens == FALSE){
      Dhat <- Dhat.store[[j]] <- deltaj
      cat("Warning: no censoring in stage", j,", censoring model will be ignored.\n")
    }else{
      Dhat <- Dhat.store[[j]] <- deltaj
    }
    
    # estimation step
    # weights
    if (is.function(weight) == FALSE) {
      w <- w.store[[j]] <- abs(A - Ahat)/(deltaj*Dhat + (1 - deltaj)*(1 - Dhat))
    } else {
      # weight function specified for pi = 1 and delta = 1, therefore
      w <- w.store[[j]] <- weight(A, Ahat, deltaj, Dhat)
    }
    # blip parameter estimates extracted via dimensions of Hbeta and Hpsi
    logY <- log(Yj[deltaj == 1])
    wj <- w[deltaj == 1]
    Hbetaj <- Hbeta[deltaj == 1, ]
    Hpsij <- Hpsi[deltaj == 1, ]
    Aj <- A[deltaj == 1]
    est <- solve(t(cbind(Hbetaj, Aj * Hpsij)) %*% cbind(wj * Hbetaj, wj * Aj * Hpsij)) %*% t(cbind(Hbetaj, Aj * Hpsij)) %*% (wj * logY)
    psi <- est[(dim(Hbeta)[2] + 1):(dim(Hbeta)[2] + dim(Hpsi)[2])]
    beta <- est[1:dim(Hbeta)[2]]
    
    # if parametric bootstrap requested
    if(var.estim == "bootstrap" & boot.opt != "standard"){
      psi.boot <- matrix(NA, nrow = B, ncol = length(psi))
      residuals <- logY - Hbetaj %*% beta - Aj * Hpsij %*% psi
      reshist <- hist(residuals, plot = FALSE)
      mean.res <- mean(residuals)
      sd.res <- sd(residuals)
      boot.data <- dataetaj
      for (k in 1:B) {
        boot.data[,which(colnames(boot.data) == all.vars(cens.mod[[1]])[1])] <- rbinom(obj$n[[j]], 1, Dhat.store[[j]])
        if(boot.opt == "empirical"){
          bins <- with(reshist, sample(length(mids), sum(boot.data[,which(colnames(boot.data) == all.vars(cens.mod[[1]])[1])]), prob = density, replace = TRUE))
          newres <- runif(obj$n[[j]], reshist$breaks[bins], reshist$breaks[bins + 1])
        }else if(boot.opt == "normal"){
          newres <- rnorm(obj$n[[j]], mean = mean.res, sd = sd.res)
        }
        boot.data[, which(colnames(boot.data) == all.vars(time[[j]]))] <- exp(Hbeta %*% beta + A*Hpsi %*% psi + newres)
        psi.boot[k,] <- DWSurv(time = list(obj$time.mod[[j]]), blip.mod = list(obj$blip.mod[[j]]), treat.mod = list(obj$treat.mod[[j]]), tf.mod = list(obj$tf.mod[[j]]), cens.mod = list(obj$cens.mod[[j]]), boot.data, var.estim="none", quiet = TRUE)$psi[[1]]
      }
      obj$covmat[[j]] <- var(psi.boot)
      obj$psi.boot[[j]] <- psi.boot
    }
    
    # use estimates to identify optimal treatments
    # depending on optimization goal
    if (optimization == "max") {
      opt <- as.numeric(Hpsi %*% psi > 0)
    }else{
      opt <- as.numeric(Hpsi %*% psi < 0)
    }
    # ensure opt is in vector form
    opt <- as.vector(opt)
    # store estimates, optimal treatments and blip parameters
    obj$beta[[j]] <- beta
    obj$psi[[j]] <- psi
    obj$opt.treat[[j]] <- opt
    # store fitted values
    obj$log.fitted[[j]] <- Hbetaj %*% beta + Hpsij %*% psi
    obj$log.obs[[j]] <- logY
    # update pseudo-Y by adding regret for everybody who enter stage j
    obj$log.regret[[j]] <- (opt - A)*(Hpsi %*% psi)
    Yj <- exp(log(Yj) + obj$log.regret[[j]])
    old.Ycumul <- Ycumul
    Ycumul[eta[[j]] == 1] <- Yj
    
    # Asymptotic variance estimator done stage-by-stage
    if (var.estim == "asymptotic" & is.function(weight) == FALSE) {
      logY <- as.matrix(log(old.Ycumul[eta[[j]] == 1]))
      theta <- as.matrix(est)
      Htheta <- Htheta.store[[j]] <- as.matrix(cbind(Hbeta, A*Hpsi))
      # residuals are the (Y-gamma-beta) residuals
      residuals <- as.vector(logY - Htheta %*% theta)
      # HD: what precedes the (Y-gamma-beta) term in the est eq.
      HD <- HD.store[[j]] <- deltaj * w * Htheta
      # estimating equation = (Y - gamma - E[Y - gamma|H])*(A - Ahat)*Hpsi
      U <- residuals * HD
      U.adj <- U.adj.store[[j]] <- U
      Ahat <- as.vector(Ahat)
      Dhat <- as.vector(Dhat)
      A <- as.vector(A)
      nj <- sum(eta[[j]])
      # term that depends on alpha treatment model
      Halpha <- as.matrix(Halpha)
      score.alpha <- score.alpha.store[[j]] <- (A - Ahat) * Halpha
      ds.dalpha <- ds.dalpha.store[[j]] <- (-1/nj) * (crossprod(Ahat * (1 - Ahat) * Halpha, Halpha))
      dU.dalpha <- (1/nj) * (crossprod(deltaj * abs(A - Ahat)/(deltaj * Dhat + (1 - deltaj) * (1 - Dhat)) * (1 - 2*A) * Ahat * (1 - Ahat) * residuals * Htheta, Halpha))
      U.adj <- U.adj.store[[j]] <- U.adj - t(dU.dalpha %*% solve(ds.dalpha) %*% t(score.alpha))
      # terms that depends on lambda censoring model
      if(nocens == FALSE){
        Hlambda <- as.matrix(Hlambda)
        score.lambda <- score.lambda.store[[j]] <- (deltaj - Dhat) * Hlambda
        ds.dlambda <- ds.dlambda.store[[j]] <- (-1/nj) * (crossprod(Dhat * (1 - Dhat) * Hlambda, Hlambda))
        dU.dlambda <- (-1/nj) * (crossprod(deltaj * abs(A - Ahat)/(deltaj * Dhat + (1 - deltaj) * (1 - Dhat))^2 * (2*deltaj - 1) * Dhat * (1 - Dhat) * residuals * Htheta, Hlambda))
        U.adj <- U.adj.store[[j]] <- U.adj - t(dU.dlambda %*% solve(ds.dlambda) %*% t(score.lambda))
      }
      # add the terms that depend on previously estimated parameters
      if(j < K){
        for(l in (j+1):K){
          Hpsitemp <- matrix(0, nrow = N, ncol = ncol(Hpsi.store[[l]]))
          Hpsitemp[eta[[l]] == 1,] <- Hpsi.store[[l]]
          psitemp <- as.vector(obj$psi[[l]])
          px <- Hpsitemp %*% psitemp
          px <- px[eta[[j]] == 1]
          y <- rep(0, N)
          y[eta[[l]] == 1] <- Y[[l]][eta[[l]] == 1]
          y <- y[eta[[j]] == 1] # survival time in interval l - ensure coded with 0 if didnt enter stage l
          optl <- rep(0, N)
          optl[eta[[l]] == 1] <- obj$opt.treat[[l]] - A.store[[l]]
          optl <- optl[eta[[j]] == 1]
          etal <- eta[[l]][eta[[j]] == 1]
          B <- (Ycumul[eta[[j]] == 1])^(-1) * etal * y * exp(px * (optl)) * (optl)
          dUj.dpsil <- (1/nj)*crossprod(as.vector(deltaj * w * B) * Hpsitemp[eta[[j]] == 1,], Htheta)
          Hpsitemp <- Hpsitemp[eta[[l]] == 1, ]
          njp1 <- obj$n[[l]]
          dUjp1.dpsijp1 <- (-1/njp1) * (crossprod(HD.store[[l]], Hpsitemp))
          # treatment component
          dUjp1.dpsijp1 <- dUjp1.dpsijp1 - t(Htheta.store[[l]]) %*% Halpha.store[[l]] %*% solve(ds.dalpha.store[[l]]) %*% t(score.alpha.store[[l]]) %*% (-1/njp1 * delta.store[[l]] * w.store[[l]] * (1 - 2 * A.store[[l]]) * Ahat.store[[l]] * (1 - Ahat.store[[l]]) * Hpsitemp)
          # censoring component if necessary
          if(nocens == FALSE){
            dUjp1.dpsijp1 <- dUjp1.dpsijp1 - t(Htheta.store[[l]]) %*% Hlambda.store[[l]] %*% solve(ds.dlambda.store[[l]]) %*% t(score.lambda.store[[l]]) %*% (1/njp1 * delta.store[[l]] * abs(A.store[[l]] - Ahat.store[[l]])/(delta.store[[l]] * Dhat.store[[l]] + (1 - delta.store[[l]]) * (1 - Dhat.store[[l]]))^2 * (2 * delta.store[[l]] - 1) * Dhat.store[[l]] * (1 - Dhat.store[[l]]) * Hpsitemp)
          }
          dUjp1.dpsijp1 <- dUjp1.dpsijp1[(length(obj$beta[[l]]) + 1):(length(obj$beta[[l]]) + length(obj$psi[[l]])),]
          U.adje <- matrix(0, ncol = length(obj$psi[[l]]), nrow = N)
          U.adje[eta[[l]] == 1,] <- U.adj.store[[l]][,(length(obj$beta[[l]]) + 1):(length(obj$beta[[l]]) + length(obj$psi[[l]]))]
          U.adje <- U.adje[eta[[j]] == 1,]
          U.adj <- U.adj - t(t(dUj.dpsil) %*% solve(dUjp1.dpsijp1) %*% t(U.adje))
        }
      }
      # derivative with respect to parameter of interest
      dUadj.dtheta <- dU.dtheta <- (-1/nj) * (crossprod(HD, Htheta))
      if(asymp.opt == "adjusted"){
        dUadj.dtheta <- dUadj.dtheta  - t(Htheta) %*% Halpha %*% solve(ds.dalpha) %*% t(score.alpha) %*% (-1/nj * deltaj * w * (1 - 2 * A) * Ahat * (1 - Ahat) * Htheta)
        if(nocens == FALSE){
          dUadj.dtheta <- dUadj.dtheta - t(Htheta) %*% Hlambda %*% solve(ds.dlambda) %*% t(score.lambda) %*% (1/nj * deltaj * abs(A - Ahat)/(deltaj * Dhat + (1 - deltaj) * (1 - Dhat))^2 * (2 * deltaj - 1) * Dhat * (1 - Dhat) * Htheta)
        }
        invdUadj.dtheta <- solve(dUadj.dtheta)
        IF <- U.adj %*% invdUadj.dtheta
        sigma.theta <- (1/nj) * var(IF)
        obj$covmat.all[[j]] <- sigma.theta
        obj$covmat[[j]] <- obj$covmat.all[[j]][(length(obj$beta[[j]]) + 1):(length(obj$beta[[j]]) + length(obj$psi[[j]])), (length(obj$beta[[j]]) + 1):(length(obj$beta[[j]]) + length(obj$psi[[j]]))]
      }else{
        invdU.dtheta <- solve(dU.dtheta)
        IF <- U %*% invdU.dtheta
        sigma.theta <- (1/nj) * var(IF)
        obj$covmat.all[[j]] <- sigma.theta
        obj$covmat[[j]] <- obj$covmat.all[[j]][(length(obj$beta[[j]]) + 1):(length(obj$beta[[j]]) + length(obj$psi[[j]])), (length(obj$beta[[j]]) + 1):(length(obj$beta[[j]]) + length(obj$psi[[j]]))]
      }
    }
  }
  
  # optimal outcome relevant only for observed events
  obj$Y.opt <- Ycumul
  
  # standard errors
  # standard bootstrap
  if (var.estim == "bootstrap" & boot.opt == "standard") {
    psi.boot <- list()
    M <- nrow(data)
    for (i in 1:B) {
      data.boot <- data[sample(1:M, replace=TRUE),]
      psi.boot[[i]] <- DWSurv(time = obj$time.mod, blip.mod = obj$blip.mod, treat.mod = obj$treat.mod, tf.mod = obj$tf.mod, cens.mod = obj$cens.mod, data.boot, var.estim="none", quiet = TRUE)$psi
    }
    psi.boot <- do.call(function(...) mapply(rbind, ..., SIMPLIFY=FALSE), psi.boot)
    covmat <- lapply(psi.boot, var)
    obj$covmat <- covmat
    obj$psi.boot <- psi.boot
  }
  
  # if variance requested, calculate non-regularity
  if (var.estim != "none" & A.bin == 1) {
    for (j in K:1) {
      Hpsi <- model.matrix(blip.mod[[j]],data)
      # lower/upper limits on psi
      psi.l <- obj$psi[[j]] - 1.96*sqrt(diag(obj$covmat[[j]]))
      psi.u <- obj$psi[[j]] + 1.96*sqrt(diag(obj$covmat[[j]]))
      # look at max/min value of blip based on sign of covariates and max/min of parameter CIs
      Hpsi.p <- apply(Hpsi,c(1,2),FUN=function(x) max(x,0))
      Hpsi.n <- apply(Hpsi,c(1,2),FUN=function(x) min(x,0))
      # lower blip is positive values * lower CI + negative values * upper CI
      blip.l <- Hpsi.p %*% psi.l + Hpsi.n %*% psi.u
      blip.u <- Hpsi.p %*% psi.u + Hpsi.n %*% psi.l
      obj$nonreg[j] <- sum(blip.l < 0 & blip.u > 0)/dim(Hpsi)[1]
    }
  }
  # return
  class(obj) <- "DTRreg"
  obj
}

###########################################################
### Generalzied dynamic weighted ordinary least squares ###
###########################################################

############################################
### G-dWOLS function for multiple stages ###
############################################

# outcome: specifies the outcome variable
# Xpsi1: tailoring variables that interact with A in the blip function
# Xpsi2: tailoring variables that interact with A^2 in the blip funciton
# blip function is such that gamma = psi1*Xpsi1*A + psi2*Xpsi2*A^2
# treat.mod: specifies the treatment model (in terms of covariates)
# treat.fam: family objects for treatment model (with link function, as specified in the standard glm function)
#            Note: this argument is only relevant if a weight function based on the continuous doses is used
#            ex: gaussian(link = "identity")
#                Gamma(link = "log")
#                inverse.gaussian(link = "1/mu^2")
#                poisson(link = "log")
# tf.mod: specifies the treatment-free model
# weight.fcn: specifies the weight function to be used
# data: specifies the data
# m: specifies the number of groups for weight functions based on binned doses
# k: specifies the number of stages
# treat.range: range of permissible treatment values.


gdwols<-function(outcome, Xpsi1, Xpsi2, treat.mod, treat.fam=NULL, tf.mod, weight.fcn, data = NULL, m, k, treat.range=NULL){ 
  
  obj <- list()
  
  # if no data specified, create data frame
  if (is.null(data)) {
    # form list of variables
    model.list <- c(Xpsi1, Xpsi2, tf.mod, treat.mod)
    # extract variables: note we drop any intercept-only terms
    data <- cbind(outcome, do.call('cbind', lapply(model.list[which(model.list != '~1')], get_all_vars)))
    # remove duplicates
    data <- data[!duplicated(lapply(data,summary))]
    Y <- outcome
  } else {
    # if data is specified, pull outcome from it
    Y <- eval(substitute(outcome),envir=data)
  }
  
  # if no treat.fam specified, use Gaussian by default
  if(is.null(treat.fam)){treat.fam<-gaussian(link = "identity")}
  
  # if no treat.range specified, use min/max observed treatment (across all stages)
  if(is.null(treat.range)){
    treat.k<-NULL
    for(j in 1:k){treat.k<-c(treat.k,model.response(model.frame(treat.mod[[j]],data)))}
    treat.range<-c(min(treat.k),max(treat.k))
    cat("Warning: no treatment range specified, [min,max] of observed treatment used.\n")
  }
  
  # start at final stage and work backwards
  for (j in k:1) {
    
    # treatment 
    A <- model.response(model.frame(treat.mod[[j]],data))
    
    # different options for determining weights
    w.vec<-NULL
    if(weight.fcn=="ipw"){w.vec<-ipw.fcn(A,treat.mod[[j]],treat.fam,data,m)}
    if(weight.fcn=="cipw"){w.vec<-ipw.cap.fcn(A,treat.mod[[j]],treat.fam,data,m)}
    if(weight.fcn=="qpom"){w.vec<-q.pom.fcn(A,treat.mod[[j]],data,m)}
    if(weight.fcn=="wo"){w.vec<-new.bin.o.fcn(A,treat.mod[[j]],data,m)}
    
    # outcome model (Y) = treatment free model + blip model
    # tailoring variables in blip function
    Hpsi1 <- model.matrix(Xpsi1[[j]], model.frame(Xpsi1[[j]], data))
    Hpsi2 <- model.matrix(Xpsi2[[j]], model.frame(Xpsi2[[j]], data))
    psi.dim<-ncol(Hpsi1)+ncol(Hpsi2)
    # Hpsi <- model.matrix(blip.mod[[j]], model.frame(blip.mod[[j]], data))
    # variables in treatment free model
    Hbeta <- model.matrix(tf.mod[[j]], model.frame(tf.mod[[j]], data))	
    beta.dim <- ncol(Hbeta)
    # outcome model
    out <- lm(Y~0+Hbeta+A:Hpsi1+I(A^2):Hpsi2,weights=w.vec)
    # treatment-free parameter estimates 
    beta <- as.numeric(out$coef[1:beta.dim])
    # blip parameter estimates 
    psi <- as.numeric(out$coef[(beta.dim+1):(beta.dim+psi.dim)])
    
    # identify optimal treatments
    # note: for simplificity, we considered a quadratic blip function of the form 
    # PSI %*% (X1*A, X2*A^2), e.g. (psi1*X1)*A + (psi2*X2)*A^2, 
    # where X1 and X2 are (subsets) of the available covariates X
    # the optimal treatment is thus simply a.opt=-(psi1*X1)/2(psi2*X2) (see paper for more details)
    opt.step<-(-1/2)*(Hpsi1 %*% psi[1:ncol(Hpsi1)])/(Hpsi2 %*% psi[(ncol(Hpsi1)+1):psi.dim])
    # NOTE: a.opt is only a maximum if (psi2*X2)<0, otherwise endpoint (only if treatment range is finite)
    if(abs(min(treat.range)) != Inf & abs(max(treat.range)) != Inf){
      Hpsi.min <- cbind((Hpsi1*min(treat.range)),(Hpsi2*(min(treat.range)^2)))
      Hpsi.max <- cbind((Hpsi1*max(treat.range)),(Hpsi2*(max(treat.range)^2)))
      endpoints.diff<-(Hpsi.min %*% psi)-(Hpsi.max %*% psi)
      opt <- (as.numeric(Hpsi2 %*% psi[(ncol(Hpsi1)+1):psi.dim] < 0)*opt.step)+
        (as.numeric(Hpsi2 %*% psi[(ncol(Hpsi1)+1):psi.dim] >= 0)*as.numeric(endpoints.diff>=0)*min(treat.range))+
        (as.numeric(Hpsi2 %*% psi[(ncol(Hpsi1)+1):psi.dim] >= 0)*as.numeric(endpoints.diff<0)*max(treat.range))
    } else {
      opt <- opt.step
      warn.vec<-which(Hpsi2 %*% psi[(ncol(Hpsi1)+1):psi.dim] > 0)
      if(length(warn.vec)>0 ){
        cat("Warning: stage",j,"optimal treatment may be a minimum for some observations. \n")
      }
      
    }
    
    # optimal dose must fall withing treatment range
    opt <- pmin(pmax(opt,min(treat.range)),max(treat.range))
    opt<-as.vector(opt)
    
    Hpsi.obs <- cbind(sweep(Hpsi1,1,A,'*'),sweep(Hpsi2,1,(A^2),'*'))
    Hpsi.opt <- cbind(sweep(Hpsi1,1,opt,'*'),sweep(Hpsi2,1,(opt^2),'*'))
    # update (pseudo-) Y by adding regret
    Y <- Y + (Hpsi.opt %*% psi) - (Hpsi.obs %*% psi)
    
    # store estimates
    obj$psi[[j]] <- psi
    obj$beta[[j]] <- beta
    obj$opt.treat[[j]] <- opt
  }
  
  return(obj)
  
}


#########################################################
### Different weight functions for two stage function ###
#########################################################

# different weight functions

# IPW 
ipw.fcn<-function(A,treat.mod,treat.fam,data,m){
  
  Xalpha<-model.matrix(treat.mod,model.frame(treat.mod,data))
  model.output<-glm(A~0+Xalpha,data=data,family=treat.fam)
  model.disp<-summary(model.output)$dispersion
  # model.output<-lm(A~0+Xalpha,data=data)
  
  if(treat.fam$family=="gaussian"){
    weight.vec<-1/dnorm(A,mean=model.output$fitted.values,sd=sqrt(model.disp))
  }
  
  if(treat.fam$family=="Gamma"){
    weight.vec<-1/dgamma(A,shape=1/model.disp,scale=model.output$fitted.values*model.disp)
  }
  
  return(weight.vec)
  
}

# capped IPW 
ipw.cap.fcn<-function(A,treat.mod,treat.fam,data,m){
  
  Xalpha<-model.matrix(treat.mod,model.frame(treat.mod,data))
  model.output<-glm(A~0+Xalpha,data=data,family=treat.fam)
  model.disp<-summary(model.output)$dispersion
  # model.output<-lm(A~0+Xalpha,data=data)
  
  if(treat.fam$family=="gaussian"){
    weights<-1/dnorm(A,mean=model.output$fitted.values,sd=sqrt(model.disp))
  }
  
  if(treat.fam$family=="Gamma"){
    weights<-1/dgamma(A,shape=1/model.disp,scale=model.output$fitted.values*model.disp)
  }
  
  cap<-quantile(weights,0.99)
  weight.vec<-pmin(weights,cap)
  
  return(weight.vec)
  
}


# quantile binning based on POM
q.pom.fcn<-function(A,treat.mod,data,m){
  
  Xalpha<-model.matrix(treat.mod,model.frame(treat.mod,data))
  
  a.binned<-ntile(A, m)
  a.cat<-as.factor(a.binned)
  
  if(ncol(Xalpha)==1){
    pom<-polr(a.cat~1)
  }
  else{
    pom<-polr(a.cat~Xalpha[,-1])
  }
  
  a.pom.prob<-pom$fitted.values[cbind(seq_along(a.binned), a.binned)]
  
  weight.vec<-(1/m)/a.pom.prob
  
  return(weight.vec)
  
}


# overlap weights
new.bin.o.fcn<-function(A,treat.mod,data,m){
  
  Xalpha<-model.matrix(treat.mod,model.frame(treat.mod,data))
  
  a.binned<-ntile(A,m)
  a.cat<-as.factor(a.binned)
  
  if(ncol(Xalpha)==1){
    pom<-polr(a.cat~1)
  }
  else{
    pom<-polr(a.cat~Xalpha[,-1])
  }
  
  a.pom.prob<-pom$fitted.values[cbind(seq_along(a.binned), a.binned)]
  
  weight.vec<-(1/a.pom.prob)/apply(1/pom$fitted.values,1,sum)
  
  return(weight.vec)
  
}

