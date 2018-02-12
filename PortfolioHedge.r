rm(list=ls());
graphics.off();
require(RMySQL);
require(fGarch);
require(xts);
source("~/kkasi/r/libxxie.r");
source("~/cake/libeix.r");

## How much of spy to long/short
## @sample: returns of [portfolio, factor]
## @lookback: 1:lookback is the section of the sample to be used for model fitting
## @money: monetary amount on the portfolio. length less than length(sample) - lookback
## @confidence level of expected shortfall
## @loss.tol The negative return that we are willing to tolerate
ptfl.hedge <- function(sample, lookback, money, confidence, loss.tol)
{
    ## Learn from the history
    models <- list();
    for (i in 1:2) {
        models[[i]] <- fit.garch(sample[1:lookback, i], max.order=c(2,2));
    }
    C <- cov(cbind(models[[1]]$inno, models[[2]]$inno));

    PnL <- function(mu, cov.mtx, alpha, Beta) {
        mean.pnl <- sum(mu * Beta);
        sd.pnl <- sqrt(sum(t(Beta) %*% cov.mtx %*% Beta));
        ## integral <- integrate(f=function(x) x * dnorm(x, mean=mean.pnl, sd=sqrt(variance.pnl)),
        ##                       lower=0, upper=alpha);
        ## upper <- qnorm(p=alpha, mean=mean.pnl, sd=sqrt(variance.pnl));
        integral <- integrate(
            f=function(x) {
                qnorm(p=x, mean=mean.pnl, sd=sd.pnl)
            },
            lower=0, upper=alpha,
            rel.tol=1.0e-4
        );
        pnl <- integral$value/alpha/sum(abs(Beta));
        return(pnl);
    }
    ## compute hedge for each day.
    hedge <- matrix(NA, nrow=length(money), ncol=2);
    for (t in 1:length(money)) {
        sig <- matrix(0, nrow=2, ncol=2);
        mu <- rep(NA, 2);
        for (i in 1:2) {
            sig[i, i] <- sum(models[[i]]$alpha * models[[i]]$r.t^2 + models[[i]]$beta * models[[i]]$sigma.t^2) +
                models[[i]]$omega;
            sig[i, i] <- sqrt(sig[i, i]);
            mu[i] <- models[[i]]$mu;
            ## Update the model parameters
            a <- length(models[[i]]$alpha);
            b <- length(models[[i]]$beta);
            models[[i]]$sigma.t <- c(tail(models[[i]]$sigma.t, n=b-1), sig[i, i]);
            models[[i]]$r.t <- c(tail(models[[i]]$r.t, n=a-1), sample[lookback + t, i]);
        }
        cov.mtx <- sig %*% C %*% t(sig);
        ## we shall solve f(x) - loss.tol = 0
        expt.stfl <- function(x) {
            Beta <- c(1, x);
            y <- PnL(mu, cov.mtx, confidence, Beta);
            return(-y);
        }
        hedge[t, 1] <- 0;
        hedge[t, 2] <- loss.tol;
        if (expt.stfl(0) > loss.tol) {
            ## We are taking more risk than we can bear
            hedge.min <- optimize(f=expt.stfl, lower=-3, upper=0);
            if (hedge.min$objective < loss.tol) {
                solution <- uniroot(f=function(x) expt.stfl(x) - loss.tol, interval=c(hedge.min$minimum, 0));
                hedge[t, 1] <- solution$root;
                hedge[t, 2] <- loss.tol;
            } else {
                hedge[t, 1] <- hedge.min$minimum;
                hedge[t, 2] <- hedge.min$objective;
            }
            hedge[t, 1] <- hedge[t, 1] * money[t];
        }
    }
    return(hedge);
}

##
## X: a matrix with the following columns, from 1st to the last.
## X$atv_stk: Returns of the portfolio
## X$spy: returns of SPY
## X$net_expo: net exposure of the portfolio
## X$closing: Prices of SPY at market close
## X$value: Worth of the portfolio at market close.
## dd.entry: drawdown level to start hedging
## dd.exit: drawdown level to exit hedging
## cal.hedge <- function(X, lookback, dd.entry, dd.exit)
cal.hedge <- function(X, lookback, confidence=0.05, loss.tol=0.05)
{
    adj.min <- 0.005;
    n <- dim(X)[1];
    hedge <- rep(0, n);
    ES <- rep(NA, 0);
    CdDD <- cummax(X$value) - X$value;
    ## Conditional profit and loss
    PnL <- function(mu, cov.mtx, alpha, Beta) {
        mean.pnl <- sum(mu * Beta);
        sd.pnl <- sqrt(sum(t(Beta) %*% cov.mtx %*% Beta));
        ## integral <- integrate(f=function(x) x * dnorm(x, mean=mean.pnl, sd=sqrt(variance.pnl)),
        ##                       lower=0, upper=alpha);
        ## upper <- qnorm(p=alpha, mean=mean.pnl, sd=sqrt(variance.pnl));
        integral <- integrate(
            f=function(x) {
                qnorm(p=x, mean=mean.pnl, sd=sd.pnl)
            },
            lower=0, upper=alpha,
            rel.tol=1.0e-4
        );
        pnl <- integral$value/alpha/sum(abs(Beta));
        return(pnl);
    }
    
    for (t in (lookback):n) {
        ## Compute the entry point to start hedging
        models <- list();
        for (i in 1:2) {
            models[[i]] <- fit.garch(X[(t-lookback+1):t, i], max.order=c(2,2));
        }
        C <- cov(cbind(models[[1]]$inno, models[[2]]$inno));
        sig <- matrix(0, nrow=2, ncol=2);
        mu <- rep(NA, 2);
        for (i in 1:2) {
            sig[i, i] <- sum(models[[i]]$alpha * models[[i]]$r.t^2) + sum(models[[i]]$beta * models[[i]]$sigma.t^2) +
                models[[i]]$omega;
            sig[i, i] <- sqrt(sig[i, i]);
            mu[i] <- models[[i]]$mu;
            ## Update the model parameters
            a <- length(models[[i]]$alpha);
            b <- length(models[[i]]$beta);
            models[[i]]$sigma.t <- c(tail(models[[i]]$sigma.t, n=b-1), sig[i, i]);
            models[[i]]$r.t <- c(tail(models[[i]]$r.t, n=a-1), X[lookback + t, i]);
        }
        cov.mtx <- sig %*% C %*% t(sig);
        ## we shall solve f(x) - loss.tol = 0
        expt.stfl <- function(x) {
            Beta <- c(1, x);
            y <- PnL(mu, cov.mtx, confidence, Beta);
            return(-y);
        }
        ES[t] <- expt.stfl(0);
        if (ES[t] > loss.tol) {
            ## We are taking more risk than we can bear
            model <- lm(X$atv_stk~X$spy, subset=(t-lookback+1):t);
            hedge[t] <- -X$net_expo[t] * coef(model)[2];
            ## hedge[t] <- max(hedge[t], -X$net_expo[t]);
            if (ES[t] > loss.tol && abs(hedge[t]/hedge[t-1] - 1) < adj.min) {
                hedge[t] <- hedge[t-1];
            } else if (ES[t] <= loss.tol) {
                hedge[t] <- 0;
            }
        }
    }
    
    V <- matrix(0, nrow=n, ncol=2);
    V[, 1] <- hedge;
    for (t in 2:n) {
        V[t, 2] <- V[t-1, 2] + V[t-1, 1]/X$closing[t-1] * X$closing[t] - V[t, 1];
    }
    H <- apply(V, MARGIN=1, FUN=sum);
    Y <- H + X$value;
    Z <- cummax(Y);
    DD <- (Z - Y)/Z;
    return(list(hit=H, DD=DD, hedge=hedge, ES=ES));
}

cal.hedge.4 <- function(X, lookback=40, confidence)
{
    n <- dim(X)[1];
    ## Number of shares to hedge.
    hedge <- rep(0, n);
    ## Number of shares held
    holding <- rep(0, n);
    for (t in (lookback):n) {
        holding[t] <- holding[t-1];
        result <- t.test(X$spy[(t-lookback+1):t]);
        if (result$p.value < confidence && result$estimate < 0) {
            ## The SPY drift is significantly negative.
            model <- lm(X$atv_stk~X$spy, subset=(t-lookback+1):t);
            holding[t] <- -X$net_expo[t] * coef(model)[2]/X$closing[t];
        } else if (result$p.value < confidence && result$estimate > 0) {
            ## The SPY drift is significantly positive.
            holding[t] <- 0;
        }
    }
    
    cash <- rep(0, n);
    for (t in 2:n) {
        cash[t] <- cash[t-1] - (holding[t] - holding[t-1])* X$closing[t];
    }
    Y <- cash + holding * X$closing + X$value;
    Z <- cummax(Y);
    DD <- (Z - Y)/Z;
    return(list(hit=1 - Y/X$value, DD=DD, hedge=hedge, Y=Y));
}


##
## X: a matrix with the following columns, from 1st to the last.
## X$atv_stk: Returns of the portfolio
## X$spy: returns of SPY
## X$net_expo: net exposure of the portfolio
## X$closing: Prices of SPY at market close
## X$value: Worth of the portfolio at market close.
cal.hedge.3 <- function(X, lookback=40, threshold=0.2)
{
    n <- dim(X)[1];
    ## Number of shares to hedge.
    hedge <- rep(0, n);
    ## Number of shares held
    holding <- rep(0, n);
    sharpe <- rep(NA, n);
    ## Number of trading days in a year
    N <- 252;
    CdDD <- cummax(X$value) - X$value;
    ## for (t in 161:n) {
    for (t in (lookback):n) {
        holding[t] <- holding[t-1];
        ## Compute the entry point to start hedging
        dist.ret <- fit.dist(X$spy[(t-lookback+1):t]);
        if (dist.ret$bic == Inf) next;
        sharpe[t] <- dist.ret$mu/dist.ret$sig;
        if (sharpe[t] < -threshold) {
            ## The SPY drift is significantly negative.
            model <- lm(X$atv_stk~X$spy, subset=(t-lookback+1):t);
            holding[t] <- -X$net_expo[t] * coef(model)[2]/X$closing[t];
        } else if (sharpe[t] > threshold) {
            holding[t] <- 0;
        }
    }
    
    cash <- rep(0, n);
    for (t in 2:n) {
        cash[t] <- cash[t-1] - (holding[t] - holding[t-1])* X$closing[t];
    }
    Y <- cash + holding * X$closing + X$value;
    Z <- cummax(Y);
    DD <- (Z - Y)/Z;
    return(list(hit=1 - Y/X$value, DD=DD, hedge=hedge, Y=Y));
}

##
## X: a matrix with the following columns, from 1st to the last.
## X$atv_stk: Returns of the portfolio
## X$spy: returns of SPY
## X$net_expo: net exposure of the portfolio
## X$closing: Prices of SPY at market close
## X$value: Worth of the portfolio at market close.
cal.hedge.2 <- function(X, lookback=40, threshold=0.1)
{
    n <- dim(X)[1];
    ## Number of shares to hedge.
    hedge <- rep(0, n);
    ## Number of shares held
    holding <- rep(0, n);
    drift <- rep(NA, n);
    ## Number of trading days in a year
    N <- 252;
    CdDD <- cummax(X$value) - X$value;
    ## for (t in 161:n) {
    for (t in (lookback):n) {
        holding[t] <- holding[t-1];
        ## Compute the entry point to start hedging
        dist.ret <- fit.dist(X$spy[(t-lookback+1):t]);
        if (dist.ret$bic == Inf) next;
        drift[t] <- dist.ret$mu;
        if (dist.ret$mu * N < -threshold) {
            ## The SPY drift is significantly negative.
            model <- lm(X$atv_stk~X$spy, subset=(t-lookback+1):t);
            holding[t] <- -X$net_expo[t] * coef(model)[2]/X$closing[t];
        } else if (dist.ret$mu * N > threshold) {
            holding[t] <- 0;
            ## drift positive.
            next;
        } else {
            ## fit the price distribution
            ## stop("unclear trend");
            ## warning(sprintf("unclear trend at %d", t));
            ## model <- lm(X$atv_stk~X$spy, subset=(t-lookback+1):t);
            ## dist.price <- fit.dist(X$closing[(t-lookback+1):t]);
            ## if (dist.price$bic == Inf) next;
            ## if (dist.price$dist == "norm") {
            ##     ub <- qnorm(p=0.95, mean=dist.price$mu, sd=dist.price$sig);
            ##     lb <- qnorm(p=0.05, mean=dist.price$mu, sd=dist.price$sig);
            ## } else if (dist.price$dist == "t") {
            ##     ub <- qt(p=0.95, df=dist.price$df) + dist.price$mu;
            ##     lb <- qt(p=0.05, df=dist.price$df) + dist.price$mu;
            ## } else {
            ##     stop("Unknown distribution");
            ## }
            ## if (X$closing[t] > ub) {
            ##     holding[t] <- -X$net_expo[t] * coef(model)[2]/X$closing[t];
            ## } else if (X$closing[t] > lb && X$closing[t] < ub) {
            ##     holding[t] <- 0;
            ## }
            next;
        }
    }
    
    cash <- rep(0, n);
    for (t in 2:n) {
        cash[t] <- cash[t-1] - (holding[t] - holding[t-1])* X$closing[t];
    }
    Y <- cash + holding * X$closing + X$value;
    Z <- cummax(Y);
    DD <- (Z - Y)/Z;
    return(list(hit=1 - Y/X$value, DD=DD, hedge=hedge, Y=Y));
}


cal.hedge.1 <- function(X, lookback, dd.entry, dd.exit)
{
    adj.min <- 0.005;
    n <- dim(X)[1];
    hedge <- rep(0, n);
    CdDD <- 1 - X$value/cummax(X$value);
    for (t in (lookback + 1):n) {
        if (hedge[t-1] == 0 && CdDD[t] <= dd.entry) {
            next;
        }
        model <- lm(X$atv_stk~X$spy, subset=(t-lookback+1):t);
        hedge[t] <- -X$net_expo[t] * coef(model)[2];
        hedge[t] <- max(hedge[t], -X$net_expo[t]);
        if (CdDD[t] <= dd.entry && CdDD[t] > 0.01 && abs(hedge[t]/hedge[t-1] - 1) < adj.min) {
            hedge[t] <- hedge[t-1];
        } else if (CdDD[t] <= dd.exit) {
            hedge[t] <- 0;
        }
    }

    V <- matrix(0, nrow=n, ncol=2);
    V[, 1] <- hedge;
    for (t in 2:n) {
        V[t, 2] <- V[t-1, 2] + V[t-1, 1]/X$closing[t-1] * X$closing[t] - V[t, 1];
    }
    H <- apply(V, MARGIN=1, FUN=sum);
    Y <- H + X$value;
    Z <- cummax(Y);
    DD <- (Z - Y)/Z;
    return(list(Y=Y, hit=-H/X$value, DD=DD, hedge=hedge));
}

database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
                     dbname='market', host="localhost");

## stmt <- paste(
##     "select S2.tm, S2.closing,",
##     "atenvia.net_expo, atenvia.value, atenvia.IntDD,",
##     "atv_stk, (S2.closing - S1.closing)/S1.closing as spy",
##     "from atv_returns join spy_daily as S1 join spy_daily as S2",
##     "join atenvia",
##     "on atv_returns.tm = S2.tm",
##     "and atenvia.tm = S2.tm",
##     "and S2.tm = (select tm from spy_daily where tm > S1.tm limit 1);"
## );
stmt <- paste(
    "select S2.tm, S2.closing,",
    "atenvia.net_expo, atenvia.value, atenvia.IntDD,",
    "atv_stk, (S2.closing - S1.closing)/S1.closing as spy",
    "from atv_returns join spy_daily as S1 join spy_daily as S2",
    "join atenvia",
    "on atv_returns.tm = S2.tm",
    "and atenvia.tm = S2.tm",
    "and S2.tm = (select tm from spy_daily where tm > S1.tm limit 1);"
);
results <- dbSendQuery(database, stmt);
data <- fetch(results, n=-1);
dbClearResult(results);
dbDisconnect(database);

X <- data[, c(6, 7, 3, 2, 4)];
n <- dim(X)[1];
lookback <- 40;
## peaks <- unlist(lapply(1:n, function(i) max(X$value[1:i])))
## draw.down <- unlist(lapply(1:n, function(i) max(X$value[1:i])))/X$value - 1;

result1 <- cal.hedge.1(X, lookback=40, dd.entry=0.05, dd.exit=0.02);
result2 <- cal.hedge.2(X, lookback=40, threshold=0.1);
result3 <- cal.hedge.3(X, lookback=40, threshold=0.3);
result4 <- cal.hedge.4(X, lookback=40, confidence=0.1);
summary <- matrix(c(
    max(result1$DD), mean(result1$DD), max(result1$hit), mean(result1$hit),
    max(result2$DD), mean(result2$DD), max(result2$hit), mean(result2$hit),
    max(result3$DD), mean(result3$DD), max(result3$hit), mean(result3$hit),
    max(result4$DD), mean(result4$DD), max(result4$hit), mean(result4$hit)),
    byrow=TRUE, nrow=4, ncol=4);
summary

pdf('/tmp/hedge.pdf');
graphics.off();
plot(xts(X$value, order.by=as.Date(data$tm)),
     main="Portfolio & Hedge",
     type="l", col="#000000");
## lines(xts(Y, order.by=as.Date(data$tm)), col="#FF0000");
lines(xts(result1$Y, order.by=as.Date(data$tm)), col="#FF0000");
lines(xts(result2$Y, order.by=as.Date(data$tm)), col="#0000FF");
lines(xts(result3$Y, order.by=as.Date(data$tm)), col="#00FF00");
lines(xts(result4$Y, order.by=as.Date(data$tm)), col="#6600CC");
dev.off();




