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
            sig[i, i] <- sum(models[[i]]$alpha * models[[i]]$r.t^2 + models[[i]]$beta * models[[i]]$sigma.t^2 ) + models[[i]]$omega;
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
cal.hedge <- function(X, lookback, IntDD, dd.entry, dd.exit)
{
    adj.min <- 0.005;
    n <- dim(X)[1];
    hedge <- rep(0, n);
    CdDD <- cumamx(X$value) - X$value;
    for (t in (lookback + 1):n) {
        ## Compute the entry point to start hedging
        ## garch.para <- vector("list", length=2);
        ## for (i in 1:2) {
        ##     garch.para[[i]] <- fit.garch(X[(t-lookback+1):t, i]);
        ## }
        ## if (garch.para[[1]].score == Inf || garch.para[[2]] == Inf) {
        ##     dd.entry 
        ## }
        if (hedge[t-1] == 0 && CdDD[t] <= dd.entry) {
            next;
        }
        model <- lm(X$atv_stk~X$spy, subset=(t-lookback+1):t);
        hedge[t] <- -X$net_expo[t] * coef(model)[2];
        hedge[t] <- max(hedge[t], -X$net_expo[t]);
        if (CdDD[t] <= dd.entry && CdDD[t] > dd.exit && abs(hedge[t]/hedge[t-1] - 1) < adj.min) {
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
    return(list(hit=H, DD=DD, hedge=hedge));
}


database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
                     dbname='market', host="localhost");

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
X <- fetch(results, n=-1);
dbClearResult(results);
dbDisconnect(database);

n <- dim(X)[1];
lookback <- 120;
## peaks <- unlist(lapply(1:n, function(i) max(X$value[1:i])))
## draw.down <- unlist(lapply(1:n, function(i) max(X$value[1:i])))/X$value - 1;
IntDD <- X$IntDD;

dd.entry <- seq(from=0.015, to=0.027, by=0.001);
dd.exit <- seq(from=0.005, to=0.017, by=0.001);
result <- array(NA, dim=c(length(dd.entry), length(dd.exit), 3));
for (i in 1:length(dd.entry)) {
    for (j in 1:length(dd.exit)) {
        hd <- cal.hedge(X, lookback, IntDD, dd.entry[i], dd.exit[j]);
        result[i, j, 1] <- max(hd$DD) * 0.7 - 0.3 * min(hd$hit);
        result[i, j, 2] <- max(hd$DD);
        result[i, j, 3] <- mean(hd$hit);
    }
}
idx <- which(result[,,1] == min(result[,,1]), arr.ind=TRUE);




plot(xts(DD, order.by=as.Date(X$tm)), type="l", col="#000000");

graphics.off();
##par(mfrow=c(2, 1));
plot(xts(X$value, order.by=as.Date(X$tm)), type="l", col="#000000");
lines(xts(Y, order.by=as.Date(X$tm)), col="#00FF00");

## plot(1:n, X$value + hedge[, 1], col="#0000FF", type="l");
## lines(1:n, X$value);

lookback <- 120;
refit <- 20;
## [$ amount of SPY, expected shortfall]
hedge <- matrix(0, nrow=n, ncol=2);
for (t in (lookback + 1):(n-refit+1)) {
    if ((t-1) %% refit == 0) {
        I <- (t - lookback):(t + refit - 1);
        hedge[t:(t+refit - 1), ] <- ptfl.hedge(cbind(X$atv_stk[I], X$spy[I]), lookback, X$net_expo[t:(t+refit - 1)], 0.05, 0.05);
    }
}
plot(xts(X$value + abs(hedge[, 1]), order.by=as.Date(X$tm)), col="#000000", lwd=2, ylim=c(0, max(X$value + abs(hedge[, 1]))));
lines(xts(x=X$value, order.by=as.Date(X$tm)), col="#00FF00", lwd=1);
abline(h=0, lwd=2, col="#FF0000");


