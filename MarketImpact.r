rm(list=ls());
require("RMySQL")
require("MASS")
require("fGarch")
require("xts")
require("stats4")
source("~/kkasi/r/libxxie.r");

database = dbConnect(MySQL(), user='sinbaski', password='q1w2e3r4',
                     dbname='market', host="localhost");

res <- dbSendQuery(database,
                   paste(
                       "select closing from spy_daily;"
                   ));
P <- fetch(res, n=-1);
dbClearResult(res);
X <- diff(log(P$closing));

garch11 <- garchFit(formula = ~garch(1,1), data=X,
                    cond.dist="std",
                    include.mean=FALSE,
                    trace=FALSE
         );
nu <- coef(garch11)[4];
    
res <- dbSendQuery(database,
                   paste(
                       "select distinct(date(tm)) as D from spy_cleaned;"
                   ));
X <- fetch(res, n=-1);
dbClearResult(res);
tm <- X$D;

res <- dbSendQuery(database,
                   paste(
                       "select count(*) as N from spy_cleaned;"
                   ));
X <- fetch(res, n=-1);
dbClearResult(res);
N <- X$N;

Y <- rep(NA, N-length(tm));
k <- 1;
for (d in 1:length(tm)) {
    res <- dbSendQuery(database, 
                       paste(
                           "select sum(V) as V from spy_cleaned where date(tm) = \"", tm[d], "\";"
                       ));
    X <- fetch(res, n=-1);
    dbClearResult(res);
    V <- X$V;
    
    res <- dbSendQuery(database, 
                       paste(
                           "select tm, P, V from spy_cleaned where date(tm) = \"", tm[d], "\";"
                       ));
    X <- fetch(res, n=-1);
    dbClearResult(res);
    ret <- diff(log(X$P));
    ratio <- X$V[-1]/V;
    sigma <- sqrt(sum(ret^2));

    Y[k:(k+length(ret)-1)] <- ret/sqrt(ratio)/sigma;
    k <- k + length(ret);
}

Z <- Y;
loglik.Y <- function(xi, p) {
    sum(log((dt(Z - xi, df=nu) + dt(Z + xi, df=nu)) * p + (1 - 2*p) * dt(Z, df=nu)))
}

## fit <- mle(function(xi, p) -loglik.Y(xi, p), start=list(xi=0.5, p=0.4), method="SANN");

## NANs produced
## fit <- mle(function(xi, p) -loglik.Y(xi, p), start=list(xi=0.5, p=0.4), method="Nelder-Mead");

## NANs produced
## fit <- mle(function(xi, p) -loglik.Y(xi, p), start=list(xi=0.5, p=0.4), method="BFGS");

## There were 14 warnings (use warnings() to see them)
fit <- mle(function(xi, p) -loglik.Y(xi, p), start=list(xi=0.5, p=0.4), upper=c(1.52, 0.5), lower=c(1.0e-6, 0.01), method="L-BFGS-B");

days <- c(
    "2017-10-05",
    "2017-09-27",
    "2017-10-04",
    "2017-09-26",
    "2017-10-03",
    "2017-09-27",
    "2017-09-27",
    "2017-09-27",
    "2017-10-02",
    "2017-09-29",
    "2017-10-06",
    "2017-10-16",
    "2017-09-01",
    "2017-09-22",
    "2017-09-18",
    "2017-10-03",
    "2017-09-28",
    "2017-10-11",
    "2017-10-03"
);

volume <- c(
    1586300,
    1203180,
    977844,
    963352,
    900809,
    535878,
    514387,
    493656,
    474800,
    406400,
    402699,
    400000,
    390427,
    387672,
    376029,
    357866,
    328921,
    312000,
    300319
);

prices <- c(
    253.899198,
    250.150000,
    252.750000,
    249.070000,
    252.480000,
    250.180803,
    250.320178,
    249.980000,
    252.120000,
    251.060000,
    254.092173,
    255.110000,
    248.130000,
    249.511926,
    249.951000,
    252.690000,
    250.260714,
    254.605000,
    252.390000
);

slippage <- rep(NA, length(days));

for (d in 1:length(days)) {
    res <- dbSendQuery(database,
                       paste(
                           "select V from daily_stat where D=\"", days[d], "\""
                       ));
    X <- fetch(res, n=-1);
    V <- X$V;
    dbClearResult(res);    

    res <- dbSendQuery(database,
                       paste(
                           "select P from spy_cleaned where date(tm)=\"", days[d], "\""
                       ));
    X <- fetch(res, n=-1);
    ret <- diff(log(X$P));
    dbClearResult(res);
    
    sigma <- sqrt(sum(ret^2));
    slippage[d] <- prices[d] * sqrt(volume[d]/V) * sigma * (qt(p=0.6, df=nu) + (if (coef(fit)[1] > 1.0e-6) coef(fit)[1] else 0));
    ## slippage <- prices[d] * sqrt(volume[d]/V) * sigma * (qt(p=seq(from=0.1, to=0.9, by=0.1), df=nu) + (if (coef(fit)[1] > 1.0e-6) coef(fit)[1] else 0));
}
dbDisconnect(database);

## garch21 <- garchFit(formula = ~garch(2,1), data=X$ret,
##                     cond.dist="std",
##                     include.mean=FALSE,
##                     trace=FALSE
##          );

## garch12 <- garchFit(formula = ~garch(1,2), data=X$ret,
##                     cond.dist="std",
##                     include.mean=FALSE,
##                     trace=FALSE
##          );

## model1 <- arima(log(X$volume), order=c(1,1,1));
## fit1 <- fitdistr(model1$residuals, "t",
##                  control=list(reltol=1.0e-2),
##                  ## method="SANN",
##                  ## method="CG",
##                  ## method="BFGS",
##                  start=list(df=4, m=0, s=sqrt(2)/sd(model1$residuals)));

## q <- quantile(X$ret, probs=seq(from=0.95, to=0.99, by=0.005));
## q <- quantile(X$ret, probs=seq(from=0.01, to=0.05, by=0.005));
## E <- rep(NA, length(q));
## for (i in 1:length(q)) {
##     I <- X$ret < q[i];
##     E[i] <- mean(X$ret[I]/garch11@sigma.t[I]);
## }
## avg <- 316.84079290044497;
## plot(q, E);
## c.m <- max(abs(E))/avg;



## x <- 1.0e-3;
## delta <- 0.01
## nu <- coef(garch11)[4];
## Q <- x^2 * X$volume / (coef(garch11)[1] + qt(1 - delta, df=nu))^2 / garch11@h.t;

## size <- 1.0e+4;
## alpha <- nu/2;
## c.m <- gamma(alpha+1/2)/gamma(alpha)*sqrt(alpha);
## q <- garch11@sigma.t * sqrt(size/X$volume);
## price <- 100;
## ## q <- q * (c.m + qt(p=seq(from=0.1, to=0.9, by=0.1), df=nu));

## slippage <- q * price * 100;
## den <- density(slippage);
## plot(den$x, den$y);
