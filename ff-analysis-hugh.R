require("DoE.base") # for half-normal plot
require("stringr")
source("interaction.plot.R")
# IMPORTANT TO ADD THIS, so that aliasing will work.
options(contrasts=c("contr.sum","contr.poly"))

# read in results for different replicates and consolidate.  
replicates <- 1:2
tmp <- NULL
for (i in replicates){
  load(paste("replicate_",i,"_results.RData",sep=""))
  tmp <- rbind(tmp,results)
}
results <- tmp
for(i in 2:7) results[,i] <- as.factor(results[,i])
rm(tmp, replicates)

# Note that if you would like to re-run the experiment, use run_cluster.R
# It takes something like 20 minutes per replicate.

# Uncomment next line if you want A ... G as factor names
#colnames(results)[2:8]  <- LETTERS[1:7]


results$response <- with(results, 1 - model.mse.test / meanSStot)
results$response <- with(results, log(response/(1-response)))

# FF design will be a crossed array of "model" by the other 6 (noise) factors.
# n=A q=B ENE=C beta.mu=D sigma=E x.cor=F model=G
# Design in A - F is a resolution 4 quarter fraction, 2^{6-2}_{IV} with generators
# ABCE = +1 (or -1)
# BCDF = +1 (or -1)

# convert data frame with actual levels to design notation with +/- levels
pm1 <- matrix(0,nrow(results),7, dimnames = list(NULL,LETTERS[1:7]))
for (i in 2:8){ 
  pm1[,i-1] <- 2*(as.numeric(results[,i])-1.5)
}

# Make the quarter fraction:
pm1 <- as.data.frame(pm1)
fraction <- with(pm1,(A*B*C*E==1) & (B*C*D*F==1))
table(fraction) # check - should have 64 out of 256 TRUE
table(with(pm1[fraction,],A*D*E*F)) # check - should be all 1
results2 <- results[fraction,]

# easier to look at aliaing if we use labels A - G
tmp <- pm1[fraction,]
tmp <- cbind(Y=runif(64),tmp)
options(max.print = 10000)

# alias pattern of all main effects.  5th order terms and higher are assumed 0
# and so not included in the model
tmp <- pm1[fraction,]
tmp <- cbind(Y=runif(64),tmp)
options(max.print = 10000)

# alias pattern of all main effects
junk <- alias(aov(Y ~ (A+B+C+D+E+F+G)^4, data = tmp))
junk[[2]][,-1]
options(max.print = 1000)

FrFac.myaov.unreplicated <- aov(response ~ (n+q+ENE+beta.mu+sigma+model+x.cor)^4, 
                                data = results2[results2$replicate==1,])
summary(FrFac.myaov.unreplicated) # Not reported in paper

cleancoef <- function(myaov){
  # Helper function to clean up the "1"s in coefficient names.
  mycoef <- coef(myaov)
  cn <- str_remove_all(names(mycoef),"1")
  names(mycoef) <- cn
  mycoef
}
# Fit a reduced version of this model to get an estimate of residual error, 
# and give SEs for effect estimates.
reduced.FrFac.myaov.unreplicated <- 
  aov(response ~ n+q+ENE+beta.mu+sigma+model+x.cor+
        sigma:model + beta.mu:model + q:x.cor+ENE:model, 
      data = results2[results2$replicate==1,])


########################################
# Figure 8 - Half normal plot for quarter fraction only
# Effect estimates are twice the estimated regression
# coefficients, since the latter correspond to +1/-1 coding.  Intercept is
# excluded because it's not an effect

pdf(file = "Figure8.pdf", width = 6, height = 6)
par(mar=c(4,4,1,1))
halfnormal(cleancoef(FrFac.myaov.unreplicated)[-1]*2, main = '',
            xlim = c(0,4), alpha = 0.25)
# Vertical line drawn at 2 standard errors of the effect.  
# Extra *2 is because regression coefficients are doubled to get effect estimates.
abline(v = 2*2*summary.lm(reduced.FrFac.myaov.unreplicated)$coefficients[2,2],
       lty = 2)
dev.off()

###############################################
# Figure 9- main effects and 2select 2fis
myheight = c(0,5.5)
pdf(file = "Figure9.pdf", width = 9, height = 4.5)
par(mfrow=c(1,3))
plot.design( response ~ n+q+ENE+beta.mu+sigma+model+x.cor, data = results2,
             ylim = myheight, main = "(a)")
text(4,-.6,labels="beta.mu",xpd=NA)
text(6,-.6,labels="model",xpd=NA)
with(results2,
     {
       interaction.plot(beta.mu, model, response, ylim = myheight, main = "(b)")
       interaction.plot(sigma, model, response, ylim = myheight, main = "(c)")
     }
)
dev.off()

# Full fraction, 2 replicates
myaov <- aov( response ~ (n+q+ENE+beta.mu+sigma+model+x.cor)^4, data = results)
# Full fraction, 1 replicate
myaov.unreplicated <- aov( response ~ (n+q+ENE+beta.mu+sigma+model+x.cor)^4, 
                           data = results[results$replicate==1,])
#Quarter fraction, 2 replicates
FrFac.myaov <- aov(formula = formula(myaov), data = results2)

###################
# Figure 15: comparision of half normal plots for all 4 possible analyses
pdf(file = "Figure15.pdf", width = 9, height = 9)
par(mfrow=c(2,2))
halfnormal(cleancoef(myaov)[-1]*2, main='Full Factorial', xlim = c(0,4.1))
halfnormal(cleancoef(myaov.unreplicated)[-1]*2, 
           main='Full Factorial, unreplicated', xlim = c(0,4.1))
halfnormal(cleancoef(FrFac.myaov)[-1]*2, 
           main = 'Quarter Fraction', xlim = c(0,4.1), alpha = 0.10)
halfnormal(cleancoef(FrFac.myaov.unreplicated)[-1]*2, 
           main = 'Quarter Fraction, unreplicated', xlim = c(0,4.1), alpha = 0.10)
dev.off()

# How do estimates of the residual standard deviation vary according to replication and fractionation?
# For the unreplicated fractional factorial we have to assume that some are small, so I delete those.
summary.lm(myaov)$sigma
summary.lm(myaov.unreplicated)$sigma
summary.lm(FrFac.myaov)$sigma
summary.lm(FrFac.myaov.unreplicated)$sigma
summary.lm(reduced.FrFac.myaov.unreplicated)$sigma




