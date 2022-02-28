# analysis of KMM data using a Binomial GLM

options(contrasts = c("contr.sum","contr.poly"))
y = scan('KrishnamoorthyTable1.txt',what=numeric(0))
x = expand.grid(tail=c('L','R','T'),sigma=c(1,2,3),method=c('GV','AN','SL','MS'),
                p0=c(.2,.3,.5,.7),n=c(20,30,50))

x$sigma = as.factor(x$sigma)
x$p0 = as.factor(x$p0)
x$n = as.factor(x$n)

mydata=data.frame(x,y)

##############################################################
# Data frame that drops runs with method = AN

mydata2 <- mydata[mydata$method!='AN',]
mydata2$method <- as.factor(as.character(mydata2$method))
mydata2$control <- as.factor(paste(
  as.character(mydata2$method),
  as.character(mydata2$tail), sep = '/'
))


########  Binary GLM model with logistic link
# According to 1st author (personal communication), the number of replicates
# varies by method.  2,500 were used for GV, 10,000 for the other 3.
reps <- c(2500, 10000, 10000, 10000)
names(reps) <- c('GV','AN','SL','MS')
trials <- rep(-1, nrow(mydata))
for (i in 1:4){
  trials[mydata$method==names(reps)[i]] <- reps[i]
}
success <- mydata$y*trials

# Due to rounding in the table, we lose some resolution.  For 10000 replicate
# cases, all success counts are multiples of 10.  For 2500 replicate cases, some
# entries have counts that end in .5.  This is fixed by randomly adding or
# subtracting .5 to each non-integer count.
not.whole <- abs(round(success)-success) > .00001
set.seed(1)
success[not.whole] <- success[not.whole] + 
  sample(c(-0.5,0.5), size = sum(not.whole), replace = TRUE)
zero <- success < .00001

fail <- trials - success

# response is a 2-column matrix of successes and failures.
logitresponse <- cbind(success, fail)
logitresponse2 <- logitresponse[mydata$method != "AN",]


logitfull4 <- with(mydata, glm(logitresponse~(method+tail+n+p0+sigma)^4, 
                               family = binomial()))

newlogitfull4 <- with(mydata2, glm(logitresponse2~(method+tail+n+p0+sigma)^4, 
                                   family = binomial()))
# Look at significance of terms.  There are many terms, and quite a few of 3rd
# and 4th order are significant.
options(max.print = 2000)
summary(logitfull4)
summary(newlogitfull4)
options(max.print = 1000)

# Look at dropping back to a 3rd or 2nd order model...
newlogitfull2 <- with(mydata2, glm(logitresponse2~(method+tail+n+p0+sigma)^2, 
                                   family = binomial()))
newlogitfull3 <- with(mydata2, glm(logitresponse2~(method+tail+n+p0+sigma)^3, 
                                   family = binomial()))
anova(newlogitfull4,newlogitfull3, test = "LRT")
anova(newlogitfull4,newlogitfull2, test = "LRT")


