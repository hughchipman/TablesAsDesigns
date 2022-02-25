#### This version is a cleaned up copy of the version used to produce figures
#### for the Oct 2021 CJS submission.

# A hack to give a bit more room for the legend in the interaction plot
source("interaction.plot.R") 

options(contrasts = c("contr.sum","contr.poly"))
y = scan('KrishnamoorthyTable1.txt',what=numeric(0))
scalefactor <- 100
target <- 0.05*scalefactor
y <- y*scalefactor # so we are talking 5%, not 0.05.  SS values are easier.
x = expand.grid(tail=c('L','R','T'),
                sigma=c(1,2,3),
                method=c('GV','AN','SL','MS'),
                p0=c(.2,.3,.5,.7),
                n=c(20,30,50)
)

x$sigma = ordered(x$sigma)
x$p0 = ordered(x$p0)
x$n = ordered(x$n)

mydata=data.frame(x,y)

########################################
# Note: Table 1 is the KMM data, which is a separate file
#       Table 2 is written in the latex document

########################################
# Table 3: ANOVA table for 432-run full-factorial experiment
full1 = aov(y~(method + tail + n + p0 + sigma),data=mydata)
full2 = aov(y~(method + tail + n + p0 + sigma)^2,data=mydata)
full3 = aov(y~(method + tail + n + p0 + sigma)^3,data=mydata)
full4 = aov(y~(method + tail + n + p0 + sigma)^4,data=mydata)
summary(full4)

##############################################################
# Primary Analysis: Drop runs corresponding to method = AN 
##############################################################

mydata2 <- mydata[mydata$method!='AN',]
# make "method" a factor without the AN level:
mydata2$method <- as.factor(as.character(mydata2$method))
# combine method and tail into a single control factor, for use with
#   interaction plots
mydata2$control <- as.factor(paste(
  as.character(mydata2$method),
  as.character(mydata2$tail), sep = '/'
))

##############
# Note: Cause-and-effect diagrams in Figures 1-2 are from powerpoint

####################################
# Table 4: ANOVA table for 324-run full-factorial experiment excluding method=AN
# fit an ANOVA to the reduced data
newfull1 = aov(y~(method + tail + n + p0 + sigma),data=mydata2)
newfull2 = aov(y~(method + tail + n + p0 + sigma)^2,data=mydata2)
newfull3 = aov(y~(method + tail + n + p0 + sigma)^3,data=mydata2)
newfull4 = aov(y~(method + tail + n + p0 + sigma)^4,data=mydata2)
newfull5 = aov(y~(method + tail + n + p0 + sigma)^5,data=mydata2)
summary(newfull4)


##############
# Figure 3: Main Effects plot and select 2fis plots
#

pdf('Figure3.pdf', width = 9, height = 4.5)
par(mfrow=c(1,3))
myylim = c(0.033,0.065)*scalefactor 
plot.design(formula(full1),data=mydata2,lwd=2, 
            ylim = myylim, main = "(a)")
text(2.1,3.1,labels="tail",xpd=NA)
abline(h=target, col='blue')
# method by tail:
with(mydata2, interaction.plot(method, tail,y,legend=TRUE,lwd=2, 
                               ylim = myylim, main="(b)", fixed = TRUE))
abline(h=target,col='blue')
# p0 by tail:
with(mydata2,interaction.plot(p0, tail,y,legend=TRUE,lwd=2, 
                              ylim = myylim, main="(c)", fixed = TRUE))
abline(h=target,col='blue')
dev.off()


#####################################
# Figure 4: Histograms of response grouped by 9 combinations of 
# control factor method:tail
#

pdf('Figure4.pdf', width = 12, height = 8)
par(mfrow=c(3,3))

mylim <- c(0,0.11)*scalefactor
for (imethod in c("GV","MS","SL")){
  for (itail in c("L","R","T")){
    with(mydata2, 
         hist(y[(method==imethod) & (tail == itail)], 
              main = paste(imethod,itail,sep="/"),
              xlim = mylim, breaks = seq(min(mylim), max(mylim), by = 0.5),
              xlab = "type I error", ylim = c(0,20))
    )
    abline(v=0.05*scalefactor, col = "blue", lwd=2)
  }
}
dev.off()



############### 
# Figure 5: 3-way interactions method:tail:p0 and method:tail:sigma
# used for TRPD.

pdf('Figure5.pdf', width = 14, height = 7)
par(mfrow=c(1,2))
myylim <- c(.01,.09)*scalefactor
with(mydata2,
     interaction.plot(
       x.factor = p0, trace.factor = control, response = y,
       lty = c(1,1,1,2,2,2,3,3,3), col = c(1:3,1:3,1:3)+1,
       lwd = 3, fixed = TRUE, ylim = myylim, trace.label = "Method/Tail", main = "(a)"
     ))
abline(h=target,col='blue')
with(mydata2,
     interaction.plot(
       x.factor = sigma, trace.factor = control, response = y,
       lty = c(1,1,1,2,2,2,3,3,3), col = c(1:3,1:3,1:3)+1,
       lwd = 3, fixed = TRUE, ylim = myylim, trace.label = "Method/Tail", main = "(b)"
     ))
abline(h=target,col='blue')
dev.off()

# Not reported in the paper - residual check
pdf('Krish_resids.pdf', width = 12, height = 12)
par(mfrow=c(2,3))
plot(full4, which = 1:2)     # all runs
boxplot(residuals(full4)~ mydata$method, xlab = 'method')
plot(newfull4, which = 1:2)  # method = AN excluded, this is an improvement.
boxplot(residuals(newfull4)~ mydata2$method, xlab = 'method')
dev.off()
# Since the residuals appear reasonable, transformations are not pursued.
# For example, although a logistic transform might be worth considering, the 
# plots suggest this need not be pursued

###################
# Note: Statistical Learning example is a separate file.  It has (or refers to)
#       Table 5 and Figures 6, 7, 8 and Supplementary Figure 2.


######################
# Supplemental Figure 1 - 4 way interaction for TRPD
pdf('SuppFigure1.pdf', width = 14, height = 7)
par(mfrow=c(1,3))
myylim <- c(.01,.09)*scalefactor
for (i in 1:3){
  with(mydata2[mydata2$sigma == i,],
       interaction.plot(
         x.factor = p0, trace.factor = control, response = y,
         lty = c(1,1,1,2,2,2,3,3,3), col = c(1:3,1:3,1:3)+1,
         lwd = 3, fixed = TRUE, ylim = myylim, trace.label = "Method/Tail"
       ))
  title(paste('sigma =',i))
  abline(h=target,col='blue')
}
dev.off()


