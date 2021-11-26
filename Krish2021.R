#### This version is a cleaned up copy of the version used to produce figures
#### for the Oct 2021 CJS submission.

#### For an earlier analysis, see "Krish.R" in folder
#### Dropbox\Work\unpack\research\TablesAsDesigns\2015_SSC_Analysis

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
table(mydata2$control, mydata2$method) # check 
table(as.numeric(mydata2$control), mydata2$control)
range(mydata$y)  # range of type I errors for all runs
range(mydata2$y) # range of type I errors without AN runs

##############
# Note: Figures 1-3 are from powerpoint, corresponding to fishbone diagrams.

##############
# Figure 4: Main Effects plot and select 2fis plots
#

pdf('Krish_ME2fi.pdf', width = 9, height = 4.5)
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



####################################
# Table 4: ANOVA table for 324-run full-factorial experiment excluding method=AN
# fit an ANOVA to the reduced data
newfull1 = aov(y~(method + tail + n + p0 + sigma),data=mydata2)
newfull2 = aov(y~(method + tail + n + p0 + sigma)^2,data=mydata2)
newfull3 = aov(y~(method + tail + n + p0 + sigma)^3,data=mydata2)
newfull4 = aov(y~(method + tail + n + p0 + sigma)^4,data=mydata2)
newfull5 = aov(y~(method + tail + n + p0 + sigma)^5,data=mydata2)
summary(newfull4)

###############################
# MSEs reported in text of paper:
summary.lm(full4)$sigma^2
summary.lm(newfull4)$sigma^2

#####################################
# Figure 5: Histograms of response grouped by 9 combinations of 
# control factor method:tail
#

pdf('Krish_Taguchi3.pdf', width = 12, height = 8)
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
# Figure 6: 3-way interactions method:tail:p0 and method:tail:sigma
# used for TRPD.

pdf('Krish_Taguchi2.pdf', width = 14, height = 7)
par(mfrow=c(1,2))
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


##############################################
# What do we lose if one of the 5 factors is completely dropped from the ANOVA?
fourfac.wo.method <- aov(y~(tail + n + p0 + sigma)^4,data=mydata2)
fourfac.wo.tail <- aov(y~(method + n + p0 + sigma)^4,data=mydata2)
fourfac.wo.n <- aov(y~(method + tail + p0 + sigma)^4,data=mydata2)
fourfac.wo.p0 <- aov(y~(method + tail + n + sigma)^4,data=mydata2)
fourfac.wo.sigma <- aov(y~(method + tail + n + p0)^4,data=mydata2)

summary.lm(fourfac.wo.method)$r.sq  # 26.88% explained
summary.lm(fourfac.wo.tail)$r.sq    # 13.47% explained
summary.lm(fourfac.wo.n)$r.sq       # 88.79% explained, reported as "about 89%" in the paper
summary.lm(fourfac.wo.p0)$r.sq      # 72.95% explained
summary.lm(fourfac.wo.sigma)$r.sq   # 74.70% explained

summary.lm(newfull4)$r.sq           # 98.09% explained

########################
# Not reported in the paper - residual check
pdf('Krish_resids.pdf', width = 12, height = 12)
par(mfrow=c(2,2))
plot(full4, which = 1:2)     # all runs
plot(newfull4, which = 1:2)  # method = AN excluded, this is an improvement.
dev.off()

###################
# "cheapo" analysis - simplify numeric factors to 2 levels (extrema)
# Table 5 in the paper

whch = (mydata$sigma != 2) & (mydata$n != 30) & (mydata$p0!=0.3) & (mydata$p0!=0.5) & (mydata$method != "AN")
mysmalldata= mydata[whch,]
small2 = aov(formula(full2),data=mysmalldata)
small3 = aov(formula(full3),data=mysmalldata)
small4 <- aov(y~(method + tail + n + p0 + sigma)^4, data = mysmalldata)
small5 <- aov(y~(method + tail + n + p0 + sigma)^5, data = mysmalldata)
summary(small4)  # Table 5

###################
# Note: Statistical Learning example is a separate file.  It has (or refers to)
#       Table 6 and Figures 7, 8, 9 and 15.


######################
# Figure 10 in supplemental material - 4 way interaction for TRPD
pdf('Krish_Taguchi.pdf', width = 14, height = 7)
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


###### 
# For interaction plots, we recast the factors without missing levels.
mysmalldata2 <- mysmalldata
mysmalldata2$sigma <- as.factor(as.character(mysmalldata$sigma))
mysmalldata2$method <- as.factor(as.character(mysmalldata$method))
mysmalldata2$p0 <- as.factor(as.character(mysmalldata$p0))
mysmalldata2$n <- as.factor(as.character(mysmalldata$n))
mysmalldata2$control <- as.factor(paste(
  as.character(mysmalldata2$method),
  as.character(mysmalldata2$tail), sep = '/'
))

#######################
# Figure 11 - main effects and 2fi plots for cheapo analysis 
pdf('Krish_cheapo_ME2fi.pdf', width = 9, height = 4.5)
par(mfrow=c(1,3))
myylim = c(0.033,0.065)*scalefactor 
plot.design(formula(full1),data=mysmalldata2,lwd=2, 
            ylim = myylim, main = "(a)")
abline(h=target, col='blue')

with(mysmalldata2, interaction.plot(method, tail,y,legend=TRUE,lwd=2, 
                               ylim = myylim, main="(b)", fixed = TRUE))
abline(h=target,col='blue')
with(mysmalldata2,interaction.plot(p0, tail,y,legend=TRUE,lwd=2, 
                              ylim = myylim, main="(c)", fixed = TRUE))
abline(h=target,col='blue')
dev.off()

###################
# Figure 12 - hist of response broken down by 9 control categories, for cheapo.
pdf('Krish_cheapo_hist.pdf', width = 12, height = 8)
par(mfrow=c(3,3))
for (imethod in c("GV","MS","SL")){
  for (itail in c("L","R","T")){
    with(mysmalldata, 
         hist(y[(method==imethod) & (tail == itail)], 
              main = paste(imethod,itail,sep="/"),
              xlim = mylim, breaks = seq(min(mylim), max(mylim), by = 0.5),
              xlab = "type I error", ylim = c(0,5))
    )
    abline(v=0.05*scalefactor, col = "blue", lwd=2)
  }
}
dev.off()

################
# Figure 13, 3-way interacions.
pdf('Krish_cheapo_3fi.pdf', width = 14, height = 7)
par(mfrow=c(1,2))
with(mysmalldata2,
     interaction.plot(
       x.factor = p0, trace.factor = control, response = y,
       lty = c(1,1,1,2,2,2,3,3,3), col = c(1:3,1:3,1:3)+1,
       lwd = 3, fixed = TRUE, ylim = myylim, trace.label = "Method/Tail",
       main = "(a)"
     ))
abline(h=target,col='blue')
with(mysmalldata2,
     interaction.plot(
       x.factor = sigma, trace.factor = control, response = y,
       lty = c(1,1,1,2,2,2,3,3,3), col = c(1:3,1:3,1:3)+1,
       lwd = 3, fixed = TRUE, ylim = myylim, trace.label = "Method/Tail",
       main = "(b)"
     ))
abline(h=target,col='blue')
dev.off()

##############
# Figure 14 - 4 factor interactions.
pdf('Krish_cheapo_4fi.pdf', width = 14, height = 7)
par(mfrow=c(1,3))
myylim <- c(.01,.09)*scalefactor
for (i in 1:3){
  with(mysmalldata2[mysmalldata2$sigma == i,],
       interaction.plot(
         x.factor = p0, trace.factor = control, response = y,
         lty = c(1,1,1,2,2,2,3,3,3), col = c(1:3,1:3,1:3)+1,
         lwd = 3, fixed = TRUE, ylim = myylim, trace.label = "Method/Tail"
       ))
  title(paste('sigma =',i))
  abline(h=target,col='blue')
}
dev.off()


################
# As noted above, the last figure in the paper is for the statistical learning
# example.