require(ggplot2)
library(dplyr)

# In order to run tshift_data, you must have a dataset with the following parameters named as such
# Rep: Replicate number 
# Gene: Gene number 
# ind1: must indicate species and individual identifier
# ind2: can be from same or different species
# distance: observed distances

# ~~~~~~~~~ DATA FILES AND MANIPULATIONS ~~~~~~~~~~~~~~~~~~~
distances     = read.table('distance.stat',sep =' ',header=T)
ds            = read.csv('sp.distance.stat',sep=" ",header=F)
di            = read.csv('ind.distance.stat',sep=" ",header=F)
same_species  = read.table('same-species.distance.stat',sep =' ',header=F) #BIG
diff_species  = read.table('two-species.distance.stat',sep =' ',header=F) #small
tau.stat      = read.table('tau.stat',sep=' ',header=F)
distances$Rep = expand.grid(seq(1,500),seq(1,92))$Var2 # distances$Rep is mislabeled (?), correct
names(ds)     = c("Rep","Gene","ind1","ind2","distance")
names(di)     = c("Rep","Gene","ind1","ind2","distance")
names(same_species)  = c("Rep","Gene","ind1","ind2","distance")
names(diff_species)  = c("Rep","Gene","ind1","ind2","distance")
names(tau.stat)      = c("Rep","-","ind1","ind2","true_tau")
tau.stat$ind1 = as.character(tau.stat$ind1)
tau.stat$ind2 = as.character(tau.stat$ind2)
tau.stat      = unique(tau.stat[c('Rep','ind1','ind2','true_tau')])
all           = rbind(diff_species,same_species)
notnormalized = read.table('est.gt.distance.stat',sep=' ',header=F)
names(notnormalized)  = c("Rep","Gene","ind1","ind2","distance")


# ~~~~~~~~ SIMULATION ~~~~~~~~~~~~~~~~ (run shift test functions first)
simulation = function(r = 100, tau = c(1,0.1,0.01,0.001,0), n = 2^(1:8), nsd = c(0.1,0.01,0.001), lambda = 1){
  r1       = expand.grid(tau,n,nsd,lambda)
  t        = do.call("rbind", replicate(r, r1, simplify = FALSE))
  names(t) = c('tau','n','nsd','lambda')
  t$p      = mapply(tshift,t$tau,t$n,t$nsd)
  
  qplot(as.factor(n),log10(p+0.000000000001),data=t,geom="boxplot")+
    facet_grid( nsd~tau)+
    geom_hline( yintercept = log10(0.05),color="red")+
    theme_bw()+
    ggtitle(expression(paste('p-values from LR Test with varying ',tau,', sample size, and noise (rejection region below red line)')))+
    ylab('p-value (log scale)')+
    xlab('sample size (n= 10, 20, 200, 1000, 10000)')
}
simulation()
# interesting: if you vary lambda, it has no effect on the shape of the results


# ~~~~~~~~ SHIFT TEST FUNCTIONS ~~~~~~~~~~~ (run before simulations section)
tshift = function (tau,n,nsd,lambda=1) {
  # Function returns the p-value of the LRT of H0: tau = 0, H1: tau = tauH
  # where tauH is the MLE of tau
  
    # 1)    Simulate an exponential distribution with gaussian noise: tau = 
    exp1 = rexp(n,lambda) + rnorm(n,tau,nsd)
    
    # 2)    Compute the MLE for lambda, tau under the alternative, MLE for lambda under the null
    #       https://math.stackexchange.com/questions/693070/shifted-exponential-distribution-and-mle
    tauH    = max(min(exp1),0)
    lambdaH = 1/mean(exp1-tauH)
    lambdaN = 1/mean(exp1)
    
    # 3)    LRT of H0: tau=0, H1: tau = tauH
    pdf0   = log( lambdaN*exp(-lambdaN*exp1) )
    pdf1   = log( lambdaH*exp(-lambdaH*(exp1-tauH)) )
    prod_0 = sum( pdf0[pdf1<=0] )
    prod_1 = sum( pdf1[pdf1<=0] )
    
    # if the logliklihood of the null model is greater than the alternative, then we know we MUST reject
    # interesting to note that this only happens once in the simulated data
    # thinking about this in terms of ratios, this would be >1, then logged, then negated
    # resulting in a negative test statistic
    if (prod_0> prod_1) {return(1)}
    
    D    =  -2*(prod_0 - prod_1)
    pval =  1 - pchisq(D,1)
    
    return(pval)
}
tshift_data = function (exp1, q=F, mom_prune=F, mle_prune=F, prune=F, nest_prune=F, truncated=T) {
  # Function returns the p-value of the LRT of H0: tau = 0, H1: tau = tauH
  
  #     1) Compute the estimates for lambdaH, tauH under the alternative, MLE for lambdaN under the null using specified method
  
  # Using MLE's for all parameters; default
  tauH    = max(min(exp1),0)
  lambdaH = 1/mean(exp1-tauH)
  lambdaN = 1/mean(exp1)
  # Using MOM of tau to remove points;
  if (mle_prune!=F){
    tauH    = max(mean(exp1) - sqrt(sum((mean(exp1)-exp1)^2)/length(exp1)),0)
    exp1    = exp1[exp1>=tauH]
    tauH    = max(min(exp1),0) # then using MLE of tauH, lambdaN, lambdaH
    lambdaH = 1/mean(exp1-tauH)
    lambdaN = 1/mean(exp1)
  }
  # Using MOM of tau to remove points; 
  if (mom_prune!=F){
    tauH    = max(mean(exp1) - sqrt(sum((mean(exp1)-exp1)^2)/length(exp1)),0)
    exp1    = exp1[exp1>=tauH]
    tauH    = max(min(exp1),0) # then using MOM of tauH, lambdaN, lambdaH
    lambdaH = sqrt(length(exp1)/sum((mean(exp1)-exp1)^2)) 
    lambdaN = 1/mean(exp1)
  }
  # quantile estimation; quant = c(q1,q2)
  if (q[1]!=F){
    d1 = sort(exp1)[q[1]*length(exp1)]
    d2 = sort(exp1)[q[2]*length(exp1)]
    lambdaQ = (log(1-q[1])-log(1-q[2]))/(d2-d1+.000001) #otherwise can be infinity = bad
    tauQ = d2+ (log(1-q[2])/lambdaQ)
    exp1 = exp1[exp1>tauQ]
    tauH    = tauQ
    lambdaH = lambdaQ
    lambdaN = lambdaQ
  }
  # prune bottom <prune> of data
  if (prune!=F){
    threshold = sort(exp1)[prune*length(exp1)]
    exp1 = exp1[exp1>threshold]
    tauH = min(exp1)
    lambdaH = 1/mean(exp1-tauH)
    lambdaN = 1/mean(exp1)
  }
  # prune bottom <prune> of data before computing MOM
  if (nest_prune!=F){
    threshold = sort(exp1)[nest_prune*length(exp1)]
    exp1 = exp1[exp1>threshold]
    tauH    = max(mean(exp1) - sqrt(sum((mean(exp1)-exp1)^2)/length(exp1)),0)
    lambdaH = sqrt(length(exp1)/sum((mean(exp1)-exp1)^2)) 
    lambdaN = 1/mean(exp1)
  }
  
  # 3)    LRT of H0: tau=0, H1: tau = tauH
  pdf0   = log( lambdaN*exp(-lambdaN*exp1) )
  pdf1   = log( lambdaH*exp(-lambdaH*(exp1-tauH)) )
  prod_0 = sum( pdf0 )
  prod_1 = sum( pdf1 )
  
  if (truncated==T){
    pdf0   = log( lambdaN*exp(-lambdaN*exp1)/(1-exp(-lambdaN*2)) )
    pdf1   = log( lambdaH*exp(-lambdaH*(exp1-tauH))/(1-exp(-lambdaH*(2-tauH))) )
    prod_0 = sum( pdf0 )
    prod_1 = sum( pdf1 )
  }
  
  if (prod_0> prod_1) {return(c(1,tauH,lambdaN))}
  if (tauH==0) {return(c(1,tauH,lambdaN))}
  
  D    =  -2*(prod_0 - prod_1)
  pval =  1 - pchisq(D,1)
  
  lambda = lambdaN
  if (pval<.05){ lambda=lambdaH }
  return(c(pval,tauH,lambda))
}
test_results = function (data,nboot=F,q=F, mom_prune=F, mle_prune=F, prune=F, nest_prune=F,truncated=F,pr_curve=F){
  
  outcomes = data.frame()
  
  n_genes      = max(data$Gene)
  observations = nrow(data)
  starts       = seq(1,observations,by=n_genes)
  ends         = seq(n_genes,observations,by=n_genes)
  n            = length(starts)
  pb           = txtProgressBar(min = 0, max = n, style = 3)
  
  for (i in 1:n){
    
    this = data[starts[i]:ends[i],]
    if (nboot!=F){ this = this[floor(runif(nboot,1,n_genes+1)),]}
    out        = tshift_data(this$distance,q, mom_prune, mle_prune, prune, nest_prune, truncated)
    run        = unique(this[c('Rep','ind1','ind2')])
    run$sp1    = strsplit(toString(run$ind1),'_')[[1]][1]
    run$sp2    = strsplit(toString(run$ind2),'_')[[1]][1]
    run$p      = out[1]
    run$tau    = out[2]
    run$lambda = out[3]
    
    # reject a true null: "false positive" 
    class = 'FP'
    if (run$sp1==run$sp2 && run$p>.05){class='TN'}
    if (run$sp1!=run$sp2 && run$p<=.05){class='TP'}
    if (run$sp1!=run$sp2 && run$p>.05){class='FN'}
    run$class=class
    
    outcomes = rbind(outcomes,run)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  outcomes = left_join(outcomes, tau.stat, by = c("sp1" = "ind1", "sp2" = "ind2","Rep"="Rep"))
  print(table(outcomes$class))
  return(outcomes)
}


# ~~~~~~~~~~~ NEW SIMULATED DATASET ~~~~~~~~~~~~~~~~~~
# poor performance on our new simulated dataset :(

ggplot(data=outcomes, aes(x=tau,y=true_tau,color=log(p+.000000001)))+geom_point(size=5) + scale_color_distiller(palette = "RdPu")
n = nrow(outcomes)
outcomes$prune = rep(pr,n)

data = outcomes1
pos = data.frame(scores=data[data$sp1!=data$sp2,]$p,labels=1)
neg = data.frame(scores=data[data$sp1==data$sp2,]$p,labels=0)
formatted = rbind(pos,neg)
data.mm = mmdata(formatted$scores,formatted$labels)
e=evalmod(data.mm)
plot(e)

pr=pr.curve(pos$p,neg$p,curve=TRUE)
plot(pr)



outcomes1 = test_results(all)
outcomes2 = test_results(all,mom_prune = T)
outcomes3 = test_results(all,mle_prune = T)
outcomes4 = test_results(all,nest_prune = .05)
outcomes5 = test_results(all,nest_prune = .1)
outcomes6 = test_results(all,nest_prune = .08)

outcomes1_n = test_results(notnormalized)
outcomes2_n = test_results(notnormalized,mom_prune = T)
outcomes3_n = test_results(notnormalized,mle_prune = T)
outcomes4_n = test_results(notnormalized,nest_prune = .05)
outcomes5_n = test_results(notnormalized,nest_prune = .1)
outcomes6_n = test_results(notnormalized,nest_prune = .08)

make_curve = function(data){ return(pr.curve(data[data$sp1!=data$sp2,]$p,data[data$sp1==data$sp2,]$p,curve=TRUE))}

plot(make_curve(outcomes1))
plot(make_curve(outcomes2))
plot(make_curve(outcomes3))
plot(make_curve(outcomes4))
plot(make_curve(outcomes5))
plot(make_curve(outcomes6))

plot(make_curve(outcomes1_n))
plot(make_curve(outcomes2_n))
plot(make_curve(outcomes3_n))
plot(make_curve(outcomes4_n))
plot(make_curve(outcomes5_n))
plot(make_curve(outcomes6_n))



outcomes1$x = rep(1,nrow(outcomes1))
outcomes2$x = rep(2,nrow(outcomes2))
outcomes3$x = rep(3,nrow(outcomes3))
outcomes4$x = rep(4,nrow(outcomes4))
outcomes5$x = rep(5,nrow(outcomes5))
outcomes6$x = rep(6,nrow(outcomes6))

d = notnormalized[notnormalized$ind1=='163_0_4' & notnormalized$ind2=='80_0_3',]

outcomes = rbind(outcomes1,outcomes2,outcomes3,outcomes4,outcomes5,outcomes6)

ggplot(data=outcomes, aes(x=as.factor(x),y=true_tau-tau,color=log(p+.000000001)))+geom_boxplot() + scale_color_distiller(palette = "RdPu")

head(new_outcomes)
View(new_outcomes)

par(mfrow= c(1,2))
plot(p~tau, data=new_outcomes, main="estimated tau vs. p")
plot(p~true_tau, data=new_outcomes,main = "true tau vs. p")
require(reshape2)
moutcome = melt(new_outcomes,measure.vars=c(tau,true_tau))
head(moutcome)
ggplot(data = moutcome, aes(x=value,p))+geom_point(alpha=0.2)+facet_grid(~variable)+theme_bw()
# what is the noise term from the simulated data (presumably high?)
# problem: estimating tau=0; will never reject if tau=0
# even when tau is estimated to be small nonzero, our test still gives TP
# problem: if we estimate tau to be .2 (lol why not) we get better results (?) BAD
# this means we should test against != 0 
# using method of moments produces way better results
#     - discard all datapoints less than tau 


#same_species is fine
outcomes = test_results(same_species)

# ~~~~~~~~~~ SAVING OUTCOMES ~~~~~~~~~~~~~
write.table(outcomes,'outcomes.txt',row.names = FALSE, col.names = TRUE)








# ~~~~~~~ CODE GRAVEYARD (rest in bits) ~~~~~~~~~~~~~~~~ 

# the lambdas are normally distributed
plot(density((1/outcomes$lambda)^2))
# but after bootstrapping we might think it is tri-modal? 
plot(density(outcomes$lambda[floor(runif(1000,1,length(outcomes$lambda)))]))
# which is just amplified noise, seen from increasing # bootstraps 
plot(density(outcomes$lambda[floor(runif(10000,1,length(outcomes$lambda)))]))
# 1/variance
sum((outcomes$lambda-mean(outcomes$lambda))^2)/length(outcomes)
# the bootstrap sample dist has the same shape
plot(density(outcomes_boot$lambda))
# 1/variance
sum((outcomes_boot$lambda-mean(outcomes_boot$lambda))^2)/length(outcomes_boot)
table(outcomes$class)
outcomes2 = test_results(same_species[1:115000,])
table(outcomes2$class)
tshift_data_aic = function (data) {
  # Function returns the p-value of the LRT of H0: tau = 0, H1: tau = tauH
  # where tauH is the MLE of tau
  
  # 1)    Simulate an exponential distribution with gaussian noise: tau = 
  exp1 = data
  
  # 2)    Compute the MLE for lambda, tau under the alternative, MLE for lambda under the null
  #       https://math.stackexchange.com/questions/693070/shifted-exponential-distribution-and-mle
  tauH    = max(min(exp1),0) #Max lik estimator
  #tauH = mean(exp1) - sqrt(sum((mean(exp1)-exp1)^2)/length(exp1)) #MOM estimator
  lambdaH = 1/mean(exp1-tauH)
  lambdaN = 1/mean(exp1)
  
  # 3)    LRT of H0: tau=0, H1: tau = tauH
  #exp1=exp1[exp1>tauH]
  pdf0   = log( lambdaN*exp(-lambdaN*exp1) )
  pdf1   = log( lambdaH*exp(-lambdaH*(exp1-tauH)) )
  prod_0 = sum( pdf0[pdf1<=0] )
  prod_1 = sum( pdf1[pdf1<=0] )
  k0 = 1 # number of parameters estimated if tau=0
  k1 = 2 # number of parameters estimated if tau=tauH
  aic0 = 2*k0 - 2*prod_0
  aic1 = 2*k1 - 2*prod_1
  
  lambda = lambdaN
  if (aic1<aic0){ lambda=lambdaH }
  return(c(as.numeric(aic1<aic0),tauH,lambda))
}
test_results_aic = function(data){
  
  outcomes = data.frame()
  
  n_genes = max(data$Gene)
  observations = nrow(data)
  starts = seq(1,observations,by=n_genes)
  ends = seq(n_genes,observations,by=n_genes)
  
  for (i in 1:length(starts)){
    
    this = data[starts[i]:ends[i],]
    out = tshift_data_aic(this$distance)
    run = unique(this[c('Rep','ind1','ind2')])
    run$sp1 = strsplit(toString(run$ind1),'_')[[1]][1]
    run$sp2 = strsplit(toString(run$ind2),'_')[[1]][1]
    run$p = out[1]
    run$tau = out[2]
    run$lambda = out[3]
    
    # reject a true null: "false positive" 
    if (run$sp1==run$sp2 && run$p==1){class='FP'}
    if (run$sp1==run$sp2 && run$p==0){class='TN'}
    if (run$sp1!=run$sp2 && run$p==1){class='TP'}
    if (run$sp1!=run$sp2 && run$p==0){class='FN'}
    run$class=class
    
    outcomes = rbind(outcomes,run)
  }
  return(outcomes)
}
outcomes = test_results(distances)
outcomes_boot = test_results_bootstrap(distances,10000)
outcomes = test_results(ds)
outcomes_boot = test_results_bootstrap(ds,10000)
outcomes = test_results(di)
outcomes_boot = test_results_bootstrap(di,10000)
outcomes_aic = test_results_aic(diff_species)
outcomes_boot = test_results_bootstrap(diff_species,200)
table(outcomes_boot$class) # view counts of all classifications
table(outcomes_aic$class) # view counts of all classifications
table(outcomes$tau!=0)
table(outcomes_boot$class) # view counts of all classifications from bootstrap
tshift_aic = function (tau,n,nsd,lambda=1) {
  # Function returns the lowest aic value from aic model selection
  
  # 1)    Simulate an exponential distribution with gaussian noise: tau = tau
  exp1 = rexp(n,lambda) + rnorm(n,tau,nsd)
  
  # 2)    Compute the MLE for lambda, tau under the alternative, MLE for lambda under the null
  #       https://math.stackexchange.com/questions/693070/shifted-exponential-distribution-and-mle
  tauH    = max(min(exp1),0)
  lambdaH = 1/mean(exp1-tauH)
  lambdaN = 1/mean(exp1)
  
  # 3)    Compute AIC for tau=0 and tau=tauH
  pdf0   = log( lambdaN*exp(-lambdaN*exp1) )
  pdf1   = log( lambdaH*exp(-lambdaH*(exp1-tauH)) )
  prod_0 = sum( pdf0[pdf1<=0] )
  prod_1 = sum( pdf1[pdf1<=0] )
  k0 = 1 # number of parameters estimated if tau=0
  k1 = 2 # number of parameters estimated if tau=tauH
  aic0 = 2*k0 - 2*prod_0
  aic1 = 2*k1 - 2*prod_1
  
  # 1: tau model is better; 0: null model is better
  return(as.numeric(aic1<aic0))
}

# AIC model selection simulation 
r      = 100 # num replicates / unique parameter combo
tau    = c(1,0.1,0.01,0.001,0)
n      = 2^(1:8) 
nsd    = c(0.1,0.01,0.001)
lambda = 1

r1       = expand.grid(tau,n,nsd,lambda)
t        = do.call("rbind", replicate(r, r1, simplify = FALSE))
names(t) = c('tau','n','nsd','lambda')
t$p      = mapply(tshift_aic,t$tau,t$n,t$nsd)
t2 = aggregate(p~tau+n+nsd+lambda, t, sum)

qplot(as.factor(n),p+0.000000000001,data=t2)+
  facet_grid( nsd~tau)+
  theme_bw()+
  ggtitle(expression(paste('Model with estimated tau is selected by AIC (%) with varying ',tau,', sample size, and noise (rejection region below red line)')))+
  ylab('% of simulations where model with tau != 0 is better')+
  xlab('sample size (n= 10, 20, 200, 1000, 10000)')
