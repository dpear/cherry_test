require(ggplot2)

# In order to run the simulation, you must have a dataset with the following parameters
# Relicate 
# Gene 
# individual 1: must indicate species and individual identifuer
# individual 2 (can be from same or different species)
# distance

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

tshift_data = function (data) {
  # Function returns the p-value of the LRT of H0: tau = 0, H1: tau = tauH
  # where tauH is the MLE of tau
  
  # 1)    Simulate an exponential distribution with gaussian noise: tau = 
  exp1 = data
  
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
  if (prod_0> prod_1) {return(c(1,tauH,lambdaN))}
  
  D    =  -2*(prod_0 - prod_1)
  pval =  1 - pchisq(D,1)
  
  lambda = lambdaN
  if (pval<.05){ lambda=lambdaH }
  return(c(pval,tauH,lambda))
}

r      = 100 # number of replicates per unique parameter combination
tau    = c(1,0.1,0.01,0.001,0)
n      = 2^(4:11)#c(10,20,200,1000,10000)
nsd    = c(0.1,0.01,0.001)
lambda = 1

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
# interesting: if you vary lambda, it has no effect on the shape of the results


# REAL DATA

distances = read.table('distance.stat',sep =' ',header=T)
unique(distances[c("ind1", "ind2")]) # 92 unique rows, so 92 replicates
distances$Rep = expand.grid(seq(1,500),seq(1,92))$Var2 # distances$Rep is mislabeled, correct

test_results = function(data){
  
  outcomes = data.frame()
  
  n_rep = length(unique(distances$Rep))
  n_genes = nrow(distances)/length(unique(distances$Rep))
  observations = nrow(distances)
  starts = seq(1,observations,by=n_genes)
  ends = seq(n_genes,observations,by=n_genes)
  
  for (i in 1:length(starts)){

    this = distances[starts[i]:ends[i],]
    out = tshift_data(this$distance)
    run = unique(this[c('Rep','ind1','ind2')])
    run$sp1 = strsplit(toString(run$ind1),'_')[[1]][1]
    run$sp2 = strsplit(toString(run$ind2),'_')[[1]][1]
    run$p = out[1]
    run$tau = out[2]
    run$lambda = out[3]
    
    # reject a true null: "false positive" 
    if (run$sp1==run$sp2 && run$p<=.05){class='FP'}
    if (run$sp1==run$sp2 && run$p>.05){class='TN'}
    if (run$sp1!=run$sp2 && run$p<=.05){class='TP'}
    if (run$sp1!=run$sp2 && run$p>.05){class='FN'}
    run$class=class
    
    outcomes = rbind(outcomes,run)
  }
  return(outcomes)
}
test_results_bootstrap = function(data,nboot){
  
  outcomes = data.frame()
  
  n_rep = length(unique(distances$Rep))
  n_genes = nrow(distances)/length(unique(distances$Rep))
  observations = nrow(distances)
  starts = seq(1,observations,by=n_genes)
  ends = seq(n_genes,observations,by=n_genes)
  
  for (i in 1:length(starts)){
    
    this = distances[starts[i]:ends[i],]
    this_boot = this[floor(runif(nboot,1,501)),]
    out = tshift_data(this_boot$distance)
    run = unique(this[c('Rep','ind1','ind2')])
    run$sp1 = strsplit(toString(run$ind1),'_')[[1]][1]
    run$sp2 = strsplit(toString(run$ind2),'_')[[1]][1]
    run$p = out[1]
    run$tau = out[2]
    run$lambda = out[3]
    
    # reject a true null: "false positive" 
    if (run$sp1==run$sp2 && run$p<=.05){class='FP'}
    if (run$sp1==run$sp2 && run$p>.05){class='TN'}
    if (run$sp1!=run$sp2 && run$p<=.05){class='TP'}
    if (run$sp1!=run$sp2 && run$p>.05){class='FN'}
    run$class=class
    
    outcomes = rbind(outcomes,run)
  }
  return(outcomes)
}

# compare: the significant p-values from bootstrap are much lower 
outcomes = test_results(distances)
outcomes_boot = test_results_bootstrap(distances,10000)

write.table(outcomes,'outcomes.txt',row.names = FALSE, col.names = TRUE)





# ? 

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



