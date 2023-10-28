#Group members: Michelle Homann, Quinn Smith, Timon Keller


library(ape)
library(phylolm)
library(lme4)
library(nlme)
library(mvtnorm)
library(phytools)
library(phangorn)
library(rr2)
library(phyr)
library(logistf)


#Question 2

# Compute the rejection rates for the 5 methods.
# Note that the same starter phylogeny is used for the simulations, since the critical values of lam and alpha used in the bootstraps around H0 could depend on the topology.

#TK: This is code Tony shared with Garrett - I got very close, but did not specify alpha to be 50 - I was struggling with finding the right value. Because this code is much cleaner than mine, I adopted it. 

# 3.5.6 Type I errors and power for tests of phylogenetic signal

# For the parametric bootstrap of the LRT under H0, it is necessary to generate the critical values of the LLR. These are placed in LLR.crit.
# Note that these bootstraps to get the critical values of LLR.lam and LLR.OU are performed for the same phylogenetic tree as used for the power simulations later. This is because the critical LLR.lam and LLR.OU could depend on the topology of the phylogeny.

n <- 30
#lam <- 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha = 50
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nboot <- 20

phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
phy.OU <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = 50))$tree
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

boot0 <- data.frame(LLR.lam=rep(NA, nboot), LLR.OU=NA)
for(i in 1:nboot){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Y.sim <- rTraitCont(phy.lam, model = "BM", sigma = 1)
  Y.sim <- rTraitCont(phy.OU, model = "BM", sigma = 1)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
  z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "OUfixedRoot")
  z.0.sim <- lm(Y.sim ~ 1)
  
  boot0$LLR.lam[i] <- 2*(z.lam.sim$logLik - logLik(z.0.sim))
  boot0$LLR.OU[i] <- 2*(z.OU.sim$logLik - logLik(z.0.sim))
}
LLR.lam.crit <- sort(boot0$LLR.lam)[floor(0.95 * nboot)]
LLR.OU.crit <- sort(boot0$LLR.OU)[floor(0.95 * nboot)]

# Compute the rejection rates for the 5 methods.
# Note that the same starter phylogeny is used for the simulations, since the critical values of lam and alpha used in the bootstraps around H0 could depend on the topology.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha.list <- c(1,2,5,10,20,50)
nsims <- 20
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reject <- data.frame(alpha=rep(NA, nsims*length(alpha.list)), lam.LRT=NA, OU.LRT=NA, lam.boot0=NA, OU.boot0=NA, K.perm=NA)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counter <- 0
for(alpha in alpha.list){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #phy.lam <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam))$tree
  phy.OU <- transf.branch.lengths(phy=phy, model="OUfixedRoot", parameters=list(alpha = alpha))$tree
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(i in 1:nsims){
    counter <- counter + 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #reject$lam[counter] <- lam       
    #Y.sim <- rTraitCont(phy.lam, model = "BM", sigma = 1)
    reject$alpha[counter] <- alpha  
    Y.sim <- rTraitCont(phy.OU, model = "BM", sigma = 1)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # LRTs
    z.lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
    z.OU.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "OUfixedRoot")
    z.0.sim <- lm(Y.sim ~ 1)
    LLR.lam <- 2*(z.lam.sim$logLik - logLik(z.0.sim)[1])
    LLR.OU <- 2*(z.OU.sim$logLik - logLik(z.0.sim)[1])
    reject$lam.LRT[counter] <- (pchisq(LLR.lam, df=1, lower.tail=F) < 0.05)
    reject$OU.LRT[counter] <- (pchisq(LLR.OU, df=1, lower.tail=F) < 0.05)
    
    # Bootstraps around H0
    reject$lam.boot0[counter] <- (LLR.lam > LLR.lam.crit)
    reject$OU.boot0[counter] <- (LLR.OU > LLR.OU.crit)
    
    # Permutation test with Blomberg's K
    z.phylosig.K <- phylosig(Y.sim, tree=phy, method="K", test=TRUE, nsim=2000)
    reject$K.perm[counter] <- (z.phylosig.K$P < 0.05)
  }
}
# This code took a long time to run, so I just saved the output and then reloaded it later. I've commented out these lines though.
# write.table(reject, file="Phylo signal power curves lam n=30.csv", sep=",", row.names=F)
# reject <- read.csv(file="Phylo signal power curves lam n=30.csv")

#Create power curves
power <- aggregate(reject[,c("lam.LRT","OU.LRT","lam.boot0","OU.boot0","K.perm")], by = list(reject$alpha), FUN = mean)
names(power)[1] <- "alpha"

# Fig. 3.8, left panel. I've left the right panel as an exercise.
par(mfrow=c(1,2))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot(lam.LRT ~ alpha, data=power, typ="l", xlim=c(0,1), ylim=c(0, 1), xlab=expression(paste("Phylogenetic signal (", lambda, ")")), ylab="Fraction rejected")
plot(lam.LRT ~ alpha, data=power, typ="l", xlim=c(0,50), ylim=c(0, 1), xlab=expression(paste("Phylogenetic signal (", alpha, ")")), ylab="Fraction rejected")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lines(OU.LRT ~ alpha, data=power, col=2)
lines(lam.boot0 ~ alpha, data=power, col=3)
lines(OU.boot0 ~ alpha, data=power, col=4)
lines(K.perm ~ alpha, data=power, col=5)

lines(c(0,50), c(.05,.05), lty=2)
legend(.0,1,c("lam.LRT","OU.LRT", "lam.boot(H0)", "OU.boot(H0)", "K.perm"), col=1:5, lty=1)


# 3. The parametric bootstrap of H0 for phylogenetic signal (subsection 3.5.3) uses the log likelihood ratio (LLR) as the test statistic. It is also possible to use the phylogenetic signal parameter (λ or α) as the test statistic. This would involve the following: (i) Fit the model to the data and calculate the value of the phylogenetic signal parameter (λ or α). (ii) Refit the model with no phylogenetic signal and use the resulting parameter values to simulate a large number of datasets. (iii) Refit the model including the phylogenetic signal parameter (λ or α) for each dataset and collect these values. (iv) The resulting distribution of λ or α estimated from the simulated datasets approximates the distribution of the estimator of λ or α, allowing P-values to be calculated. Perform this bootstrap and compare the results to the bootstrap of H0 using the LLR.

#TK: I tried writing my own code so this might just all be terribly wrong.

n <- 30
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
plot(phy)

# Simulate data with Pagel's lambda transform.
lam.sig <- .6
lam.nosig <- 0

phy.sig <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.sig))$tree
phy.nosig <- transf.branch.lengths(phy=phy, model="lambda", parameters=list(lambda = lam.nosig))$tree

Y.sig <- rTraitCont(phy.sig, model = "BM", sigma = 1)
Y.nosig <- rTraitCont(phy.nosig, model = "BM", sigma = 1)

mod.sig <- phylolm(Y.sig ~ 1, phy=phy, model = "lambda")
mod.nosig <- phylolm(Y.nosig ~ 1, phy=phy, model = "lambda")

lm.nosig <- lm(Y.nosig ~ 1)

# Perform the bootstrap using the parameters estimated under H0. The function simulate() is applied to the model z.0. The simulated datasets are fit with Pagel's lambda
nboot <- 2000
Y.sim <- Y.nosig
boot.lam <- data.frame(lambda=rep(NA,nboot))
for(i in 1:nboot){
  Y.sim <- simulate(lm.nosig)[[1]]
  names(Y.sim) <- names(Y.nosig)
  lam.sim <- phylolm(Y.sim ~ 1, phy=phy, model = "lambda")
  boot.lam$lambda[i] <- lam.sim$optpar

}

par(mfrow=c(1,1), mai=c(1,1,.4,.4))
hist(boot.lam$lambda)
lines(c(mod.sig$optpar, mod.sig$optpar), c(0,1500), col="red")
#TK: How do I calculate P-values? Uncertain if this is correct.
pvalue <- mean(boot.lam$lambda > mod.sig$optpar)
text(0.4,800,paste("P =", pvalue), cex=1.5)

# 4. When investigating mixed models (chapters 1 and 2), we focused on testing hypotheses concerning the fixed effects (the slope). Knowing what you know about testing for phylogenetic signal (sections 3.4 and 3.5), what methods could you use to test hypotheses regarding the random effects in a mixed model? You don't need to do this in R: I'm only asking for ideas. But you will impress me if you can take the models from chapter 1 and perform a test of a null hypothesis concerning a random effect. 

#Within the mixed model, you can use Likelihood ratio tests and Wald tests to test hypotheses regarding the random effects. You can also make two models - one where yoy specify that there is no interaction term in your simulations and one where there is an interaction term and compare the two using anova(). OR you can try the following? 

#Stuff from CRAN - attempt at looking at random effects 
library(pez)
library(ape)

nspp <- 15
nsite <- 10
env <- 1:nsite
env <- as.numeric(scale(env))

phy <- rcoal(n = nspp)
Vphy <- vcv(phy)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

iD <- t(chol(Vphy))
intercept <- iD %*% rnorm(nspp)
slope <- iD %*% rnorm(nspp)

prob <- rep(intercept, each = nsite)
prob <- prob + rep(slope, each = nsite) *
  rep(env, nspp)
prob <- prob + rnorm(nspp * nsite)
pres <- rbinom(length(prob), size = 1, prob = exp(prob)/(1 + exp(prob)))

site <- factor(rep(1:nsite, nspp))
species <- factor(rep(1:nspp, each = nsite))
env <- rep(env, nspp)

r.intercept.spp.indep <- list(1, sp = species,
                              covar = diag(nspp))
r.intercept.spp.phy <- list(1, sp = species,
                            covar = Vphy)
r.slope.spp.indep <- list(env, sp = species,
                          covar = diag(nspp))
r.slope.spp.phy <- list(env, sp = species,
                        covar = Vphy)
r.site <- list(1, site = site, covar = diag(nsite))
rnd.effects <- list(r.intercept.spp.indep,
                    r.intercept.spp.phy, r.slope.spp.indep,
                    r.slope.spp.phy, r.site)
model <- communityPGLMM(pres ~ env, family = "binomial",
                        sp = species, site = site, random.effects = rnd.effects,
                        REML = TRUE, verbose = FALSE)

communityPGLMM.binary.LRT(model, re.number = 1)
communityPGLMM.binary.LRT(model, re.number = 2)

#Now you can specify the random effect 


