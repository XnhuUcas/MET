library(mvtnorm)

## Computing the p-value of the MET based on the case-control design;
## r: the number of cases; s: the number of controls; grr1,grr2: the genotype relative risks;
## k: the disease prevalence; alpha: the nominal significance level; c: the threshold for HWDTT;
## maf: minor allele frequency; 

MET <- function(r, s, grr1, grr2, k, alpha, c, maf)
{
  # population
  n <- r+s
  
  # scores of the recessive model, the additive model, and the dominant model  
  score <- c(0,1/2,1)
  
  # penetrance
  f_0 <- k/((1-maf)^2+grr1*2*maf*(1-maf)+grr2*maf^2)
  
  # genotype frequencies in cases
  p_0 <- f_0*(1-maf)^2/k
  p_1 <- grr1*f_0*2*maf*(1-maf)/k
  p_2 <- grr2*f_0*maf^2/k
  
  # genotype frequencies in controls
  q_0 <- (1-f_0)*(1-maf)^2/(1-k)
  q_1 <- (1-grr1*f_0)*2*maf*(1-maf)/(1-k)
  q_2 <- (1-grr2*f_0)*maf^2/(1-k)
  
  met.pvalue <- c()
    
  # observed genotype counts in cases and controls
  case <- sample(c(0,1,2),r,replace=T,prob=c(p_0,p_1,p_2))
  control <- sample(c(0,1,2),s,replace=T,prob=c(q_0,q_1,q_2))
    
  r_0 <- length(which(case==0))
  r_1 <- length(which(case==1))
  r_2 <- length(which(case==2))
  s_0 <- length(which(control==0))
  s_1 <- length(which(control==1))
  s_2 <- length(which(control==2)) 
  n_0 <- r_0+s_0
  n_1 <- r_1+s_1
  n_2 <- r_2+s_2
    
  if(n_2!=0){
  
    # estimators of genotype frequencies under the null hypothesis
    ep_1 <- eq_1 <- n_1/n
    ep_2 <- eq_2 <- n_2/n
    
    # estimators of the variances and covariances under three genetic models
    cvar_rec <- (1/r)*(ep_2-3*ep_2^2+2*ep_2^3+2*ep_1*ep_2^2+(1/2)*ep_1^2*ep_2-ep_1*ep_2)
    cvar_add <- (1/r)*(ep_2-2*ep_1*ep_2-3*ep_2^2+2*ep_2^3+3*ep_1*ep_2^2+(3/2)*ep_1^2*ep_2-(1/4)*ep_1^2+(1/4)*ep_1^3)
    cvar_dom <- (1/r)*(ep_2-3*ep_2^2-3*ep_1*ep_2+2*ep_2^3+4*ep_1*ep_2^2+(5/2)*ep_1^2*ep_2-(1/2)*ep_1^2+(1/2)*ep_1^3)
    var_rec <- (ep_2-ep_2^2)/r+(eq_2-eq_2^2)/s
    var_add <- ((ep_2+(1/4)*ep_1)-(ep_2+(1/2)*ep_1)^2)/r+((eq_2+(1/4)*eq_1)-(eq_2+(1/2)*eq_1)^2)/s
    var_dom <- ((ep_2+ep_1)-(ep_2+ep_1)^2)/r+((eq_2+eq_1)-(eq_2+eq_1)^2)/s
    vardel <- (1/r)*(ep_2-5*ep_2^2+8*ep_2^3-4*ep_2^4+(1/4)*ep_1^3-(1/4)*ep_1^4-2*ep_1*ep_2+3*ep_1^2*ep_2+9*ep_1*ep_2^2-6*ep_1^2*ep_2^2-8*ep_1*ep_2^3-2*ep_1^3*ep_2)
    
    # CATT under three genetic models
    Z_rec <- n^(1/2)*(score[1]*(s*r_1-r*s_1)+(s*r_2-r*s_2))/((r*s*(n*(score[1]^2*n_1+n_2)-(score[1]*n_1+n_2)^2))^(1/2))
    Z_add <- n^(1/2)*(score[2]*(s*r_1-r*s_1)+(s*r_2-r*s_2))/((r*s*(n*(score[2]^2*n_1+n_2)-(score[2]*n_1+n_2)^2))^(1/2))
    Z_dom <- n^(1/2)*(score[3]*(s*r_1-r*s_1)+(s*r_2-r*s_2))/((r*s*(n*(score[3]^2*n_1+n_2)-(score[3]*n_1+n_2)^2))^(1/2))
      
    # HWDTT in only case
    delta_p <- r_2/r-(r_2/r+r_1/(2*r))^2
    Z_C <- delta_p/sqrt(vardel)
            
    # expectation of the MET
    delta_pt <- p_2-(p_2+p_1/2)^2
    mu_deltap <- delta_pt/sqrt(vardel)
    mean_met <- c(0,mu_deltap)
    
    # variance and covariance matrix of the MET    
    ecvar <- c(cvar_rec/sqrt(var_rec*vardel),cvar_add/sqrt(var_add*vardel),cvar_dom/sqrt(var_dom*vardel))
    covr_met <- matrix(c(1,ecvar[1],ecvar[1],1),ncol=2)
    cova_met <- matrix(c(1,ecvar[2],ecvar[2],1),ncol=2)
    covd_met <- matrix(c(1,ecvar[3],ecvar[3],1),ncol=2)
	  
	# p-value of the MET
    if(Z_C>c){
	    lower1 <- c(abs(Z_rec),c)
        upper1 <- rep(Inf,2)
        lower2 <- c(abs(Z_rec),-c)
        upper2 <- c(Inf,c)
        lower3 <- c(abs(Z_rec),-Inf)
        upper3 <- c(Inf,-c)
        lower4 <- c(-Inf,c)
        upper4 <- c(-abs(Z_rec),Inf)
        lower5 <- c(-Inf,-c)
        upper5 <- c(-abs(Z_rec),c)
        lower6 <- c(-Inf,-Inf)
        upper6 <- c(-abs(Z_rec),-c)
        
        met.pvalue <- pmvnorm(lower1,upper1,mean_met,sigma=covr_met)[[1]]+pmvnorm(lower2,upper2,mean_met,sigma=cova_met)[[1]]+pmvnorm(lower3,upper3,mean_met,sigma=covd_met)[[1]]+pmvnorm(lower4,upper4,mean_met,sigma=covr_met)[[1]]+pmvnorm(lower5,upper5,mean_met,sigma=cova_met)[[1]]+pmvnorm(lower6,upper6,mean_met,sigma=covd_met)[[1]]
	}else if(Z_C<(-c)){
        lower1 <- c(abs(Z_dom),c)
        upper1 <- rep(Inf,2)
        lower2 <- c(abs(Z_dom),-c)
        upper2 <- c(Inf,c)
        lower3 <- c(abs(Z_dom),-Inf)
        upper3 <- c(Inf,-c)
        lower4 <- c(-Inf,c)
        upper4 <- c(-abs(Z_dom),Inf)
        lower5 <- c(-Inf,-c)
        upper5 <- c(-abs(Z_dom),c)
        lower6 <- c(-Inf,-Inf)
        upper6 <- c(-abs(Z_dom),-c)
        
        met.pvalue <- pmvnorm(lower1,upper1,mean_met,sigma=covr_met)[[1]]+pmvnorm(lower2,upper2,mean_met,sigma=cova_met)[[1]]+pmvnorm(lower3,upper3,mean_met,sigma=covd_met)[[1]]+pmvnorm(lower4,upper4,mean_met,sigma=covr_met)[[1]]+pmvnorm(lower5,upper5,mean_met,sigma=cova_met)[[1]]+pmvnorm(lower6,upper6,mean_met,sigma=covd_met)[[1]]
      }else{
        lower1 <- c(abs(Z_add),c)
        upper1 <- rep(Inf,2)
        lower2 <- c(abs(Z_add),-c)
        upper2 <- c(Inf,c)
        lower3 <- c(abs(Z_add),-Inf)
        upper3 <- c(Inf,-c)
        lower4 <- c(-Inf,c)
        upper4 <- c(-abs(Z_add),Inf)
        lower5 <- c(-Inf,-c)
        upper5 <- c(-abs(Z_add),c)
        lower6 <- c(-Inf,-Inf)
        upper6 <- c(-abs(Z_add),-c)
        
        met.pvalue <- pmvnorm(lower1,upper1,mean_met,sigma=covr_met)[[1]]+pmvnorm(lower2,upper2,mean_met,sigma=cova_met)[[1]]+pmvnorm(lower3,upper3,mean_met,sigma=covd_met)[[1]]+pmvnorm(lower4,upper4,mean_met,sigma=covr_met)[[1]]+pmvnorm(lower5,upper5,mean_met,sigma=cova_met)[[1]]+pmvnorm(lower6,upper6,mean_met,sigma=covd_met)[[1]]
      }
      
    }

list(met.pvalue)
  
}

