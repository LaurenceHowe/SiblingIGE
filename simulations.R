#Load required packages
require(MASS)

# Parameters
#Values of these parameters do not affect the non-singleton, singleton comparison model.
n <- 100000 # Sample size
d <- 0.2 # Direct genetic effect
p <- 0.1 # Indirect genetic effect of mean parental genotype
c <- 0 # Confounding between genetic variant and phenotype - e.g. population stratification, selection bias

#Values for these parameters will affect the non-singleton, singleton comparison model.
s <- 0.1 # Sibling indirect genetic effect
d_diff <- 0 # Difference between direct genetic effects for singletons/non-singletons
p_diff <- 0 # Difference between parental indirect genetic effects for singletons/non-singletons
c_diff <- 0 # Difference between confounding for singletons/non-singletons

#Simulate genotypes under normal distribution
#Model 1) Random mating - Correlation of 0.5 between sibling genotypes

model1_rm <- function(n, d, p, c, s, d_diff, p_diff, c_diff) 
{

output <- NULL

#Number of simulations
for (i in 1:100) {
  
#Simulate singleton and parental genotypes under random mating
mean_s <- c(0, 0, 0) #Parent 1, Parent 2, index individual                        
vcov_s_rm <- matrix(c(1, 0, 0.5, 
                       0, 1, 0.5,
                       0.5, 0.5, 1),                
                       ncol = 3)
data_s <- data.frame(mvrnorm(n = n, mu = mean_s, Sigma = vcov_s_rm))
names(data_s) <- c("G_P1", "G_P2", "G_I")

#Simulate singleton phenotypes
data_s$Phen_I <- (d * data_s$G_I) +
                 (0.5* p * (data_s$G_P1 + data_s$G_P2)) +
                 (c * data_s$G_I) + 
                 rnorm(1000, 0, 1)

#Extract OLS regression estimate of individual's genotype on their phenotype
ols_s <- lm(Phen_I ~ G_I, data = data_s)
ols_s_beta <- summary(ols_s)$coefficients[2, 1]


#Simulate non-singleton (sibling pair) genotypes
mean_ns <- c(0, 0, 0, 0) #Parent 1, Parent 2, Sibling 1, Sibling 2                      
vcov_ns_rm <- matrix(c(1, 0, 0.5, 0.5, 
                      0, 1, 0.5, 0.5,
                      0.5, 0.5, 1, 0.5,
                      0.5, 0.5, 0.5, 1),                
                    ncol = 4)

data_ns <- data.frame(mvrnorm(n = n, mu = mean_ns, Sigma = vcov_ns_rm))
names(data_ns) <- c("G_P1", "G_P2", "G_S1", "G_S2")
  
#Simulate Sibling 1 phenotypes with IGE from Sibling 2
#Add difference terms (e.g. d_diff) to simulate differences between singletons and non-singletons.
data_ns$Phen_S1 <- ((d - d_diff) * data_ns$G_S1) + 
                        (0.5 * (p - p_diff) * (data_ns$G_P1 + data_ns$G_P2)) +
                        (s * data_ns$G_S2) + 
                        ((c - c_diff) * data_ns$G_S1) + 
                        rnorm(1000, 0, 1)

#Extract OLS regression estimate of individual's genotype on their phenotype
ols_ns <- lm(Phen_S1 ~ G_S1, data = data_ns)
ols_ns_beta <- summary(ols_ns)$coefficients[2, 1]
  
  
  
#Extract the phenotypic variances in singleton and non-singleton samples
var_s <- var(data_s$Phen_I)
var_ns <- var(data_ns$Phen_S1)
  
temp <- cbind(ols_s_beta, ols_ns_beta, var_s, var_ns)
  
output <- rbind(output, temp)

}

mean_ols_beta_s <- mean(output[, 1])
mean_ols_beta_ns <- mean(output[, 2])
mean_var_s <- mean(output[, 3])
mean_var_ns <- mean(output[, 4])

output_f <- cbind(mean_ols_beta_s, mean_ols_beta_ns, mean_var_s, mean_var_ns)
return(output_f)
}

model1_rm(n, d, p, c, s, d_diff, p_diff, c_diff)


model1_rm(n, d, p, c, s = 0, d_diff = -0.05, p_diff, c_diff)

#Model 2) Assortative mating - Correlation of 0.6 between sibling genotypes

model2_am <- function(n, d, p, c, s, d_diff, p_diff, c_diff) 
{
  
  output <- NULL
  
  #Number of simulations
  for (i in 1:100) {
    
    #Simulate singleton and parental genotypes under random mating
    mean_s <- c(0, 0, 0) #Parent 1, Parent 2, index individual                        
    vcov_s_am <- matrix(c(1, 0.2, 0.6, 
                          0.2, 1, 0.6,
                          0.6, 0.6, 1),                
                        ncol = 3)
    data_s <- data.frame(mvrnorm(n = n, mu = mean_s, Sigma = vcov_s_am))
    names(data_s) <- c("G_P1", "G_P2", "G_I")
    
    #Simulate singleton phenotypes
    data_s$Phen_I <- (d * data_s$G_I) +
      (0.5* p * (data_s$G_P1 + data_s$G_P2)) +
      (c * data_s$G_I) + 
      rnorm(1000, 0, 1)
    
    #Extract OLS regression estimate of individual's genotype on their phenotype
    ols_s <- lm(Phen_I ~ G_I, data = data_s)
    ols_s_beta <- summary(ols_s)$coefficients[2, 1]
    
    
    #Simulate non-singleton (sibling pair) genotypes
    mean_ns <- c(0, 0, 0, 0) #Parent 1, Parent 2, Sibling 1, Sibling 2                      
    vcov_ns_am <- matrix(c(1, 0.2, 0.6, 0.6, 
                           0.2, 1, 0.6, 0.6,
                           0.6, 0.6, 1, 0.6,
                           0.6, 0.6, 0.6, 1),                
                         ncol = 4)
    
    data_ns <- data.frame(mvrnorm(n = n, mu = mean_ns, Sigma = vcov_ns_am))
    names(data_ns) <- c("G_P1", "G_P2", "G_S1", "G_S2")
    
    #Simulate Sibling 1 phenotypes with IGE from Sibling 2
    #Add difference terms (e.g. d_diff) to simulate differences between singletons and non-singletons.
    data_ns$Phen_S1 <- ((d - d_diff) * data_ns$G_S1) + 
      (0.5 * (p - p_diff) * (data_ns$G_P1 + data_ns$G_P2)) +
      (s * data_ns$G_S2) + 
      ((c - c_diff) * data_ns$G_S1) + 
      rnorm(1000, 0, 1)
    
    #Extract OLS regression estimate of individual's genotype on their phenotype
    ols_ns <- lm(Phen_S1 ~ G_S1, data = data_ns)
    ols_ns_beta <- summary(ols_ns)$coefficients[2, 1]
    
    
    
    #Extract the phenotypic variances in singleton and non-singleton samples
    var_s <- var(data_s$Phen_I)
    var_ns <- var(data_ns$Phen_S1)
    
    temp <- cbind(ols_s_beta, ols_ns_beta, var_s, var_ns)
    
    output <- rbind(output, temp)
    
  }
  
  mean_ols_beta_s <- mean(output[, 1])
  mean_ols_beta_ns <- mean(output[, 2])
  mean_var_s <- mean(output[, 3])
  mean_var_ns <- mean(output[, 4])
  
  output_f <- cbind(mean_ols_beta_s, mean_ols_beta_ns, mean_var_s, mean_var_ns)
  return(output_f)
}

model2_am(n, d, p, c, s, d_diff, p_diff, c_diff)
