# ***README***
# The following script is used to fit the marked spatial point pattern data to the model 
# proposed in the submitted manuscript titled "A Bayesian Hidden Potts Mixture Model for 
# Analyzing Lung Cancer Pathology Images."

# Before running the following code, please first load the data using "data_loader.R". The
# necessary inputs should be x, y, and z, which denote the x, y coordinates, and marks of
# each point. Note that the code is only designed for the case of Q <= 4.
# ***END***

# Load library
library(Rcpp)
library(coda)

setwd('C:/Users/Diana/Documents/Bayes/Project_Lungs/plan B/BayesHiddenPottsMixture-master')

# Load function
source("functions.R");
Rcpp::sourceCpp('functions_x.cpp');
ef_arr<- read.table("ef_params_priors.txt", header = TRUE)

setwd('C:/Users/Diana/Documents/Bayes/Project_Lungs/plan B/BayesHiddenPottsMixture-master/100381/cell_info_txt')


patient_ID<-100381 #name of folder
fl<-list.files();
#because the fl stores the files as 1 10 2 3 ....
#needs rearranging
fl[11]<-fl[2]
fl<-fl[-(2)]


num_to_proc<-c(1,4,5,6,8)#order number of 5 most mixed images
ef_percents<-c(30,30,30,30,30)#% of AOI possile values 1,2,5,10,15,20,25,30,35,40


# ========================================================================================
# ========================================================================================
# Load algorithm setting
# ========================================================================================
# ========================================================================================
Q <- 3;           # Number of marks
W <- 50;          # Lattice width
L <- 50;          # Lattice length


# ========================================================================================
# ========================================================================================
# Load hyperparameters
# ========================================================================================
# ========================================================================================
# Prior for theta
mu <- 0.5;
sigma <- 1;

# Prior for theta_0
mu_0 <- -0.5;
sigma_0 <- 1;



# Prior for d
a_d <- 0.001;
b_d <- 0.001;


# ========================================================================================
# ========================================================================================
# Load initial configuration
# ========================================================================================
# ========================================================================================
d_start <- runif(1, 1, 10);
theta_start <- rnorm(Q*(Q - 1)/2, mu, sigma);
Theta_start <- array2matrix_r(theta_start); #from the functions.r
theta_0_start <- rnorm(Q*(Q - 1)/2, mu_0, sigma_0);
Theta_0_start <- array2matrix_r(theta_0_start); #from the functions.r
Delta_start <- matrix(0, nrow = W, ncol = L);
W_s <- sample(2:(W - 6), 1);
W_s <- c(W_s, W_s + 5);
L_s <- sample(2:(L - 6), 1);
L_s <- c(L_s, L_s + 5);
Delta_start[W_s[1]:W_s[2], L_s[1]:L_s[2]] <- 1;
P_start <- matrix(sample(1:Q, W*L, replace = TRUE), ncol = L, nrow = W);






# ========================================================================================
# ========================================================================================
# CYCLE OF 5 IMAGES
# ========================================================================================
# ========================================================================================
for(j in 1:5){

filename1<-fl[num_to_proc[j]];  
data <- read.table(filename1, header = TRUE)
# Prior for delta
e <- ef_arr[ef_percents[j],1];
f <- ef_arr[ef_percents[j],2]; #prior for delta
  
# ========================================================================================
# Preprocess data (DO NOT EDIT THE CODE IN THIS BLOCK)
# ========================================================================================

HP <- data;
x_c <- seq(floor(min(HP$x)), floor(max(HP$x)) + 1, length.out = L);#sequence generation
y_c <- seq(floor(min(HP$y)), floor(max(HP$y)) + 1, length.out = W);#seq(from = 1, to = 1, by = ((to - from)/(length.out - 1))
HP$l <- NA;
  HP$w <- NA;
  for (i in 1:dim(HP)[1]) {
    temp <- 1;
    while(HP$x[i] > x_c[temp]) {
      temp <- temp + 1;
    }
    HP$l[i] <- temp - 1;
    temp <- 1;
    while(HP$y[i] > y_c[temp]) {
      temp <- temp + 1;
    }
    HP$w[i] <- temp - 1;
  }
  cell_info <- as.matrix(HP[, c(4, 5, 3)]);
  
  data_potts <- data.frame(cbind(rep(1:W, each = L), rep(1:L, times = W)));
  names(data_potts) <- c("w", "l");
  data_potts$nb <- NA;
  max_count <- rep(0, dim(data_potts)[1]);
  for (i in 1:dim(data_potts)[1]) {
    temp <- NULL;
    temp <- c(temp, which((HP$w == data_potts$w[i] & HP$l == data_potts$l[i])));
    temp <- c(temp, which((HP$w == data_potts$w[i] - 1 & HP$l == data_potts$l[i])));
    temp <- c(temp, which((HP$w == data_potts$w[i] & HP$l == data_potts$l[i] - 1)));
    temp <- c(temp, which((HP$w == data_potts$w[i] - 1 & HP$l == data_potts$l[i] - 1)));
    if (length(temp) > 0) {
      max_count[i] <- length(temp);
      data_potts$nb[i] <- paste0(temp, collapse = "", sep = ",");
    }
  }
  
  potts_info <- matrix(0, nrow = dim(data_potts)[1], ncol = 2 + Q);
  potts_info[, 1] <- data_potts$w;
  potts_info[, 2] <- data_potts$l;
  for (i in 1:dim(data_potts)[1]) {
    potts_info[i, 3:(2 + Q)] <- table(factor(HP$z[which(HP$l == data_potts$l[i] & HP$w == data_potts$w[i])], levels = 1:Q));
  }
  colnames(potts_info) <- c("w", "l", paste0("class ", 1:Q));
  
    
  
  
  
  



iter <- 20000;     # Number of iterations
burn <- iter/2;   # Number of burn in

start_time <- proc.time();
Y <- mcmc_hidden_2(W, L, Q, potts_info, Theta_start, Theta_0_start, Delta_start, P_start, 
                   d_start, iter, burn, mu, mu_0, sigma, sigma_0, e, f, a_d, b_d);
end_time <- proc.time();
time <- end_time - start_time;

saveRDS(Y, paste(filename1, ".rds",sep=""))
#---------------------------1000--------------------------------------------

c_iter<-1000;#iter for coda
c_burn <- c_iter/2;
c_thin<- 1;

for_coda_1000<-cbind(Y$theta[1:1000,],Y$theta_0,Y$d[1:1000,])

Y_mc_1000<- mcmc(data =for_coda_1000, start = c_burn + 1, end = c_iter, thin = c_thin)
s<-summary(Y_mc_1000)
theta_1_1000<-s$statistics[1,1]
theta_2_1000<-s$statistics[2,1]
theta_3_1000<-s$statistics[3,1]
theta_01_1000<-s$statistics[4,1]
theta_02_1000<-s$statistics[5,1]
theta_03_1000<-s$statistics[6,1]
d_1000<-s$statistics[7,1]

ESS_1000<-effectiveSize(Y_mc_1000)

Y_for_test_1000<- mcmc(data =for_coda_1000, start =  1, end =c_iter, thin = 1)

Hd_1000<-heidel.diag(Y_for_test_1000, eps=0.1, pvalue=0.05) #2 additional variables for output

#rm(for_coda_1000,Y_mc_1000)

#------------------------5000----------------------------------------------

c_iter<-5000;#iter for coda
c_burn <- c_iter/2;
c_thin<-1

for_coda_5000<-cbind(Y$theta[1:5000,],Y$theta_0[1:5000,],Y$d[1:5000,])

Y_mc_5000<- mcmc(data =for_coda_5000, start = c_burn + 1, end = c_iter, thin = c_thin)
s<-summary(Y_mc_5000)
theta_1_5000<-s$statistics[1,1]
theta_2_5000<-s$statistics[2,1]
theta_3_5000<-s$statistics[3,1]
theta_01_5000<-s$statistics[4,1]
theta_02_5000<-s$statistics[5,1]
theta_03_5000<-s$statistics[6,1]
d_5000<-s$statistics[7,1]


ESS_5000<-effectiveSize(Y_mc_5000)

Y_for_test_5000<- mcmc(data =for_coda_5000, start =  1, end =c_iter, thin = 1)

Hd_5000<-heidel.diag(Y_for_test_5000, eps=0.1, pvalue=0.05) #2 additional variables

#rm(for_coda_5000,Y_mc_5000)


#-----------------------10000-------------------------------------------------
c_iter<-10000;#iter for coda
c_burn <- c_iter/2;
c_thin<-1

for_coda_10000<-cbind(Y$theta[1:10000,],Y$theta_0[1:10000,],Y$d[1:10000,])

Y_mc_10000<- mcmc(data =for_coda_10000, start = c_burn + 1, end = c_iter, thin = c_thin)
s<-summary(Y_mc_10000)
theta_1_10000<-s$statistics[1,1]
theta_2_10000<-s$statistics[2,1]
theta_3_10000<-s$statistics[3,1]
theta_01_10000<-s$statistics[4,1]
theta_02_10000<-s$statistics[5,1]
theta_03_10000<-s$statistics[6,1]
d_10000<-s$statistics[7,1]

ESS_10000<-effectiveSize(Y_mc_10000)

Hd_10000<-heidel.diag(Y_for_test_10000, eps=0.1, pvalue=0.05) #2 additional variables
#rm(for_coda_10000,Y_mc_10000)

#------------------------20000------------------------------------------------
c_iter<-20000;#iter for coda
c_burn <- c_iter/2;
c_thin<-1

for_coda_all_20000<-cbind(Y$theta[1:20000,],Y$theta_0[1:20000,],Y$d[1:20000,])

Y_mc_20000<- mcmc(data =for_coda_all_20000, start = c_burn + 1, end = c_iter, thin = c_thin)
s<-summary(Y_mc_20000)
theta_1_20000<-s$statistics[1,1]
theta_2_20000<-s$statistics[2,1]
theta_3_20000<-s$statistics[3,1]
theta_01_20000<-s$statistics[4,1]
theta_02_20000<-s$statistics[5,1]
theta_03_20000<-s$statistics[6,1]
d_20000<-s$statistics[7,1]

ESS_20000<-effectiveSize(Y_mc_20000)

Hd_20000<-heidel.diag(Y_for_test_20000, eps=0.1, pvalue=0.05) #2 additional variables
#rm(for_coda_20000,Y_mc_20000)

}




# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm (DO NOT EDIT THE CODE IN THIS BLOCK)
# ========================================================================================
# ========================================================================================
#start_time <- proc.time();
#Y <- mcmc_hidden_2(W, L, Q, potts_info, Theta_start, Theta_0_start, Delta_start, P_start, 
#                   d_start, iter, burn, mu, mu_0, sigma, sigma_0, e, f, a_d, b_d);
#end_time <- proc.time();
#time <- end_time - start_time;
# The MCMC outputs are all stored in Y
# $P_1:       the number of iterations that the spins equal to mark 1
# $P_2:       the number of iterations that the spins equal to mark 2
# $P_3:       the number of iterations that the spins equal to mark 3
# $P_4:       the number of iterations that the spins equal to mark 4
# $theta:     the value of theta for each iteration
# $d:         the value of d for each iteration
# $theta_0:   the value of theta_0 for each iteration
# $Delta_ppi: the marginal posterior probability of inclusion (PPI) of Delta
# $Delta_sum: the number of spins in the AOI for each iteration
# $accept:    the acceptance rate for the parameter theta, theta_0, Delta, and d


# Obtain the MP estimate for the spins P
P_MP <- matrix(0, nrow = W, ncol = L);
for (h in 1:W) {
  for (l in 1:L) {
    if (max(Y$P_1[h, l], Y$P_2[h, l], Y$P_3[h, l]) != 0) {
      P_MP[h, l] <- which.max(c(Y$P_1[h, l], Y$P_2[h, l], Y$P_3[h, l]));
    }
  }
}


# ========================================================================================
# ========================================================================================
# Summarize result
# ========================================================================================
# ========================================================================================
print(paste0("Runtime = ", round(time[3], 1), "s"));
# print(paste0("Acceptance rate (theta): ", round(Y$accept[1], 3)));
# print(paste0("Acceptance rate (theta_0): ", round(Y$accept[2], 3)));
# print(paste0("Acceptance rate (Delta): ", round(Y$accept[3], 3)));
# print(paste0("Acceptance rate (d): ", round(Y$accept[4], 3)));
print(paste0("Estimated parameters (theta) = ", paste0(round(colMeans(Y$theta[burn:iter,]), 3), collapse = "", sep = ", ")));
print(paste0("Estimated parameters (theta_0) = ", paste0(round(colMeans(Y$theta_0[burn:iter,]), 3), collapse = "", sep = ", ")));
print(paste0("Estimated parameters (d) = ", paste0(round(mean(Y$d[burn:iter]), 3), collapse = "", sep = ", ")));

# Plot the marginal posterior of probabilities of Delta
par(pty = "s", mar = c(4.5, 2.5, 1, 1));
plot(rep(1:L, each = W), rep(1:W, times = L), asp = 1, cex.main = 2, 
     col = gray(rev(0:100/100))[as.numeric(cut(Y$Delta_ppi, breaks = 100))], pch = 15, 
     xlab = "L", ylab = "W", main = "", cex = 1, cex.lab = 1.5, cex.axis = 1.5);

# Plot the marginal probabilities of P
par(pty = "s", mar = c(4.5, 2.5, 1, 1));
plot(rep(1:L, each = W), rep(1:W, times = L), col = as.numeric(P_MP), asp = 1, cex.main = 2, 
     pch = 15, xlab = "L", ylab = "W", main = "", cex = 1, cex.lab = 1.5, cex.axis = 1.5);


