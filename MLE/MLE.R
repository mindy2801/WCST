#Code for MLE
#Wanted parameter=theta(proportion)
n <- 0 # Sample size
t <- c() # Time interval
y <- c() # Observed proportion
x <- n*y # Number of correct responses

#Example (number of subjects=6)
n <- 100
t <- c(1, 3, 6, 9, 12, 18)
y <- c(.94, .77, .40, .26, .24, .16)
x

  
lik_fun <- function(n, y, w){
  factorial(n)/(factorial(y)*factorial(n-y))*w^y*(1-w)^(n-y)
}

#Power model
power_x <- w*t