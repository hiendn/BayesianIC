library(pracma)  # Load the 'pracma' library, which is used for numerical integration.

STORE <- matrix(NA, 10 * 5, 3)  # Initialize a matrix to store results, with 50 rows and 3 columns.
COUNT <- 0  # Initialize a counter to track the row index in STORE.

# Outer loop over nn (ranging from 1 to 5), corresponding to powers of 10 for n.
for (nn in 1:5) {
  # Inner loop to perform 10 repetitions for each value of n.
  for (rep in 1:10) {
    COUNT <- COUNT + 1  # Increment the counter to track the current iteration.
    n <- 10^nn  # Define n as 10 raised to the power of nn.
    mu <- 0  # Set the mean parameter for the prior.
    beta <- 1 / (n * log(n))  # Define beta, inversely proportional to n * log(n).
    
    X <- rnorm(n, 1, 1)  # Generate n random samples from a normal distribution with mean=1 and sd=1.
    MEAN <- mean(X)  # Compute the sample mean of X.
    
    # Compute the posterior mean and variance using conjugate prior properties.
    P_mean <- (n * beta * MEAN + mu) / (n * beta + 1)  
    P_var <- 1 / (n * beta + 1)
    
    # Define the log-likelihood function for the normal distribution.
    LogLike <- function(x) {
      sum(dnorm(X, x, 1, log = TRUE))  # Sum of log-likelihood values for X given mean x.
    }
    
    # Define the integrand function for WBIC calculation.
    Integrand <- function(x) {
      LogLike(x) * dnorm(x, P_mean, sqrt(P_var))  # Weighted log-likelihood with posterior density.
    }
    
    # Store the results: n (sample size), repetition number, and WBIC value.
    STORE[COUNT, 1] <- n
    STORE[COUNT, 2] <- rep
    STORE[COUNT, 3] <- -(2 / n) * integral(
      Integrand,
      P_mean - 15 * sqrt(P_var),
      P_mean + 15 * sqrt(P_var),
      method = 'Simpson'
    )  # Perform numerical integration using Simpson's method.
    
    # Print the current results and a reference value for comparison.
    print(STORE[COUNT, ])
    print(log(2 * pi) + 1)
  }
}

# Plot the WBIC values against log10(n).
plot(
  STORE[, 3] ~ log(STORE[, 1], 10),
  ylab = 'WBIC',
  xlab = 'log10(n)',
  main = '(f)'
)
grid()  # Add a grid to the plot.
abline(h = log(2 * pi) + 1, col = 'red', lwd = 2)


## A = 1/log(n)
## B = 1/loglog(n)
## C = 1
## D = 1/sqrt(n)
## E = 1/n
## F = 1/nlog(n)
