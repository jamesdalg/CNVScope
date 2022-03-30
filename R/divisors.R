#' List of Divisors
#'
#' Generates a list of divisors of an integer number.
#' Identical to the same function within the numbers package.
#' The code has been modified from the numbers package,
#'  following GPL 3.0 guidelines on 3/30/2022, section 5.
#'  Reference for GPL v3.0 LICENSE:
#'  https://www.gnu.org/licenses/gpl-3.0.en.html.
#' 
#'
#' @name divisors
#' @keywords divisors numbers
#' @param n   an integer whose divisors will be generated.
#' @return Returns a vector integers.
#' @examples 
#' divisors(1)          # 1
#' divisors(2)          # 1 2
#' divisors(3)          # 1 2 3
#' divisors(2^5)        # 1  2  4  8 16 32
#' divisors(1000)       # 1  2  4  5  8 10 ... 100 125 200 250 500 1000
#' divisors(1001)       # 1  7 11 13 77 91 143 1001
#' @seealso [numbers::divisors()]
#' @export

divisors <- function(n) {
  if (n != floor(n) || n <= 0)
    stop("Argument 'n' must be a nonnegative integer.")
  if (n == 1) {
    return(1)
  } else if (n <= 1000) {
    return( (1:n)[(n %% 1:n) == 0] )
  } else {
    pfs <- rle(primeFactors(n))
    pfs_len <- pfs$length
    pfs_val <- pfs$values
    
    m <- length(pfs_len)
    D <- pfs_val[1]^c(0:pfs_len[1])
    if (m == 1) return(D)
    
    for (k in 2:m) {
      D <- c( outer(D, pfs_val[k]^c(0:pfs_len[k])) )
    }
    return( sort(D) )
  }
}
primeSieve <- function(n) {
  if (!is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 1)
    stop("Argument 'n' must be an integer number greater or equal 1.")
  if (n > 2^53 - 1)
    stop("Argument 'n' must be smaller than 2^53 - 1.")
  
  if (n < 2) return(c())
  p <- seq(1, n, by=2)
  q <- length(p)
  p[1] <- 2
  if (n >= 9) {
    for (k in seq(3, sqrt(n), by=2)) {
      if (p[(k+1)/2] != 0)
        p[seq((k*k+1)/2, q, by=k)] <- 0
    }    
  }
  p[p > 0]
}


Primes <- function(n1 = 1, n2 = NULL) {
  if (is.null(n2))
    return(primeSieve(n1))
  
  if (!is.numeric(n1) || length(n1) != 1 || floor(n1) != ceiling(n1) || n1 <= 0 ||
      !is.numeric(n2) || length(n2) != 1 || floor(n2) != ceiling(n2) || n2 <= 0 )
    stop("Arguments 'n1' and 'n2' must be integers.")
  
  if (n2 > 2^53 - 1)  stop("Upper bound 'n2' must be smaller than 2^53-1.")
  if (n1 > n2)        stop("Upper bound must be greater than lower bound.")
  
  if (n2 <= 1000) {
    P <- primeSieve(n2)
    return(P[P >= n1])
  }
  
  myPrimes <- primeSieve(floor(sqrt(n2)))
  
  N <- seq.int(n1, n2)
  n <- length(N)
  A <- numeric(n)
  if (n1 == 1) A[1] <- -1
  
  for (p in myPrimes) {
    r <- n1 %% p                                    # rest modulo p
    if (r == 0) { i <- 1 } else { i <- p - r + 1 }  # find next divisible by p
    if (i <= n && N[i] == p) { i <- i + p }         # if it is p itself, skip
    while (i <= n) { A[i] <- 1; i <- i + p }        # mark those divisible by p
  }
  return(N[A == 0])
}

primeFactors <- function(n) {
  if (!is.numeric(n) || length(n) != 1 || n != round(n) || n < 1) {
    warning("Argument 'n' must be a nonnegative integer.")
    return(NULL)
  }
  if (n < 4) return(n) 
  
  if (n <= 2^53 - 1) {
    f <- c()
    p <- Primes(floor(sqrt(n)))
    d <- which(n %% p == 0)
    if (length(d) == 0) return(n)  # n is prime
    
    for (q in p[d]) {
      while (n %% q == 0) {
        f <- c(f, q)
        n <- n/q
      }
    }
    if (n > 1) f <- c(f, n)
    
  } else {
    warning("Argument 'n' too big: use 'gmp::factorize()' for this.")
    f <- NA
  }
  
  return(f)
}
