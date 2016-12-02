
permute_test <- function(x1, x2, n, stat=mean, tail="both") {
  n1 <- length(x1)
  n2 <- length(x2)
  x <- c(x1, x2)
  n <- min(n, choose(n1+n2, n1))
  diffs <- numeric(n)

  if (tail == "both") {
    diff_fun <- function(a,b) abs(a - b)
  } else if (tail == "left") {
    diff_fun <- function(a,b) a - b
  } else {
    diff_fun <- function(a,b) b - a
  }

  for (i in seq_along(diffs)) {
    i1 <- sample.int(n1+n2, size=n1)
    s1 <- stat(x[i1])
    s2 <- stat(x[-i1])
    diffs[i] <- diff_fun(s1, s2)
  }
  pval <- sum(diff_fun(stat(x1), stat(x2)) > diffs) / n

  return(list(diffs=diffs, pval=pval))
}

x1 <- c(6.52559065,
        3.968626967,
        2.020725942,
        9.733961167,
        4.509249753,
        19.27649692)

x2 <- c(11.01514109,
        7.686568372,
        2.362907813,
        16.653328,
        7.505553499,
        12.76714533)

test <- permute_test(x1, x2, n=1000, tail="left")
