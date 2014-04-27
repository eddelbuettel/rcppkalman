

library(RcppKalman)

## example is from the vignette in package expm
M <- matrix(c(4, 1, 1, 2, 4, 1, 0, 1, 4), 3, 3)

## expected output
expM <- matrix(c(147.8666, 127.7811, 127.7811, 183.7651, 183.7651, 163.6796, 71.79703, 91.88257, 111.96811), 3, 3)

## we only have the expected result to about six digits
all.equal(expm(M), expM, tolerance=1.0e-6)


