# Linear Model Notes
# April 11th, 2017

# https://ms.mcmaster.ca/~bolker/emdbook/book.pdf

# p. 402: analysis of covariance (ANCOVA)

lm(Y ~ f * x)

where *f* is a factor and *x* is a covariate:  
    - the formula *Y~f+x* would specify parallel slopes  
    - *Y~f* would specify zero slopes but different intercepts  
    - *Y~x* would specify a single slope