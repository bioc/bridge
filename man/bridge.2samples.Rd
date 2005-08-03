\name{bridge.2samples}
\alias{bridge.2samples}
\title{Bayesian Robust Inference for Differential Gene Expression
  (BRIDGE) with two Samples}
\description{
  Test for differentially expressed genes in a two sample set-up. This code can be used with both cDNA microarrays or Affymetrix chips.  
}
\usage{
bridge.2samples(sample1,sample2,B=1000,min.iter=0,batch=10,mcmc.obj=NULL,all.out=TRUE,affy=FALSE,verbose=FALSE,log=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\details{
  This code fits a robust Bayesian hierarchical model for 
  testing for differential expression.
  Outliers are modeled explicitly using a $t$-distribution. 
  The model includes an exchangeable prior for the variances which allow
  different variances for the genes but still shrink extreme empirical variances. More details can be found in the
  references below. 
}
\arguments{
\item{sample1}{
    The matrix of intensity from the sample 1. Each row corresponds to a different gene.}
  \item{sample2}{
    The matrix of intensity from the sample 2. Each row corresponds to a different gene.}
  \item{B}{
    The number of iteration used in the MCMC algorithm.}
  \item{min.iter}{
    The length of the burn-in period in the MCMC algorithm.\code{min.iter} should be less than
    B.
  }
  \item{batch}{The thinning value to be used in the MCMC. Only every \code{batch}-th iteration will be stored.}
  \item{mcmc.obj}{An object of type bridge2, as returned by \code{bridge.2samples}. \code{mcmc.obj} is used to initialized the MCMC.
    If no \code{mcmc.obj}, the MCMC is initialized to the least squares estimates.}
  \item{all.out}{A logical value indicating if all the parameters
    should be output. If \code{all.out} is FALSE, only the posterior
    mean is output. This could be used to save memory. }
  \item{affy}{A logical value indicating if the data correspond to affy
    data or cDNA microarray data. If affy=FALSE, a bivariate
    distribution is used.}
  \item{verbose}{A logical value indicating if
    the current MCMC iteration number should be printed out.}
  \item{log}{A logical value indicating if
    the data are log transformed.}
  
}
\value{
  An object of type \code{bridge2} containing the sampled values from the
  posterior distribution.
  \item{gamma1}{A matrix, each row contains the sampled values from the
    corresponding gene effect in sample 1.}
  \item{gamma2}{A matrix, each row contains the sampled values from the
    corresponding gene effect in sample 2.}
  \item{lambda.gamma1}{
    A vector containing the sampled values for the precision of the gene
    effect prior in sample 1.}
  \item{lambda.gamma2}{
    A vector containing the sampled values for the precision of the gene
    effect prior in sample 2.}
  \item{rho}{A vector containing the sampled values from between sample
    correlation coefficient \code{rho}. If affy=TRUE, \code{rho} is
    identically zero.}
  \item{lambda_eps1}{A matrix, each row contains the sampled values from the
    corresponding gene specific error precision in sample 1.}
  \item{lambda_eps2}{A matrix, each row contains the sampled values from the
    corresponding gene specific error precision in sample 2.}
  \item{a.eps1}{A vector containing the sampled values for the mean of
    the prior of the error precisions in sample 1.}
  \item{b.eps1}{A vector containing the sampled values for the variance of
    the prior of the error precisions in sample 1.}
  \item{a.eps2}{A vector containing the sampled values for the mean of
    the prior of the error precisions in sample 2.}
  \item{b.eps2}{A vector containing the sampled values for the variance of
    the prior of the error precisions in sample 2.}  
  \item{w1}{A matrix, each element (i,j) correspond to the posterior mean
    of the sampled weights of replicate j in gene i and sample 1. To save memory, we
    only store the posterior mean of the weigths.}
  \item{nu1}{A matrix containing the sampled degrees of freedom in
    sample 1.}
  \item{w2}{A matrix, each element (i,j) correspond to the posterior mean
    of the sampled weights of replicate j in gene i and sample 2. To save memory, we
    only store the posterior means of the weigths. If affy=FALSE, this
    is identically equal to w2 as we have a bivariate t-distribution.}
  \item{nu2}{A matrix containing the sampled degrees of freedom in
    sample 2. If affy=FALSE, this
    is identically equal to nu2 as we have a bivariate t-distribution.}  
  \item{p}{The mixing proportion in the two component mixture.}
  \item{post.p}{The posterior probabilities of differential expression.}
  \item{move}{The proportion of moves between components. This should be
    used as a diagnostic tool.}
}
\references{
  Robust Estimation of cDNA Microarray Intensities with Replicates
  Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
  Department of Statistics, University of Washington, Box 354322,
  Seattle, WA 98195-4322}

\author{Raphael Gottardo}

\seealso{\code{\link{bridge.3samples}}}
\examples{
data(hiv)

bridge.hiv<-bridge.2samples(hiv[1:10,c(1:4)],hiv[1:10,c(5:8)],B=2000,min.iter=0,batch=1,mcmc.obj=NULL,verbose=FALSE)
}
\keyword{models}
\keyword{robust}
