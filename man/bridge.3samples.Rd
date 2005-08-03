\name{bridge.3samples}
\alias{bridge.3samples}

\title{Bayesian Robust Inference for Differential Gene Expression
  (BRIDGE) with three Samples }
\description{
  Test for differentially expressed genes in a three sample set-up. This
  code can be used with both cDNA microarrays or Affymetrix chip.
}
\usage{
bridge.3samples(sample1,sample2,sample3,B=1000,min.iter=0,batch=10,mcmc.obj=NULL,all.out=TRUE,verbose=FALSE,log=FALSE)
}

\details{
  This code fits a robust Bayesian hierarchical model for 
  testing for differential expression.
  Outliers are modeled explicitly using a $t$-distribution. 
  The model includes an exchangeable prior for the variances which allow different variances for the genes but still
  shrink extreme empirical variances. This function DO NOT perform
  normalization.
  The data should be normalized before hands such as
  centering the mean expression of each experiment.
  %The data should be
  %on the raw scale, i.e. untransformed.
  More details can be found in the
  references below. 
}
\arguments{
\item{sample1}{
    The matrix of intensity from the sample 1. Each row corresponds to a different gene.}
  \item{sample2}{
    The matrix of intensity from the sample 2. Each row corresponds to a
    different gene.}
  \item{sample3}{
    The matrix of intensity from the sample 3. Each row corresponds to a
    different gene.}  
  \item{B}{
    The number of iteration used in the MCMC algorithm.}
  \item{min.iter}{
    The length of the burn-in period in the MCMC algorithm.\code{min.iter} should be less than
    B.
  }
  \item{batch}{The thinning value to be used in the MCMC. Only every \code{batch}-th iteration will be stored.}
  \item{mcmc.obj}{An object of type bridge2, as returned by
    \code{bridge.2samples}. \code{mcmc.obj} can be used to initialized the MCMC.
    If no \code{mcmc.obj}, the MCMC is initialized to the least squares estimates.}
  \item{all.out}{A logical value indicating if all the parameters
    should be output. If \code{all.out} is FALSE, only the posterior
    mean is output. This could be used to save memory. }
  \item{verbose}{A logical value indicating if
    the current MCMC iteration number should be printed out.}
  \item{log}{A logical value indicating if
    the data are log transformed.}
}
\value{
  An object of type \code{bridge3} containing the sampled values from the
  posterior distribution.
  \item{gamma1}{A matrix, each row contains the sampled values from the
    corresponding gene effect in sample 1.}
  \item{gamma2}{A matrix, each row contains the sampled values from the
    corresponding gene effect in sample 2.}
  \item{gamma3}{A matrix, each row contains the sampled values from the
    corresponding gene effect in sample 3.}  
  \item{lambda.gamma1}{
    A vector containing the sampled values for the precision of the gene
    effect prior for the component corresponding to sample 1.}
  \item{lambda.gamma2}{
    A vector containing the sampled values for the precision of the gene
    effect prior for the component corresponding to sample 2.}
  \item{lambda.gamma3}{
    A vector containing the sampled values for the precision of the gene
    effect prior for the component corresponding to sample 3.}
  \item{lambda.gamma12}{
    A vector containing the sampled values for the precision of the gene
    effect prior for the component where sample 1 and sample 2 are combined.}
  \item{lambda.gamma13}{
    A vector containing the sampled values for the precision of the gene
    effect prior for the component where sample 1 and sample 3 are combined.}
  \item{lambda.gamma23}{
    A vector containing the sampled values for the precision of the gene
    effect prior for the component where sample 2 and sample 3 are combined.}
  \item{lambda.gamma123}{
    A vector containing the sampled values for the precision of the gene
    effect prior for the component where all the samples are combined.}
  \item{lambda_eps1}{A matrix, each row contains the sampled values from the
    corresponding gene precision in sample 1.}
  \item{lambda_eps2}{A matrix, each row contains the sampled values from the
    corresponding gene precision in sample 2.}
  \item{lambda_eps3}{A matrix, each row contains the sampled values from the
    corresponding gene precision in sample 3.}
  \item{a.eps1}{A vector containing the sampled values for the mean of
    the prior of the genes precision in sample 1.}
  \item{b.eps1}{A vector containing the sampled values for the variance of
    the prior of the genes precision in sample 1.}
    \item{a.eps2}{A vector containing the sampled values for the mean of
    the prior of the genes precision in sample 2.}
  \item{b.eps2}{A vector containing the sampled values for the variance of
    the prior of the genes precision in sample 2.}
  \item{a.eps3}{A vector containing the sampled values for the mean of
    the prior of the genes precision in sample 3.}
  \item{b.eps3}{A vector containing the sampled values for the variance of
    the prior of the genes precision in sample 3.}  
  \item{w1}{A matrix, each element (i,j) correspond to the posterior mean
    of the sampled weights of replicate j in gene i and sample 1. To save memory, we
    only store the posterior means of the weigths.}
  \item{w2}{A matrix, each element (i,j) correspond to the posterior mean
    of the sampled weights of replicate j in gene i and sample 2. To save memory, we
    only store the posterior means of the weigths.}
  \item{w3}{A matrix, each element (i,j) correspond to the posterior mean
    of the sampled weights of replicate j in gene i and sample 3. To save memory, we
    only store the posterior means of the weigths.}
  \item{w1}{A matrix, each element (i,j) correspond to the posterior mean
    of the sampled weights of replicate j in gene i and sample 1. To save memory, we
    only store the posterior means of the weigths.}
  \item{nu1}{A matrix containing the sampled degrees of freedom in sample
    1.}
  \item{nu2}{A matrix containing the sampled degrees of freedom in sample
    2.}
  \item{nu3}{A matrix containing the sampled degrees of freedom in sample
    3.}
  \item{w.mix}{The posterior mixing proportions in the mixture
    component. }
  \item{prop.model}{The posterior proportions of each component for each
    each gene.}
  \item{move}{The proportion of moves between components. This should be
    used as a diagnostic tool.}
}
\references{Robust Estimation of cDNA Microarray Intensities with Replicates
  Raphael Gottardo, Adrian E. Raftery, Ka Yee Yeung, and Roger Bumgarner
  Department of Statistics, University of Washington, Box 354322, Seattle, WA 98195-4322}
\author{Raphael Gottardo}


\seealso{\code{\link{bridge.2samples}}}
\examples{
sample1<-matrix(exp(rnorm(150)),50,3)
sample2<-matrix(exp(rnorm(200)),50,4)
sample3<-matrix(exp(rnorm(150)),50,3)

mcmc.bridge3<-bridge.3samples(sample1,sample2,sample3,B=10,min.iter=0,batch=10,mcmc.obj=NULL,all.out=TRUE,verbose=FALSE)
}
\keyword{models}
\keyword{robust}
