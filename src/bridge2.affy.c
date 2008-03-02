#include "bridge2_util.h" 


void all_gibbs2(double **y1, int *R1, double **y2, int *R2, int *G, double *mu1, double *mu2, int *model, double *post_p, double *lambda_mu1, double *lambda_mu2, double *lambda_mu12, double *lambda1, double *lambda2, double **weight1, double **weight2, double *w, double *nu1, double *nu2, double *nu_choice, int *nb_nu, double *a_eps1, double *b_eps1, double *a_eps2, double *b_eps2, double *move, int *robust);

void up_date_weight_nu2(double **y, int R, int G, double *lambda, double *mu,
		       double **weight, double *nu_choice, int nb_nu, double *nu);
void up_date_lambda_mu2(int G, double *mu1, double *mu2, int *model, double *lambda_mu1, double *lambda_mu2, 
		       double *lambda_mu12);
void dirichlet2(double *gamma, int n, double *w);
void up_date_w2(int G, int *model, double *w);
void gibbs_lambda2(double *y, int R, int G, double mu, double *weight,  double *lambda, double a_eps, double b_eps);
void up_date_a_b2(double *lambda1, int G, double *a_eps, double *b_eps);
void gibbs_mu2(double *y1, int R1, double *y2, int R2, double *mu1, double *mu2, int *model, double *post_p,
	      double lambda_mu1, double lambda_mu2, double lambda_mu12, double lambda1, double lambda2, 
	      double *weight1, double *weight2, double *w);


void gene_express_2s(double *y1_vec, int *R1, double *y2_vec, int *R2, int *G, double *mu1, 
		     double *mu2, double *mu1_out, double *mu2_out, int *model, 
		     int *model_out, double *post_p, double *lambda_mu1, double *lambda_mu1_out, 
		     double *lambda_mu2, double *lambda_mu2_out,
		     double *lambda_mu12, double *lambda_mu12_out,   
		     double *lambda1, double *lambda1_out, double *lambda2, double *lambda2_out, 
		     double *a_eps1, double *a_eps1_out,
		     double *b_eps1, double *b_eps1_out,
		     double *a_eps2, double *a_eps2_out,
		     double *b_eps2, double *b_eps2_out,
		     double *move,
		     double *weight1_vec, double *weight2_vec,  
		     double *nu1, double *nu1_out, double *nu2, double *nu2_out, 
		     double *nu_choice, int *nb_nu,
		     double *w, double *w_out, int *min_iter, 
		     int *batch, int *B, int *all_out, int *verbose, int *robust)
{
  int old_model;
  int count_model=0;
  int k;
  int count=0,count2=0;
  double **y1, **y2;
  double **weight1,**weight2;
  int i,g;

  y1=dmatrix(*G, *R1);
  y2=dmatrix(*G, *R2);
  weight1=dmatrix(*G, *R1);
  weight2=dmatrix(*G, *R2);
  
  vec_mat(y1_vec,G,R1,y1);
  vec_mat(y2_vec,G,R2,y2);
  vec_mat(weight1_vec,G,R1,weight1);
  vec_mat(weight2_vec,G,R2,weight2);
  

  GetRNGstate();

 
  for(k=0;k<*B;k++)
    { 
      if((((k+1)*100)%(10**B)==0) & (*verbose==1))
	{
	  printf("%d percent completed \n",(((k+1)*100)/(*B)));
	}
      all_gibbs2(y1, R1, y2, R2, G, mu1, mu2, model, post_p, lambda_mu1, lambda_mu2, lambda_mu12, lambda1, lambda2, weight1, weight2, w, nu1, nu2, nu_choice, nb_nu, a_eps1, b_eps1, a_eps2, b_eps2, move, robust);

      
      if(k>=*min_iter)
	{
	  if(((count2+1)%(*batch))==0) /** Batch sampling **/
	    {
	      if(*all_out==1)
		{
		  for(g=0;g<*G;g++)
		    {
		      mu1_out[count**G+g]=mu1[g];
		      mu2_out[count**G+g]=mu2[g];
		      model_out[count**G+g]=model[g];
		      
		      lambda1_out[count**G+g]=lambda1[g];	      	      
		      lambda2_out[count**G+g]=lambda2[g];
		      
		      for(i=0;i<*R1;i++)
			weight1_vec[i**G+g]+=weight1[g][i]/((*B-*min_iter)/(*batch));
		      for(i=0;i<*R2;i++)
			weight2_vec[i**G+g]+=weight2[g][i]/((*B-*min_iter)/(*batch));
		    }
		  	      
		  for(i=0;i<2;i++)
		    w_out[count*2+i]=w[i];
		  
		  for(i=0;i<*R1;i++)
		    nu1_out[count**R1+i]=nu1[i];
		  for(i=0;i<*R2;i++)
		    nu2_out[count**R2+i]=nu2[i];
		  
		  a_eps1_out[count]=*a_eps1;
		  b_eps1_out[count]=*b_eps1;
		  a_eps2_out[count]=*a_eps2;
		  b_eps2_out[count]=*b_eps2;
		  
		  lambda_mu1_out[count]=*lambda_mu1;
		  lambda_mu2_out[count]=*lambda_mu2;
		  lambda_mu12_out[count]=*lambda_mu12;
		}
	      else /** Posterior mean **/
		{
		  for(g=0;g<*G;g++)
		    {
		      mu1_out[g]+=mu1[g]/((*B-*min_iter)/(*batch));
		      mu2_out[g]+=mu2[g]/((*B-*min_iter)/(*batch));
		      model_out[g]=model[g];
		      
		      lambda1_out[g]+=lambda1[g]/((*B-*min_iter)/(*batch));	      	      
		      lambda2_out[g]+=lambda2[g]/((*B-*min_iter)/(*batch));
		      
		      for(i=0;i<*R1;i++)
			weight1_vec[i**G+g]+=weight1[g][i]/((*B-*min_iter)/(*batch)+1.);
		      for(i=0;i<*R2;i++)
			weight2_vec[i**G+g]+=weight2[g][i]/((*B-*min_iter)/(*batch)+1.);
		    }
		  	      
		  for(i=0;i<2;i++)
		    w_out[i]+=w[i]/((*B-*min_iter)/(*batch));
	      
		  for(i=0;i<*R1;i++)
		    nu1_out[i]+=nu1[i]/((*B-*min_iter)/(*batch));
		  for(i=0;i<*R2;i++)
		    nu2_out[i]+=nu2[i]/((*B-*min_iter)/(*batch));
		  
		  *a_eps1_out+=*a_eps1/((*B-*min_iter)/(*batch));
		  *b_eps1_out+=*b_eps1/((*B-*min_iter)/(*batch));
		  *a_eps2_out+=*a_eps2/((*B-*min_iter)/(*batch));
		  *b_eps2_out+=*b_eps2/((*B-*min_iter)/(*batch));
		  
		  *lambda_mu1_out+=*lambda_mu1/((*B-*min_iter)/(*batch));
		  *lambda_mu2_out+=*lambda_mu2/((*B-*min_iter)/(*batch));
		  *lambda_mu12_out+=*lambda_mu12/((*B-*min_iter)/(*batch));
		}
	      count++;
	    }
	  count2++;
	}
    }
  
  PutRNGstate();

  free_dmatrix(y1, *G);
  free_dmatrix(y2, *G);
  free_dmatrix(weight1, *G);
  free_dmatrix(weight2, *G);
  
}


void all_gibbs2(double **y1, int *R1, double **y2, int *R2, int *G, double *mu1, double *mu2, int *model, double *post_p, double *lambda_mu1, double *lambda_mu2, double *lambda_mu12, double *lambda1, double *lambda2, double **weight1, double **weight2, double *w, double *nu1, double *nu2, double *nu_choice, int *nb_nu, double *a_eps1, double *b_eps1, double *a_eps2, double *b_eps2, double *move, int *robust)
{
  int g;
  int old_model;
  
  for(g=0;g<*G;g++)
    {
      old_model=model[g];
      gibbs_mu2(y1[g], *R1, y2[g], *R2, mu1+g, mu2+g, model+g, post_p+g, 
		*lambda_mu1, *lambda_mu2, *lambda_mu12,lambda1[g], lambda2[g], weight1[g], weight2[g], w);
      /** Count the number of moves between components **/
      if(model[g]!=old_model)
	move[g]+=1;
      gibbs_lambda2(y1[g], *R1, *G, mu1[g], weight1[g], lambda1+g, *a_eps1, *b_eps1);
      gibbs_lambda2(y2[g], *R2, *G, mu2[g], weight2[g], lambda2+g, *a_eps2, *b_eps2);
    }

  up_date_lambda_mu2(*G, mu1, mu2, model, lambda_mu1, lambda_mu2, lambda_mu12);
  if(*robust==1)
  {
    up_date_weight_nu2(y1, *R1, *G, lambda1, mu1, weight1, nu_choice, *nb_nu, nu1);
    up_date_weight_nu2(y2, *R2, *G, lambda2, mu2, weight2, nu_choice, *nb_nu, nu2);
  }
  up_date_w2(*G, model, w);
  up_date_a_b2(lambda1, *G, a_eps1, b_eps1);
  up_date_a_b2(lambda2, *G, a_eps2, b_eps2);
  
}
     

void up_date_a_b2(double *lambda1, int G, double *a_eps, double *b_eps)
{
  double sum_lambda=0;
  double sum_log_lambda=0;
  int i;

  for(i=0;i<G;i++)
    {
      sum_lambda+=lambda1[i];
      sum_log_lambda+=log(lambda1[i]);
    }

  
  *a_eps=slice_sampling_a(*a_eps, 2, 100, sum_log_lambda, sum_lambda, *b_eps, G); 
  *b_eps=slice_sampling_b(*b_eps, 2, 100, sum_log_lambda, sum_lambda, *a_eps, G); 
}
     

void gibbs_mu2(double *y1, int R1, double *y2, int R2, double *mu1, double *mu2, int *model, double *post_p, 
	      double lambda_mu1, double lambda_mu2, double lambda_mu12, double lambda1, double lambda2, double *weight1, double *weight2, double *w)
{
  int i;
  double sum_y1=0,sum_y2=0;
  double sum_w1=0,sum_w2=0;
  double K[2],sum_K=0;
  double w_star[2];
  double u=runif(0,1);
  double sum_y1_squared=0;
  double sum_y2_squared=0;

  /** Compute the sum of the y **/
  
  for(i=0;i<R1;i++)
    {
      sum_y1+=weight1[i]*y1[i];
      sum_y1_squared+=weight1[i]*y1[i]*y1[i];
      sum_w1+=weight1[i];
    }
  for(i=0;i<R2;i++)
    {
      sum_y2+=weight2[i]*y2[i];
      sum_y2_squared+=weight2[i]*y2[i]*y2[i];
      sum_w2+=weight2[i];
    }  

  
  //  K[0]=w[0]*sqrt(lambda_mu12)/sqrt(sum_w1*lambda1+sum_w2*lambda2+lambda_mu12)*exp(fmin2(0.5*(sum_y1*lambda1+sum_y2*lambda2)*(sum_y1*lambda1+sum_y2*lambda2)/(sum_w1*lambda1+sum_w2*lambda2+lambda_mu12)-0.5*lambda1*sum_y1_squared-0.5*lambda2*sum_y2_squared,100));
  
  //  K[1]=w[1]*sqrt(lambda_mu1*lambda_mu2)/sqrt((sum_w1*lambda1+lambda_mu1)*(sum_w2*lambda2+lambda_mu2))*exp(fmin2(0.5*(lambda1*sum_y1)*(lambda1*sum_y1)/(sum_w1*lambda1+lambda_mu1)-0.5*lambda1*sum_y1_squared+0.5*(lambda2*sum_y2)*(lambda2*sum_y2)/(sum_w2*lambda2+lambda_mu2)-0.5*lambda2*sum_y2_squared,100));


  K[0]=w[0]*sqrt(lambda_mu12)/sqrt(sum_w1*lambda1+sum_w2*lambda2+lambda_mu12);

  K[1]=w[1]*sqrt(lambda_mu1*lambda_mu2)/sqrt((sum_w1*lambda1+lambda_mu1)*(sum_w2*lambda2+lambda_mu2))*exp(0.5*(lambda1*sum_y1)*(lambda1*sum_y1)/(sum_w1*lambda1+lambda_mu1)-0.5*lambda1*sum_y1_squared+0.5*(lambda2*sum_y2)*(lambda2*sum_y2)/(sum_w2*lambda2+lambda_mu2)-0.5*lambda2*sum_y2_squared-0.5*(sum_y1*lambda1+sum_y2*lambda2)*(sum_y1*lambda1+sum_y2*lambda2)/(sum_w1*lambda1+sum_w2*lambda2+lambda_mu12)+0.5*lambda1*sum_y1_squared+0.5*lambda2*sum_y2_squared);


  sum_K=K[0]+K[1];
  w_star[0]=K[0]/(sum_K);
  w_star[1]=(K[0]+K[1])/(sum_K);
  

  if(u<w_star[0]) /** They are all equal **/ 
    {
      *mu1=rnorm((lambda1*sum_y1+lambda2*sum_y2)/(sum_w1*lambda1+sum_w2*lambda2+lambda_mu12),1./sqrt(sum_w1*lambda1+sum_w2*lambda2+lambda_mu12));
      *mu2=*mu1;
      *model=0;
    }
  else 
    {
      *mu1=rnorm(lambda1*sum_y1/(sum_w1*lambda1+lambda_mu1),1./sqrt(sum_w1*lambda1+lambda_mu1));
      *mu2=rnorm(lambda2*sum_y2/(sum_w2*lambda2+lambda_mu2),1./sqrt(sum_w2*lambda2+lambda_mu2));
      *model=1;
    }
  /** Rao-Blackwellization **/
  (*post_p)=(*post_p)+(1.-w_star[0]);
  
}

void gibbs_lambda2(double *y, int R, int G, double mu, double *weight,  double *lambda, double a_eps, double b_eps)
{
  double sum_res_y=0;
  int i;

  /** Compute the sum of the y **/
  
  sum_res_y=0;
  for(i=0;i<R;i++)
    sum_res_y+=weight[i]*(y[i]-mu)*(y[i]-mu);

  *lambda=rgamma(a_eps*a_eps/b_eps+R/2.,1./(sum_res_y/2.+a_eps/b_eps));
  
}

void up_date_w2(int G, int *model, double *w)
{
  
  int n[]={0,0};
  double w_new[2];
  double gamma[2];
  double gamma0=1;
  int i;

  for(i=0;i<G;i++)
    {
      if(model[i]==0) /** All the same **/
	{
	  n[0]++;
	}
      else if(model[i]==1)
	{
	  n[1]++;
	}
    }
  gamma[0]=gamma0+n[0];
  gamma[1]=gamma0+n[1];
  
  /** Up-date the weight **/
  dirichlet(gamma, 2, w);
  
}

void dirichlet2(double *gamma, int n, double *w)
{
  int i;
  double sum_w=0;
  
  for(i=0;i<n;i++)
    {
      w[i]=rgamma(gamma[i],1);
      sum_w+=w[i];
    }
  for(i=0;i<n;i++)
    w[i]=w[i]/sum_w;
  
}

void up_date_lambda_mu2(int G, double *mu1, double *mu2, int *model, double *lambda_mu1, double *lambda_mu2,  
		       double *lambda_mu12)
{
  
  double sum1=0, sum2=0, sum12=0;
  int n1=0, n2=0, n12=0;
  int i;

  for(i=0;i<G;i++)
    {
      if(model[i]==0) /** All the same **/
	{
	  sum12+=mu1[i]*mu1[i];
	  n12++;
	}
      else
	{
	  sum1+=mu1[i]*mu1[i];
	  n1++;
	  sum2+=mu2[i]*mu2[i];
	  n2++;	  
	}
    }
  *lambda_mu1=rgamma(n1/2.+1.,1./(sum1/2.+0.005));
  *lambda_mu2=rgamma(n2/2.+1.,1./(sum2/2.+0.005));
  *lambda_mu12=rgamma(n12/2.+1.,1./(sum12/2.+0.005));
}

void up_date_weight_nu2(double **y, int R, int G, double *lambda, double *mu,
		       double **weight, double *nu_choice, int nb_nu, double *nu)
{
  
  int i,j;
  double sum_l_nu=0, sum_l_nu_new=0;
  int nu_new=0;
  
  /** Up-date the weight and the df for the t-distribution **/
  /** A block up-date is used **/
  /** First the parameter w is integrated out and df is up-dated **/
  /** Then the w are up-dated by Gibbs sampling **/
  

  for(j=0;j<R;j++)
    {
      sum_l_nu_new=0.;
      sum_l_nu=0.;
      nu_new=nu_choice[(int)(runif(0.,1.)*(nb_nu))];  

      for(i=0;i<G;i++)
	{ 
	  sum_l_nu_new+=lgammafn((nu_new+1.)/2.)-lgammafn(nu_new/2.)+log(2./nu_new)/2.
	    -((nu_new+1.)/2.)*log(1.+lambda[i]*(y[i][j]-mu[i])*(y[i][j]-mu[i])/nu_new);
	  //sum_l_nu_new+=dgamma(weight[i][j],nu_new/2.,2./nu_new,1);

	  sum_l_nu+=lgammafn((nu[j]+1.)/2.)-lgammafn(nu[j]/2.)+log(2./nu[j])/2.
	  -((nu[j]+1.)/2.)*log(1.+lambda[i]*(y[i][j]-mu[i])*(y[i][j]-mu[i])/nu[j]);
	  //sum_l_nu+=dgamma(weight[i][j],nu[j]/2.,2./nu[j],1);
	}
      
      if((sum_l_nu_new-sum_l_nu)>log(runif(0.,1.)))
	{
	  nu[j]=nu_new;
	}
    }
  for(i=0;i<G;i++)
    for(j=0;j<R;j++)
      {
	weight[i][j]=rgamma((nu[j]+1.)/2.,1./(lambda[i]/2.*(y[i][j]-mu[i])*(y[i][j]-mu[i])+nu[j]/2.));
	//weight[i][j]=1.;
	//weight[i][j]=rgamma((nu[j]+1.)/2.,1./(nu[j]/2.));
      }
}



