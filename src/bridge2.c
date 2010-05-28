#include "bridge2_util.h" 

void mcmc(double **data1,double **data2, int *n1, int *n2, int *nb_col1,
  int *B, int *dye_swap,
  double *gamma_1, double *gamma1_p, 
  double *gamma_2, double *gamma2_p, 
  double *p,
  double *p_p,
  double *mu, 
  double *beta2, 
  double *alpha2,
  double *delta22,
  double *eta, 
  double *lambda_eps1, double *lambda_eps1_p,
  double *lambda_eps2, double *lambda_eps2_p,
  double *a_eps1, double *b_eps1,
  double *a_eps1_p, double *b_eps1_p,
  double *a_eps2, double *b_eps2,
  double *a_eps2_p, double *b_eps2_p,	  
  double *w, double *df_choice, int *nb_df,
  double *df, double *df_p, double *w_p,
  double *rho, double *rho_p,
  double *shift, 
  double *lambda_gamma1, double *lambda_gamma2,
  double *lambda_gamma1_p, double *lambda_gamma2_p,
  double *lambda_gamma, double *lambda_gamma_p,
  double *post_prob,
  double *move,
  int *min_iter, int *batch, int *all_out, int *verbose, int *robust);


void ex_R_link_bridge2(double *data_vec1, double *data_vec2, 
  int *n1, int *n2, int *nb_col1, int *B, int *dye_swap,
  double *gamma_1, double *gamma1_p, 
  double *gamma_2, double *gamma2_p, 
  double *p,
  double *p_p,
  double *mu, 
  double *beta2, 
  double *alpha2, 
  double *delta22,
  double *eta, 
  double *lambda_eps1, double *lambda_eps1_p,
  double *lambda_eps2, double *lambda_eps2_p,
  double *a_eps1, double *b_eps1,
  double *a_eps1_p, double *b_eps1_p,
  double *a_eps2, double *b_eps2,
  double *a_eps2_p, double *b_eps2_p,
  double *w, double *df_choice, int *nb_df,
  double *df, double *df_p, double *w_p,
  double *rho, double *rho_p,
  double *shift,
  double *lambda_gamma1, double *lambda_gamma2,
  double *lambda_gamma1_p, double *lambda_gamma2_p,
  double *lambda_gamma, double *lambda_gamma_p,
  double *post_prob,
  double *move,
  int *min_iter, int *batch, int *all_out,int *verbose, int *robust)
{

  double **data1;
  double **data2;

  GetRNGstate();
  data1=dmatrix(*n1, *n2);  
  vec_mat(data_vec1,n1,n2,data1);
  data2=dmatrix(*n1, *n2);  
  vec_mat(data_vec2,n1,n2,data2);


  mcmc(data1, data2, n1, n2, nb_col1,
    B, dye_swap,
    gamma_1, gamma1_p, gamma_2, gamma2_p,
    p,
    p_p,
    mu, 
    beta2, 
    alpha2, 
    delta22,
    eta, 
    lambda_eps1, lambda_eps1_p,
    lambda_eps2, lambda_eps2_p,
    a_eps1, b_eps1,
    a_eps1_p, b_eps1_p,
    a_eps2, b_eps2,
    a_eps2_p, b_eps2_p,
    w, df_choice, nb_df,
    df, df_p, w_p,
    rho, rho_p, 
    shift,
    lambda_gamma1, lambda_gamma2,
    lambda_gamma1_p, lambda_gamma2_p,
    lambda_gamma, lambda_gamma_p,
    post_prob,
    move,
    min_iter, batch, all_out, verbose, robust);

  PutRNGstate();

  free_dmatrix(data1, *n1);
  free_dmatrix(data2, *n1);


}



void mcmc(double **data1,double **data2, int *n1, int *n2, int *nb_col1,
  int *B, int *dye_swap, double *gamma_1, double *gamma1_p, 
  double *gamma_2, double *gamma2_p, 
  double *p,
  double *p_p,
  double *mu, 
  double *beta2, 
  double *alpha2, 
  double *delta22,
  double *eta, 
  double *lambda_eps1, double *lambda_eps1_p,
  double *lambda_eps2, double *lambda_eps2_p,
  double *a_eps1, double *b_eps1,
  double *a_eps1_p, double *b_eps1_p,
  double *a_eps2, double *b_eps2,
  double *a_eps2_p, double *b_eps2_p,
  double *w, double *df_choice, int *nb_df,
  double *df, double *df_p, double *w_p,
  double *rho, double *rho_p,
  double *shift,
  double *lambda_gamma1, double *lambda_gamma2,
  double *lambda_gamma1_p, double *lambda_gamma2_p,
  double *lambda_gamma, double *lambda_gamma_p,
  double *post_prob,
  double *move,
  int *min_iter, int *batch, int *all_out, int *verbose, int *robust)
{
  int i,j,k;
  int count=0,count2=0;

  double rho_new;
  double lambda_eps1_new, lambda_eps2_new;
  double l_epsilon=0,l_epsilon_new=0;
  double dens_epsilon=0,dens_epsilon_new=0;
  double sum_rho=0,sum_rho_new=0;
  double dens_rho,dens_rho_new;
  double SSR1=0,SSR2=0,SS12=0;
  double pi=3.14;
  double Sgamma1=0;
  double Sgamma2=0;
  double Sgamma=0;
  double SSgamma1=0, SSgamma1sd=0, SSgamma2=0, SSgamma2sd=0;
  double SSgamma=0, SSgammasd=0;
  double SSeps_sd=0.;
  double SSeps_sd_new=0.;
  double *sum_l_w, *sum_l_w_new;
  double dens_t=0,dens_t_new=0;
  double *df_new;


  /** Parameters for the full conditional of mu **/
  double post_mean_mu=0,post_prec_mu=0;
  double post_mean_beta=0,post_prec_beta=0;
  double post_mean_alpha=0,post_prec_alpha=0;
  double post_mean_delta=0,post_prec_delta=0;
  double post_mean_eta=0,post_prec_eta=0;
  double mean_gamma1=0,mean_gamma2=0;
  double l_ab=0,l_ab_new=0;
  double a_eps_new=0,b_eps_new=0;


  /** Count the number of acceptance with Metropolis **/
  int counteps=0;
  int countrho=0;
  int countdf=0;
  int counta=0;
  int counta_eps=0;
  int countb_eps=0;
  int countgamma=0;
  int countp=0;

  /** Parameter used in the prior of the gamma's **/
  double *mu_gamma1, *mu_gamma2, *mu_gamma;
  double norm_const_gamma=0, norm_const_gamma1=0, norm_const_gamma2=0;
  double gamma_1_new=0, gamma_2_new=0;
  int nb_gamma_same=0,nb_gamma_diff=0;
  double p_star;
  double SSR1_gamma_new, SSR2_gamma_new, SS12_gamma_new;
  double dens_gamma=0, dens_gamma_new=0;
  double p_new=0;
  double likelihood_p=0., likelihood_p_new=0.;

  /** Parameter used to up-date lambda_eps **/
  double sum_res1;
  double sum_res2;
  /** Parameter used in the slice sampling **/
  double sum_lambda1=0;
  double sum_log_lambda1=0;
  double sum_lambda2=0;
  double sum_log_lambda2=0;  
  double old_b=0,width_b=5;
  double old_a=0,width_a=5;



  /** use for the t-distribution **/
  df_new=dvector(*n2,1);
  sum_l_w=dvector(*n2,0);
  sum_l_w_new=dvector(*n2,0);
  /** Set the prior means to zero **/
  mu_gamma1=dvector(1,0);
  mu_gamma2=dvector(1,0);
  mu_gamma=dvector(1,0);


  for(k=0;k<*B;k++)
  {

    if((((k+1)*100)%(10**B)==0) & (*verbose==1))
      printf("%d percent completed \n",(((k+1)*100)/(*B)));





      /** Update the gamma's **/
      /** Active genes **/
    Sgamma1=0.;Sgamma2=0.;
    Sgamma=0;
    mean_gamma1=0;
    mean_gamma2=0;
    nb_gamma_diff=0.;
    nb_gamma_same=0.;
    likelihood_p=0;
    likelihood_p_new=0;


    for(i=0;i<*n1;i++)
    {

      SSgamma1=0;SSgamma2=0;
      SSgamma1sd=0;SSgamma2sd=0;
      SSgamma=0;
      SSgammasd=0;
      norm_const_gamma=0.;
      norm_const_gamma1=0.;
      norm_const_gamma2=0.;

    /** First part of the loop before dye swap **/
      for(j=0;j<(*nb_col1);j++)
      {

        SSgamma1+=w[j*(*n1)+i]*(lambda_eps1[i]*(data1[i][j]-*mu-eta[j]));  
        SSgamma1sd+=w[j*(*n1)+i]*(lambda_eps1[i]);

        norm_const_gamma1+=w[j*(*n1)+i]*lambda_eps1[i]*(data1[i][j]-*mu-eta[j])*(data1[i][j]-*mu-eta[j]);
        norm_const_gamma2+=w[j*(*n1)+i]*lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-eta[j])*(data2[i][j]-*mu-*alpha2-eta[j]);


        SSgamma2sd+=w[j*(*n1)+i]*(lambda_eps2[i]);
        SSgamma2+=w[j*(*n1)+i]*(lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-eta[j]));

      }

    /** Second part of the loop after dye swap **/
      for(j=(*nb_col1);j<*n2;j++)
      {

        SSgamma1+=w[j*(*n1)+i]*(lambda_eps1[i]*(data1[i][j]-*mu-*beta2-eta[j]));
        SSgamma1sd+=w[j*(*n1)+i]*(lambda_eps1[i]);

        SSgamma2sd+=w[j*(*n1)+i]*(lambda_eps2[i]);
        SSgamma2+=w[j*(*n1)+i]*(lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-eta[j]));

        norm_const_gamma1+=w[j*(*n1)+i]*lambda_eps1[i]*(data1[i][j]-*mu-*beta2-eta[j])*(data1[i][j]-*mu-*beta2-eta[j]);

        norm_const_gamma2+=w[j*(*n1)+i]*lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-eta[j]);

      }


      norm_const_gamma=sqrt(*lambda_gamma/(SSgamma1sd+SSgamma2sd+*lambda_gamma))
      //*exp(0.5*(SSgamma1+SSgamma2)*(SSgamma1+SSgamma2)/(SSgamma1sd+SSgamma2sd+*lambda_gamma));
        *exp(-0.5*(norm_const_gamma1+norm_const_gamma2-(SSgamma1+SSgamma2)*(SSgamma1+SSgamma2)/(SSgamma1sd+SSgamma2sd+*lambda_gamma)));

      norm_const_gamma1=sqrt(*lambda_gamma1/(SSgamma1sd+*lambda_gamma1))
      //*exp(0.5*(SSgamma1)*(SSgamma1)/(SSgamma1sd+*lambda_gamma1));
        *exp(-0.5*(norm_const_gamma1-(SSgamma1)*(SSgamma1)/(SSgamma1sd+*lambda_gamma1)));
      norm_const_gamma2=sqrt(*lambda_gamma2/(SSgamma2sd+*lambda_gamma2))
      //*exp(0.5*(SSgamma2)*(SSgamma2)/(SSgamma2sd+*lambda_gamma2));
        *exp(-0.5*(norm_const_gamma2-(SSgamma2)*(SSgamma2)/(SSgamma2sd+*lambda_gamma2)));

    //*p=.5;
      p_star=(1.-(1.-*p)*norm_const_gamma/((1.-*p)*norm_const_gamma+*p*norm_const_gamma1*norm_const_gamma2));

    //printf("p_star:%g\n",p_star);

      if((1.-p_star)>runif(0,1)) 
      {
        gamma_1_new=rnorm((SSgamma1+SSgamma2)/(*lambda_gamma+SSgamma1sd+SSgamma2sd),1./sqrt(*lambda_gamma+SSgamma1sd+SSgamma2sd));
        gamma_2_new=gamma_1_new;

      }
      else
      {
        gamma_1_new=rnorm((SSgamma1)/(*lambda_gamma1+SSgamma1sd),1./sqrt(*lambda_gamma1+SSgamma1sd));

        gamma_2_new=rnorm((SSgamma2)/(*lambda_gamma2+SSgamma2sd),1./sqrt(*lambda_gamma2+SSgamma2sd));
      }


      SSR1=0;SSR2=0;SS12=0;
      SSR1_gamma_new=0;SSR2_gamma_new=0;SS12_gamma_new=0;
    /** First part of the loop before dye swap **/
      for(j=0;j<(*nb_col1);j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j]);
        SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
        SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);

        SSR1_gamma_new+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1_new-eta[j])*(data1[i][j]-*mu-gamma_1_new-eta[j]);
        SSR2_gamma_new+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2_new-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2_new-eta[j]);
        SS12_gamma_new+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1_new-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2_new-eta[j]);

      }

    /** Second part of the loop after dye swap **/
      for(j=(*nb_col1);j<*n2;j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]);
        SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
        SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);

        SSR1_gamma_new+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1_new-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1_new-eta[j]);
        SSR2_gamma_new+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2_new-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2_new-eta[j]);
        SS12_gamma_new+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1_new-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2_new-eta[j]);


      }
    /** Likelihood **/
      dens_gamma=-1./(2*(1-*rho*(*rho)))*(lambda_eps1[i]*SSR1-2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*SS12+lambda_eps2[i]*SSR2);

      dens_gamma_new=-1./(2*(1-*rho*(*rho)))*(lambda_eps1[i]*SSR1_gamma_new-2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*SS12_gamma_new+lambda_eps2[i]*SSR2_gamma_new);

      if(gamma_1[i]==gamma_2[i])
      {
        /** Add the prior **/
        dens_gamma+=log(1-*p)+dnorm(gamma_1[i],0,1/sqrt(*lambda_gamma),1);
        /** Add the proposal **/
        dens_gamma-=log(1-p_star)+dnorm(gamma_1[i],(*lambda_gamma*(*mu_gamma)+SSgamma1+SSgamma2)/(*lambda_gamma+SSgamma1sd+SSgamma2sd),1./sqrt(*lambda_gamma+SSgamma1sd+SSgamma2sd),1);

      }
      else
      {
        dens_gamma+=log(*p)+dnorm(gamma_1[i],0,1./sqrt(*lambda_gamma1),1)+dnorm(gamma_2[i],0,1./sqrt(*lambda_gamma2),1);
        dens_gamma-=log(p_star)+dnorm(gamma_1[i],(*lambda_gamma1*(*mu_gamma1)+SSgamma1)/(*lambda_gamma1+SSgamma1sd),1./sqrt(*lambda_gamma1+SSgamma1sd),1)+dnorm(gamma_2[i],(*lambda_gamma2*(*mu_gamma2)+SSgamma2)/(*lambda_gamma2+SSgamma2sd),1./sqrt(*lambda_gamma2+SSgamma2sd),1.);
      }

      if(gamma_1_new==gamma_2_new)
      {
        /** Add the prior **/
        dens_gamma_new+=log(1.-*p)+dnorm(gamma_1_new,0,1./sqrt(*lambda_gamma),1);
        /** Add the proposal **/
        dens_gamma_new-=log(1.-p_star)+dnorm(gamma_1_new,(*lambda_gamma*(*mu_gamma)+SSgamma1+SSgamma2)/(*lambda_gamma+SSgamma1sd+SSgamma2sd),1./sqrt(*lambda_gamma+SSgamma1sd+SSgamma2sd),1);

      }
      else
      {
        dens_gamma_new+=log(*p)+dnorm(gamma_1_new,0,1./sqrt(*lambda_gamma1),1)+dnorm(gamma_2_new,0,1./sqrt(*lambda_gamma2),1);
        dens_gamma_new-=log(p_star)+dnorm(gamma_1_new,(*lambda_gamma1*(*mu_gamma1)+SSgamma1)/(*lambda_gamma1+SSgamma1sd),1./sqrt(*lambda_gamma1+SSgamma1sd),1)+dnorm(gamma_2_new,(*lambda_gamma2*(*mu_gamma2)+SSgamma2)/(*lambda_gamma2+SSgamma2sd),1./sqrt(*lambda_gamma2+SSgamma2sd),1.);
      }


    /** Accept the up-dated variables **/
      if((dens_gamma_new-dens_gamma)>log(runif(0,1)))
      {
        /** count the # of times we jump from one component to the other **/
        if((gamma_1[i]==gamma_2[i]) & (gamma_1_new!=gamma_2_new))
          move[i]+=1./(*B);
        if((gamma_1[i]!=gamma_2[i]) & (gamma_1_new==gamma_2_new))
          move[i]+=1./(*B);


        gamma_1[i]=gamma_1_new;
        gamma_2[i]=gamma_2_new;
        countgamma++;
      }

    /** Variable to up-date lambda_gamma1 and lambda_gamma2 **/
      if(gamma_1[i]==gamma_2[i])
      {
        likelihood_p+=log(1-*p)+dnorm(gamma_1[i],0,1./sqrt(*lambda_gamma),1.);
        likelihood_p_new+=log(1-p_new)+dnorm(gamma_1[i],0,1./sqrt(*lambda_gamma),1.);

        Sgamma+=(gamma_1[i]-*mu_gamma)*(gamma_1[i]-*mu_gamma);
        nb_gamma_same++;

      }
      else
      {
        likelihood_p+=log(*p)+dnorm(gamma_1[i],0,1./sqrt(*lambda_gamma1),1.)+dnorm(gamma_2[i],0,1./sqrt(*lambda_gamma2),1.);
        likelihood_p_new+=log(p_new)+dnorm(gamma_1[i],0,1./sqrt(*lambda_gamma1),1.)+dnorm(gamma_2[i],0,1./sqrt(*lambda_gamma2),1.);


        Sgamma1+=(gamma_1[i]-*mu_gamma1)*(gamma_1[i]-*mu_gamma1);
        Sgamma2+=(gamma_2[i]-*mu_gamma2)*(gamma_2[i]-*mu_gamma2);
        nb_gamma_diff++;
        post_prob[i]+=1./(*B);
      }
    }

    *p=slice_sampling_p(*p, 0.05, 100, nb_gamma_same, nb_gamma_diff);



      /** Up-date lambda_gamma1 and lambda_gamma2 **/
      /** These are the original priors **/
    *lambda_gamma1=rgamma(1.+nb_gamma_diff/2.,1./(0.005+Sgamma1/2.));
    *lambda_gamma2=rgamma(1.+nb_gamma_diff/2.,1./(0.005+Sgamma2/2.));
    *lambda_gamma=rgamma(1.+nb_gamma_same/2.,1./(0.005+Sgamma/2.));

      /** Test the sensibilty to the priors **/
      //*lambda_gamma1=rgamma(1.+nb_gamma_diff/2.,1./(0.5+Sgamma1/2.));
      //*lambda_gamma2=rgamma(1.+nb_gamma_diff/2.,1./(0.5+Sgamma2/2.));
      //*lambda_gamma=rgamma(1.+nb_gamma_same/2.,1./(0.5+Sgamma/2.));




      /** Up-date the epsilon **/
    for(i=0;i<*n1;i++)
    {


      SSeps_sd=0.;
      SSeps_sd_new=0.;

      SSR1=0;SSR2=0;SS12=0;


    /** First part of the loop before dye swap **/
      for(j=0;j<(*nb_col1);j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j]);
        SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
        SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);


      }


    /** Second part of the loop after dye swap **/
      for(j=(*nb_col1);j<*n2;j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]);
        SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
        SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);


      }

    /** Gibbs sampling if lambda_eps' are independent **/
      lambda_eps1_new=rgamma(*a_eps1*(*a_eps1)/(*b_eps1)+*n2/2.,1./(*a_eps1/(*b_eps1)+SSR1/2.));

      lambda_eps2_new=rgamma(*a_eps2*(*a_eps2)/(*b_eps2)+*n2/2.,1./(*a_eps2/(*b_eps2)+SSR2/2.));


    /** up-date lambda_epsilon **/
      l_epsilon=0.,l_epsilon_new=0.;
      SSeps_sd=0.;
      SSeps_sd_new=0.;

      SSR1=0;SSR2=0;SS12=0;


    /** First part of the loop before dye swap **/
      for(j=0;j<(*nb_col1);j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j]);
        SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
        SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);

        SSeps_sd+=log(w[j*(*n1)+i])+log(lambda_eps1[i])/2.+log(lambda_eps2[i])/2.;
        SSeps_sd_new+=log(w[j*(*n1)+i])+log(lambda_eps1_new)/2.+log(lambda_eps2_new)/2.;

      }


    /** Second part of the loop after dye swap **/
      for(j=(*nb_col1);j<*n2;j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]);
        SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
        SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);

        SSeps_sd+=log(w[j*(*n1)+i])+log(lambda_eps1[i])/2.+log(lambda_eps2[i])/2.;
        SSeps_sd_new+=log(w[j*(*n1)+i])+log(lambda_eps1_new)/2.+log(lambda_eps2_new)/2.;

      }
    /** Store another sum (likelihood) to do Metropolis **/
      l_epsilon+=lambda_eps1[i]*SSR1-2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*SS12+lambda_eps2[i]*SSR2;
      l_epsilon_new+=lambda_eps1_new*SSR1-2*sqrt((lambda_eps1_new)*(lambda_eps2_new))*(*rho)*SS12+lambda_eps2_new*SSR2;



    /** Up-date epsilon1 and epsilon2 by metropolis **/
      dens_epsilon=-1./2*l_epsilon/(1-*rho*(*rho))+SSeps_sd
        +dgamma(lambda_eps1[i],*a_eps1*(*a_eps1)/(*b_eps1),(*b_eps1)/(*a_eps1),1.)+dgamma(lambda_eps2[i],*a_eps2*(*a_eps2)/(*b_eps2),(*b_eps2)/(*a_eps2),1.);
      dens_epsilon_new=-1./2*l_epsilon_new/(1-*rho*(*rho))+SSeps_sd_new
        +dgamma(lambda_eps1_new,*a_eps1*(*a_eps1)/(*b_eps1),(*b_eps1)/(*a_eps1),1.)+dgamma(lambda_eps2_new,*a_eps2*(*a_eps2)/(*b_eps2),(*b_eps2)/(*a_eps2),1.);


    /** Accept the up-dated variables **/
      if((dens_epsilon_new-dens_epsilon
        +dgamma(lambda_eps1[i],*a_eps1*(*a_eps1)/(*b_eps1)+*n2/2.,1./(*a_eps1/(*b_eps1)+SSR1/2.),1.)
        +dgamma(lambda_eps2[i],*a_eps2*(*a_eps2)/(*b_eps2)+*n2/2.,1./(*a_eps2/(*b_eps2)+SSR2/2.),1.)
        -dgamma(lambda_eps1_new,*a_eps1*(*a_eps1)/(*b_eps1)+*n2/2.,1./(*a_eps1/(*b_eps1)+SSR1/2.),1.)
        -dgamma(lambda_eps2_new,*a_eps2*(*a_eps2)/(*b_eps2)+*n2/2.,1./(*a_eps2/(*b_eps2)+SSR2/2.),1.))
        >log(runif(0,1)))
      {
        lambda_eps1[i]=lambda_eps1_new;
        lambda_eps2[i]=lambda_eps2_new;

        counteps++;
      }
    }


    sum_lambda1=0;
    sum_log_lambda1=0;
    sum_lambda2=0;
    sum_log_lambda2=0;

    for(i=0;i<*n1;i++)
    {
      sum_lambda1+=lambda_eps1[i];
      sum_lambda2+=lambda_eps2[i];
      sum_log_lambda1+=log(lambda_eps1[i]);
      sum_log_lambda2+=log(lambda_eps2[i]);
    }


    *a_eps1=slice_sampling_a(*a_eps1, 2., 10, sum_log_lambda1, sum_lambda1, *b_eps1, *n1);       
    *b_eps1=slice_sampling_b(*b_eps1, 5., 10, sum_log_lambda1, sum_lambda1, *a_eps1, *n1); 
    *a_eps2=slice_sampling_a(*a_eps2, 2., 10, sum_log_lambda2, sum_lambda2, *b_eps2, *n1);       
    *b_eps2=slice_sampling_b(*b_eps2, 5., 10, sum_log_lambda2, sum_lambda2, *a_eps2, *n1); 



    SSR1=0;SSR2=0;SS12=0;	  
      /** Up-date rho **/
    for(i=0;i<*n1;i++)
    {
    /** First part of the loop before dye swap **/
      for(j=0;j<(*nb_col1);j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=lambda_eps1[i]*w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j]);
        SSR2+=lambda_eps2[i]*w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
        SS12+=sqrt(lambda_eps1[i]*lambda_eps2[i])*w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);

      }

    /** Second part of the loop after dye swap **/
      for(j=(*nb_col1);j<*n2;j++)
      {
        /** General sum of squares used in the likelihood **/
        SSR1+=lambda_eps1[i]*w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]);
        SSR2+=lambda_eps2[i]*w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
        SS12+=sqrt(lambda_eps1[i]*lambda_eps2[i])*w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);

      }

    }
    *rho=slice_sampling_rho(*rho, 0.1, 10, SSR1, SSR2, SS12, *n1**n2);           
      //*rho=0;  
      /** Up-date the weight for the t-distribution **/

    for(j=0;j<*n2;j++)
    {

      /** Up-date the degrees of freedom **/
      df_new[j]=df_choice[(int)(runif(0.,1.)*(*nb_df))];  
      sum_l_w[j]=0.;
      sum_l_w_new[j]=0.;
    }



    for(i=0;i<(*n1);i++)
    { 
      for(j=0;j<(*nb_col1);j++)
      {

        sum_l_w[j]+=log(2./df[j])-(df[j]/2.+1)*log(1+1/((1-*rho*(*rho)))
          *(lambda_eps1[i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j])
          -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)
          *(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])
          +lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]))/df[j])
          +lgammafn(df[j]/2.+1)-lgammafn(df[j]/2.);

        sum_l_w_new[j]+=log(2./df_new[j])-(df_new[j]/2.+1)*log(1+1/((1-*rho*(*rho)))*
          (lambda_eps1[i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j])
          -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-gamma_1[i]-eta[j])
          *(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])
          +lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]))/df_new[j])
          +lgammafn(df_new[j]/2.+1)-lgammafn(df_new[j]/2.);

      }

      for(j=(*nb_col1);j<*n2;j++)
      {

        sum_l_w[j]+=log(2./df[j])-(df[j]/2.+1)*log(1+1/((1-*rho*(*rho)))*
          (lambda_eps1[i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])
          -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])
          +lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]))/df[j])+lgammafn(df[j]/2.+1)-lgammafn(df[j]/2.);


        sum_l_w_new[j]+=log(2./df_new[j])-(df_new[j]/2.+1)*log(1+1/((1-*rho*(*rho)))*
          (lambda_eps1[i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])
          -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])
          +lambda_eps2[i]*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j]))/df_new[j])+lgammafn(df_new[j]/2.+1)-lgammafn(df_new[j]/2.);
      }
    }


  /** Up-date the degrees of freedom for the t-distribution **/
    if(*robust==1)
    {
      for(j=0;j<(*nb_col1);j++)
      {
        if((sum_l_w_new[j]-sum_l_w[j])>log(runif(0,1)))
        {
          df[j]=df_new[j];
          countdf++;
        }
        for(i=0;i<*n1;i++)
          w[j*(*n1)+i]=rgamma(df[j]/2.+1,1./(df[j]/2.+1/(2*(1-*rho*(*rho)))*
          (lambda_eps1[i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j])
          -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])
          +lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]))));
      }

      for(j=(*nb_col1);j<*n2;j++)
      {
        if((sum_l_w_new[j]-sum_l_w[j])>log(runif(0,1)))
        {
          df[j]=df_new[j];
          countdf++;
        }
        for(i=0;i<*n1;i++)
          w[j*(*n1)+i]=rgamma(df[j]/2.+1,1./(df[j]/2.+1/(2*(1-*rho*(*rho)))*
          (lambda_eps1[i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])
          -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])
          +lambda_eps2[i]*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j]))));
      }
    }

    if(k>=*min_iter)
    {
      if(((count2+1)%(*batch))==0) /** Batch sampling **/
      {

        if(*all_out==1)
        {
          for(i=0;i<*n1;i++)
          {
            for(j=0;j<*n2;j++)
            {
          /** Compute the mean of the weights **/
              w_p[j*(*n1)+i]+=w[j*(*n1)+i]/((*B-*min_iter)/(*batch));		  
              df_p[count*(*n2)+j]=df[j];
            }
            gamma1_p[count*(*n1)+i]=gamma_1[i];
            gamma2_p[count*(*n1)+i]=gamma_2[i];
            lambda_eps1_p[count*(*n1)+i]=lambda_eps1[i];
            lambda_eps2_p[count*(*n1)+i]=lambda_eps2[i];
          }
        }
        else /* posterior mean **/
        {
          for(i=0;i<*n1;i++)          
          {
            for(j=0;j<*n2;j++)
            {
              /** Compute the mean of the weights **/
              w_p[j*(*n1)+i]+=w[j*(*n1)+i]/((*B-*min_iter)/(*batch));		  
              df_p[j]+=df[j]/((*B-*min_iter)/(*batch));
            }
            lambda_eps1_p[i]+=lambda_eps1[i]/((*B-*min_iter)/(*batch));
            lambda_eps2_p[i]+=lambda_eps2[i]/((*B-*min_iter)/(*batch));

            gamma1_p[i]+=(gamma_1[i])/((*B-*min_iter)/(*batch));
            gamma2_p[i]+=(gamma_2[i])/((*B-*min_iter)/(*batch));
          }
        }


        if(*all_out==1)
        {
          a_eps1_p[count]=*a_eps1;
          b_eps1_p[count]=*b_eps1;
          a_eps2_p[count]=*a_eps2;
          b_eps2_p[count]=*b_eps2;

          lambda_gamma1_p[count]=*lambda_gamma1;
          lambda_gamma2_p[count]=*lambda_gamma2;
          lambda_gamma_p[count]=*lambda_gamma;
          p_p[count]=*p;
          rho_p[count]=*rho;
        }
        else /** posterior mean **/
        {
          *a_eps1_p+=*a_eps1/((*B-*min_iter)/(*batch));
          *b_eps1_p+=*b_eps1/((*B-*min_iter)/(*batch));
          *a_eps2_p+=*a_eps2/((*B-*min_iter)/(*batch));
          *b_eps2_p+=*b_eps2/((*B-*min_iter)/(*batch));

          *lambda_gamma1_p+=*lambda_gamma1/((*B-*min_iter)/(*batch));
          *lambda_gamma2_p+=*lambda_gamma2/((*B-*min_iter)/(*batch));
          *lambda_gamma_p=+*lambda_gamma/((*B-*min_iter)/(*batch));
          *p_p+=*p/((*B-*min_iter)/(*batch));
          *rho_p=*rho/((*B-*min_iter)/(*batch));
        }


        count++;
      }
      count2++;
    }



  }
  //printf("\nSummary on the up-dates\n");
  //printf("\nThe number of acceptance for df is %g\n",(double)countdf/((double)((*B)*(*n2))));
  //printf("\nThe number of acceptance for gamma is %g\n",(double)countgamma/((double)((*B)*(*n1))));

  Free(df_new);
  Free(sum_l_w);
  Free(sum_l_w_new);

  Free(mu_gamma1);
  Free(mu_gamma2);
  Free(mu_gamma);
}

