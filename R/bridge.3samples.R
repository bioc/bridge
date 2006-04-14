"bridge.3samples" <-
  function(sample1,sample2,sample3,B=1000,min.iter=0,batch=10,mcmc.obj=NULL,all.out=TRUE,verbose=FALSE,log=FALSE)
{
###  Only take the finite observations
  if(!log)
    {
      sample1<-log(sample1)
      sample2<-log(sample2)
      sample3<-log(sample3)
    }
  
  G<-dim(sample1)[1]
  R1<-dim(sample1)[2]
  R2<-dim(sample2)[2]
  R3<-dim(sample3)[2]


  if(length(mcmc.obj)>0)
    {
      if(class(mcmc.obj)!="bridge3")
        stop("'mcmc.obj' should be of type 'bridge3'" , call. = TRUE)
      
      n.iter<-length(mcmc.obj$a.eps)
      gamma1<-mcmc.obj$gamma1[,n.iter]
      gamma2<-mcmc.obj$gamma2[,n.iter]
      gamma3<-mcmc.obj$gamma3[,n.iter]
      
      lambda1<-mcmc.obj$lambda1[,n.iter]
      lambda2<-mcmc.obj$lambda2[,n.iter]
      lambda3<-mcmc.obj$lambda3[,n.iter]
      a.eps1<-mcmc.obj$a.eps1[n.iter]
      b.eps1<-mcmc.obj$b.eps1[n.iter]
      a.eps2<-mcmc.obj$a.eps2[n.iter]
      b.eps2<-mcmc.obj$b.eps2[n.iter]
      a.eps3<-mcmc.obj$a.eps3[n.iter]
      b.eps3<-mcmc.obj$b.eps3[n.iter]
      
      
      lambda.gamma1<-mcmc.obj$lambda.gamma1[n.iter]
      lambda.gamma2<-mcmc.obj$lambda.gamma2[n.iter]
      lambda.gamma3<-mcmc.obj$lambda.gamma3[n.iter]
      lambda.gamma12<-mcmc.obj$lambda.gamma12[n.iter]
      lambda.gamma23<-mcmc.obj$lambda.gamma23[n.iter]
      lambda.gamma13<-mcmc.obj$lambda.gamma13[n.iter]
      lambda.gamma123<-mcmc.obj$lambda.gamma123[n.iter]
      
      w.mix<-mcmc.obj$w.mix[,n.iter]
      model<-mcmc.obj$model[,n.iter]
      
    }
  else # Least squares 
    {
      gamma1<-mat.mean(sample1)[,1]
      gamma2<-mat.mean(sample2)[,1]
      gamma3<-mat.mean(sample3)[,1]
      
      lambda1<-1./(mat.mean(sample1)+0.001)[,2]^2
      lambda2<-1./(mat.mean(sample2)+0.001)[,2]^2
      lambda3<-1./(mat.mean(sample3)+0.001)[,2]^2
      a.eps1<-median(lambda1)
      b.eps1<-mad(lambda1)
      a.eps2<-median(lambda2)
      b.eps2<-mad(lambda2)
      a.eps3<-median(lambda3)
      b.eps3<-mad(lambda3)
      
      
      lambda.gamma1<-1./10
      lambda.gamma2<-1./10
      lambda.gamma3<-1./10
      lambda.gamma12<-1./10
      lambda.gamma23<-1./10
      lambda.gamma13<-1./10
      lambda.gamma123<-1./10
      
      w.mix<-rep(.2,5)
      model<-rep(4,G)
      
    }
  
  nu.choice<-c(1:10,seq(20,100,10))
  nb.nu<-length(nu.choice)
  vec1<-as.double(t(sample1))
  vec1[is.finite(vec1)==FALSE]<- -9999999
  vec2<-as.double(t(sample2))
  vec2[is.finite(vec2)==FALSE]<- -9999999
  vec3<-as.double(t(sample3))
  vec3[is.finite(vec3)==FALSE]<- -9999999
  

### Main code linked to a c function
  if(all.out==TRUE)
    length<-(B-min.iter)/batch
  else
    length<-1
  
  
  obj<-.C("gene_express_3s",
          as.double(t(sample1)),
          as.integer(R1),
          as.double(t(sample2)),
          as.integer(R2),
          as.double(t(sample3)),
          as.integer(R3),           
          as.integer(G),
          as.double(gamma1),
          as.double(gamma2),
          as.double(gamma3),
          gamma1=double(G*length),
          gamma2=double(G*length),
          gamma3=double(G*length),
          as.integer(model),
          model=integer(G*(B-min.iter)/batch),
          as.double(lambda.gamma1),
          lambda.gamma1=double(length),
          as.double(lambda.gamma2),
          lambda.gamma2=double(length),
          as.double(lambda.gamma3),
          lambda.gamma3=double(length),
          as.double(lambda.gamma12),
          lambda.gamma12=double(length),
          as.double(lambda.gamma23),
          lambda.gamma23=double(length),
          as.double(lambda.gamma13),
          lambda.gamma13=double(length),
          as.double(lambda.gamma123),
          lambda.gamma123=double(length),
          as.double(lambda1),
          lambda1=double(G*length),
          as.double(lambda2),
          lambda2=double(G*length),
          as.double(lambda3),
          lambda3=double(G*length),
          as.double(a.eps1),
          a.eps1=double(length),
          as.double(b.eps1),
          b.eps1=double(length),
          as.double(a.eps2),
          a.eps2=double(length),
          as.double(b.eps2),
          b.eps2=double(length),
          as.double(a.eps3),
          a.eps3=double(length),
          as.double(b.eps3),
          b.eps3=double(length),
          move=double(G),
          w1=as.double(rep(0.001,R1*G)),
          w2=as.double(rep(0.001,R2*G)),
          w3=as.double(rep(0.001,R3*G)),
          as.double(rep(100,R1)),
          nu1=double(R1*length),
          as.double(rep(100,R2)),
          nu2=double(R2*length),
          as.double(rep(100,R3)),
          nu3=double(R3*length),
          as.double(nu.choice),
          as.double(w.mix),
          w.mix=double(5*length),
          as.integer(min.iter),
          as.integer(batch),
          as.integer(B),
          as.integer(all.out),
          as.integer(verbose),
          PACKAGE="bridge")
  
  prop.model<-matrix(0,G,5)
  model<-t(matrix(obj$model,((B-min.iter)/batch),G,byrow=TRUE))
  gamma1<-t(matrix(obj$gamma1,(length),G,byrow=TRUE))
  gamma2<-t(matrix(obj$gamma2,(length),G,byrow=TRUE))
  gamma3<-t(matrix(obj$gamma3,(length),G,byrow=TRUE))
  lambda1<-t(matrix(obj$lambda1,(length),G,byrow=TRUE))
  lambda2<-t(matrix(obj$lambda2,(length),G,byrow=TRUE))
  lambda3<-t(matrix(obj$lambda3,(length),G,byrow=TRUE))
  w.mix<-matrix(obj$w.mix,(length),5,byrow=TRUE)
  w1<-t(matrix(obj$w1,R1,G,byrow=TRUE))
  w2<-t(matrix(obj$w2,R2,G,byrow=TRUE))
  w3<-t(matrix(obj$w3,R3,G,byrow=TRUE))
  nu1<-t(matrix(obj$nu1,(length),R1,byrow=TRUE))
  nu2<-t(matrix(obj$nu2,(length),R2,byrow=TRUE))
  nu3<-t(matrix(obj$nu3,(length),R3,byrow=TRUE))
  
  for(g in 1:G)
    {
      prop.model[g,]<-c(sum(model[g,]==0)/((B-min.iter)/batch),sum(model[g,]==1)/((B-min.iter)/batch),
                        sum(model[g,]==2)/((B-min.iter)/batch),sum(model[g,]==3)/((B-min.iter)/batch),
                        sum(model[g,]==4)/((B-min.iter)/batch))
    }
  
  new.mcmc<-list(gamma1=gamma1, gamma2=gamma2, gamma3=gamma3, model=model, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3,
                 a.eps1=obj$a.eps1, b.eps1=obj$b.eps1,
                 a.eps2=obj$a.eps2, b.eps2=obj$b.eps2,
                 a.eps3=obj$a.eps3, b.eps3=obj$b.eps3,
                 w.mix=w.mix,
                 w1=w1, w2=w2, w3=w3,
                 nu1=nu1,nu2=nu2,nu3=nu3,
                 lambda.gamma1=obj$lambda.gamma1, lambda.gamma2=obj$lambda.gamma2, lambda.gamma3=obj$lambda.gamma3, lambda.gamma12=obj$lambda.gamma12,
                 lambda.gamma13=obj$lambda.gamma13, lambda.gamma23=obj$lambda.gamma23, lambda.gamma123=obj$lambda.gamma123, 
                 prop.model=prop.model,move=obj$move/B)
  
  class(new.mcmc)<-"bridge3"
  
  return(new.mcmc)
  
}
