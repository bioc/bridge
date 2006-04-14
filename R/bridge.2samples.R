"bridge.2samples" <-
  function(sample1,sample2,B=1000,min.iter=0,batch=10,mcmc.obj=NULL,all.out=TRUE,affy=FALSE,verbose=FALSE,log=FALSE)
{
###  Only take the finite observations
  n<-dim(sample1) 

  if(!log)
    {
      sample1<-log(sample1)
      sample2<-log(sample2)
    }
  
  vec1<-as.double(t(sample1))
  vec1[is.finite(vec1)==FALSE]<- -9999999
  vec2<-as.double(t(sample2))
  vec2[is.finite(vec2)==FALSE]<- -9999999
    
  nu.choice<-c(1:10,seq(20,100,10))
  
  w.out<-rep(0,n[1]*n[2])
  nu.out<-double(n[2]*(B-min.iter)/batch)           
  a.gamma<-1
  b.gamma<-0.005


  
  if(length(mcmc.obj)>0)
    {
      if(class(mcmc.obj)!="bridge2")
        stop("'mcmc.obj' should be of type 'bridge2'" , call. = TRUE)
      
      
      n.iter<-length(mcmc.obj$gamma1[1,])
      
                                        
      lambda.eps1<-mcmc.obj$lambda.eps1[,n.iter]
      lambda.eps2<-mcmc.obj$lambda.eps2[,n.iter]
      lambda.gamma1<-mcmc.obj$lambda.gamma1[n.iter]
      lambda.gamma2<-mcmc.obj$lambda.gamma2[n.iter]
      lambda.gamma<-mcmc.obj$lambda.gamma[n.iter]
      p<-mcmc.obj$p[n.iter]
      
      a.eps1<-mcmc.obj$a.eps1[n.iter]
      b.eps1<-mcmc.obj$b.eps1[n.iter]
      a.eps2<-mcmc.obj$a.eps2[n.iter]
      b.eps2<-mcmc.obj$b.eps2[n.iter]
      
      
      rho<-mcmc.obj$rho[n.iter]
      
      nu.in<-mcmc.obj$nu[,n.iter]
      w.in<-matrix(mcmc.obj$w)
      w.in<-as.double(t(w.in))
      
    }
  else
    {
      
      gamma1<-mat.mean(sample1)[,1]
      gamma2<-mat.mean(sample2)[,1]

      lambda.eps1<-1/(mat.mean(sample1)[,2]+0.0001)^2
      lambda.eps2<-1/(mat.mean(sample2)[,2]+0.0001)^2
      
      lambda.gamma1<-1/var(gamma1)
      lambda.gamma2<-1/var(gamma2)
      lambda.gamma<-0.5
      p<-0.05
      
      a.eps1<-median(lambda.eps1)
      b.eps1<-mad(lambda.eps1)
      a.eps2<-median(lambda.eps2)
      b.eps2<-mad(lambda.eps2)
      nu.in<-rep(10,n[2])
      w.in<-rep(1,n[1]*n[2])

      rho<-0
      
    }
  
  
### Main code linked to a c function
  if(all.out==TRUE)
    length<-(B-min.iter)/batch
  else
    length<-1

  if(affy)
    {
      G<-dim(sample1)[1]
      R1<-dim(sample1)[2]
      R2<-dim(sample2)[2]
      model<-rep(1,G)
      w.mix=c(.5,.5)
      
      obj<-.C("gene_express_2s",
              as.double(vec1),
              as.integer(R1),
              as.double(vec2),
              as.integer(R2),
              as.integer(G),
              as.double(gamma1),
              as.double(gamma2),
              gamma1=double(G*length),
              gamma2=double(G*length),
              as.integer(model),
              model=integer(G*(B-min.iter)/batch),
              post.p=as.double(rep(0,G)),
              as.double(lambda.gamma1),
              lambda.gamma1=double(length),
              as.double(lambda.gamma2),
              lambda.gamma2=double(length),
              as.double(lambda.gamma),
              lambda.gamma=double(length),
              as.double(lambda.eps1),
              lambda.eps1=double(G*length),
              as.double(lambda.eps2),
              lambda.eps2=double(G*length),
              as.double(a.eps1),
              a.eps1=double(length),
              as.double(b.eps1),
              b.eps1=double(length),
              as.double(a.eps2),
              a.eps2=double(length),
              as.double(b.eps2),
              b.eps2=double(length),
              move=double(G),
              w1=as.double(rep(0.001,R1*G)),
              w2=as.double(rep(0.001,R2*G)),
              as.double(rep(100,R1)),
              nu1=double(R1*length),
              as.double(rep(100,R2)),
              nu2=double(R2*length),
              as.double(nu.choice),
              as.integer(length(nu.choice)),
              as.double(w.mix),
              w.mix=double(2*length),
              as.integer(min.iter),
              as.integer(batch),
              as.integer(B),
              as.integer(all.out),
              as.integer(verbose),
              PACKAGE="bridge")
    }
  else
    obj<-.C("ex_R_link_bridge2",
            as.double(vec1),
            as.double(vec2),
            as.integer(n[1]),
            as.integer(n[2]),
            as.integer(n[2]/2),
            as.integer(B),
            as.integer(0),
            as.double(gamma1),
            gamma1=double(length*n[1]),
            as.double(gamma2),
            gamma2=double(length*n[1]),
            as.double(p),
            p=double(length),
            as.double(0),
            as.double(0),
            as.double(0),
            as.double(0),
            as.double(rep(0,n[2])),
            as.double(lambda.eps1),
            lambda.eps1=double(length*n[1]),
            as.double(lambda.eps2),
            lambda.eps2=double(length*n[1]),
            as.double(a.eps1),
            as.double(b.eps1),
            a.eps1=double(length),
            b.eps1=double(length),
            as.double(a.eps2),
            as.double(b.eps2),
            a.eps2=double(length),
            b.eps2=double(length),
            as.double(w.in),
            as.double(nu.choice),
            as.integer(length(nu.choice)),
            as.double(nu.in),
            nu=double(n[2]*length),
            w=as.double(w.out),
            as.double(rho),
            rho=double(length),
            as.double(0),
            as.double(lambda.gamma1),
            as.double(lambda.gamma2),
            lambda.gamma1=double(length),
            lambda.gamma2=double(length),
            as.double(lambda.gamma),
            lambda.gamma=double(length),
            post.p=double(n[1]),
            move=double(n[1]), 
            as.integer(min.iter),
            as.integer(batch), as.integer(all.out),
            as.integer(verbose),
            PACKAGE="bridge")
  
### Create a new object
  if(affy)
    {
      gamma1<-t(matrix(obj$gamma1,(length),G,byrow=TRUE))
      gamma2<-t(matrix(obj$gamma2,(length),G,byrow=TRUE))
      lambda.eps1<-t(matrix(obj$lambda.eps1,(length),G,byrow=TRUE))
      lambda.eps2<-t(matrix(obj$lambda.eps2,(length),G,byrow=TRUE))
      w.mix<-t(matrix(obj$w.mix,(length),2,byrow=TRUE))
      w1<-t(matrix(obj$w1,R1,G,byrow=TRUE))
      w2<-t(matrix(obj$w2,R2,G,byrow=TRUE))
      nu1<-t(matrix(obj$nu1,(length),R1,byrow=TRUE))
      nu2<-t(matrix(obj$nu2,(length),R2,byrow=TRUE))
      
      new.mcmc<-list(gamma1=gamma1, gamma2=gamma2, model=model, lambda.eps1=lambda.eps1, lambda.eps2=lambda.eps2,
                     a.eps1=obj$a.eps1, b.eps1=obj$b.eps1,a.eps2=obj$a.eps2, b.eps2=obj$b.eps2,
                     w.mix=w.mix,
                     w1=w1, w2=w2, 
                     nu1=nu1, nu2=nu2, 
                     lambda.gamma1=obj$lambda.gamma1,
                     lambda.gamma2=obj$lambda.gamma2,
                     lambda.gamma=obj$lambda.gamma,
                     move=obj$move/B,post.p=obj$post.p/B)
      
    }  
  else
    {
      new.mcmc<-list(gamma1=t(matrix(obj$gamma1,(B-min.iter)/batch,n[1],byrow=TRUE)),
                     gamma2=t(matrix(obj$gamma2,(B-min.iter)/batch,n[1],byrow=TRUE)),
                     lambda.eps1=t(matrix(obj$lambda.eps1,(B-min.iter)/batch,n[1],byrow=TRUE)),
                     lambda.eps2=t(matrix(obj$lambda.eps2,(B-min.iter)/batch,n[1],byrow=TRUE)),
                     lambda.gamma1=obj$lambda.gamma1,lambda.gamma2=obj$lambda.gamma2,
                     rho=obj$rho, lambda.gamma=obj$lambda.gamma, 
                     w1=matrix(obj$w,n[1],n[2]),
                     w2=obj$w,
                     nu1=t(matrix(obj$nu,(B-min.iter)/batch,n[2],byrow=TRUE)),
                     nu2=obj$nu,
                     a.eps1=obj$a.eps1,b.eps1=obj$b.eps1,a.eps2=obj$a.eps2,
                     b.eps2=obj$b.eps2, p=obj$p,post.p=obj$post.p,move=obj$move)
    }
  
### Give it the right class
  
  class(new.mcmc)<-"bridge2"
  
  return(new.mcmc)
}
