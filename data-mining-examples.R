
userDA=function(){
  
  data=readline("Enter the data file name : " )  
  cat("Select the data coding format(1 = 'a b c' or 2 = 'a,b,c'):  ") 
  fm = scan(n=1, quiet=TRUE)    
  if(fm==1) {form = ""} else {form = ","}    
  cat("Does the data have variable name?(1 = No or 2 = Yes):  ")
  aw = scan(n=1, quiet=TRUE)    
  if(aw==1) {head = FALSE} else {head = TRUE} 
  data=read.table(data, header=head, sep=form)
  
  # location of class variable and categorical variables 
  
  cat("where is the class variable ? (write down the order):  ")
  resps = scan(n=1,quiet=TRUE)
  data=cbind(y=data[,resps],data[,-c(resps)] ) ;  K=nlevels(factor(data[,1]))
  n=nrow(data)  ; p=ncol(data)
  cat("where is the categorical variable ? (write down the order with comma if many or just pass:  ")  
  catgv= scan(,quiet=TRUE)
  if(length(catgv)==0)  {data=data} else {for (i in 1:length(catgv))
    data[,catgv[i]]=factor(data[,catgv[i]],labels=1:nlevels(data[,catgv[i]]))}
  
  
  method=readline("Which method do you apply for classification ?((LDA,QDA,CLDA, Logistic, Tree, BaggingEsb) :")
  etool=readline("Which evaluation tool do you apply ?(CV or testset) :")
  if(etool=="CV"){nfold=readline("#-fold CV- enter #: ") } else
  {testset1=readline("Enter the dataset name :   ")
  testset=read.table(testset1,header=head, sep=form)
  testset=cbind(testset[,resps],testset[,-resps]) }
  
  ######### Run classification #########
  #LDA, QDA, CLDA, Logit, Tree, BaggingEsb
  
  k=c()
  for(i in 1:K){k=c(k,nrow(data[data[,1]==i,]))}
  
  
  da2=function(fitdat,testdat,kclass){
    mu=c(); sig=c() ; num=c()
    for(i in 1:K){
      mu1=c()
      for(j in 2:p){
        mu1=c(mu1,mean(fitdat[fitdat[,1]==i,][,j]))     }
      mu=cbind(mu,mu1)
      sig=rbind(sig,cov(fitdat[fitdat[,1]==i,][,-1]))
      num=c(num,nrow(fitdat[fitdat[,1]==i,]))  }  
    return(list(mu,sig,num))}
  
  ###### LDA method#######
  if(method=="LDA"){
    da=function(fitdat,testdat,kclass) {
      dda=da2(fitdat,testdat,kclass)
      mu=dda[[1]] ; sig=dda[[2]] ; num=dda[[3]]
      sig.s=0
      for(i in 1:K){
        sig.s=sig.s+(num[i]-1)*sig[((p-1)*(i-1)+1):((p-1)*i),]}
      sp=1/(sum(num)-K)*sig.s
      resul.lda=c()
      for(i in 1:K){
        mu1=as.matrix(mu[,i])
        lda.d=t(mu1)%*%solve(sp)%*%t(testdat[,-1])-1/2*as.numeric(t(mu1)%*%solve(sp)%*%mu1)
        resul.lda=cbind(resul.lda,t(lda.d))  }
      fit.lda=c()
      for(i in 1:nrow(testdat)){
        fit.lda=rbind(fit.lda,resul.lda[i,]==max(resul.lda[i,]))    }
      cbind(kclass,fit.lda%*%matrix(1:K,K,1))}  }            
  
  ###### QDA method# #######
  if(method=="QDA"){
    da=function(fitdat,testdat,kclass){
      dda=da2(fitdat,testdat,kclass)
      mu=dda[[1]] ; sig=dda[[2]] ; num=dda[[3]]
      resul.qda=c()
      for(i in 1:K){
        mu1=mu[,i]
        for(j in 2:nrow(testdat)){
          mu1=rbind(mu1,mu[,i])    }
        sig1=as.matrix(sig[((p-1)*(i-1)+1):((p-1)*i),])
        qda.d=-1/2*log(det(sig1))-0.5*as.matrix((testdat[,-1]-mu1))%*%solve(sig1)%*%t(as.matrix(testdat[,-1]-mu1))
        resul.qda=cbind(resul.qda,diag(qda.d))}
      fit.qda=c()
      for(i in 1:nrow(testdat)){
        fit.qda=rbind(fit.qda,resul.qda[i,]==max(resul.qda[i,])) }
      cbind(kclass,fit.qda%*%matrix(1:K,K,1))    }  }   
  
  ###### Canonical LDA #######
  if(method=="CLDA"){ 
    da=function(fitdat,testdat,kclass){
      dda=da2(fitdat,testdat,kclass)
      mu=dda[[1]] ; sig=dda[[2]] ; num=dda[[3]]
      sig.s=0
      for(i in 1:K){
        sig.s=sig.s+(num[i]-1)*sig[((p-1)*(i-1)+1):((p-1)*i),]   }
      sp=1/(sum(num)-K)*sig.s
      
      ev=eigen(solve(sig.s)%*%(cov(t(mu))*3))$values
      ev=order(Re(as.matrix(Im(ev)==0)*as.matrix(ev)),decreasing=T)[1:min(K-1,19)]
      e1=eigen(solve(sig.s)%*%(cov(t(mu))*3))$vectors
      for(i in 1:min(K-1,19)){
        e1[,ev[i]]=1/sqrt(t(e1[,ev[i]]%*%sp%*%e1[,ev[i]]))*e1[,ev[i]]     }
      resul.clda=c()
      for(i in 1:K){
        mu1=mu[,i]
        for(j in 2:nrow(testdat)){
          mu1=rbind(mu1,mu[,i])    }
        clda.d=c(0)
        for(j in 1:min(K-1,19)){
          clda.d=clda.d+Re((as.matrix((testdat[,-1]-mu1))%*%as.matrix(e1[,j]))^2)    }
        resul.clda=cbind(resul.clda,clda.d) }
      fit.clda=c()
      for(i in 1:nrow(testdat)){
        fit.clda=rbind(fit.clda,resul.clda[i,]==min(resul.clda[i,]))  }
      cbind(kclass,fit.clda%*%matrix(1:K,K,1))   }   } 
  
  ###### Logit regression #######
  if(method=="Logistic") {
    da=function(fitdat,testdat,kclass){		
      y=as.matrix(fitdat[,1]-1) ; K=nlevels(fitdat[,1])
      x=as.matrix(cbind(1,fitdat[,-1]))
      betac=matrix(1,nrow=dim(x)[2],K-1) ;ebeta=matrix(1,nrow=dim(testdat)[1],K-1)
      for ( i in 1:K-1) {
        yn=as.matrix(y[y==0|y==i,])
        xn=rbind(x[y[,1]==0,],x[y[,1]==i,])
        loglik=function(beta){		
          -(t(yn)%*%xn%*%beta+sum(log(1/(1+exp(xn%*%beta)))) ) }
        logit=function(beta){
          sum(abs(t(x)%*%as.matrix(1/(1+exp(-x%*%beta)))-t(x)%*%y))   }
        beta2=optim(rep(0.1,p),loglik,method ="BFGS")$par  
        betac[,i]=optim(beta2,logit,method ="BFGS")$par
        ebeta[,i]=exp(as.matrix(cbind(1,testdat[,-1]))%*%betac[,i]) } 
      theta1=matrix(1,nrow(testdat),K)
      for (i in 1:nrow(testdat)) {
        for ( j in 1:(K-1) ) {
          theta1[i,j]=ebeta[i,j]/(1+sum(ebeta[i,])) }}
      for (i in 1:nrow(testdat)) {
        theta1[i,K]=1/(1+sum(ebeta[i,])) }  
      class=1:K	
      cbind(kclass,as.matrix(class[apply(theta1,1,which.max)])) }    } 
  
  ##### 1-level tree#####
  if(method=="Tree"){
    da=function(fitdat,testdat,kclass){
      impt=1-(sum(fitdat[,1]==1)/n)^2-(1-sum(fitdat[,1]==1)/n)^2  # impurity measure 
      dvn=c()
      for(i in 1:(p-1)){dvn=c(dvn,nlevels(factor(fitdat[,i+1])))}
      gns=c()  # goodness of split
      gns2=c()
      for(i in 1:(p-1)){
        gns1=cbind(1:(dvn[i]-1),0)   # divide var by dvn[i]-1
        for(j in 1:(dvn[i]-1)){
          t1=fitdat[fitdat[,i+1]<=as.numeric(levels(factor(fitdat[,i+1])))[j],1]
          t12=sum(t1==2)/length(t1) ; t11=sum(t1==1)/length(t1)
          impt1=1-t12^2-t11^2
          t2=fitdat[fitdat[,i+1]>as.numeric(levels(factor(fitdat[,i+1])))[j],1]
          t22=sum(t2==2)/length(t2) ; t21=sum(t2==1)/length(t2)
          impt2=1-t22^2-t21^2
          impt3=(length(t1)/(length(t1)+length(t2)))*impt1+(length(t2)/(length(t1)+length(t2)))*impt2
          gns1[j,2]=impt-impt3      }
        gns=c(gns,as.numeric(levels(factor(fitdat[,i+1])))[gns1[gns1[,2]==max(gns1[,2]),1]])
        gns2=c(gns2,max(gns1[,2]))   }
      
      cbind(kclass,(testdat[,vs+1]<=gns[vs])+s1)  }         }
  
  #########BaggingEsb############
  if (method=="BaggingEsb"){
    da  =function( fitdat, testdat, kclass) {
      qda=function(fitdat,testdat,kclass){
        dda=da2(fitdat,testdat,kclass)
        mu=dda[[1]] ; sig=dda[[2]] ; num=dda[[3]]
        resul.qda=c()
        for(i in 1:K){
          mu1=mu[,i]
          for(j in 2:nrow(testdat)){
            mu1=rbind(mu1,mu[,i])    }
          sig1=as.matrix(sig[((p-1)*(i-1)+1):((p-1)*i),])
          qda.d=-1/2*log(det(sig1))-0.5*as.matrix((testdat[,-1]-mu1))%*%solve(sig1)%*%t(as.matrix(testdat[,-1]-mu1))
          resul.qda=cbind(resul.qda,diag(qda.d))}
        fit.qda=c()
        for(i in 1:nrow(testdat)){
          fit.qda=rbind(fit.qda,resul.qda[i,]==max(resul.qda[i,]))}
        cbind(kclass,fit.qda%*%matrix(1:K,K,1)) }  
      test.n=dim(testdat)[1]
      result <- matrix(,nrow=test.n,ncol=51)    
      for( i in 1:51 ){                                      
        train <- sample(x=1:n, size=n, replace=TRUE)           # Bootstrap
        data.train <- fitdat[train,]
        result[,i] <- qda(data.train,testdat,kclass)[,2]  }
      pred.y <- c()
      for(i in 1:test.n){
        y.hat <-  which.max(tabulate( result[i,] ))     
        if( length(y.hat) >= 2 ){      
          y.hat[i] <- sample(x=y.hat, size = 1 )      
        } else { pred.y[i] <- y.hat  } }
      pred.y <- as.numeric(pred.y)     
      cbind(kclass, as.matrix(pred.y)) }  
  } 
  
  #########etool= CV #########
  if(etool=="CV"){
    nfold=as.numeric(nfold) 
    num1=round((k/n)*(n/nfold))
    redat = replicate(nfold,c())  
    data1=cbind(data,1:n)         
    for(i in 1:K){    # Stratified n-fold cv data split to nfold
      kclass=c(1:k[i])
      redat1=cbind(kclass,data1[data1[,1]==i,])
      ind=as.matrix(sample(kclass,replace=F))
      colnames(ind)=c("kclass")
      redat1=merge(ind,redat1,by.x="kclass",sort = F)[,-1]
      for(j in 1:(nfold-1)){
        redat[[j]]=rbind(redat[[j]],as.matrix(redat1[(1+(j-1)*num1[i]):(j*num1[i]),]))       }
      redat[[nfold]]=rbind(redat[[nfold]],as.matrix(redat1[(1+(nfold-1)*num1[i]):nrow(redat1),])) }  
    fit.result=c()
    for(j in 1:nfold){      
      testdat=redat[[j]]
      fitdat=c()
      for(i in 1:nfold){
        if(i!=j)  fitdat=rbind(fitdat,redat[[i]])   }
      fit.result=rbind(fit.result, da(fitdat[,-(p+1)],testdat[,-(p+1)],testdat[,(p+1)] )) 	}  
    fitdat=rbind(fitdat,testdat)   
    ind=matrix(1:n,n,1)  ; colnames(ind)=c("kclass")
    fit.res=merge(ind,da(fitdat[,-(p+1)],fitdat[,-(p+1)],fitdat[,(p+1)]))[,2]  # resub
    fit.test=merge(ind,fit.result)[,2]   # testset
    testset=data
  }
  #########etool= testdata #########
  if(etool=="testset"){
    fit.res1=da(as.matrix(data),as.matrix(data),1:nrow(data))
    fit.test1=da(as.matrix(data),as.matrix(testset),1:nrow(testset))
    fit.res=fit.res1[,2]
    fit.test=fit.test1[,2]
  }
  
  
  ###### data output #####
  class.resul=table(data[,1],fit.res)  
  class.resul2=table(testset[,1],fit.test) 
  res=c() ; res2=c()
  if(K>2){
    for(i in 3:K){
      res=c(res,paste("           ",i,paste("      ",class.resul[i,1:K],collapse= ""),"\n"))
      res2=c(res2,paste("           ",i,paste("      ",class.resul2[i,1:K],collapse= ""),"\n"))
    }}
  result1=c(method,"\n","  ID, Actual class, Resub pred ","\n",
            " ----------------------------- ","\n",
            paste("",1:n,",",data[1:n,1],",",fit.res[1:n],"\n"))
  result2=c( "\n"," Confusion Matrix (Resubstitution)","\n",
             " ---------------------------------- ","\n",
             "                   Predicted Class","\n","             ", paste("      ",1:K,collapse= ""),"\n",
             " Actual    ",paste(1,paste("      ",class.resul[1,1:K],collapse= ""),"\n"),
             " Class     ",paste(2,paste("      ",class.resul[2,1:K],collapse= ""),"\n"),res,
             "\n"," Model Summary (Resubstitution)","\n"," ------------------------------","\n",
             " Overall accuracy =", round(sum(diag(class.resul))/sum(class.resul),digits=3),"\n")
  
  result3=c("\n"," ID, Actual class, Test pred ","\n",
            " ----------------------------- ","\n",
            paste("",1:nrow(testset),",",testset[1:nrow(testset),1],",",fit.test[1:nrow(testset)],"\n"))
  result4=c( "\n"," Confusion Matrix (",etool,")","\n",
             " ---------------------------------- ","\n",
             "                   Predicted Class","\n","             ", paste("      ",1:K,collapse= ""),"\n",
             " Actual    ",paste(1,paste("      ",class.resul2[1,1:K],collapse= ""),"\n"),
             " Class     ",paste(2,paste("      ",class.resul2[2,1:K],collapse= ""),"\n"),res2,
             "\n"," Model Summary (",etool,")","\n"," ------------------------------","\n",
             " Overall accuracy =", round(sum(diag(class.resul2))/sum(class.resul2),digits=3),"\n")
  
  
  
  name=readline("Enter the output file name :  ")
  if(method=="BaggingEsb") cat(result3, result4, file=name) else
    cat(result1,result2,result3, result4, file=name)
  
  cat("Want to check the result in R console, too ? (1=No or 2=Yes): ")
  aw1 = scan(n=1, quiet=TRUE)
  if (aw1==1)  {
    cat("Check the result at your outpur file. \n")  		
  } else {cat(result1,result2,result3, result4, file="")} 
  
} # userDA() 


