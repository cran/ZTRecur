// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// TOTALM number of subjects in the study
// TOTALM2 number of visits in the study
// TMPvisitage age at each of the visits
// TMPcount number of visits for each subject at each age unit 
// tmpY the 0/1 element matrix for whether in censoring interval 
// tmpYs0 the 0/1 element matrix for whether in censoring interval & in stratum 1

// [[Rcpp::export]]
  
List getmoredata(arma::vec birth, 
                 int WINDOWL, int WINDOWR, arma::vec dates,
                 arma::vec p1stvt, 
                 int TOTALM, int TOTALM2, int ONUM,
                 arma::vec numvst, int AGEUNIT){
  arma::vec tmpbirth=arma::zeros<arma::vec>(TOTALM2); 
  arma::mat TMPcount=arma::zeros<arma::mat>(TOTALM,ONUM);
  arma::mat tmpY=arma::zeros<arma::mat>(TOTALM,ONUM),tmpYs0=arma::zeros<arma::mat>(TOTALM,ONUM);
  arma::vec TMPvisitage,TMPcensoringl,TMPcensoringr;

  int i,j,k;

  //update birthdate for each event
  tmpbirth=tmpbirth.zeros();

  for(i=0; i<TOTALM; i++)
    {
      for(k=0;k<numvst(i); k++)
        {
         j=p1stvt(i)-1+k;
          tmpbirth(j)=birth(i);
        }
    }
      
  //get age at visits (in age unit)
  TMPvisitage=arma::ceil((dates-tmpbirth)/365.25*AGEUNIT);//AC: change floor to ceil  ############1
  
  //c. get censoring times (in age unit) ------------------
  TMPcensoringl=arma::ceil((WINDOWL-birth)/365.25*AGEUNIT);//AC: change floor to ceil
  TMPcensoringr=arma::ceil((WINDOWR-birth)/365.25*AGEUNIT);//AC: change floor to ceil
  //std::cout<<"before i loop over TOTALM\n"<<std::endl;      
  for(i=0; i<TOTALM; i++)
    {  
    if(TMPcensoringl(i)<1)
    TMPcensoringl(i)=1;
    if(TMPcensoringr(i)>(ONUM-1))
    TMPcensoringr(i)=ONUM;
    }
  //std::cout<<"before part d\n"<<std::endl;  
  //d. get counts matrix (in age unit) ------------------
  TMPcount=TMPcount.zeros();
    for(i=0; i<TOTALM; i++)
      { for(k=0; k<numvst(i); k++)
        {
          j=TMPvisitage(p1stvt(i)-1+k);
          if(j<1){
            Rprintf("\n\ni=%d j=%d k=%d\n",i,j,k);//fflush(stdout);
            Rprintf("p1stvt(%d)=%.0f, p1stvt(i)-1+k=%.0f\n",i,p1stvt(i),p1stvt(i)-1+k);//fflush(stdout);
            Rprintf("TMPvisitage(%.0f)=%.0f\n",p1stvt(i)-1+k,TMPvisitage(p1stvt(i)-1+k));//fflush(stdout);
          }else{
            if(j<=ONUM){
              TMPcount(i,j-1)=TMPcount(i,j-1)+1;
            } 
          }       
        }
      }
        //std::cout<<"before part d 2"<<std::endl;

  tmpY=tmpY.zeros();
	tmpYs0=tmpYs0.zeros();//-AC
        //std::cout<<"before part e\n"<<std::endl;  
  for(i=0; i<TOTALM; i++)
    {
      for(j=0; j<ONUM; j++)
        {
         if(TMPcensoringl(i)<=(j+1)&&TMPcensoringr(i)>=(j+1)){
            tmpY(i,j)=1.0;
  		      if((j+1)<=TMPvisitage(p1stvt(i)-1)){
		          tmpYs0(i,j)=1.0;//-AC
		          }
            }
        }
    }

    return List::create(TMPcount,tmpY,TMPvisitage,tmpYs0);//-AC
}

//The following c codes for generating matrices for projects 
// with missing birthdates.

// [[Rcpp::export]]

List getmoredataComb2(arma::vec birthl, arma::vec birthr,
                      int WINDOWL, int WINDOWR,
                 arma::vec dates,
                 arma::vec p1stvt, int reptw,
                 int TOTALM, int TOTALM2, int ONUM,
                 arma::vec numvst, int AGEUNIT){
  arma::vec tmpbirth=arma::zeros<arma::vec>(TOTALM2),TMPbirthdate;
  arma::mat TMPcount=arma::zeros<arma::mat>(TOTALM,ONUM),
    TMPcountnew=arma::zeros<arma::mat>(TOTALM,ONUM),
    TMPcountnewf=arma::zeros<arma::mat>(TOTALM,ONUM);
  arma::mat tmpY=arma::zeros<arma::mat>(TOTALM,ONUM),
    TMPYnew=arma::zeros<arma::mat>(TOTALM,ONUM),
    TMPYnewf=arma::zeros<arma::mat>(TOTALM,ONUM);
  arma::mat tmpYs0=arma::zeros<arma::mat>(TOTALM,ONUM),
    TMPYs0new=arma::zeros<arma::mat>(TOTALM,ONUM),
    TMPYs0newf=arma::zeros<arma::mat>(TOTALM,ONUM);
  arma::mat dNStarnew=arma::zeros<arma::mat>(TOTALM,ONUM),
    dNStarnewf=arma::zeros<arma::mat>(TOTALM,ONUM),
    dNStars0new=arma::zeros<arma::mat>(TOTALM,ONUM),
    dNStars0newf=arma::zeros<arma::mat>(TOTALM,ONUM);
  arma::vec TMPvisitage,TMPcensoringl,TMPcensoringr;
  arma::vec TMPvisitagetot=arma::zeros<arma::vec>(TOTALM2);
  arma::vec TMPvisitagef=arma::zeros<arma::vec>(TOTALM2);
  
  int i,j,k;
  
  for (int it=0;it<reptw;it++)        
  {
    //generate birthdates from a uniform distribution
    TMPbirthdate=arma::floor(birthl+arma::randu<arma::vec>(TOTALM)/(1+0.00001)%(birthr-birthl));
    tmpbirth=tmpbirth.zeros();
  
  for(i=0; i<TOTALM; i++)
  {
    for(k=0;k<numvst(i); k++)
    {
      j=p1stvt(i)-1+k;
      tmpbirth(j)=TMPbirthdate(i);// make birthdate has the same length with visit
    }
  }
  
  //get age at visits (in age unit)
  TMPvisitage=arma::ceil((dates-tmpbirth)/365.25*AGEUNIT);
  TMPvisitagetot=TMPvisitagetot+TMPvisitage; 
  
  //c. get censoring times (in age unit) ------------------
  TMPcensoringl=arma::ceil((WINDOWL-TMPbirthdate)/365.25*AGEUNIT);//AC: change floor to ceil
  TMPcensoringr=arma::ceil((WINDOWR-TMPbirthdate)/365.25*AGEUNIT);//AC: change floor to ceil
     
  for(i=0; i<TOTALM; i++)
  {  
    if(TMPcensoringl(i)<1)
      TMPcensoringl(i)=1;
    if(TMPcensoringr(i)>(ONUM-1))
      TMPcensoringr(i)=ONUM;
  }
  
  TMPvisitage.elem(arma::find(TMPvisitage > ONUM)).fill(ONUM);
  
  //std::cout<<"before part d\n"<<std::endl;  
  //d. get counts matrix (in age unit) ------------------
  TMPcount.zeros();
  for(i=0; i<TOTALM; i++)
  { for(k=0; k<numvst(i); k++)
  {
    j=TMPvisitage(p1stvt(i)-1+k);
    if(j<1){
      Rprintf("\n\ni=%d j=%d k=%d\n",i,j,k);//fflush(stdout);
      Rprintf("p1stvt(%d)=%.0f, p1stvt(i)-1+k=%.0f\n",i,p1stvt(i),p1stvt(i)-1+k);//fflush(stdout);
      Rprintf("TMPvisitage(%.0f)=%.0f\n",p1stvt(i)-1+k,TMPvisitage(p1stvt(i)-1+k));//fflush(stdout);
    }else{//-AC
      //if(j<ONUM){AC: Change for ceil function##############2
      if(j<=ONUM){
        //TMPcount(i,j)=TMPcount(i,j)+1;AC: Change for ceil function#################3
        TMPcount(i,j-1)=TMPcount(i,j-1)+1;
      } 
    }       
  }
  }
  //std::cout<<"before part d 2"<<std::endl;
  TMPcountnew=TMPcountnew+TMPcount;
  
  tmpY.zeros();
  tmpYs0.zeros();//-AC
  //std::cout<<"before part e\n"<<std::endl;  
  for(i=0; i<TOTALM; i++)
  {
    for(j=0; j<ONUM; j++)
    {
      if(TMPcensoringl(i)<=(j+1)&&TMPcensoringr(i)>=(j+1)){
        tmpY(i,j)=1.0;
        if((j+1)<=TMPvisitage(p1stvt(i)-1)){
          tmpYs0(i,j)=1.0;//-AC
        }
      }
    }
  }
  
  TMPYnew=TMPYnew+tmpY;
  TMPYs0new=TMPYs0new+tmpYs0;
  dNStarnew=dNStarnew+tmpY%TMPcount;
  dNStars0new=dNStars0new+tmpYs0%TMPcount;
  }
  
  TMPcountnewf=(TMPcountnew/reptw);
  TMPYnewf=(TMPYnew/reptw);
  TMPYs0newf=(TMPYs0new/reptw);
  dNStarnewf=(dNStarnew/reptw);
  dNStars0newf=(dNStars0new/reptw);
  TMPvisitagef=arma::round((TMPvisitagetot/reptw));
  return List::create(TMPcountnewf,TMPYnewf,TMPvisitagef,TMPYs0newf,
                      dNStarnewf,dNStars0newf);//-AC
}


// time-indpt effect and cov
struct AGreturn{
  arma::mat m1; //PscoreDerivInv
  arma::mat m2; //S1tByk
  arma::vec v1; //beta
  arma::vec v2; //Pscore
  arma::vec v3; //S0t
  int n;
} AGreturnv,AGreturnva,AGreturnva1,AGreturnva2;

AGreturn AGcoh(arma::vec betaini, arma::mat covat0z, 
               arma::mat YcYsz,arma::vec YN,arma::mat YZN,
               double assumedtiny, int reptmax,double etapara,int timeUN)
{
  // betaDiff: the average absolute difference between two successive estimates of beta
  // reptmax: the max number of iterations. 
  // m: counter of iterations.
  arma::vec beta=betaini, beta0;
  arma::vec betaZ,expbetaZ,Sbar,S0t,Pscore;
  arma::mat PscoreDeriv,S1tByk,PscoreDerivInv;
  double betaDiff=1000;//  betaRelDiff
  int m=0;// initial value of iteration counter
  //int nSubj=covat0.n_rows;
  //int nSubjz=covat0z.n_rows;
  int nCov=covat0z.n_cols;
  
  while((betaDiff>etapara)&&(m<reptmax))
  {
    beta0=beta;
    
    /*2.calculate the needed quantities*/
    /*to get expbetaZ*/
    //printf("check: getting tmpexpbetaZn");
    betaZ=covat0z*beta0;//dim: # of covariate combination * 1
    //std::cout<<betaZ(1)<<std::endl;
    expbetaZ=arma::exp(betaZ);//dim: # of covariate combination * 1
    //std::cout<<expbetaZ(100)<<std::endl;
    
    /*to get S0t*/
    //printf("check: getting S0tn");
    S0t=YcYsz*expbetaZ;//dim: # of time unit*1
    //std::cout<<S0t(1)<<std::endl;
    
    /*to get S1tByk*/  
    arma::mat tmpcov=covat0z.each_col()%expbetaZ;
    S1tByk=YcYsz*tmpcov;//dim: # of time units * # of covaraites

    
    /*to get partial score*/
    arma::vec invS0 = YN / (S0t + assumedtiny);
    Sbar=arma::trans(S1tByk)*invS0;//dim: # of covariate *1
    //Pscore: value of partial score function for given beta's. 
    Pscore=YZN-Sbar;
    
    arma::mat TMPPG02=arma::zeros<arma::mat>(nCov,nCov);
    arma::mat TMPPG01=(S1tByk.t() * (S1tByk.each_col() % (invS0/(S0t + assumedtiny))));
    //tmpS2=arma::zeros<arma::mat>(nCov,nCov);
    
    //calculate the second derivative of partial log-likelihood function
    // or negative observed fisher information
    for (int j = 0; j < timeUN; j++)
    {
      arma::rowvec Yj = YcYsz.row(j);               // 1 × nSubjz
      arma::vec wj = Yj.t() % expbetaZ;             // nSubjz × 1
      
      // S2(j)
      TMPPG02 += (covat0z.t() * (covat0z.each_col() % wj)) * (invS0(j));
    }
    
    
    //PscoreDeriv: the second derivative of the partial log-likelihood function
    //TMPPG01=(S1/S2)^2
    //TMPPG02=S2/S0
    PscoreDeriv=TMPPG01-TMPPG02;
    
    //std::cout<<PscoreDeriv<<std::endl;
    /*to get inverse of the 2nd derivative*/
    PscoreDerivInv= arma::pinv(PscoreDeriv,1.0e-8);//1.0e-8

    beta=beta0-PscoreDerivInv*Pscore;
    
    /*4. is the new beta good enough*/
    betaDiff=arma::sum(arma::abs((beta0-beta)/beta0))/nCov;
    //betaRelDiff=sum(abs((beta0-beta)/beta0))/nCov;
    m++;
  }
  
  if(m>=reptmax&&m>2)
  {Rcpp::Rcout << "getting AG beta does not converge at betaDiff="<<betaDiff<<", m="<< m << std::endl;}
  
  
  /*5. finalize beta estm*/
  //printf("AG1 betaestm: m=%dn",m);
  beta0=(beta+beta0)/2.0;
  AGreturnv.m1=PscoreDerivInv;
  AGreturnv.m2=S1tByk;
  AGreturnv.v1=beta0;
  AGreturnv.v2=Pscore;//score function
  AGreturnv.v3=S0t;	
  AGreturnv.n=m;	
  return AGreturnv;
}


//Estimate beta.1 and beta.2 separately or no strata
// [[Rcpp::export]]

List getAG1coh (arma::vec betaini, arma::mat covat0,arma::mat covat0z, 
                arma::mat YcYs,arma::mat YcYsz, double assumedtiny, 
                int reptmax, double etapara, int timeUN, 
                int calcvar,int nSubj0,arma::mat YNmat){
  // AC: I added the YcYsz, which is the Y_c*Y_s in the ratio, Z_bar.
  arma::vec S0t=arma::zeros<arma::vec>(timeUN);
  int nSubj=covat0.n_rows;
  //int nSubjz=covat0z.n_rows;
  int nCov=covat0z.n_cols;
  
  //a. get Y(j)*dN(j): total obs visits at age j ---
  //add matrix TMPcount for each row (dim=1)
  arma::vec YN=arma::sum(arma::trans(YNmat),1);  //timeUN*1
  arma::vec YNnsub=arma::sum(YNmat,1);//nSubj*1
  
  //c. get Y(.)Z(k)N(.): weighted sum of num obs visits ---
  // arma::vec YZN=arma::trans(covat0)*numvst;
  arma::vec YZN=arma::zeros<arma::vec>(nCov);
  arma::mat temp1=arma::zeros<arma::mat>(nSubj,nCov);
  
  temp1=covat0.each_col()%YNnsub;//nSubj*nCov
  YZN=arma::sum(arma::trans(temp1),1);
  
  //fit until converge or reach max iteration
  AGreturnva=AGcoh(betaini, covat0z, YcYsz,YN,YZN,
                   assumedtiny,reptmax,etapara, timeUN);
  // extract the estimates of betas
  arma::vec beta=AGreturnva.v1;
  int m=AGreturnva.n;
  
  // Initialize the covariance matrix of betas
  arma::mat betavar=arma::zeros<arma::mat>(nCov,nCov);
  
  //fit again with initial beta equal the average of last two beta result to find variance
  AGreturnva=AGcoh(beta,covat0z, YcYsz, YN,YZN,assumedtiny,1,
                   etapara, timeUN);
  
  //beta=AGreturnva.v1;//estimates of beta
  arma::mat PscoreDerivInv=AGreturnva.m1;
  arma::mat S1tByk=AGreturnva.m2;
  arma::vec Pscore=AGreturnva.v2;//partial score function; dim: # of covariates*1
  S0t=AGreturnva.v3;
  
  //std::cout <<"PscoreDerivInv" <<PscoreDerivInv<< std::endl;
  
  arma::vec Lambda=arma::zeros<arma::vec>(timeUN);
  arma::mat temp0=arma::zeros<arma::mat>(nSubj,nCov);
  arma::mat PscoreMat=arma::zeros<arma::mat>(nSubj,nCov),TMPGvar;
  
  if(calcvar==1){ //calculate the variance
    arma::mat Sbar=S1tByk.each_col()/(S0t+assumedtiny);//timeUN*nCov
    
    temp0=YNmat*Sbar;//nSubj*timeUN*timeUN*nCov
    
    PscoreMat=temp1-temp0;// score function NSubj*nCov
    
    TMPGvar=arma::trans(PscoreMat)*PscoreMat;// dim: # of covariate * # of covariates
    //std::cout <<"arma::trans(PscoreMat)*PscoreMat" <<TMPGvar<< std::endl;
    TMPGvar=(TMPGvar*nSubj0-Pscore*arma::trans(Pscore))/(nSubj0-1);//check: it's the right formula
    //std::cout <<"Pscore*arma::trans(Pscore)/nSubj0" <<Pscore*arma::trans(Pscore)/nSubj0<< std::endl;
    //std::cout <<"TMPGvar" <<TMPGvar<< std::endl;
    
    betavar=PscoreDerivInv*TMPGvar*PscoreDerivInv;
    
  } //end of calculating the variance
  
  
  /*6. get Lambda0 estm*/
  // The previous S0all is for the last iteration of betas.
  // Since we consider the average of the last two iteration of betas, it's
  // better to recompute S0 again just in case. I believe the differences are small
  S0t=YcYsz*arma::exp(covat0z*beta);
  Lambda=arma::cumsum((YN/(S0t+assumedtiny)));
  
  return List::create(beta,betavar,Lambda,S0t,S1tByk,PscoreDerivInv,m,
                      Pscore);
}

//Estimate beta.1 and beta.2 jointly
// [[Rcpp::export]]

List getAG1cohS2 (arma::vec betaini1,arma::vec betaini2,
                  arma::vec lambdaini1,arma::vec lambdaini2,
                  arma::mat covat0,arma::mat covat0z, 
                  arma::mat Yc,arma::mat Yc1,
                  arma::mat Nstar, double assumedtiny, 
                  int reptmax, double etapara, int timeUN, 
                  arma::mat count, arma::vec FlagBirth, 
                  arma::vec a1star,
                  arma::vec CL, //arma::vec CR,
                  int calcvar,int nSubj0, int comCoef,
                  int comLamb){//, int tauL2
  // AC: I added the YcYsz, which is the Y_c*Y_s in the ratio, Z_bar.
  arma::vec S0t1=arma::zeros<arma::vec>(timeUN);
  arma::vec S0t2=arma::zeros<arma::vec>(timeUN);
  int nSubj=covat0.n_rows;
  int nSubjz=covat0z.n_rows;
  int nCov=covat0z.n_cols;
  
  arma::mat Yc1pz,Yc2pz,Yc1p,Yc2p;
  arma::vec YN1, YN2;   
  
  arma::vec YZN1=arma::zeros<arma::vec>(nCov);
  arma::vec YZN2=arma::zeros<arma::vec>(nCov);
  
  arma::mat EIs1z=arma::zeros<arma::mat>(nSubjz,timeUN);
  arma::mat EIs2z=arma::zeros<arma::mat>(nSubjz,timeUN);
  
  arma::vec beta1,beta2,beta01, beta02;
  arma::vec lambda1=lambdaini1;
  arma::vec lambda2=lambdaini2;
  arma::mat hr1,hr2,hr1z;//,hr2z
  
  arma::mat betavar1=arma::zeros<arma::mat>(nCov,nCov);
  arma::mat betavar2=arma::zeros<arma::mat>(nCov,nCov);
  
  double betaRelDiff1=1000,betaRelDiff2=1000;
  int m=0;

  while(((betaRelDiff1>(2*etapara))||(betaRelDiff2>(2*etapara)))&&(m<reptmax))
  {
    if(m==0){
      beta01=betaini1;
      beta02=betaini2;
    }else{
      beta01=beta1;
      beta02=beta2;
    }
    
    //a. Generate harzard functions for strata 1 & 2
    hr1=arma::exp(covat0*beta01)*arma::trans(lambda1);//dim: n*time unit
    hr2=arma::exp(covat0*beta02)*arma::trans(lambda2);
    hr1z=arma::exp(covat0z*beta01)*arma::trans(lambda1);//dim: covariate combination*time unit
    //hr2z=arma::exp(covat0z*beta02)*arma::trans(lambda2);
    
    
    arma::vec Nr=arma::zeros<arma::vec>(nSubj);
    arma::vec Dr=arma::zeros<arma::vec>(nSubj);
    Yc1p=Yc1;
    
    //Rcout<<"Flag 1\n";
    for(int i=0;i<nSubj;i++){
      if(FlagBirth(i)==0){
        //born before the study
        
        //a1 is not missing, which means that we observed at least one event.
        double tmpAC=0.0,tmpAC1=0.0,tmpAC2=0.0;
        
        for(int j=0; j<a1star(i); j++){
          tmpAC=tmpAC-hr1(i,j);
          if(j<CL(i)){
            tmpAC1=tmpAC1-hr1(i,j);
          }else{
            tmpAC2=tmpAC2-hr2(i,j);
          }
        }
        Nr(i)=hr1(i,a1star(i)-1)*std::exp(tmpAC);
        Dr(i)=Nr(i)+(1-std::exp(tmpAC1))*hr2(i,a1star(i)-1)*std::exp(tmpAC2);
        
        Yc1p.row(i)=Yc1.row(i)*(Nr(i)/Dr(i));
      }
      
    }
    //Rcout<<"Flag 2\n";
    Yc2p=Yc-Yc1p;
    
    EIs1z=arma::zeros<arma::mat>(nSubjz,timeUN);
    EIs2z=arma::zeros<arma::mat>(nSubjz,timeUN);
    
    for(int i=0;i<nSubjz;i++){
      double tmpACz=0.0;
      EIs1z(i,0)=1.0;
      for(int j=1; j<timeUN; j++){
        tmpACz=tmpACz-hr1z(i,j-1);
        EIs1z(i,j)=std::exp(tmpACz);
      }
    }
    EIs2z=1.0-EIs1z;
    
    Yc1pz=arma::trans(EIs1z)%Nstar;// time unit* covariate combs
    Yc2pz=arma::trans(EIs2z)%Nstar;
    
    //a. get Y(j)*dN(j): total obs visits at age j ---
    //add matrix TMPcount for each row (dim=1)
    YN1=arma::sum(arma::trans(count%Yc1p),1);  
    YN2=arma::sum(arma::trans(count%Yc2p),1);   
    
   
    YZN1=arma::sum(covat0.t()*(count%Yc1p),1);//nCov*timeUN
    YZN2=arma::sum(covat0.t()*(count%Yc2p),1);//nCov*timeUN
    
    //fit until converge or reach max iteration
    //if(betaDiff1>(500*etapara)){
    if(comCoef==0){
      AGreturnva1=AGcoh(beta01, covat0z, Yc1pz,YN1,YZN1,
                        assumedtiny,reptmax,etapara, timeUN);
      // extract the estimates of betas
      beta1=AGreturnva1.v1;
    
      AGreturnva2=AGcoh(beta02, covat0z, Yc2pz,YN2,YZN2,assumedtiny,
                        reptmax,etapara, timeUN);
      // extract the estimates of betas
      beta2=AGreturnva2.v1;
    }else{
      AGreturnva1=AGcoh(beta01, covat0z,(Yc1pz+Yc2pz),
                        (YN1+YN2),(YZN1+YZN2),assumedtiny,reptmax,etapara, 
                        timeUN);
      // extract the estimates of betas
      beta1=AGreturnva1.v1;
      beta2=beta1;
    }
    
    lambda1=arma::zeros<arma::vec>(timeUN);
    lambda2=arma::zeros<arma::vec>(timeUN);
    
    
    
    /*6. get Lambda0 estm*/
    // The previous S0all is for the last iteration of betas.
    // Since we consider the average of the last two iteration of betas, it's
    // better to recompute S0 again just in case. I believe the differences are small
    //Rcout<<"beta1 before S0t1"<<beta1<<"\n";
    //Rcout<<"beta2 before S0t2"<<beta2<<"\n";
    //Yc1pz dim: time unit* covariate combs
    S0t1=Yc1pz*arma::exp(covat0z*beta1);//dim: time *1
    S0t2=Yc2pz*arma::exp(covat0z*beta2);
    
    if(comLamb==0){
      lambda1=YN1/(S0t1+assumedtiny);
      lambda2=YN2/(S0t2+assumedtiny);
      //the denominator is 0 so the lambda value is extremely large.
      lambda2(0)=0;
    }else{
      lambda1=(YN1+YN2)/(S0t1+S0t2+assumedtiny);
      lambda2=lambda1;
    }
    
    betaRelDiff1=arma::sum(arma::abs((beta01-beta1)/beta01))/nCov;
    betaRelDiff2=arma::sum(arma::abs((beta02-beta2)/beta02))/nCov;
    
    m++;
  }
  
  if(m>=reptmax&&m>2)
  {Rcpp::Rcout << "getting AG beta does not converge at m="<< m << std::endl;}
  
  arma::vec Lambda1=arma::cumsum(lambda1);
  arma::vec Lambda2=arma::cumsum(lambda2);
  
  if(comCoef==0){
    //fit again with initial beta equal the average of last two beta result to find variance
    AGreturnva1=AGcoh(beta1, covat0z, Yc1pz,YN1,YZN1,
                      assumedtiny,1,etapara, timeUN);
    AGreturnva2=AGcoh(beta2, covat0z, Yc2pz,YN2,YZN2,
                      assumedtiny,1,etapara, timeUN);
  }else{
    AGreturnva1=AGcoh(beta1, covat0z,(Yc1pz+Yc2pz),
                      (YN1+YN2),(YZN1+YZN2),assumedtiny,1,etapara, 
                      timeUN);
    AGreturnva2=AGreturnva1;
  }
  
  
  arma::mat PscoreDerivInv1=AGreturnva1.m1;
  arma::mat PscoreDerivInv2=AGreturnva2.m1;
  arma::vec Pscore1=AGreturnva1.v2;//partial score function; dim: # of covariates*1
  arma::vec Pscore2=AGreturnva2.v2;//partial score function; dim: # of covariates*1
  
  arma::mat S1tByk1=AGreturnva1.m2;
  arma::mat S1tByk2=AGreturnva2.m2;
  
  S0t1=AGreturnva1.v3;
  S0t2=AGreturnva2.v3;
  
  arma::vec temp10=arma::zeros<arma::vec>(nCov);
  arma::vec temp20=arma::zeros<arma::vec>(nCov);
  arma::mat PscoreMat1=arma::zeros<arma::mat>(nSubj,nCov);
  arma::mat PscoreMat2=arma::zeros<arma::mat>(nSubj,nCov);
  arma::mat temp11=arma::zeros<arma::mat>(nSubj,nCov);
  arma::mat temp21=arma::zeros<arma::mat>(nSubj,nCov);
  
  if(calcvar==1){ //calculate the variance
    /*to get partial scores variance estm*/
    for(int i=0; i<nSubj; i++) { 
      for(int k=0; k<nCov;k++) {       
        temp10(k)=0.0;
        temp11(i,k)=0.0;
        
        temp20(k)=0.0;
        temp21(i,k)=0.0;
        /* YX: revised codes to get variance estm*/
        for(int j=0;j<timeUN;j++){
          if(count(i,j)>0){
            if(comCoef==0){
              temp10(k)+=S1tByk1(j,k)/(S0t1(j)+assumedtiny)*Yc1p(i,j)*count(i,j);
              temp11(i,k)+=Yc1p(i,j)*covat0(i,k)*count(i,j);
              
              temp20(k)+=S1tByk2(j,k)/(S0t2(j)+assumedtiny)*Yc2p(i,j)*count(i,j);
              temp21(i,k)+=Yc2p(i,j)*covat0(i,k)*count(i,j);
            }else{
              temp10(k)+=S1tByk1(j,k)/(S0t1(j)+assumedtiny)*(Yc1p(i,j)+Yc2p(i,j))*count(i,j);
              temp11(i,k)+=(Yc1p(i,j)+Yc2p(i,j))*covat0(i,k)*count(i,j);
            }
            
          }
        }
        PscoreMat1(i,k)=temp11(i,k)-temp10(k);// score function
        if(comCoef==0){
          PscoreMat2(i,k)=temp21(i,k)-temp20(k);// score function
        }else{
          PscoreMat2(i,k)= PscoreMat1(i,k);
        }
      }
    }
    
    arma::mat TMPGvar1=arma::trans(PscoreMat1)*PscoreMat1;// dim: # of covariate * # of covariates
    TMPGvar1=(TMPGvar1*nSubj0-Pscore1*arma::trans(Pscore1))/(nSubj0-1);//check: it's the right formula
    
    arma::mat TMPGvar2=arma::trans(PscoreMat2)*PscoreMat2;// dim: # of covariate * # of covariates
    TMPGvar2=(TMPGvar2*nSubj0-Pscore2*arma::trans(Pscore2))/(nSubj0-1);//check: it's the right formula
    
    // Initialize the covariance matrix of betas
    //Rcout<<"Flag 11\n";
    
    for(int k=0; k<nCov; k++) {
      for(int k1=0; k1<nCov; k1++) {
        betavar1(k,k1)=0.0;
        betavar2(k,k1)=0.0;
        for(int m1=0;m1<nCov;m1++) {
          for(int m2=0;m2<nCov;m2++) {
            // sandwitch form/ robust variance. 
            betavar1(k,k1)+=PscoreDerivInv1(k,m1)*TMPGvar1(m1,m2)*PscoreDerivInv1(m2,k1);
            betavar2(k,k1)+=PscoreDerivInv2(k,m1)*TMPGvar2(m1,m2)*PscoreDerivInv2(m2,k1);
          }	
        }
      }
    }
    
  } //end of calculating the variance
  
  return List::create(beta1,beta2,Lambda1,Lambda2,betavar1,betavar2,
                      Pscore1,Pscore2,Yc1p,Yc2p,PscoreMat1,PscoreMat2,
                      PscoreDerivInv1,PscoreDerivInv2,m,S0t1,S0t2,
                      YN1,YN2);
}

// local constant estimates 
struct LCreturn{
  arma::mat m1; //PscoreDerivInv
  arma::mat m2; //S1tByk
  arma::vec v1; //beta
  arma::vec v2; //Pscore//score function
  arma::vec v3; //S0t
} LCreturnv,LCreturnva,LCreturnva1,LCreturnva2;


LCreturn AGcoha( arma::vec betaini, arma::mat covat0,arma::mat covat0z, 
                arma::mat YcYs,arma::mat YcYsz, arma::mat YN,arma::mat YZNW,
                double assumedtiny, int reptmax,double etapara,int timeUN,
                int iter,int windowKH0,int nSubj,int nCov)
{
  // betaDiff: the average absolute difference between two successive estimates of beta
  // reptmax: the max number of iterations. 
  // m: counter of iterations.
  arma::vec beta, beta0,betaZ,expbetaZ,Sbar,temp0,S0t;
  arma::mat PscoreDeriv,TMPPG01,TMPPG02,S1tByk,tmpS2,Pscore,PscoreDerivInv;
  double betaDiff;//TMPPL, TMPPL0, tmpp,
  int m=0;// initial value of iteration counter

  betaDiff=1000;// initial value of difference in estimates
  Sbar=arma::zeros<arma::vec>(nCov);
  TMPPG01=arma::zeros<arma::mat>(nCov,nCov);
  TMPPG02=arma::zeros<arma::mat>(nCov,nCov);
  tmpS2=arma::zeros<arma::mat>(nCov,nCov);
  //int nSubj=covat0.n_rows;
  //int nCov=covat0.n_cols;
  int nSubjz=covat0z.n_rows;
  
  /*1.set initial of beta*/
  beta=betaini;
  
  while((betaDiff>etapara)&&(m<reptmax))
  {
    beta0=beta;
    
    /*2.calculate the needed quantities*/
    /*to get expbetaZ*/
    betaZ=covat0z*beta0;//dim: # of covariate combination * 1
    expbetaZ=arma::exp(betaZ);//dim: # of covariate combination * 1
    
    /*to get S0t*/
    S0t=YcYsz*expbetaZ;//dim: # of time unit*1
    
    /*to get S1tByk*/  
    arma::mat tmpcov=arma::zeros<arma::mat>(nSubjz,nCov);
    tmpcov.each_col()+=expbetaZ;

    S1tByk=YcYsz*(covat0z%tmpcov);//dim: # of time units * # of covaraites
    
    /*to get partial score*/
    Sbar.zeros();//dim: # of covariate *1
    for(int k=0; k<nCov; k++)
    {
      for(int j=0; j<timeUN; j++)
      {
        if((j<=(iter+windowKH0))&&(j>=(iter-windowKH0)))
          Sbar(k)+=0.75/windowKH0*(1-(j-iter)*(j-iter)/windowKH0/windowKH0)*S1tByk(j,k)/(S0t(j)+assumedtiny)*YN(j);
      }
    }
    //YZNW: # of time unit * # covariate
    Pscore=arma::trans(YZNW.row(iter))-Sbar;
    
    //calculate the second derivative of partial log-likelihood function
    // or negative observed fisher information
    TMPPG01.zeros();
    TMPPG02.zeros();
    //S1tByk: dim: # time * # covariate
    for(int k=0; k<nCov; k++)
    {
      for(int k1=0; k1<nCov; k1++)
      {
        for(int j=0; j<timeUN; j++)
        {
          if((j<=(iter+windowKH0))&&(j>=(iter-windowKH0)))
          {   
            TMPPG01(k,k1)+=0.75/windowKH0*(1-(j-iter)*(j-iter)/windowKH0/windowKH0)*S1tByk(j,k)*
              S1tByk(j,k1)/pow((S0t(j)+assumedtiny),2.0)*YN(j);
            
            tmpS2(k,k1)=0.0;
            
            for(int i=0; i<nSubjz; i++)
            {
              tmpS2(k,k1)+=YcYsz(j,i)*covat0z(i,k)*covat0z(i,k1)*expbetaZ(i);
            }		
            TMPPG02(k,k1)+=0.75/windowKH0*(1-(j-iter)*(j-iter)/windowKH0/windowKH0)*
              tmpS2(k,k1)/(S0t(j)+assumedtiny)*YN(j);
            //TMPPG02(k,k1)+=0.5/windowKH0*tmpS2(k,k1)/(S0t(j)+assumedtiny)*YN(j);
            
          }
        }
        
      }
      
    }
    //PscoreDeriv: the second derivative of the partial log-likelihood function
    PscoreDeriv=TMPPG01-TMPPG02;
    
    //std::cout<<"PscoreDeriv"<<PscoreDeriv<<std::endl;
    /*to get inverse of the 2nd derivative*/
    PscoreDerivInv= pinv(PscoreDeriv, 1.0e-4);
    //std::cout<<"PscoreDerivInv"<<PscoreDerivInv<<std::endl;
    
    /*3.update the current beta
     Pscore:score function
     PscoreDerivInv: inverse of the 2nd derivative of log-likelihood
     Use Newton's method:
     */
    //std::cout<<"Old beta"<<beta0<<std::endl;
    beta=beta0-PscoreDerivInv*Pscore;
    //std::cout<<"new beta"<<beta<<std::endl;
    
    /*4. is the new beta good enough*/
    //betaDiff=sum(abs(beta0-beta))/nCov;
    betaDiff=sum(abs((beta0-beta)/beta0))/nCov;
    
    m++;
    //std::cout <<"Difference in beta"<< betaDiff<< std::endl;
  }
  
  if(m>=reptmax&&m>2)
  {Rcpp::Rcout << "getting AG beta does not converge at betaDiff="<<betaDiff<<", m="<< m << std::endl;}
  
  
  /*5. finalize beta estm*/
  //printf("AG1 betaestm: m=%dn",m);
  beta0=(beta+beta0)/2.0;
  LCreturnv.m1=PscoreDerivInv;
  LCreturnv.m2=S1tByk;
  LCreturnv.v1=beta0;//-AC
  LCreturnv.v2=Pscore;//score function
  LCreturnv.v3=S0t;	
  return LCreturnv;
}

//Estimate beta.1 and beta.2 seperately
// [[Rcpp::export]]

List getLCcoh(arma::mat betaMat, arma::mat covat0,arma::mat covat0z,
                 arma::mat YcYs,arma::mat YcYsz,
                 double assumedtiny, int reptmax, double etapara, int timeUN,
                 arma::vec windowKH, int tauL, int tauR, int calcvar,int nSubj0, 
                 arma::mat YNmat){ 
  // AC: I added the YcYsz, which is the Y_c*Y_s in the ratio, Z_bar.
  int nSubj=covat0.n_rows;
  //int nSubjz=covat0z.n_rows;
  //int windowKH=0;-AC
  int nCov=covat0z.n_cols;
  
  arma::mat EstBetaMat=arma::zeros<arma::mat>(nCov,timeUN);
  arma::mat BetaVarMat=arma::zeros<arma::mat>(nCov,timeUN);
  arma::mat betaZ,expbetaZ;
  arma::vec S0t;
  arma::vec Lambda=arma::zeros<arma::vec>(timeUN);
  //betaDiffLambdaNA=arma::zeros<arma::vec>(timeUN);temp0=arma::zeros<arma::vec>(nCov)-AC
  arma::mat PscoreMat=arma::zeros<arma::mat>(nSubj,nCov),TMPGvar;
  int j, k, l;
  
  
  //a. get Y(j)*dN(j): total obs visits at age j ---
  //add matrix betaDiffcount for each row (dim=1)
  arma::vec YN=arma::sum(arma::trans(YNmat),1); //corrected: AC 
  arma::mat YZN=covat0.t()*YNmat;//nCov*timeUN
  
  
  //b2. get Y(.)Z(k)dN(.)W(l,j):weighted total obs visits at age j around l
  //windowKH=windowKH1;-AC
  arma::mat YZNW=arma::zeros<arma::mat>(timeUN,nCov);//dim: # time * # covariate
  int windowKH00;
  for(l=(tauL-1); l<tauR; l++)
  {
    windowKH00=windowKH(l);
    for(j=0;j<timeUN;j++)
    {
      if((j<=(l+windowKH00))&&(j>=(l-windowKH00)))
      {
        YZNW.row(l)+=0.75/windowKH00*(1-(j-l)*(j-l)/windowKH00/windowKH00)*YZN.col(j).t();//corrected:-AC
      }
    }
  }
  
  arma::vec betaini;
  int windowKH0;
  
  //fit until converge or reach max iteration
  // iter is the target time point
  for (int iter=tauL-1;iter<tauR;iter++)
  {
    betaini=arma::trans(betaMat.row(iter));
    windowKH0=windowKH(iter);
    //std::cout << iter<<std::endl;
    LCreturnva=AGcoha(betaini, covat0,covat0z, YcYs,YcYsz,YN,
                     YZNW,assumedtiny,reptmax,etapara, timeUN,iter,
                     windowKH0,nSubj,nCov);
    // extract the estimates of betas
    arma::vec beta=LCreturnva.v1;
    
    //fit again with initial beta equal the average of last two beta result to find variance
    LCreturnva=AGcoha(beta,covat0,covat0z, YcYs,YcYsz, YN,
                     YZNW,assumedtiny,1,etapara, timeUN,iter,windowKH0,
                     nSubj,nCov);
    
    arma::mat PscoreDerivInv=LCreturnva.m1; 
    arma::mat S1tByk=LCreturnva.m2;
    arma::vec Pscore=LCreturnva.v2;//partial score function
    S0t=LCreturnva.v3;
    //arma::vec Lambda=arma::zeros<arma::vec>(timeUN);
    
    arma::mat betavar=arma::zeros<arma::mat>(nCov,nCov);
    
    if(calcvar==1){ //calculate the variance
      arma::mat PscoreMat=arma::zeros<arma::mat>(nSubj,nCov),TMPGvar;
      
      /*to get partial scores variance estm*/
      for(int i=0; i<nSubj; i++)
      { 
        for(int j=0; j<timeUN; j++)
        {     
          if((j<=(iter+windowKH0))&&(j>=(iter-windowKH0))){
            if(YNmat(i,j)>0){
              double weight=0.75/windowKH0*(1.0-(j-iter)*(j-iter)/windowKH0/windowKH0);
              arma::rowvec Sj=S1tByk.row(j)/(S0t(j)+assumedtiny);
              PscoreMat.row(i)+=weight*(covat0.row(i)-Sj)*YNmat(i,j);
            }
          }
        }
      }
      
      TMPGvar=arma::trans(PscoreMat)*PscoreMat;
      TMPGvar=(TMPGvar*nSubj0-Pscore*arma::trans(Pscore))/(nSubj0-1); // should be equivalent to Sigma in the paper
      
      betavar.zeros();
      betavar=PscoreDerivInv*TMPGvar*PscoreDerivInv;
      //update beta and betavar matirx: collect the variance
      BetaVarMat.col(iter)=betavar.diag();//dim: # covariate * # time
    }
    
    EstBetaMat.col(iter)=beta;//dim: # covariate * # time
  }
  
  //Take its transpose  
  EstBetaMat= EstBetaMat.t();//dim: # time * # covariate
  BetaVarMat=BetaVarMat.t();//dim: # time * # covariate
  
  //Extend the estimates of betas to [0,timeUN-1]
  for(l=0;l<tauL-1;l++)
    for(k=0; k<(nCov); k++)
      EstBetaMat(l,k)=EstBetaMat(tauL-1,k);
  
  for(l=tauR;l<timeUN;l++)//corrected: -AC
    for(k=0; k<(nCov); k++)
      EstBetaMat(l,k)=EstBetaMat(tauR-1,k);
  
  
  /*6. get Lambda0 estm*/
  // The previous S0all is for estimated beta at each time unit.
  // We need to recalculate it for the final estimation results
  S0t=arma::sum(YcYsz%arma::exp(EstBetaMat*covat0z.t()),1);
  
  //std::cout << "pass"<<std::endl;
  Lambda=arma::cumsum((YN/(S0t+assumedtiny)));
  
  //double AIC=0.0;
  //AIC=2*rrnewDFV-2*rrnewloglikeV;
  return List::create(EstBetaMat,BetaVarMat,Lambda);//,AIC,rrnewDFV,rrnewloglikeV
}



//Estimate beta.1 and beta.2 jointly
// [[Rcpp::export]]

List getLCcohS2 (arma::mat betaMat1,arma::mat betaMat2,
                    arma::vec lambdaini1,arma::vec lambdaini2,
                    arma::mat covat0,arma::mat covat0z,
                    arma::mat Yc,arma::mat Yc1,
                    arma::mat Nstar,double assumedtiny, 
                    int reptmax, double etapara, int timeUN,
                    arma::mat count,arma::vec FlagBirth, 
                    arma::vec a1star,arma::vec CL, //arma::vec CR,
                    arma::vec windowKH, //-AC
                    int tauL1, int tauR1,
                    int tauL2, int tauR2, int comCoef,
                    int comLamb){ //, int calcvar,int nSubj0
  // AC: I added the YcYsz, which is the Y_c*Y_s in the ratio, Z_bar.
  int nSubj=covat0.n_rows;
  int nSubjz=covat0z.n_rows;
  //int windowKH=0;-AC
  int nCov=covat0z.n_cols;
  
  arma::vec S0t1=arma::zeros<arma::vec>(timeUN);
  arma::vec S0t2=arma::zeros<arma::vec>(timeUN);
  
  arma::mat Yc1pz,Yc2pz,Yc1p,Yc2p;
  arma::vec YN1, YN2;   
  
  arma::vec YZN1=arma::zeros<arma::vec>(nCov);
  arma::vec YZN2=arma::zeros<arma::vec>(nCov);
  
  arma::mat YZNW1=arma::zeros<arma::mat>(timeUN,nCov);//dim: # time * # covariate
  arma::mat YZNW2=arma::zeros<arma::mat>(timeUN,nCov);
  
  arma::mat EIs1z=arma::zeros<arma::mat>(nSubjz,timeUN);
  arma::mat EIs2z=arma::zeros<arma::mat>(nSubjz,timeUN);
  
  arma::mat EstBetaMat01=arma::zeros<arma::mat>(nCov,timeUN);
  arma::mat EstBetaMat02=arma::zeros<arma::mat>(nCov,timeUN);
  arma::mat EstBetaMat1=arma::zeros<arma::mat>(nCov,timeUN);
  arma::mat EstBetaMat2=arma::zeros<arma::mat>(nCov,timeUN);
  
  //arma::mat BetaVarMat=arma::zeros<arma::mat>(nCov,timeUN);
  arma::vec lambda1=lambdaini1;
  arma::vec lambda2=lambdaini2;
  arma::mat hr1,hr2,hr1z;//,hr2z
  
  double betaRelDiff1,betaRelDiff2;
  betaRelDiff1=1000;
  betaRelDiff2=1000;
  int m=0;
  int j,l;
  
  while(((betaRelDiff1>(50*etapara))||(betaRelDiff2>(50*etapara)))&&(m<reptmax))
  {
    if(m==0){
      EstBetaMat01=arma::trans(betaMat1);
      EstBetaMat02=arma::trans(betaMat2);
    }else{
      EstBetaMat01=EstBetaMat1;
      EstBetaMat02=EstBetaMat2;
    }
    
    //a. Generate harzard functions for strata 1 & 2
    arma::mat tmplambda1=arma::zeros<arma::mat>(nSubj,timeUN);
    tmplambda1.each_row()+=arma::trans(lambda1);
    arma::mat tmplambda2=arma::zeros<arma::mat>(nSubj,timeUN);
    tmplambda2.each_row()+=arma::trans(lambda2);
    arma::mat tmplambda1z=arma::zeros<arma::mat>(nSubjz,timeUN);
    tmplambda1z.each_row()+=arma::trans(lambda1);
    
    hr1=arma::exp(covat0*EstBetaMat01)%tmplambda1;//dim: n*time unit
    hr2=arma::exp(covat0*EstBetaMat02)%tmplambda2;
    hr1z=arma::exp(covat0z*EstBetaMat01)%tmplambda1z;//dim: n*time unit
    //hr2z=arma::exp(covat0z*beta02)*arma::trans(lambda2);
    
    
    arma::vec Nr=arma::zeros<arma::vec>(nSubj);
    arma::vec Dr=arma::zeros<arma::vec>(nSubj);
    Yc1p=Yc1;
    
    //Rcout<<"Flag 1\n";
    for(int i=0;i<nSubj;i++){
      if(FlagBirth(i)==0){
        //born before the study
        
        //a1 is not missing, which means that we observed at least one event.
        double tmpAC=0.0,tmpAC1=0.0,tmpAC2=0.0;
        
        for(int j=0; j<a1star(i); j++){
          tmpAC=tmpAC-hr1(i,j);
          if(j<CL(i)){
            tmpAC1=tmpAC1-hr1(i,j);
          }else{
            tmpAC2=tmpAC2-hr2(i,j);
          }
        }
        Nr(i)=hr1(i,a1star(i)-1)*std::exp(tmpAC);
        Dr(i)=Nr(i)+(1-std::exp(tmpAC1))*hr2(i,a1star(i)-1)*std::exp(tmpAC2);
        
        Yc1p.row(i)=Yc1.row(i)*(Nr(i)/Dr(i));
      }
    }
    //Rcout<<"Flag 2\n";
    Yc2p=Yc-Yc1p;
    
    EIs1z=arma::zeros<arma::mat>(nSubjz,timeUN);
    EIs2z=arma::zeros<arma::mat>(nSubjz,timeUN);
    
    for(int i=0;i<nSubjz;i++){
      double tmpACz=0.0;
      EIs1z(i,0)=1.0;
      for(int j=1; j<timeUN; j++){
        tmpACz=tmpACz-hr1z(i,j-1);
        EIs1z(i,j)=std::exp(tmpACz);
      }
    }
    EIs2z=1.0-EIs1z;
    
    Yc1pz=arma::trans(EIs1z)%Nstar;
    Yc2pz=arma::trans(EIs2z)%Nstar;
    
    //a. get Y(j)*dN(j): total obs visits at age j ---
    //add matrix betaDiffcount for each row (dim=1)
    YN1=arma::sum(arma::trans(count%Yc1p),1);  
    YN2=arma::sum(arma::trans(count%Yc2p),1); 
    
    arma::mat YZN1=covat0.t()*(count%Yc1p);//nCov*timeUN
    arma::mat YZN2=covat0.t()*(count%Yc2p);//nCov*timeUN
    
    //b2. get Y(.)Z(k)dN(.)W(l,j):weighted total obs visits at age j around l
    YZNW1=arma::zeros<arma::mat>(timeUN,nCov);//dim: # time * # covariate
    YZNW2=arma::zeros<arma::mat>(timeUN,nCov);
    
    int windowKH00;
    for(l=(tauL1-1); l<tauR1; l++)
    {
      windowKH00=windowKH(l);
      for(j=0;j<timeUN;j++)
      {
        if((j<=(l+windowKH00))&&(j>=(l-windowKH00)))
        {
          YZNW1.row(l)+=0.75/windowKH00*(1-(j-l)*(j-l)/windowKH00/windowKH00)*YZN1.col(j).t();
          
        }
      }
    }
    
    for(l=(tauL2-1); l<tauR2; l++)
    {
      windowKH00=windowKH(l);
      
      for(j=0;j<timeUN;j++)
      {
        if((j<=(l+windowKH00))&&(j>=(l-windowKH00)))
        {
          YZNW2.row(l)+=0.75/windowKH00*(1-(j-l)*(j-l)/windowKH00/windowKH00)*YZN2.col(j).t();
        }
      }
    }
    
    arma::vec betaini1,betaini2;
    int windowKH0;
    
    if(comCoef==0){
      //fit until converge or reach max iteration
      // iter is the target time point
      for (int iter=tauL1-1;iter<tauR1;iter++)
      {
        //Rcout<<"iter:"<<iter<<"\n";
        betaini1=EstBetaMat01.col(iter);
        windowKH0=windowKH(iter);
        //std::cout << iter<<std::endl;
        LCreturnva1=AGcoha(betaini1, covat0,covat0z, Yc1p,Yc1pz,YN1,
                          YZNW1,assumedtiny,reptmax,etapara, timeUN,iter,
                          windowKH0,nSubj,nCov);
        arma::vec beta1=LCreturnva1.v1;
         
         EstBetaMat1.col(iter)=beta1;//dim: # covariate * # time
      }
      
      for (int iter=tauL2-1;iter<tauR2;iter++)
      {
        //Rcout<<"iter_s2:"<<iter<<"\n";
        betaini2=EstBetaMat02.col(iter);
        windowKH0=windowKH(iter);
        //std::cout << iter<<std::endl;
        
        LCreturnva2=AGcoha(betaini2, covat0,covat0z, Yc2p,Yc2pz,YN2,
                          YZNW2,assumedtiny,reptmax,etapara, timeUN,iter,
                          windowKH0,nSubj,nCov);
        arma::vec beta2=LCreturnva2.v1;
         
         EstBetaMat2.col(iter)=beta2;
      }
    }else{
      for (int iter=(std::min(tauL1,tauL2)-1);iter<(std::max(tauR1,tauR2));
      iter++)
      {
        //Rcout<<"iter:"<<iter<<"\n";
        betaini1=EstBetaMat01.col(iter);
        windowKH0=windowKH(iter);
        //std::cout << iter<<std::endl;
        LCreturnva1=AGcoha(betaini1, covat0,covat0z, (Yc1p+Yc2p),(Yc1pz+Yc2pz),
                          (YN1+YN2),(YZNW1+YZNW2),assumedtiny,reptmax,etapara, 
                          timeUN,iter,windowKH0,nSubj,nCov);
        arma::vec beta1=LCreturnva1.v1;
        //arma::vec beta2=beta1;
        
        
        EstBetaMat1.col(iter)=beta1;//dim: # covariate * # time
        //EstBetaMat2.col(iter)=beta2;//dim: # covariate * # time
      }
    }
    
    
    //Extend the estimates of betas to [0,timeUN-1]
    if(comCoef==0){
      for(l=0;l<tauL1-1;l++){
        EstBetaMat1.col(l)=EstBetaMat1.col(tauL1-1);
      }
      
      for(l=0;l<tauL2-1;l++){
        EstBetaMat2.col(l)=EstBetaMat2.col(tauL2-1);
      }
      
      for(l=tauR1;l<timeUN;l++){
        EstBetaMat1.col(l)=EstBetaMat1.col(tauR1-1);
      }
      
      for(l=tauR2;l<timeUN;l++){
        EstBetaMat2.col(l)=EstBetaMat2.col(tauR2-1);
      }
    }else{
      for(l=0;l<(std::min(tauL1,tauL2)-1);l++){
        EstBetaMat1.col(l)=EstBetaMat1.col(std::min(tauL1,tauL2)-1);
      }
      
      for(l=(std::max(tauR1,tauR2));l<timeUN;l++){
        EstBetaMat1.col(l)=EstBetaMat1.col(std::max(tauR1,tauR2)-1);
      }
      EstBetaMat2=EstBetaMat1;
    }
    
    
    lambda1=arma::zeros<arma::vec>(timeUN);
    lambda2=arma::zeros<arma::vec>(timeUN);
    
    /*6. get Lambda0 estm*/
    
    S0t1=arma::sum(Yc1pz%arma::exp(arma::trans(covat0z*EstBetaMat1)),1);//dim: time *1
    S0t2=arma::sum(Yc2pz%arma::exp(arma::trans(covat0z*EstBetaMat2)),1);
    // The previous S0all is for estimated beta at each time unit.
    // We need to recalculate it for the final estimation results
    
    if(comLamb==0){
      lambda1=YN1/(S0t1+assumedtiny);
      lambda2=YN2/(S0t2+assumedtiny);
      //the denominator is 0 so the lambda value is extremely large.
      lambda2(0)=0;
    }else{
      lambda1=(YN1+YN2)/(S0t1+S0t2+assumedtiny);
      lambda2=lambda1;
    }
    
    betaRelDiff1=arma::accu(abs((EstBetaMat01-EstBetaMat1)/EstBetaMat01))/(nCov*timeUN);
    betaRelDiff2=arma::accu(abs((EstBetaMat02-EstBetaMat2)/EstBetaMat02))/(nCov*timeUN);
    
    m++;
    
  }
  //arma::mat betaZ,expbetaZ;
  
  if(m>=reptmax&&m>2)
  {Rcpp::Rcout << "getting AG beta does not converge at m="<< m << std::endl;}
  
  arma::vec Lambda1=arma::cumsum(lambda1);
  arma::vec Lambda2=arma::cumsum(lambda2);

  //Take its transpose  
  EstBetaMat1= EstBetaMat1.t();//dim: # time * # covariate
  EstBetaMat2= EstBetaMat2.t();

  return List::create(EstBetaMat1,EstBetaMat2,Lambda1,Lambda2,
                      Yc1p,Yc2p,Yc1pz,Yc2pz,S0t1,S0t2);//,AIC,rrnewDFV,rrnewloglikeV
}
