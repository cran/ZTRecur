# Internal function: Hidden from end users
find_CL=function(x){
  return(which(x>0)[1])
}


find_CR=function(x){
  ind=sum(x>0)
  return(which(x>0)[ind])
}


cumHZ=function(hr,Lbound,Rbound){
  #dimension of hr is #subjects * # time units
  res=rep(NA,nrow(hr))
  for(i in 1:nrow(hr)){
    if(Lbound[i]<=Rbound[i]){
      res[i]=sum(hr[i,Lbound[i]:Rbound[i]]) 
    }else{
      res[i]=0
    }
  }
  return(as.vector(res))
}


gen.lambda=function(Lambda,ONUM,AGEUNIT,mat=FALSE,lambda.var=FALSE){
  if(!mat){
    lambda=c(Lambda[1],(Lambda[2:ONUM]-Lambda[1:(ONUM-1)]))*AGEUNIT
  }else{
    lambda=cbind(Lambda[,1],(Lambda[,2:ONUM]-Lambda[,1:(ONUM-1)]))*AGEUNIT
  }
  
  if(lambda.var){
    res=list(lambda=lambda,lambda.var=apply(lambda,2,var))
  }else{
    res=lambda
  }
  return(res)
}


##################################
#' Stratified Regression Analysis for Zero-Truncated Recurrent Event Data
#'
#' Fits a stratified Cox-type model with time-constant or time-varying 
#' regression coefficients using zero-truncated recurrent event data integrated 
#' with population census information. The method accounts for coarsened follow-up 
#' information and providing bootstrap or multiplier resampling variance estimates.
#'
#' @param base Matrix of subject-level information. Each row corresponds to one
#' subject and contains demographic information (covariates and either the exact
#' birthdate (\code{birthJulian}) or birthdate interval (\code{birthl} and
#' \code{birthr})), the number of observed events (\code{nmvst}), and the row 
#' index of the subject's first observed event in the \code{visit} matrix
#' (\code{p1stvt}). The variables must be provided in the following order:
#' covariates, \code{nmvst}, \code{p1stvt}, and either \code{birthJulian} or 
#' \code{birthl} and \code{birthr}.
#' 
#' @param visit Matrix of event-level information. Each row corresponds to one
#' observed event and contains the subject ID (\code{patid}), an indicator of
#' whether the event is a subsequent visit (\code{reptind}; 0 for the first
#' observed visit and 1 for subsequent visits for the same subject), the Julian
#' date of the event (\code{date}), and the integer age at the time of the event
#' (\code{age}). Subject IDs must be integers from 1 through the number of rows
#' in the \code{base} matrix. If birthdate information is unavailable, the
#' integer age (\code{age}) must be provided for each event.
#' 
#' @param AGEUNIT Integer; the number of time units per year (e.g., 6 for 
#' bimonthly data).
#' @param ONUM Integer; the maximum number of time units under observation 
#' (e.g., 18 * AGEUNIT for <18 years).
#' @param WINDOWL Numeric; Julian date representing the left boundary of the 
#' data extraction window.
#' @param WINDOWR Numeric; Julian date representing the right boundary of the 
#' data extraction window.
#' @param etapara Numeric; tolerance parameter controlling the convergence
#' criterion when solving the estimating equations.
#' @param reptmax Integer; the maximum number of iterations allowed when solving 
#' estimating equations.
#' @param windowKH Numeric; bandwidth size of the smoothing kernel, specified in 
#' time units.
#' @param tauL Numeric (in time units); left boundary used for local constant 
#' estimation to avoid boundary issues.
#' @param tauR Numeric (in time units); right boundary used for local constant 
#' estimation to avoid boundary issues.
#' @param method Character string specifying the estimation method. Options are
#' estimation under the Andersen-Gill model (\code{"L"}), local constant estimation 
#' under a Cox-type model with time-varying coefficients (\code{"LC"}), or both
#' estimation procedures (\code{"all"}).
#' @param adjust.LC Logical; if \code{TRUE}, applies boundary adjustment using 
#' \code{tauR.1} and \code{tauL.2} for stratum-specific data sparsity. Default 
#' is \code{FALSE}.
#' @param tauR.1 Numeric (in time units); right boundary parameter for stratum 1 
#' when \code{adjust.LC = TRUE}.
#' @param tauL.2 Numeric (in time units); left boundary parameter for stratum 2 
#' when \code{adjust.LC = TRUE}.
#' @param N_star Data frame containing aggregated population census data, 
#' including integer age (\code{AGE}),mycovars, and number of subjects in that 
#' subgroup (\code{pop}).
#' @param n Integer; the approximate size of $O$, the random sample drawn from the 
#' target population.
#' @param B Integer or \code{NA}; number of bootstrap or multiplier-resampling
#' replicates used for variance estimation. Defaults to \code{NA}, in which case
#' no variance estimator is calculated.
#' @param seeds Numeric; seed value for reproducible sampling. Default is 
#' \code{1234}.
#' @param com.coef Binary; \code{0} (default) allows stratum-specific regression 
#' coefficients, while \code{1} constrains coefficients to be common across strata.
#' @param com.lamb Binary; \code{0} (default) allows stratum-specific baseline 
#' functions, while \code{1} constrains baselines to be common across strata.
#' @param knownBirth Logical; indicates whether exact birthdates are known. Defaults
#' to \code{FALSE}. If \code{FALSE}, exact birthdates are unknown, and
#' \code{birth.interval} can be used to calculate the corresponding birthdate
#' intervals in Julian dates.
#' @param reptw Integer; number of birthdate sets generated for the missing 
#' birthdate procedure A (Hu & Rosychuk, 2016). Default is \code{100}.
#' @param resampling Logical; if \code{TRUE}, applies the multiplier resampling 
#' method described in Chen, Hu, & Rosychuk (2026). If \code{FALSE} (default), 
#' standard bootstrapping is used.
#' @param r.distn Character string specifying the multiplier distribution used for
#' resampling when \code{resampling = TRUE}. Options are the standard normal
#' distribution (\code{"std Normal"}; default) or a Poisson distribution with
#' variance 1 (\code{"Poisson"}).
#' @param mycovars Character vector specifying the variable names of covariates 
#' to define subgroups.
#' @param save.res Logical; if \code{TRUE}, writes output files to disk. 
#' Default is \code{FALSE}.
#' @param resultpath Character string indicating the file path destination 
#' for saved results, such as "C:\\project\\". Required when \code{save.res = TRUE}.
#' @param output.lambda Logical; if \code{TRUE}, includes baseline
#' function estimates in the returned output. Default is \code{FALSE}.
#' 
#' @return A list containing estimated stratified regression coefficients, 
#' standard error estimates (bootstrap- or resampling-based), estimated 
#' stratified cumulative baseline functions, and estimated stratified baseline 
#' functions when \code{output.lambda = TRUE}.
#'   
#' @details
#'  The detailed estimation procedures are described in Chen, Hu & Rosychuk 
#'  (2026a) for time-constant regression coefficients and in Chen et al. (2026b) 
#'  for time-varying regression coefficients. The approach accounts for missing 
#'  birthdates using Procedure A of Hu & Rosychuk (2016).
#' 
#' @references
#' Chen, A. A., Hu, X. J., & Rosychuk, R. J. (2026a). Stratified regression 
#' analysis of zero-truncated recurrent event data. \emph{The Annals of Applied 
#' Statistics}, 20(1), 434-451.
#' 
#' Chen, A. A., Hu, X. J., Rosychuk, R. J., & Zeng, L. (2026b). Analyzing 
#' zero‐truncated recurrent event data by stratified regression with time‐varying 
#' coefficients. \emph{Canadian Journal of Statistics}, e70064.
#' 
#' Hu, X. J., & Rosychuk, R. J. (2016). Imputation procedures for missing 
#' birthdates in administrative health data. \emph{Statistics in Medicine}, 
#' 35(14), 2341-2355.
#' 
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' \donttest{
#' # When birthdates are unknown, obtain resampling-based variance estimates
#' # using standard normal multipliers..
#' strati.fit.uk.R <- strati.reg(base = my.base[,c(-1,-9)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 0, com.lamb = 0, knownBirth = FALSE, 
#'   reptw = 100, resampling = TRUE,r.distn = "std Normal",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#'   # When birthdates are known, obtain resampling-based variance estimates
#'   # using Poisson multipliers with variance 1.
#'   strati.fit.k.R.p <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 0, com.lamb = 0, knownBirth = TRUE, 
#'   reptw = 100, resampling = TRUE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#'   # When birthdates are known, obtain bootstrap-based variance estimates.
#'   strati.fit.k.B <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 0, com.lamb = 0, knownBirth = TRUE, 
#'   reptw = 100, resampling = FALSE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#'   # Fit models with common regression coefficients across strata.
#'   strati.fit.k.R.p.combreg <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 1, com.lamb = 0, knownBirth = TRUE, 
#'   reptw = 100, resampling = TRUE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#'   # Fit models with a common baseline function across strata.
#'   strati.fit.k.R.p.combB <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 0, com.lamb = 1, knownBirth = TRUE, 
#'   reptw = 100, resampling = TRUE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#' }
#' 
strati.reg=function(base, visit,
                    AGEUNIT,ONUM,WINDOWL,WINDOWR,etapara,reptmax,
                    windowKH,tauL,tauR,method,adjust.LC=FALSE,
                    tauR.1,tauL.2,N_star,n,B=NA,seeds,
                    com.coef=0, com.lamb=0,knownBirth=FALSE,reptw=100,resampling=FALSE,
                    r.distn="std Normal",mycovars,save.res=FALSE,resultpath=NULL,
                    output.lambda=FALSE)
{ 
  if(is.matrix(base)==FALSE || is.matrix(visit)==FALSE){
    stop("Base and/or visit are not passed as matrices. Stopping.\n");
  }
  
  TOTALM=dim(base)[1]
  TOTALM2=dim(visit)[1]
  assumedtiny=1.0e-10   # may be taken as zero
  
  # number of covariates Z
  if(knownBirth){
    dimK=dim(base)[2]-3 # number of covariates
  }else{
    dimK=dim(base)[2]-4 # number of covariates
  }
  
  message("The input covariate matrix has dimension ",TOTALM," by ",dimK,"\n")
  message("The input event matrix has dimension ",TOTALM2," by ",dim(visit)[2],"\n")
  
  covat0=base[,1:dimK]
  
  nmvst=base[,dimK+1]
  p1stvt=base[,dimK+2]
  
  if(knownBirth){
    birth=base[,dimK+3]
  }else{
    birthl=base[,dimK+3]
    birthr=base[,dimK+4]
  }
  
  patid=visit[,1]
  reptind=visit[,2]
  dates=visit[,3]
 
  #########################################
  #method="L","LC"
  #method="all" calls each of the 2 methods; ="L" corresponds to the estimation under AG model;
  # ="LC" corresponds to the local constant estimation under the model with age-varying coefficients.
  t1=Sys.time()
  if(knownBirth){
    a=getmoredata(birth,WINDOWL,WINDOWR,dates,p1stvt,
                  TOTALM,TOTALM2,ONUM,nmvst,AGEUNIT)
  }else{
    a=getmoredataComb2(birthl,birthr,WINDOWL,WINDOWR,dates,p1stvt,
                       reptw,TOTALM,TOTALM2,ONUM,nmvst,AGEUNIT)
  }
  t2=Sys.time()
  t=t2-t1
  #str(t)
  message(paste0("Prepare data for modelling, which takes ",round(t,3)," ",attr(t, "units"),".\n"))
  #rm(t)

  TMPcount=a[[1]]
  tmpY=a[[2]]
  #summary(apply(tmpY,1,sum))
  # Indicator for each stratum
  tmpYs1=a[[4]]
  tmpYs2=tmpY-tmpYs1
  
  # we should remove TMPvisitage=108, which is greater than 18 yrs old
  TMPvisitage=a[[3]]
  # the first observed event time in time units
  a1_star=TMPvisitage[p1stvt,]
  #table(a1_star,useNA = "ifany")
  # left censoring time in time units
  CL=apply(tmpY,1,find_CL)
  # right censoring time in time units
  CR=apply(tmpY,1,find_CR)
  #all(CR>=CL)
  #summary(CR-CL)
  
  if(knownBirth){
    flag.Born.in.study=as.numeric(tmpY[,1]==1)
  }else{
    flag.Born.in.study=as.numeric(tmpY[,1]>0.5)
  }
  #table(flag.Born.in.study,useNA = "ifany")
  
  if(!knownBirth){
    dN=a[[5]]
    #range(dN)
    #all(TMPcount==dN)
  }else{
    dN=TMPcount*tmpY
    #all(TMPcount*tmpY==TMPcount)
  }
  
  if(is.element(method,c("L","LC","all"))==FALSE){
    warning("ERROR: Need to have one of L, LC, or all for the method.\n\n")
    return(list(AG1=NA,AG1LC=NA))
  }
  
  N_star.c=N_star%>%dplyr::group_by(AGE,across(all_of(mycovars)))%>%dplyr::summarise(pop=sum(pop),.groups="drop")
  N_star.c$comb=apply(N_star.c[,mycovars],1,paste0,collapse="")
  N_star.wide=reshape(data=as.data.frame(N_star.c[,c("AGE","comb","pop")]),idvar = "AGE",timevar = "comb",
                      direction = "wide")
  covat0z=as.matrix(unique(N_star.c[,mycovars]))
  N_star.wide=N_star.wide[order(N_star.wide$AGE),]
  rownames(N_star.wide)=N_star.wide$AGE
  N_star.wide$AGE=NULL
  N_star.long=apply(N_star.wide,2,rep,each=6)
  rm(N_star.c,N_star.wide)
  
  #initialize the estimators
  betaini=rep(0,ncol(covat0))
    
  if(method=="all"|| method=="L")
  {
    message("Estimation under an AG model\n")
    if(!exists("getAG1coh"))
      stop("Compilation failed. The function getAG1coh does not exist")
    
    if(!exists("getAG1cohS2"))
      stop("Compilation failed. The function getAG1cohS2 does not exist")
    
      b=getAG1coh(betaini,covat0,covat0z,tmpY,N_star.long,assumedtiny,
                reptmax,etapara,ONUM,0,n,dN)
    
    
      beta.1=beta.2=b[[1]]
      Lambda.1=b[[3]]
    
    # No need to *AGEUNIT. it is dLambada actually.
      lambda.1=as.vector(c(Lambda.1[1,1],(Lambda.1[2:nrow(Lambda.1),1]-
                                          Lambda.1[1:(nrow(Lambda.1)-1),1])))
    
      lambda.2=lambda.1
  
      b=getAG1cohS2(beta.1,beta.2,lambda.1,lambda.2,
                  covat0,covat0z,tmpY,tmpYs1,
                  N_star.long,assumedtiny,
                  reptmax,etapara,ONUM,TMPcount,
                  flag.Born.in.study,a1_star,
                  CL,0,n,com.coef,com.lamb)#CR,
       
    # etapara*1000<0.1
    if(mean(abs(b[[7]]/n))>etapara){
      stop("The estimating function for stratum 1 is not close to 0!")
    }
    
    if(mean(abs(b[[8]]/n))>etapara){
      stop("The estimating function for stratum 2 is not close to 0!")
    }
    
    if(b[[15]]>=reptmax){
      stop("It doesn't converge!")
    }

    
    if(!is.na(B)){
      t1=Sys.time()
      set.seed(seeds)
      beta1_mat=beta2_mat=matrix(NA,nrow=B,ncol=length(b[[1]]))
      Lambda1_mat=Lambda2_mat=matrix(NA,nrow=B,ncol=ONUM)
      
      bt=1
      while(bt<=B){
        message(paste0("AG: Boot B=", bt,"\n"))
        if(resampling){
          if(r.distn=="std Normal"){
            lst_w=rnorm(n=nrow(covat0))+1
          }else if(r.distn=="Poisson"){
            lst_w=rpois(n=nrow(covat0),lambda = 1)
          }
          #mean(lst_w)
          #var(lst_w)
          mat_w=matrix(rep(lst_w,ONUM),ncol=ONUM)
        }else{
          lst.sub=sample(x=c(1:nrow(covat0)),size=nrow(covat0),replace = TRUE)
        }
        
        tryCatch({
          if(resampling){
            boot=getAG1cohS2(beta.1,beta.2,lambda.1,lambda.2,
                             covat0,covat0z,tmpY,tmpYs1,
                             N_star.long,assumedtiny,
                             reptmax,etapara,ONUM,TMPcount*mat_w,
                             flag.Born.in.study,a1_star,
                             CL,0,n,com.coef,com.lamb)
          }else{
            boot=getAG1cohS2(beta.1,beta.2,lambda.1,lambda.2,
                             covat0[lst.sub,],covat0z,tmpY[lst.sub,],tmpYs1[lst.sub,],
                             N_star.long,assumedtiny,
                             reptmax,etapara,ONUM,TMPcount[lst.sub,],
                             flag.Born.in.study[lst.sub],a1_star[lst.sub],
                             CL[lst.sub],0,n,com.coef,com.lamb)
          }
          
          if(mean(abs(boot[[7]]/n))>etapara){
            stop("The estimating function for stratum 1 is not close to 0!")
          }
          
          if(mean(abs(boot[[8]]/n))>etapara){
            stop("The estimating function for stratum 2 is not close to 0!")
          }
          
          if(boot[[15]]>=reptmax){
            stop("It doesn't converge!")
          }
          
          beta1_mat[bt,]=boot[[1]]
          beta2_mat[bt,]=boot[[2]]
          Lambda1_mat[bt,]=boot[[3]]
          Lambda2_mat[bt,]=boot[[4]]
          
          res.AG=list(b=b,beta1_mat=beta1_mat, beta2_mat=beta2_mat,
                   Lambda1_mat=Lambda1_mat, Lambda2_mat=Lambda2_mat)
          
          if(save.res){
            if (is.null(resultpath) || length(resultpath) != 1L ||
                is.na(resultpath) || !nzchar(resultpath)) {
              stop("Please provide an explicit output directory in 'resultpath'.")
            }
            
            if (!dir.exists(resultpath)) {
              stop("The output directory does not exist: ", resultpath)
            }
            
            save(res.AG, file=paste(resultpath,"Stratified_Result_constant_coef",
                                    ifelse(resampling,"_R","_B"),B,".RData",sep=""))
          }
          bt=bt+1
          },
        error=function(e){warning("ERROR: ",conditionMessage(e),"\n")})
      }
      t2=Sys.time()
      t=t2-t1
      message(paste0("Resampling method for estimating variances (constant coefs) is done, which takes ",
                 round(t,3)," ",attr(t, "units"),".\n"))
      rm(t)
      
      if(output.lambda){
        AG1_s1=list(beta=b[[1]],betavar=cov(beta1_mat),Lambda=b[[3]],
                    Lambdavar=apply(Lambda1_mat,2,var),
                    lambda=gen.lambda(Lambda=b[[3]],ONUM,AGEUNIT,mat=FALSE,lambda.var=FALSE),
                    lambdavar=gen.lambda(Lambda=Lambda1_mat,ONUM,AGEUNIT,mat=TRUE,lambda.var=TRUE)$lambda.var,
                    Pscore=b[[7]],
                    beta_mat=beta1_mat,Lambda_mat=Lambda1_mat)  
        # cov(beta2_mat)=sqrt(diag(cov(beta2_mat)))
        AG1_s2=list(beta=b[[2]],betavar=cov(beta2_mat),Lambda=b[[4]],
                    Lambdavar=apply(Lambda2_mat,2,var),
                    lambda=gen.lambda(Lambda=b[[4]],ONUM,AGEUNIT,mat=FALSE,lambda.var=FALSE),
                    lambdavar=gen.lambda(Lambda=Lambda2_mat,ONUM,AGEUNIT,mat=TRUE,lambda.var=TRUE)$lambda.var,
                    Pscore=b[[8]],
                    beta_mat=beta2_mat,Lambda_mat=Lambda2_mat)
      }else{
        AG1_s1=list(beta=b[[1]],betavar=cov(beta1_mat),Lambda=b[[3]],
                    Lambdavar=apply(Lambda1_mat,2,var),
                    Pscore=b[[7]],
                    beta_mat=beta1_mat,Lambda_mat=Lambda1_mat)  
        
        AG1_s2=list(beta=b[[2]],betavar=cov(beta2_mat),Lambda=b[[4]],
                    Lambdavar=apply(Lambda2_mat,2,var),
                    Pscore=b[[8]],
                    beta_mat=beta2_mat,Lambda_mat=Lambda2_mat)
      }
      
    }else{
      AG1_s1=list(beta=b[[1]],betavar=NA,Lambda=b[[3]],Lambdavar=NA,
                  Pscore=b[[7]])  

      AG1_s2=list(beta=b[[2]],betavar=NA,Lambda=b[[4]],Lambdavar=NA,
                  Pscore=b[[8]])
    }
  }
  

  if(method=="all"|| method=="LC")
  {
    if(!exists("getLCcoh"))
      stop("Compilation failed. The function getLCcoh does not exist")
    
    if(!exists("getLCcohS2"))
      stop("Compilation failed. The function getLCcohS2 does not exist")
    
    if(tauR>ONUM)
      stop("tauR should be less than ONUM: the largest age units")
    
    message("Local Constant Estimation: AG1LC\n")
    #initialize the estimators
    betaMat1=betaMat2=matrix(rep(betaini,ONUM),ncol = ncol(covat0z),byrow = TRUE)
    
    blc=getLCcoh(betaMat1,covat0,covat0z,tmpY,N_star.long,assumedtiny,
                    reptmax,etapara,ONUM,windowKH,tauL,tauR,0,n,dN)
    betaMat1=betaMat2=blc[[1]]
    # cumulative baseline
    Lambda.1=blc[[3]]
    # baseline HR
    lambda.1=as.vector(c(Lambda.1[1,1],(Lambda.1[2:nrow(Lambda.1),1]-
                                          Lambda.1[1:(nrow(Lambda.1)-1),1])))
    
    lambda.2=lambda.1

      
    blc=getLCcohS2(betaMat1,betaMat2,lambda.1,lambda.2,
                        covat0,covat0z,tmpY,tmpYs1,
                        N_star.long,assumedtiny,
                        reptmax,etapara,ONUM,TMPcount,
                        flag.Born.in.study,a1_star,
                        CL,windowKH,tauL,ifelse(adjust.LC,tauR.1,tauR),
                        ifelse(adjust.LC,tauL.2,tauL),tauR,
                        com.coef,com.lamb)
      
      if(!is.na(B)){
        t1=Sys.time()
        set.seed(seeds)
        
        beta1_mat=beta2_mat=array(NA,dim = c(dim(blc[[1]]),B))
        Lambda1_mat=Lambda2_mat=matrix(NA,nrow=B,ncol=ONUM)
        
        bt=1
        while(bt<=B){
          message(paste0("LC: Boot B=", bt,"\n"))
          
          if(resampling){
            if(r.distn=="std Normal"){
              lst_w=rnorm(n=nrow(covat0))+1
            }else if(r.distn=="Poisson"){
              lst_w=rpois(n=nrow(covat0),lambda = 1)
            }
            #mean(lst_w)
            #var(lst_w)
            mat_w=matrix(rep(lst_w,ONUM),ncol=ONUM)
          }else{
            lst.sub=sample(x=c(1:nrow(covat0)),size=nrow(covat0),replace = TRUE)
          }
          
          tryCatch({
            if(resampling){
              bootlc=getLCcohS2(betaMat1,betaMat2,lambda.1,lambda.2,
                                  covat0,covat0z,tmpY,tmpYs1,
                                  N_star.long,assumedtiny,
                                  reptmax,etapara,ONUM,TMPcount*mat_w,
                                  flag.Born.in.study,a1_star,
                                  CL,windowKH,tauL,ifelse(adjust.LC,tauR.1,tauR),
                                  ifelse(adjust.LC,tauL.2,tauL),tauR,
                                  com.coef,com.lamb)
            }else{
              bootlc=getLCcohS2(betaMat1,betaMat2,lambda.1,lambda.2,
                                   covat0[lst.sub,],covat0z,tmpY[lst.sub,],
                                   tmpYs1[lst.sub,],
                                   N_star.long,assumedtiny,
                                   reptmax,etapara,ONUM,TMPcount[lst.sub,],
                                   flag.Born.in.study[lst.sub],a1_star[lst.sub],
                                   CL[lst.sub],windowKH,tauL,ifelse(adjust.LC,tauR.1,tauR),
                                   ifelse(adjust.LC,tauL.2,tauL),tauR,
                                   com.coef,com.lamb) 
            }
            
            beta1_mat[,,bt]=bootlc[[1]]
            Lambda1_mat[bt,]=bootlc[[3]]
            beta2_mat[,,bt]=bootlc[[2]]
            Lambda2_mat[bt,]=bootlc[[4]]
            
            res.LC=list(blc=blc,beta1_mat=beta1_mat, beta2_mat=beta2_mat,
                        Lambda1_mat=Lambda1_mat, Lambda2_mat=Lambda2_mat)
            
            if(save.res){
              if (is.null(resultpath) || length(resultpath) != 1L ||
                  is.na(resultpath) || !nzchar(resultpath)) {
                stop("Please provide an explicit output directory in 'resultpath'.")
              }
              
              if (!dir.exists(resultpath)) {
                stop("The output directory does not exist: ", resultpath)
              }
              
              save(res.LC, file=paste(resultpath,"Stratified_Result_age_varying_coef",
                                      ifelse(resampling,"_R","_B"),B,".RData",sep=""))
            }
            
            bt=bt+1
          },
          error=function(e){warning("ERROR: ",conditionMessage(e),"\n")})
        }
        t2=Sys.time()
        t=t2-t1
        message(paste0("Resampling method for estimating variances (age-varying coefs) is done, which takes ",round(t,3),
                   " ",attr(t, "units"),".\n"))
        rm(t)
        
        if(output.lambda){
          AG1LC_s1=list(beta=blc[[1]],betavar=apply(simplify2array(beta1_mat),c(1,2),var),
                        Lambda=blc[[3]],
                        Lambdavar=apply(Lambda1_mat,2,var),
                        lambda=gen.lambda(Lambda=blc[[3]],ONUM,AGEUNIT,mat=FALSE,lambda.var=FALSE),
                        lambdavar=gen.lambda(Lambda=Lambda1_mat,ONUM,AGEUNIT,mat=TRUE,lambda.var=TRUE)$lambda.var,
                        beta1_mat=beta1_mat,Lambda1_mat=Lambda1_mat)  
          #betavar=b[[6]]
          AG1LC_s2=list(beta=blc[[2]],betavar=apply(simplify2array(beta2_mat),c(1,2),var),
                        Lambda=blc[[4]],
                        Lambdavar=apply(Lambda2_mat,2,var),
                        lambda=gen.lambda(Lambda=blc[[4]],ONUM,AGEUNIT,mat=FALSE,lambda.var=FALSE),
                        lambdavar=gen.lambda(Lambda=Lambda2_mat,ONUM,AGEUNIT,mat=TRUE,lambda.var=TRUE)$lambda.var,
                        beta2_mat=beta2_mat,Lambda2_mat=Lambda2_mat)
        }else{
          AG1LC_s1=list(beta=blc[[1]],betavar=apply(simplify2array(beta1_mat),c(1,2),var),
                        Lambda=blc[[3]],
                        Lambdavar=apply(Lambda1_mat,2,var),
                        beta1_mat=beta1_mat,Lambda1_mat=Lambda1_mat)  
          #betavar=b[[6]]
          AG1LC_s2=list(beta=blc[[2]],betavar=apply(simplify2array(beta2_mat),c(1,2),var),
                        Lambda=blc[[4]],
                        Lambdavar=apply(Lambda2_mat,2,var),
                        beta2_mat=beta2_mat,Lambda2_mat=Lambda2_mat)
        }
        
      }else{
      AG1LC_s1=list(beta=blc[[1]],betavar=NA,Lambda=blc[[3]],Lambdavar=NA)
      #print(AG1LC_s1)
      AG1LC_s2=list(beta=blc[[2]],betavar=NA,Lambda=blc[[4]],Lambdavar=NA)
      }
 
  }

  if(!exists("AG1_s1"))
  {AG1_s1=NA}
  
  if(!exists("AG1_s2"))
  {AG1_s2=NA}
  
  if(!exists("AG1LC_s1"))
  {AG1LC_s1=NA}
  
  if(!exists("AG1LC_s2"))
  {AG1LC_s2=NA}
  
  RESULTA=list(AG1_s1=AG1_s1,AG1_s2=AG1_s2,AG1LC_s1=AG1LC_s1,
               AG1LC_s2=AG1LC_s2)
  return(RESULTA)
}



#' Plot Estimated Regression Coefficients from Stratified Cox-Type Models
#'
#' Displays estimates and pointwise confidence intervals for time-varying 
#' and time-constant regression coefficients from fitted stratified models.
#'
#' @param fit Object containing output from a stratified Cox-type model fit 
#'   (e.g., from \code{\link{strati.reg}}).
#' @param cov_name Character vector of covariate names matching the model 
#' specification.
#' @param i Integer; index of the covariate in \code{cov_name} to be plotted. 
#' Default is \code{1}.
#' @param age Numeric vector representing the age or time points corresponding to 
#' the estimates.
#' @param tauL Numeric (in time units); left boundary used for local constant 
#' estimation to avoid boundary issues.
#' @param tauR Numeric (in time units); right boundary used for local constant 
#' estimation to avoid boundary issues.
#' @param yRange Numeric vector of length 2 specifying the y-axis limits 
#' \code{c(ymin, ymax)}. If \code{NA} (default), limits are calculated automatically.
#' @param adjust.LC Logical; if \code{TRUE}, applies boundary adjustment using 
#' \code{tauR.1} and \code{tauL.2} for stratum-specific data sparsity. Default 
#' is \code{FALSE}.
#' @param tauR.1 Numeric (in time units); right boundary parameter for stratum 1 
#' when \code{adjust.LC = TRUE}.
#' @param tauL.2 Numeric (in time units); left boundary parameter for stratum 2 
#' when \code{adjust.LC = TRUE}.
#' @param comCoef Logical; indicates whether common regression coefficients across strata were 
#'   assumed in the model. Default is \code{FALSE}.
#' @param model.v Character string; legend label for the model with time-varying 
#' coefficients. Default is \code{"Age-varying coef"}.
#' @param model.c Character string; legend label for the model with time-constant
#'  coefficients. Default is \code{"Age-constant coef"}.
#' @param location Character string specifying the position of the legend (e.g., 
#' \code{"topright"}).
#' @param cex.axis Numeric magnification factor for axis annotation relative to 
#' current setting. Default is \code{1.2}.
#' @param cex.lab Numeric magnification factor for axis labels relative to current 
#' setting. Default is \code{1.3}.
#' @param cex.legend Numeric magnification factor for legend text. Default is 
#' \code{1.2}.
#' @param c Numeric; critical value used for constructing pointwise confidence 
#' intervals. Default is \code{1.96} (corresponding to a 95\% confidence interval).
#' @param xlab Character string specifying the x-axis label. Defaults to
#' \code{"Age (in years)"}.
#'
#' @return No return value; generates a stratum-specific coefficient plot for 
#' the specified covariate.
#' 
#' @details
#'  The function assumes that \code{method = "all"} was used in
#'  \code{\link{strati.reg}}.
#'  
#' @references
#' Chen, A. A., Hu, X. J., & Rosychuk, R. J. (2026a). Stratified regression 
#' analysis of zero-truncated recurrent event data. \emph{The Annals of Applied 
#' Statistics}, 20(1), 434-451.
#' 
#' Chen, A. A., Hu, X. J., Rosychuk, R. J., & Zeng, L. (2026b). Analyzing 
#' zero‐truncated recurrent event data by stratified regression with time‐varying 
#' coefficients. \emph{Canadian Journal of Statistics}, e70064.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The 'strati.fit.k.R.p' is the output from strati.reg().
#' strati.fit.k.R.p <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 0, com.lamb = 0, knownBirth = TRUE, 
#'   reptw = 100, resampling = TRUE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#' # Plot age-varying regression coefficient estimates for each stratum.
#' # The age range is adjusted for the local constant estimates in both strata.
#' plot_stratified_coef(fit = strati.fit.k.R.p, 
#' cov_name = c("MaleInd","CityAInd","CityBInd"), 
#'   i = 1, age = c(1:108)/6, tauL = 9, tauR = 105, yRange=c(-3.5,1),
#'   adjust.LC = TRUE, tauR.1=105,tauL.2=18,comCoef = FALSE)
#'   
#' plot_stratified_coef(fit = strati.fit.k.R.p, 
#' cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 2, age = c(1:108)/6, tauL = 9, tauR = 105, yRange=c(-2,1),
#' adjust.LC = TRUE, tauR.1=105,tauL.2=18,comCoef = FALSE, c = 1.96)
#' 
#' plot_stratified_coef(fit = strati.fit.k.R.p, 
#' cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 3, age = c(1:108)/6, tauL = 9, tauR = 105, yRange=c(-2,1),
#' adjust.LC = TRUE, tauR.1=105,tauL.2=18,comCoef = FALSE,c = 1.96)
#' 
#' # Plot common regression coefficient estimates across strata.
#' strati.fit.k.R.p.combreg <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 1, com.lamb = 0, knownBirth = TRUE, 
#'   reptw = 100, resampling = TRUE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#' plot_stratified_coef(fit = strati.fit.k.R.p.combreg, 
#' cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 1, age = c(1:108)/6, tauL = 9, tauR = 105, yRange=c(-2.5,-1),
#' adjust.LC = TRUE, tauR.1=105,tauL.2=18,comCoef = TRUE,c = 1.96)
#' 
#' plot_stratified_coef(fit = strati.fit.k.R.p.combreg, 
#' cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 2, age = c(1:108)/6, tauL = 9, tauR = 105, yRange=c(-1.4,-0.5),
#' adjust.LC = TRUE, tauR.1=105,tauL.2=18,comCoef = TRUE,c = 1.96)
#' 
#' plot_stratified_coef(fit = strati.fit.k.R.p.combreg, 
#' cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 3, age = c(1:108)/6, tauL = 9, tauR = 105, yRange=c(-2,-1),
#' adjust.LC = TRUE, tauR.1=105,tauL.2=18,comCoef = TRUE,c = 1.96)
#' }
plot_stratified_coef=function(fit,cov_name,i=1,age,tauL,tauR,yRange=NA,adjust.LC=FALSE,
                              tauR.1=NA,tauL.2=NA,comCoef=FALSE,
                              model.v="Age-varying coef",model.c="Age-constant coef",
                              location="topright",cex.axis=1.2,
                              cex.lab=1.3,cex.legend=1.2,c=1.96,
                              xlab="Age (in years)"){
  
  est.AG.s1=fit$AG1_s1$beta[i]
  est.AG.s2=fit$AG1_s2$beta[i]
  
  est.ext.s1=fit$AG1LC_s1$beta[,i]
  est.ext.s2=fit$AG1LC_s2$beta[,i]
  
  
  # Bootstrapping-based/resampling-based SE
  se.AG.s1=sqrt(diag(fit$AG1_s1$betavar))[i]
  se.AG.s2=sqrt(diag(fit$AG1_s2$betavar))[i]
  
  se.ext.s1=sqrt(fit$AG1LC_s1$betavar[,i])
  se.ext.s2=sqrt(fit$AG1LC_s2$betavar[,i])
  
  c.AG.s1=c.AG.s2=c.ext.s1=c.ext.s2=c
  
  if(any(is.na(yRange))){
    maxy2.s1=max(abs(est.ext.s1[tauL:tauR]+c.ext.s1*se.ext.s1[tauL:tauR]),
                 abs(est.ext.s1[tauL:tauR]-c.ext.s1*se.ext.s1[tauL:tauR]))
    maxy2.s2=max(abs(est.ext.s2[tauL:tauR]+c.ext.s2*se.ext.s2[tauL:tauR]),
                 abs(est.ext.s2[tauL:tauR]-c.ext.s2*se.ext.s2[tauL:tauR]))
    maxy=max(maxy2.s1,maxy2.s2)*1.1
    yRange=c(-maxy,maxy)
  }  

  
  # stratum 1
  ## AG1 VS AG1LC
  ## AG1##########################################################
  if(!comCoef){
    plot(c(min(age),max(age)),c(-ceiling(max(est.ext.s1[tauL:tauR])),
                                ceiling(max(est.ext.s1[tauL:tauR]))),
         type="n", ylim=yRange,xaxt="n",xlab=xlab,ylab="Coefficient estimate",
         sub=paste0(#ifelse(Bspline,"(Basis Function) ","(LC) "),
           "(a) ", cov_name[i],ifelse(comCoef," (Common Coefficient)"," (stratum 1)")),
         cex.lab=cex.lab,cex.axis=cex.axis,las=1)
    axis(side=1,at=seq(0,18,2),cex.axis=cex.axis)
    lines(c(min(age),max(age)),rep(0,2), col=2, lty=4, lwd=2)
    
    cord.x <- c(age,rev(age)) 
    cord.y <- c(rep(est.AG.s1-c.AG.s1*se.AG.s1,length(age)),
                rep(est.AG.s1+c.AG.s1*se.AG.s1,length(age)))
    
    
    polygon(cord.x,cord.y,col=adjustcolor('green',alpha.f = 0.3),border=NA,lty=2)
    lines(c(min(age),max(age)), rep(est.AG.s1,2),lty=1,lwd=2,col="darkgreen")
    
    lines(c(min(age),max(age)), rep(est.AG.s1-c.AG.s1*se.AG.s1,2),lty=2,lwd=1,col="darkgreen")
    lines(c(min(age),max(age)), rep(est.AG.s1+c.AG.s1*se.AG.s1,2),lty=2,lwd=1,col="darkgreen")
    
    
    ##AG1LC##########################################################
    cord.x<-c(age[tauL:ifelse(adjust.LC,tauR.1,tauR)],age[ifelse(adjust.LC,tauR.1,tauR):tauL])
    cord.y<-c(est.ext.s1[tauL:ifelse(adjust.LC,tauR.1,tauR)]-
                c.ext.s1*se.ext.s1[tauL:ifelse(adjust.LC,tauR.1,tauR)],
              est.ext.s1[ifelse(adjust.LC,tauR.1,tauR):tauL]+
                c.ext.s1*se.ext.s1[ifelse(adjust.LC,tauR.1,tauR):tauL])
    
    
    polygon(cord.x,cord.y,col=adjustcolor('skyblue',alpha.f = 0.5),
            border = NA)
    lines(age[tauL:ifelse(adjust.LC,tauR.1,tauR)],est.ext.s1[tauL:ifelse(adjust.LC,tauR.1,tauR)],
          lty=1,lwd=2,col="blue")
    
    lines(age[tauL:ifelse(adjust.LC,tauR.1,tauR)],est.ext.s1[tauL:ifelse(adjust.LC,tauR.1,tauR)]+
            c.ext.s1*se.ext.s1[tauL:ifelse(adjust.LC,tauR.1,tauR)],
          lty=2,lwd=1,col="blue")
    lines(age[tauL:ifelse(adjust.LC,tauR.1,tauR)],est.ext.s1[tauL:ifelse(adjust.LC,tauR.1,tauR)]-
            c.ext.s1*se.ext.s1[tauL:ifelse(adjust.LC,tauR.1,tauR)],
          lty=2,lwd=1,col="blue")
    
    legend(location,legend=c(model.v, model.c,"Line zero"),bty="n",
           lty=c(1,1,4),col=c("blue","darkgreen","red"),cex=cex.legend,lwd=rep(2,3))
  } 
  
  # stratum 2
  ## AG1 VS AG1LC
  ## AG1##########################################################
  plot(c(min(age),max(age)),c(-ceiling(max(est.ext.s2[tauL:tauR])),
                              ceiling(max(est.ext.s2[tauL:tauR]))),
       xlab=xlab, type="n", ylab="Coefficient estimate",
       sub=paste0(#ifelse(Bspline,"(Basis Function) ","(LC) "),
         ifelse(comCoef,"" , "(b) "), cov_name[i],ifelse(comCoef,
        " (Common Coefficient)"," (stratum 2)")),ylim=yRange,
       xaxt="n",cex.lab=cex.lab,cex.axis=cex.axis,las=1)
  axis(side=1,at=seq(0,18,2),cex.axis=cex.axis)
  lines(c(min(age),max(age)),rep(0,2), col=2, lty=4, lwd=2)
  
  
  cord.x <- c(age,rev(age)) 
  cord.y <- c(rep(est.AG.s2-c.AG.s2*se.AG.s2,length(age)),
              rep(est.AG.s2+c.AG.s2*se.AG.s2,length(age)))
  
  
  polygon(cord.x,cord.y,col=adjustcolor('green',alpha.f = 0.3),border=NA,lty=2)
  lines(c(min(age),max(age)), rep(est.AG.s2,2),lty=1,lwd=2,col="darkgreen")
  
  lines(c(min(age),max(age)), rep(est.AG.s2-c.AG.s2*se.AG.s2,2),lty=2,lwd=1,col="darkgreen")
  lines(c(min(age),max(age)), rep(est.AG.s2+c.AG.s2*se.AG.s2,2),lty=2,lwd=1,col="darkgreen")
  
  
  ##AG1LC##########################################################
  #if(AG.only==FALSE){
  cord.x<-c(age[ifelse(adjust.LC,tauL.2,tauL):tauR],
            age[tauR:ifelse(adjust.LC,tauL.2,tauL)])
  cord.y<-c(est.ext.s2[ifelse(adjust.LC,tauL.2,tauL):tauR]-
              c.ext.s2*se.ext.s2[ifelse(adjust.LC,tauL.2,tauL):tauR],
            est.ext.s2[tauR:ifelse(adjust.LC,tauL.2,tauL)]+
              c.ext.s2*se.ext.s2[tauR:ifelse(adjust.LC,tauL.2,tauL)])
  
  
  polygon(cord.x,cord.y,col=adjustcolor('skyblue',alpha.f = 0.5),
          border = NA)
  lines(age[ifelse(adjust.LC,tauL.2,tauL):tauR],
        est.ext.s2[ifelse(adjust.LC,tauL.2,tauL):tauR],lty=1,lwd=2,col="blue")
  
  lines(age[ifelse(adjust.LC,tauL.2,tauL):tauR],
        est.ext.s2[ifelse(adjust.LC,tauL.2,tauL):tauR]-
          c.ext.s2*se.ext.s2[ifelse(adjust.LC,tauL.2,tauL):tauR],lty=2,lwd=1,col="blue")
  lines(age[ifelse(adjust.LC,tauL.2,tauL):tauR],
        est.ext.s2[ifelse(adjust.LC,tauL.2,tauL):tauR]+
          c.ext.s2*se.ext.s2[ifelse(adjust.LC,tauL.2,tauL):tauR],lty=2,lwd=1,col="blue")
  
  
  legend(location,legend=c(model.v,model.c, "Line zero"),bty="n",
         lty=c(1,1,4),col=c("blue","darkgreen","red"),cex=cex.legend,lwd=rep(2,3))

}




#' Plot (Cumulative) Baseline Functions from Stratified Cox-Type Models
#'
#' Displays estimated (cumulative) baseline functions along with pointwise 
#' confidence intervals for stratified models with both time-varying and 
#' time-constant regression coefficients.
#'
#' @param fit Object containing output from a stratified Cox-type model fit 
#'   (e.g., from \code{\link{strati.reg}}).
#' @param ONUM Integer; maximum number of time units (i.e., \code{AGEUNIT}* number 
#' of years).
#' @param AGEUNIT Integer; number of time units per year (e.g., 6 for bimonthly 
#' data).
#' @param subtitle Character string specifying plot subtitle. Default is 
#' \code{"stratified"}.
#' @param leg.loc Character string indicating legend placement 
#' (e.g., \code{"topright"}).
#' @param c Numeric; critical value used for constructing pointwise confidence 
#' intervals. Default is \code{1.96} (corresponding to a 95\% confidence interval).
#' @param leg.cex Numeric magnification factor for legend text size. Default is 
#' \code{1.2}.
#' @param cumulative Logical; if \code{TRUE} (default), plots estimated cumulative 
#' baseline functions; if \code{FALSE}, plots estimated baseline functions.
#' @param ylim Numeric vector of length 2 specifying y-axis limits 
#' \code{c(ymin, ymax)}. If \code{NA} (default), limits are computed automatically.
#' @param comLamb Logical; indicates whether common baseline functions 
#'   across strata were assumed in the model fit. Default is \code{FALSE}.
#' @param cex.lab Numeric magnification factor for axis labels relative to current 
#' setting. Default is \code{1.3}.
#' @param cex.axis Numeric magnification factor for axis tick labels relative to 
#' current setting. Default is \code{1.3}.
#' @param model.v Character string; legend label for the model with time-varying 
#' coefficients. Default is \code{"Age-varying coef"}.
#' @param model.c Character string; legend label for the model with time-constant
#'  coefficients. Default is \code{"Age-constant coef"}.
#' @param xlab Character string specifying the x-axis label. Defaults to
#' \code{"Age (in years)"}.
#'
#' @return No return value; generates a plot of stratum-specific baseline functions.
#'
#' @details
#'  The function assumes that \code{method = "all"} and \code{B!=NA} was used in
#'  \code{\link{strati.reg}}.
#'  
#' @references
#' Chen, A. A., Hu, X. J., & Rosychuk, R. J. (2026a). Stratified regression 
#' analysis of zero-truncated recurrent event data. \emph{The Annals of Applied 
#' Statistics}, 20(1), 434-451.
#' 
#' Chen, A. A., Hu, X. J., Rosychuk, R. J., & Zeng, L. (2026b). Analyzing 
#' zero‐truncated recurrent event data by stratified regression with time‐varying 
#' coefficients. \emph{Canadian Journal of Statistics}, e70064.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The 'strati.fit.k.R.p' is the output from strati.reg().
#' strati.fit.k.R.p <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 0, com.lamb = 0, knownBirth = TRUE, 
#'   reptw = 100, resampling = TRUE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#' # Plot stratified cumulative baseline function estimates:
#'   plot_stratified_baseline(fit = strati.fit.k.R.p, ONUM = 108, AGEUNIT = 6, 
#'   subtitle = "stratified analysis", leg.loc="topleft",cumulative = TRUE, 
#'   comLamb = FALSE,c = 1.96)
#'   
#'   # Plot stratified baseline function estimates:
#'   plot_stratified_baseline(fit = strati.fit.k.R.p, ONUM = 108, AGEUNIT = 6, 
#'   subtitle = "stratified analysis", leg.loc="topleft",cumulative = FALSE, 
#'   comLamb = FALSE,c = 1.96)
#'   
#'   # The object 'strati.fit.k.R.p.combB' is the output from strati.reg().
#'   strati.fit.k.R.p.combB <- strati.reg(base = my.base[,c(-1,-7,-8)], visit = my.visit, 
#' AGEUNIT = 6,ONUM = 108, WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), etapara = 1e-4, reptmax = 50,
#'   windowKH = rep(9,108), tauL = 9, tauR = 105, method = "all",
#'   N_star =  my.agge.info, n=50000,
#'   B = 200, seeds = 12345, com.coef = 0, com.lamb = 1, knownBirth = TRUE, 
#'   reptw = 100, resampling = TRUE,r.distn = "Poisson",
#'   mycovars = c("MaleInd","CityAInd","CityBInd"), save.res = FALSE, 
#'   output.lambda = TRUE)
#'   
#'   # Plot the common cumulative baseline function estimates.
#'   plot_stratified_baseline(fit = strati.fit.k.R.p.combB, ONUM = 108, AGEUNIT = 6, 
#'   subtitle = "stratified analysis", leg.loc="topleft",cumulative = TRUE, 
#'   comLamb = TRUE,c = 1.96)
#'   
#'   
#'   # Plot the common baseline function estimates.
#'   plot_stratified_baseline(fit = strati.fit.k.R.p.combB, ONUM = 108, AGEUNIT = 6, 
#'   subtitle = "stratified analysis", leg.loc="topleft",cumulative = FALSE, 
#'   comLamb = TRUE,c = 1.96)
#' }
plot_stratified_baseline=function(fit,ONUM,AGEUNIT,
                                  subtitle="stratified",leg.loc="topright",c=1.96,
                                 leg.cex=1.2,cumulative=TRUE,ylim=NA,
                                 comLamb=FALSE,cex.lab=1.3, cex.axis=1.3,
                                 model.v="Age-varying coef", 
                                 model.c="Age-constant coef",
                                 xlab="Age (in years)"){
  
  agev=1:(ONUM)/AGEUNIT
  
  if(cumulative){
    est.AG.s1=fit$AG1_s1$Lambda
    est.AG.s2=fit$AG1_s2$Lambda
    
    est.LC.s1=fit$AG1LC_s1$Lambda
    est.LC.s2=fit$AG1LC_s2$Lambda
    
    se.AG.s1=sqrt(fit$AG1_s1$Lambdavar)
    se.AG.s2=sqrt(fit$AG1_s2$Lambdavar)
    
    se.LC.s1=sqrt(fit$AG1LC_s1$Lambdavar)
    se.LC.s2=sqrt(fit$AG1LC_s2$Lambdavar)
  }else{
    est.AG.s1=fit$AG1_s1$lambda
    est.AG.s2=fit$AG1_s2$lambda
    
    est.LC.s1=fit$AG1LC_s1$lambda
    est.LC.s2=fit$AG1LC_s2$lambda
    
    se.AG.s1=sqrt(fit$AG1_s1$lambdavar)
    se.AG.s2=sqrt(fit$AG1_s2$lambdavar)
    
    se.LC.s1=sqrt(fit$AG1LC_s1$lambdavar)
    se.LC.s2=sqrt(fit$AG1LC_s2$lambdavar)
  }
  
  if(any(is.na(ylim))){
    yLB=min(est.AG.s1-c*se.AG.s1,est.AG.s2-c*se.AG.s2,
            est.LC.s1-c*se.LC.s1,est.LC.s2-c*se.LC.s2)
    yUB=max(est.AG.s1+c*se.AG.s1,est.AG.s2+c*se.AG.s2,
            est.LC.s1+c*se.LC.s1,est.LC.s2+c*se.LC.s2)
    
  }else{
    yLB=ylim[1]
    yUB=ylim[2]
  }
  
  
  plot(c(min(agev),max(agev)),c(yLB,yUB),
       xlab=xlab, type="n", main=paste0(ifelse(cumulative,"Estimated Cumulative Baseline","Estimated Baseline")," \n(",subtitle,")"),
       ylab=ifelse(cumulative,"Estimated Cumulative Baseline","Estimated Baseline"),
       ylim=c(yLB,yUB),xaxt="n",cex.lab=cex.lab, cex.axis=cex.axis,las=1)
  axis(side=1,at=seq(0,18,2),cex.axis=cex.axis)
  
  lines(agev, est.AG.s1, lty=1, col="black", lwd=2)
  lines(agev, est.AG.s1+c*se.AG.s1, lty=2, col="black", lwd=1)
  lines(agev, est.AG.s1-c*se.AG.s1, lty=2, col="black", lwd=1)
  
  lines(agev, est.LC.s1, lty=1, col="red", lwd=2)
  lines(agev, est.LC.s1+c*se.LC.s1, lty=2, col="red", lwd=1)
  lines(agev, est.LC.s1-c*se.LC.s1, lty=2, col="red", lwd=1)
  
  
  if(!comLamb){
    lines(agev, est.AG.s2, lty=1, col="green", lwd=2)
    lines(agev, est.AG.s2+c*se.AG.s2, lty=2, col="green", lwd=1)
    lines(agev, est.AG.s2-c*se.AG.s2, lty=2, col="green", lwd=1)
    
    lines(agev, est.LC.s2, lty=1, col="blue", lwd=2)
    lines(agev, est.LC.s2+c*se.LC.s2, lty=2, col="blue", lwd=1)
    lines(agev, est.LC.s2-c*se.LC.s2, lty=2, col="blue", lwd=1)
    
    legend(leg.loc,legend=c(paste0(model.c," (s=1)"),paste0(model.c," (s=2)"),
                            paste0(model.v," (s=1)"),paste0(model.v," (s=2)")),
           lty=rep(1,4), col=c("black","green","red","blue"),lwd=rep(2,4),
           bty="n",cex=leg.cex)
  }else{
    legend(leg.loc,legend=c(model.c,model.v),
           lty=rep(1,2), col=c("black","red"),lwd=rep(2,2),
           bty="n",cex=leg.cex)
  }
  
  
  
  
}
