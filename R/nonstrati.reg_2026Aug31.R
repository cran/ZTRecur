# Internal function: Hidden from end users
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
#' Non-Stratified Regression Analysis for Zero-Truncated Recurrent Event Data
#' 
#' Fits a non-stratified Cox-type model with time-constant or time-varying 
#' regression coefficients using zero-truncated recurrent event data integrated 
#' with population census information. The method accounts for coarsened follow-up 
#' information and yields sandwich-form variance estimates along with empirical 
#' bootstrap or multiplier resampling alternatives.
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
#' @param N_star Data frame containing aggregated population census data, 
#' including integer age (\code{AGE}),mycovars, and number of subjects in that 
#' subgroup (\code{pop}).
#' @param n Integer; the approximate size of $O$, the random sample drawn from the 
#' target population.
#' @param knownBirth Logical; indicates whether exact birthdates are known. Defaults
#' to \code{FALSE}. If \code{FALSE}, exact birthdates are unknown, and
#' \code{birth.interval} can be used to calculate the corresponding birthdate
#' intervals in Julian dates.
#' @param reptw Integer; number of birthdate sets generated for the missing 
#' birthdate procedure A (Hu & Rosychuk, 2016). Default is \code{100}.
#' @param B Integer or \code{NA}; number of bootstrap or resampling replicates for 
#' variance estimation. Set to \code{NA} (default) to compute only the robust 
#' sandwich-form variance estimator.
#' @param seeds Numeric; seed value for reproducible sampling. Default is 
#' \code{1234}.
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
#' @return A list containing estimates of the regression coefficients, 
#' standard error estimates (sandwich and/or bootstrap- or resampling-based), 
#' estimated cumulative baseline functions, and estimated baseline functions when
#' \code{output.lambda = TRUE}.
#'   
#' @details
#'  The estimation procedure is similar to that described in Xiong, Hu & Rosychuk
#'  (2026). The approach accounts for missing birthdates using Procedure A of
#'  Hu & Rosychuk (2016) and obtains variance estimates using the multiplier
#'  resampling procedure described in Chen, Hu & Rosychuk (2026).
#' 
#' @references
#' Chen, A. A., Hu, X. J., & Rosychuk, R. J. (2026). Stratified regression 
#' analysis of zero-truncated recurrent event data. \emph{The Annals of Applied 
#' Statistics}, 20(1), 434-451.
#' 
#' Hu, X. J., & Rosychuk, R. J. (2016). Imputation procedures for missing 
#' birthdates in administrative health data. \emph{Statistics in Medicine}, 
#' 35(14), 2341-2355.
#' 
#' Xiong, Y., Hu, X. J., & Rosychuk, R. J. (2026). Learning differences between two 
#' decades of mental health related emergency department visits by youth via 
#' recurrent events data analyses. \emph{Journal of the Royal Statistical Society Series 
#' C: Applied Statistics}, qlag001.
#'
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' \donttest{
#' # Example non-stratified regression model with unknown birthdates.
#' # Sandwich-form variance estimates.
#' nonstrati.fit.uk <- nonstrati.reg(base = my.base[,c(-1,-9)], visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = FALSE, reptw = 100, B = NA, seeds = 12345, resampling = TRUE, 
#'   r.distn = "std Normal",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#'   
#'   
#' # Example non-stratified regression model with known birthdates.
#' # Sandwich-form variance estimates.
#' nonstrati.fit.k <- nonstrati.reg(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = TRUE, reptw = 100, B = NA, seeds = 12345, resampling = TRUE, 
#'   r.distn = "std Normal",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#'   
#' # Example with known birthdates and bootstrap-based variance estimates.
#' nonstrati.fit.k.B <- nonstrati.reg(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = TRUE, reptw = 100, B = 200, seeds = 12345, resampling = FALSE, 
#'   r.distn = "std Normal",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#'   
#' # Example with known birthdates and multiplier resampling-based variance
#' # estimates using standard normal multipliers.
#' nonstrati.fit.k.R <- nonstrati.reg(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = TRUE, reptw = 100, B = 200, seeds = 12345, resampling = TRUE, 
#'   r.distn = "std Normal",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#'   
#' # Example with known birthdates and multiplier resampling-based variance
#' # estimates using Poisson multipliers with variance 1.
#' nonstrati.fit.k.R.p <- nonstrati.reg(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = TRUE, reptw = 100, B = 200, seeds = 12345, resampling = TRUE, 
#'   r.distn = "Poisson",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#' }
#' 
nonstrati.reg=function(base, visit,
                    AGEUNIT,ONUM,WINDOWL,WINDOWR,etapara,reptmax,
                    windowKH,tauL,tauR,method,N_star,
                    n,knownBirth=FALSE,reptw=100,B=NA,seeds=1234,resampling=FALSE,
                    r.distn="std Normal",
                    mycovars,save.res=FALSE,resultpath=NULL,output.lambda=FALSE){
  # that will cause the code to fail
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
  #method="all" calls each of the 2 methods; ="L" corresponds to the estimation under an AG model;
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
  TMPvisitage=a[[3]]
  
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
    
    b=getAG1coh(betaini,covat0,covat0z,tmpY,N_star.long,assumedtiny,
                reptmax,etapara,ONUM,1,n,dN)
    
    if(mean(abs(b[[8]]/n))>etapara){
      stop("The estimating function is not close to 0!")
    }
    
    if(b[[7]]>=reptmax){
      stop("It doesn't converge!")
    }
    # resampling method is not really needed. We have a close form for variance.
    if(!is.na(B)){
      t1=Sys.time()
      set.seed(seeds)
      beta_mat=matrix(NA,nrow=B,ncol=length(b[[1]]))
      Lambda_mat=matrix(NA,nrow=B,ncol=ONUM)
      
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
          #any(duplicated(lst.sub))
        }
        
        tryCatch({
          if(resampling){
            boot=getAG1coh(betaini,covat0,covat0z,tmpY,N_star.long,assumedtiny,
                        reptmax,etapara,ONUM,0,n,dN*mat_w)
            
          }else{
            boot=getAG1coh(betaini,covat0[lst.sub,],covat0z,tmpY[lst.sub,],N_star.long,assumedtiny,
                           reptmax,etapara,ONUM,0,n,dN[lst.sub,])
            
          }
          
          if(mean(abs(boot[[8]]/n))>etapara){
            stop("The estimating function is not close to 0!")
          }
          
          if(boot[[7]]>=reptmax){
            stop("It doesn't converge!")
          }
          
          beta_mat[bt,]=boot[[1]]
          Lambda_mat[bt,]=boot[[3]]
          
          res.AG=list(b=b,beta_mat=beta_mat, Lambda_mat=Lambda_mat)
          
          if(save.res){
            if (is.null(resultpath) || length(resultpath) != 1L ||
                is.na(resultpath) || !nzchar(resultpath)) {
              stop("Please provide an explicit output directory in 'resultpath'.")
            }
            
            if (!dir.exists(resultpath)) {
              stop("The output directory does not exist: ", resultpath)
            }
            
            save(res.AG, file=paste(resultpath,"Non-stratified_Result_constant_coef",
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
        AG1=list(beta=b[[1]],betavar=b[[2]],Lambda=b[[3]],Lambdavar=apply(Lambda_mat,2,var),
                 lambda=gen.lambda(Lambda=b[[3]],ONUM,AGEUNIT,mat=FALSE,lambda.var=FALSE),
                 lambdavar=gen.lambda(Lambda=Lambda_mat,ONUM,AGEUNIT,mat=TRUE,lambda.var=TRUE)$lambda.var,
                 Pscore=b[[8]],beta_mat=beta_mat,Lambda_mat=Lambda_mat)
      }else{
        AG1=list(beta=b[[1]],betavar=b[[2]],Lambda=b[[3]],Lambdavar=apply(Lambda_mat,2,var),
                 Pscore=b[[8]],beta_mat=beta_mat,Lambda_mat=Lambda_mat)
      }
    }else{
      AG1=list(beta=b[[1]],betavar=b[[2]],Lambda=b[[3]],Lambdavar=NA,
               Pscore=b[[8]])
    }
  }  

  if(method=="all"|| method=="LC")
  {
    if(!exists("getLCcoh"))
      stop("Compilation failed. The function getLCcoh does not exist")
    
    if(tauR>ONUM)
      stop("tauR should be less than ONUM: the largest age units")
    
    message("Local Constant Estimation: AG1LC\n")
    betaMat.0=matrix(rep(betaini,ONUM),ncol = ncol(covat0z),byrow = TRUE)
    
    blc=getLCcoh(betaMat.0,covat0,covat0z,tmpY,N_star.long,assumedtiny,
                    reptmax,etapara,ONUM,windowKH,tauL,tauR,1,n,dN)
  
    if(!is.na(B)){
      t1=Sys.time()
      set.seed(seeds)
      
      beta_mat=array(NA,dim = c(dim(blc[[1]]),B))
      Lambda_mat=matrix(NA,nrow=B,ncol=ONUM)
      
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
          #any(duplicated(lst.sub))
        }
        
        tryCatch({
          if(resampling){
            bootlc=getLCcoh(betaMat.0,covat0,covat0z,tmpY,N_star.long,assumedtiny,
                           reptmax,etapara,ONUM,windowKH,tauL,tauR,0,n,dN*mat_w)
            
          }else{
            bootlc=getLCcoh(betaMat.0,covat0[lst.sub,],covat0z,tmpY[lst.sub,],N_star.long,assumedtiny,
                           reptmax,etapara,ONUM,windowKH,tauL,tauR,0,n,dN[lst.sub,])
            
          }
          
          beta_mat[,,bt]=bootlc[[1]]
          Lambda_mat[bt,]=bootlc[[3]]
          
          res.LC=list(blc=blc,beta_mat=beta_mat, Lambda_mat=Lambda_mat)
          
          if(save.res){
            if (is.null(resultpath) || length(resultpath) != 1L ||
                is.na(resultpath) || !nzchar(resultpath)) {
              stop("Please provide an explicit output directory in 'resultpath'.")
            }
            
            if (!dir.exists(resultpath)) {
              stop("The output directory does not exist: ", resultpath)
            }
            
            save(res.LC, file=paste(resultpath,"Non-stratified_Result_age_varying_coef",
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
        AG1LC=list(beta=blc[[1]],betavar=blc[[2]],Lambda=blc[[3]],
                   Lambdavar=apply(Lambda_mat,2,var),
                 lambda=gen.lambda(Lambda=blc[[3]],ONUM,AGEUNIT,mat=FALSE,lambda.var=FALSE),
                 lambdavar=gen.lambda(Lambda=Lambda_mat,ONUM,AGEUNIT,mat=TRUE,lambda.var=TRUE)$lambda.var,
                 beta_mat=beta_mat,Lambda_mat=Lambda_mat)
      }else{
        AG1LC=list(beta=blc[[1]],betavar=blc[[2]],Lambda=blc[[3]],
                   Lambdavar=apply(Lambda_mat,2,var),
                 beta_mat=beta_mat,Lambda_mat=Lambda_mat)
      }
    }else{
      AG1LC=list(beta=blc[[1]],betavar=blc[[2]],Lambda=blc[[3]])
    }
  }
  
  
  if(!exists("AG1"))
  {AG1=NA}
  
  if(!exists("AG1LC"))
  {AG1LC=NA}

  RESULTA=list(AG1=AG1,AG1LC=AG1LC)
  return(RESULTA)
}


#' Plot Estimated Regression Coefficients from Non-Stratified Cox-Type Models
#'
#' Displays estimates and pointwise confidence intervals for time-varying 
#' and time-constant regression coefficients from fitted non-stratified models.
#'
#' @param fit.non.s Object containing the output from a non-stratified model fit 
#'   (e.g., from \code{\link{nonstrati.reg}}).
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
#' @param resampling.se Logical; if \code{TRUE} (default), uses standard errors 
#' derived from bootstrap/resampling. If \code{FALSE}, uses sandwich standard errors.
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
#' @return No return value; generates a coefficient plot for the specified covariate.
#'
#' @details
#'  The function assumes that \code{method = "all"} was used in
#'  \code{\link{nonstrati.reg}}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The 'nonstrati.fit.uk' and 'nonstrati.fit.k' are outputs from nonstrati.reg()
#' nonstrati.fit.uk <- nonstrati.reg(base = my.base[,c(-1,-9)], visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = FALSE, reptw = 100, B = NA, seeds = 12345, resampling = TRUE, 
#'   r.distn = "std Normal",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#'   
#' # 95% pointwise confidence intervals using sandwich-form standard errors.
#' plot_nonstrati_coef(fit.non.s = nonstrati.fit.uk, cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 1, age = c(1:108)/6, tauL = 9, tauR = 105,resampling.se = FALSE, c = 1.96)
#' 
#' plot_nonstrati_coef(fit.non.s = nonstrati.fit.uk, cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 2, age = c(1:108)/6, tauL = 9, tauR = 105,resampling.se = FALSE, c = 1.96)
#' 
#' plot_nonstrati_coef(fit.non.s = nonstrati.fit.uk, cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 3, age = c(1:108)/6, tauL = 9, tauR = 105,resampling.se = FALSE, c = 1.96)
#' 
#' # 95% pointwise confidence intervals using bootstrap- or
#' # multiplier-resampling-based standard errors.
#' nonstrati.fit.k <- nonstrati.reg(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = TRUE, reptw = 100, B = 100, seeds = 12345, resampling = TRUE, 
#'   r.distn = "std Normal",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#'   
#' plot_nonstrati_coef(fit.non.s = nonstrati.fit.k, cov_name = c("MaleInd","CityAInd","CityBInd"), 
#' i = 1, age = c(1:108)/6, tauL = 9, tauR = 105,resampling.se = TRUE, c = 1.96)
#' }
plot_nonstrati_coef=function(fit.non.s,cov_name,i=1,age,tauL,tauR,yRange=NA,
                             resampling.se=TRUE,model.v="Age-varying coef", 
                             model.c="Age-constant coef",
                             location="topright",cex.axis=1.2,
                             cex.lab=1.3,cex.legend=1.2,c=1.96,
                             xlab="Age (in years)"){
  
  est.AG=fit.non.s$AG1$beta[i]
  est.ext=fit.non.s$AG1LC$beta[,i]
  
  if(!resampling.se){
    # sandwich-form SE
    se.AG=sqrt(diag(fit.non.s$AG1$betavar))[i]
    se.ext=sqrt(fit.non.s$AG1LC$betavar[,i])
  }else{
    # Bootstrapping-based/resampling-based SE
    se.AG=sd(fit.non.s$AG1$beta_mat[,i])
    se.ext=apply(fit.non.s$AG1LC$beta_mat[,i,],1,sd)
  }
  c.AG=c.ext=c
  
  if(any(is.na(yRange))){
    maxy1=max(abs(est.ext[tauL:tauR]+c.ext*se.ext[tauL:tauR]),
              abs(est.ext[tauL:tauR]-c.ext*se.ext[tauL:tauR]))
    maxy=1.1*maxy1
    yRange=c(-maxy,maxy)
  }  
  
  
  ###############non-stratified################## AG1#################
  plot(c(min(age),max(age)),c(-ceiling(max(est.ext[tauL:tauR])),ceiling(max(est.ext[tauL:tauR]))),
       xlab=xlab, type="n", ylab="Coefficient estimate",
       sub=paste0(cov_name[i]," (non-stratified)"),ylim=yRange,xaxt="n",cex.lab=cex.lab,
       cex.axis=cex.axis,las=1)
  axis(side=1,at=seq(0,18,2),cex.axis=cex.axis)
  lines(c(min(age),max(age)),rep(0,2), col=2, lty=4, lwd=2)
  
  cord.x <- c(age,rev(age))
  cord.y <- c(rep(est.AG-c.AG*se.AG,length(age)),
              rep(est.AG+c.AG*se.AG,length(age)))
  
  
  polygon(cord.x,cord.y,col=adjustcolor('green',alpha.f = 0.3),border=NA,lty=2)
  
  lines(c(min(age),max(age)), rep(est.AG,2),lty=1,lwd=2,col="darkgreen")
  lines(c(min(age),max(age)), rep(est.AG-c.AG*se.AG,2),lty=2, lwd=1, col="darkgreen")
  lines(c(min(age),max(age)), rep(est.AG+c.AG*se.AG,2),lty=2, lwd=1, col="darkgreen")
  
  
  ##AG1LC##########################################################
  # if(AG.only==FALSE){
  cord.x<-c(age[tauL:tauR],age[tauR:tauL])
  cord.y<-c(est.ext[tauL:tauR]-c.ext*se.ext[tauL:tauR],
            est.ext[tauR:tauL]+c.ext*se.ext[tauR:tauL])
  
  
  polygon(cord.x,cord.y,col=adjustcolor('skyblue',alpha.f = 0.5),
          border = NA)
  lines(age[tauL:tauR],est.ext[tauL:tauR],lty=1,lwd=2,col="blue")
  
  lines(age[tauL:tauR],est.ext[tauL:tauR]-c.ext*se.ext[tauL:tauR],lty=2,lwd=1,col="blue")
  lines(age[tauL:tauR],est.ext[tauL:tauR]+c.ext*se.ext[tauL:tauR],lty=2,lwd=1,col="blue")
  
  
  leg=c(model.v, model.c, "Line zero")
  
  legend(location,legend=leg,bty="n",
         lty=c(1,1,4),col=c("blue","darkgreen","red"),cex=cex.legend,lwd=rep(2,3))
  
}



#' Plot (Cumulative) Baseline Functions from Non-Stratified Cox-Type Models
#'
#' Displays estimated (cumulative) baseline functions along with pointwise 
#' confidence intervals for non-stratified models with both time-varying and 
#' time-constant regression coefficients.
#'
#' @param fit.non.s Object containing the output from a non-stratified model fit 
#'   (e.g., from \code{\link{nonstrati.reg}}).
#' @param ONUM Integer; maximum number of time units (i.e., \code{AGEUNIT}* number 
#' of years).
#' @param AGEUNIT Integer; number of time units per year (e.g., 6 for bimonthly 
#' data).
#' @param subtitle Character string specifying plot subtitle. Default is 
#' \code{"non-stratified"}.
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
#' @return No return value; generates a (cumulative) baseline function plot.
#' 
#' @details
#'  The function assumes that \code{method = "all"} and \code{B!=NA} was used in
#'  \code{\link{nonstrati.reg}}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The 'nonstrati.fit.k.B' is the output from nonstrati.reg()
#' 
#' # Note: 95% pointwise confidence intervals for the baseline functions can
#' # be plotted only when bootstrap- or multiplier-resampling-based variance
#' # estimates are available from nonstrati.reg().
#' nonstrati.fit.k.B <- nonstrati.reg(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, 
#' AGEUNIT = 6, ONUM = 108, 
#'   WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), 
#'   etapara = 1e-4, reptmax = 50, windowKH = rep(9,108), 
#'   tauL = 9, tauR = 105, method = "all", N_star = my.agge.info, n=50000, 
#'   knownBirth = TRUE, reptw = 100, B = 200, seeds = 12345, resampling = FALSE, 
#'   r.distn = "std Normal",mycovars = c("MaleInd","CityAInd","CityBInd"), 
#'   save.res = FALSE, output.lambda = TRUE)
#'
#' # Plot the estimated cumulative baseline functions:
#' plot_nonstrati_baseline(fit.non.s = nonstrati.fit.k.B, ONUM = 108, AGEUNIT = 6, 
#'   subtitle = "non-stratified analysis", cumulative = TRUE, c = 1.96,
#'   leg.loc="topleft")
#'   
#' # Plot the estimated baseline functions:
#' plot_nonstrati_baseline(fit.non.s = nonstrati.fit.k.B, ONUM = 108, AGEUNIT = 6, 
#'   subtitle = "non-stratified analysis", cumulative = FALSE, c = 1.96,
#'   leg.loc="topleft",ylim=c(0.05,0.17))
#' }
plot_nonstrati_baseline=function(fit.non.s,ONUM,AGEUNIT,
                                 subtitle="non-stratified",leg.loc="topright",c=1.96,
                                 leg.cex=1.2,cumulative=TRUE,ylim=NA,
                                 cex.lab=1.3, cex.axis=1.3,model.v="Age-varying coef", 
                                 model.c="Age-constant coef",xlab="Age (in years)"){
  
  agev=1:(ONUM)/AGEUNIT
  
  if(cumulative){
    est.AG=fit.non.s$AG1$Lambda
    est.LC=fit.non.s$AG1LC$Lambda
    
    se.AG=sqrt(fit.non.s$AG1$Lambdavar)
    se.LC=sqrt(fit.non.s$AG1LC$Lambdavar)
  }else{
    est.AG=fit.non.s$AG1$lambda
    est.LC=fit.non.s$AG1LC$lambda
    
    se.AG=sqrt(fit.non.s$AG1$lambdavar)
    se.LC=sqrt(fit.non.s$AG1LC$lambdavar)
  }
  
  
  if(any(is.na(ylim))){
    yLB=min(est.AG-c*se.AG,est.LC-c*se.LC)
    yUB=max(est.AG+c*se.AG,est.LC+c*se.LC)
  }else{
    yLB=ylim[1]
    yUB=ylim[2]
  }
  
  plot(c(min(agev),max(agev)),c(yLB,yUB),
       xlab=xlab, type="n", main=paste0(ifelse(cumulative,
                "Estimated Cumulative Baseline","Estimated Baseline")," \n(",subtitle,")"),
       ylab=ifelse(cumulative,"Estimated Cumulative Baseline","Estimated Baseline"),
       ylim=c(yLB,yUB),xaxt="n",cex.lab=cex.lab, cex.axis=cex.axis,las=1)
  axis(side=1,at=seq(0,18,2),cex.axis=cex.axis)
  
  lines(agev, est.AG, lty=1, col="black", lwd=2)
  lines(agev, est.AG+c*se.AG, lty=2, col="black", lwd=1)
  lines(agev, est.AG-c*se.AG, lty=2, col="black", lwd=1)
  
  lines(agev, est.LC, lty=1, col="red", lwd=2)
  lines(agev, est.LC+c*se.LC, lty=2, col="red", lwd=1)
  lines(agev, est.LC-c*se.LC, lty=2, col="red", lwd=1)
  
  legend(leg.loc,legend=c(model.c,model.v),
         lty=c(1,1), col=c("black","red"),lwd=rep(2,2),
         bty="n",cex=leg.cex)
}
