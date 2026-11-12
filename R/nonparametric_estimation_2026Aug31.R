# Internal function: Hidden from end users
dmu.z=function(mycovars,value.ls,dNStarnewf,covat0,tmpY.z){
  ls=which(apply(covat0,1,paste0,collapse="")==value.ls)
  ls.2=which(apply(tmpY.z[,mycovars],1,paste0,collapse="")==value.ls)
    
  dNY=apply(dNStarnewf[ls,],2,sum)
  den=tmpY.z[ls.2,"pop"]
  dmu=dNY/den
  res=list(dNY=dNY,den=den,dmu=dmu)
  return(res)
}

utils::globalVariables(c("%>%", "patid", "birthL", "birthR", "AGE", "pop"))

#' Generate Birthdate Intervals Based on Integer Ages at Event Times
#'
#' @param visit Matrix of event-level information. It contains subject ID 
#' (\code{patid}), Julian date of the event (\code{date}), and integer age at 
#' that event date (\code{age}).
#' @param MAXAGE Maximum age considered. The default value is 18.
#' @param WINDOWL The start date of data extraction window, such as 
#' as.Date("2010-04-01").
#' @param WINDOWR The end date of data extraction window, such as 
#' as.Date("2025-03-31").
#' @returns A subject-level data frame with patient ID (\code{patid}), 
#' left (\code{birthL}) and right (\code{birthR}) endpoints of birthdate 
#' intervals in Julian date.
#' 
#' @details
#' The calculation of birthdate interval is described in Hu & Rosychuk (2016).
#' 
#'@references
#' Hu, X. J., & Rosychuk, R. J. (2016). Marginal regression analysis of recurrent 
#' events with coarsened censoring times. \emph{Biometrics}, 72(4), 1113-1122.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' # Note: The \code{visit} must be a data frame.
#' birth.b <- birth.interval(visit=as.data.frame(my.visit),
#' WINDOWL=as.Date("2010-04-01"), WINDOWR=as.Date("2025-03-31"))
#' 
birth.interval=function(visit,MAXAGE=18,WINDOWL,WINDOWR){
  visit$birthL=visit$date-(visit$age+1)*365.25
  visit$birthL2=as.numeric(WINDOWL)-(MAXAGE*365.25)
  visit$birthL=pmax(visit$birthL,visit$birthL2)
  visit$birthL2=NULL
  
  visit$birthR=visit$date-visit$age*365.25
  visit$birthR2=as.numeric(WINDOWR)
  visit$birthR=pmin(visit$birthR,visit$birthR2)
  visit$birthR2=NULL
  # range(visit$birthL)
  # range(visit$birthR)
  #all(visit$birthL<visit$birthR)
  #summary(visit$birthR-visit$birthL)
  birth.b=visit%>%group_by(patid)%>%summarise(birthL=max(birthL),
                                              birthR=min(birthR))
  return(birth.b)
}

##################################
#' Nonparametric Estimation of Marginal Mean and Rate Functions
#'
#' Performs nonparametric estimation of marginal mean and rate functions across 
#' covariate groups using zero-truncated recurrent event data integrated with 
#' supplementary population information. The procedure accounts for coarsened 
#' individual follow-up information and provides variance estimates via bootstrap 
#' or multiplier resampling techniques.
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
#' @param mycovars Character vector specifying the variable names of covariates 
#' to define subgroups.
#' @param AGEUNIT Integer; the number of time units per year (e.g., 6 for 
#' bimonthly data).
#' @param ONUM Integer; the maximum number of time units under observation 
#' (e.g., 18 * AGEUNIT for <18 years).
#' @param WINDOWL Numeric; Julian date representing the left boundary of the 
#' data extraction window.
#' @param WINDOWR Numeric; Julian date representing the right boundary of the 
#' data extraction window.
#' @param agge.info Data frame containing aggregated population census data, 
#' including integer age (\code{AGE}),mycovars, and number of subjects in that 
#' subgroup (\code{pop}).
#' @param knownBirth Logical; indicates whether exact birthdates are known. Defaults
#' to \code{FALSE}. If \code{FALSE}, exact birthdates are unknown, and
#' \code{birth.interval} can be used to calculate the corresponding birthdate
#' intervals in Julian dates.
#' @param reptw Integer; number of birthdate sets generated for the missing 
#' birthdate procedure A (Hu & Rosychuk, 2016). Default is \code{100}.
#' @param B Integer; number of bootstrap samples or resampling replicates for 
#' variance estimation. Default is \code{200}.
#' @param seeds Numeric; seed value for reproducible sampling. Default is 
#' \code{1234}.
#' @param resampling Logical; if \code{TRUE}, applies the multiplier resampling 
#' method described in Chen, Hu, & Rosychuk (2026). If \code{FALSE} (default), 
#' standard bootstrapping is used.
#' @param r.distn Character string specifying the multiplier distribution used for
#' resampling when \code{resampling = TRUE}. Options are the standard normal
#' distribution (\code{"std Normal"}; default) or a Poisson distribution with
#' variance 1 (\code{"Poisson"}).
#' @param save.res Logical; if \code{TRUE}, writes output files to disk. 
#' Default is \code{FALSE}.
#' @param resultpath Character string indicating the file path destination 
#' for saved results, such as "C:\\project\\". Required when \code{save.res = TRUE}.
#'
#' @return A list containing the estimated marginal mean functions and rate functions, 
#'   their corresponding variance estimates, and numerical results for each bootstrap 
#'   or resampling sample.
#'
#' @details
#'  The detailed estimation procedure is described in Chen, Hu & Rosychuk (2026a).
#'  The approach accounts for missing birthdates using Procedure A of Hu &
#'  Rosychuk (2016) and calculates variance estimates using the multiplier
#'  resampling procedure described in Chen, Hu & Rosychuk (2026b).
#' 
#'@references
#' Chen, A. A., Hu, X. J., & Rosychuk, R. J. (2026a). Statistical Learning of 
#' Pediatric Mental Health-Related Emergency Department Visits Across COVID-19 
#' Pandemic Periods. \emph{arXiv preprint arXiv:2607.26210}.
#'
#' Chen, A. A., Hu, X. J., & Rosychuk, R. J. (2026b). Stratified regression 
#' analysis of zero-truncated recurrent event data. \emph{The Annals of Applied 
#' Statistics}, 20(1), 434-451.
#'
#' Hu, X. J., & Rosychuk, R. J. (2016). Imputation procedures for missing 
#' birthdates in administrative health data. \emph{Statistics in Medicine}, 
#' 35(14), 2341-2355.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # When birthdates are unknown, bootstrap-based variance estimates are used.
#' nonpara.res.uk <- nonpara_esti(base = my.base[,c(-1,-9)], 
#' visit = my.visit, mycovars = c("MaleInd","CityAInd","CityBInd"),
#'   AGEUNIT = 6, ONUM = (18*6), WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), agge.info = my.agge.info,
#'   seeds = 1234)
#'   
#' # When birthdates are known, bootstrap-based variance estimates are used.
#' nonpara.res.k <- nonpara_esti(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, mycovars = c("MaleInd","CityAInd","CityBInd"),
#'   AGEUNIT = 6, ONUM = (18*6), WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), agge.info = my.agge.info,
#'   knownBirth=TRUE,seeds = 1234) 
#'   
#' # When birthdates are known, multiplier resampling-based variance estimates
#' # are obtained using standard normal multipliers.
#' nonpara.res.k.r.std <- nonpara_esti(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, mycovars = c("MaleInd","CityAInd","CityBInd"),
#'   AGEUNIT = 6, ONUM = (18*6), WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), agge.info = my.agge.info,
#'   knownBirth=TRUE,seeds = 1234,resampling=TRUE) 
#'   
#' # When birthdates are known, multiplier resampling-based variance estimates
#' # are obtained using Poisson multipliers.
#' nonpara.res.k.r.poi <- nonpara_esti(base = my.base[,c(-1,-7,-8)], 
#' visit = my.visit, mycovars = c("MaleInd","CityAInd","CityBInd"),
#'   AGEUNIT = 6, ONUM = (18*6), WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), agge.info = my.agge.info,
#'   knownBirth=TRUE,seeds = 1234,resampling=TRUE,r.distn="Poisson")
#' }
#' 
nonpara_esti=function(base, visit, mycovars,AGEUNIT,ONUM,WINDOWL,WINDOWR,
                      agge.info,knownBirth=FALSE,reptw=100,B=200,seeds=1234,
                      resampling=FALSE,r.distn="std Normal",save.res=FALSE,resultpath=NULL)
{
  #do a few checks; should do many more here so that people can't pass things 
  # that will cause the code to fail
  if(is.matrix(base)==FALSE || is.matrix(visit)==FALSE){
    stop("Base or visit are not passed as matrices. Stopping.\n");
  }
  
  TOTALM=dim(base)[1]
  TOTALM2=dim(visit)[1]
  
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
  t1=Sys.time()
  if(knownBirth){
    a=getmoredata(birth,WINDOWL,WINDOWR,dates,p1stvt,
                  TOTALM,TOTALM2,ONUM,nmvst,AGEUNIT)
    TMPcount=a[[1]]
    tmpY=a[[2]]
    dN=TMPcount*tmpY
  }else{
    a=getmoredataComb2(birthl,birthr,WINDOWL,WINDOWR,dates,p1stvt,
                       reptw,TOTALM,TOTALM2,ONUM,nmvst,AGEUNIT)
    dN=a[[5]]
  }
  t2=Sys.time()
  t=t2-t1
  #str(t)
  message(paste0("Prepare data for modelling, which takes ",round(t,3)," ",attr(t, "units"),".\n"))
  #rm(t)
  
  tmpY.z=apply(agge.info,2,rep,each=6)
  
  value.ls=sort(unique(apply(covat0,1,paste0,collapse="")))
  
  res.zMat=lapply(value.ls,dmu.z,mycovars=mycovars,
                 dNStarnewf=dN,covat0=covat0,tmpY.z=tmpY.z)
  
  # dNY.zMat=do.call(cbind,lapply(res.zMat,'[[',"dNY"))
  # colnames(dNY.zMat)=value.ls
  # den.zMat=do.call(cbind,lapply(res.zMat,'[[',"den"))
  # colnames(den.zMat)=value.ls
  dmu.zMat=do.call(cbind,lapply(res.zMat,'[[',"dmu"))
  colnames(dmu.zMat)=value.ls
  
  mu.zMat=apply(dmu.zMat,2,cumsum)
  rate.zMat=dmu.zMat*AGEUNIT
  
  if(!is.na(B)){
    t1=Sys.time()
    set.seed(seeds)
    mu.BList=replicate(B,matrix(NA,nrow=ONUM,
                      ncol=ncol(dmu.zMat)),simplify = FALSE)
    rate.BList=replicate(B,matrix(NA,nrow=ONUM,
                      ncol=ncol(dmu.zMat)),simplify = FALSE)
    
    bt=1
    while(bt<=B){
      message(paste0("Boot B=", bt,"\n"))
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
      
      if(resampling){
        res.zMat.b=lapply(value.ls,dmu.z,mycovars=mycovars,
                          dNStarnewf=dN*mat_w,
                          covat0=covat0,tmpY.z=tmpY.z)
      }else{
        res.zMat.b=lapply(value.ls,dmu.z,mycovars=mycovars,
                          dNStarnewf=dN[lst.sub,],
                          covat0=covat0[lst.sub,],tmpY.z=tmpY.z)
      }
      
      
      dmu.zMat.b=do.call(cbind,lapply(res.zMat.b,'[[',"dmu"))
      colnames(dmu.zMat.b)=value.ls
      mu.BList[[bt]]=apply(dmu.zMat.b,2,cumsum)
      rate.BList[[bt]]=dmu.zMat.b*AGEUNIT
      
      res.nonpara=list(mu.zMat=mu.zMat,rate.zMat=rate.zMat,
                  mu.BList=mu.BList,rate.BList=rate.BList)
      
      if(save.res){
        if (is.null(resultpath) || length(resultpath) != 1L ||
            is.na(resultpath) || !nzchar(resultpath)) {
          stop("Please provide an explicit output directory in 'resultpath'.")
        }
        
        if (!dir.exists(resultpath)) {
          stop("The output directory does not exist: ", resultpath)
        }
        
        save(res.nonpara, file=paste(resultpath,"Nonparametric_Result",
                                ifelse(resampling,"_R","_B"),B,".RData",sep=""))
      }
      
      bt=bt+1
    }
    t2=Sys.time()
    t=t2-t1
    message(paste0("Resampling method for estimating variances is done, which takes ",
               round(t,3)," ",attr(t, "units"),".\n"))
    rm(t)
    
    arr.mu=simplify2array(mu.BList)
    var.mu.mat=apply(arr.mu,c(1,2),var)
    
    arr.rate=simplify2array(rate.BList)
    var.rate.mat=apply(arr.rate,c(1,2),var)
  }    
  
  if(!exists("mu.BList"))
  {mu.BList=NA}
  
  if(!exists("rate.BList"))
  {rate.BList=NA}
  
  if(!exists("var.mu.mat"))
  {var.mu.mat=NA}
  
  if(!exists("var.rate.mat"))
  {var.rate.mat=NA}
  
  RESULTA=list(mean=mu.zMat,mean.var=var.mu.mat,rate=rate.zMat,rate.var=var.rate.mat,
               mean.Boot=mu.BList,rate.Boot=rate.BList)
  
  return(RESULTA)
}


#' Plot Nonparametric Marginal Mean or Rate Estimates
#'
#' Displays estimated marginal mean or rate functions along with pointwise 
#' confidence intervals for different covariate groups based on output from 
#' \code{\link{nonpara_esti}}.
#'
#' @param res.nonpar An object returned by \code{\link{nonpara_esti}} containing 
#' nonparametric estimates.
#' @param rate.p Logical; if \code{TRUE} (default), plots marginal rate estimates. 
#' If \code{FALSE}, plots marginal mean estimates.
#' @param c Numeric; critical value used for constructing pointwise confidence 
#' intervals. Default is \code{1.96} (corresponding to a 95\% confidence interval).
#' @param max.L.m Numeric; upper limit for the y-axis. If \code{NA} (default), 
#' the maximum value is automatically calculated based on the data.
#' @param ONUM Integer; maximum number of time units (i.e., \code{AGEUNIT}* number 
#' of years).
#' @param AGEUNIT Integer; number of time units per year (e.g., 6 for bimonthly 
#' data).
#' @param cex Numeric; magnification factor for legend text size. Default is 
#' \code{1}.
#' @param two.plot Logical; if \code{TRUE}, splits display into two separate 
#' panels, useful when visualizing a large number of covariate combinations. 
#' Default is \code{FALSE}.
#' @param xlab Character string specifying the x-axis label. Defaults to
#' \code{"Age (in years)"}.
#'
#' @return No return value; generates a plot of marginal mean or rate functions.
#'
#'@references
#' Chen, A. A., Hu, X. J., & Rosychuk, R. J. (2026). Statistical Learning of 
#' Pediatric Mental Health-Related Emergency Department Visits Across COVID-19 
#' Pandemic Periods. \emph{arXiv preprint arXiv:2607.26210}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # The object 'nonpara.res.uk' is the output from nonpara_esti().
#' nonpara.res.uk <- nonpara_esti(base = my.base[,c(-1,-9)], 
#' visit = my.visit, mycovars = c("MaleInd","CityAInd","CityBInd"),
#'   AGEUNIT = 6, ONUM = (18*6), WINDOWL = as.numeric(as.Date("2010-04-01")), 
#'   WINDOWR = as.numeric(as.Date("2025-03-31")), agge.info = my.agge.info,
#'   seeds = 1234)
#' 
#' # Plot estimated marginal means across covariate subgroups.
#' plot_nonpar(res.nonpar = nonpara.res.uk, rate.p = FALSE, c = 1.96, ONUM = 108, 
#'   AGEUNIT = 6, two.plot = FALSE)
#'   
#' # Plot estimated marginal rates across covariate subgroups.
#' plot_nonpar(res.nonpar = nonpara.res.uk, rate.p = TRUE, c = 1.96, ONUM = 108, 
#'   AGEUNIT = 6, two.plot = FALSE,max.L.m=0.20)
#' }
#' 
plot_nonpar=function(res.nonpar,rate.p=TRUE,c=1.96,max.L.m=NA,ONUM,AGEUNIT,cex=1,
                     two.plot=FALSE,xlab="Age (in years)"){
  if(rate.p){
    est=res.nonpar$rate
    LB=est-c*sqrt(res.nonpar$rate.var)
    UB=est+c*sqrt(res.nonpar$rate.var)
  }else{
    est=res.nonpar$mean
    LB=est-c*sqrt(res.nonpar$mean.var)
    UB=est+c*sqrt(res.nonpar$mean.var)
  }
  
  if(ncol(est)<=12){
    col=c("black","red","blue","green","orange","purple",
          "brown4","gray","darkolivegreen1","aquamarine",
          "deeppink","deepskyblue")[1:ncol(est)]
  }else{
    col=rainbow(ncol(est))
  }
  
  if(is.na(max.L.m)){
    max.L.m=max(UB)
  }
  
  agev=c(1:ONUM)/AGEUNIT
  leg=colnames(est)
  
  if(!two.plot){
    matplot(agev,est,type="l",lty=1,lwd=2,
            col=col,xlab=xlab,
            ylab=ifelse(rate.p,"Rate function","Mean function"),
            main=ifelse(rate.p,"Marginal rate function","Marginal mean function"),
            ylim = c(0,max.L.m))
    matlines(agev,LB,lty=2,lwd=1,col=col)
    matlines(agev,UB,lty=2,lwd=1,col=col)
    
    legend("topleft",legend=leg,lty=rep(1,ncol(est)),
           lwd=rep(2,ncol(est)),
           col=col,cex=cex)
  }else{
    h=ncol(est)/2
    matplot(agev,est[,1:h],type="l",lty=1,lwd=2,
            col=col[1:h],xlab=xlab,
            ylab=ifelse(rate.p,"Rate function","Mean function"),
            main=ifelse(rate.p,"Marginal rate function","Marginal mean function"),
            ylim = c(0,max.L.m))
    matlines(agev,LB[,1:h],lty=2,lwd=1,col=col[1:h])
    matlines(agev,UB[,1:h],lty=2,lwd=1,col=col[1:h])
    
    legend("topleft",legend=leg[1:h],lty=rep(1,h),
           lwd=rep(2,h),
           col=col[1:h],cex=cex)
    
    
    matplot(agev,est[,(h+1):(2*h)],type="l",lty=1,lwd=2,
            col=col[(h+1):(2*h)],xlab=xlab,
            ylab=ifelse(rate.p,"Rate function","Mean function"),
            main=ifelse(rate.p,"Marginal rate function","Marginal mean function"),
            ylim = c(0,max.L.m))
    matlines(agev,LB[,(h+1):(2*h)],lty=2,lwd=1,col=col[(h+1):(2*h)])
    matlines(agev,UB[,(h+1):(2*h)],lty=2,lwd=1,col=col[(h+1):(2*h)])
    
    legend("topleft",legend=leg[(h+1):(2*h)],lty=rep(1,h),
           lwd=rep(2,h),
           col=col[(h+1):(2*h)],cex=cex)
  }
}


