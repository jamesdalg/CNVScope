#'  Find the negative log p-value of a pair of vectors.
#'
#' Finds the negative log p-value of a matrix, if it exists.
#' Checks first to see if there is a p-value to return.
#' @keywords lm linear regression 
#' @param x a vector that is regressed in the fashion y~x.
#' @param y a vector that is regressed in the fashion y~x.
#' @param repval the replacement value if the regression cannot be performed, default 300 (the vectors are identical if this is used).
#' @param lowrepval The low replacement value in the case that a regression p-value is undefined.
#' @param signed change the sign of the negative log p-value based on the sign of beta?
#'  e.g. if the line has a negative slope, so will the returned value.
#'   If there is a positive slope, there will be a positive negative log p-value.
#'   if this option is disabled, then no sign changes will happen based on the sign of the slope. 
#' @return The negative log p-value or replacement value.
#' @examples
#' #small example
#' xval<-c(1,1,1,1,1)
#' yval<-c(1,2,3,4,5)
#' a<-c(3,4,5,6,7)
#' extractNegLogPval(x=xval,y=yval) #no possible p-value if one vector is constant.
#' #Some edge cases this may not be correct (if the data lies near a constant),
#' # but the indiviual sample data should reveal true trends.
#' suppressWarnings(cor(xval,yval)) #you can't get a correlation value either.
#' cor(a,a) #gives correlation of 1.
#' extractNegLogPval(a,a) 
#' #gives replacement value.
#' suppressWarnings(extractNegLogPval(x=a,y=yval))
#' #gives 107.3909 and warns about a nearly perfect fit.
#' @export
extractNegLogPval<-function(x,y,repval=300,lowrepval=0,signed=F)
{
  if(all.equal(x,y)==T){neglogpval<-repval} else{
    coef<-summary(lm(y~x))$coefficients
    if(nrow(coef)>=2&ncol(coef)>=4) {
      neglogpval<-(-log(coef[2,4]))
      if(signed){neglogpval<-sign(coef[2,1])*neglogpval}
    } else {neglogpval<-lowrepval}
    
  }
  return(neglogpval)
}