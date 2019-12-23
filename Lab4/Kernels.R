
# Square Exp Kernel  ------------------------------------------------------

SquaredExpKernel <- function(x1,x2,sigmaF,l){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}


# -------------------------------------------------------------------------

# Square Exponential Kernel for Kernel matrix -----------------------------

SEkernel=function(sigmaf = 1,ell=1){
  rval <- function(x, y = NULL) {
    if (!is(x, "vector")) 
      stop("x must be a vector")
    if (!is(y, "vector") && !is.null(y)) 
      stop("y must a vector")
    if (is(x, "vector") && is.null(y)) {
      return(1)
    }
    if (is(x, "vector") && is(y, "vector")) {
      if (!length(x) == length(y)) 
        stop("number of dimension must be the same on both data points")
      
      #res=(sigmaf^2)*exp(ell * (2 * crossprod(x, y) - crossprod(x) -crossprod(y)))
      r2=crossprod(x-y)
      res=(sigmaf^2)*exp(-(r2/(2*ell^2)))
      
      return(res)
    
    }
  }
  return(new("kernel", .Data = rval, kpar = list(sigmaf = sigmaf,ell=ell)))
}

# -------------------------------------------------------------------------

# Periodic Exp Kernel -----------------------------------------------------


PerExpKernel <- function(x1,x2,sigmaF,l1,l2,d){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    r1=-2*sin( (pi*abs(x1-x2[i])) /d)^2
    K[,i] <- sigmaF^2*exp(r1/l1^2)*exp(-0.5*( (x1-x2[i])/l2^2)^2 )
  }
  return(K)
}

# -------------------------------------------------------------------------

# Periodic Exp Kernel for Kernel Matrix -----------------------------------

perKernel=function(sigmaf,ell1,ell2,d){
  rval <- function(x, y = NULL) {
    if (!is(x, "vector")) 
      stop("x must be a vector")
    if (!is(y, "vector") && !is.null(y)) 
      stop("y must a vector")
    if (is(x, "vector") && is.null(y)) {
      return(1)
    }
    if (is(x, "vector") && is(y, "vector")) {
      if (!length(x) == length(y)) 
        stop("number of dimension must be the same on both data points")
      r2=crossprod(x-y)
      p1=(-2)*sin((pi*abs(x-y))/d)^2
      res=(sigmaf^2)*exp(p1/ell1^2)*exp(-r2/(2*ell2^2))
      return(res)
    }
  }
  return(new("kernel", .Data = rval, kpar = list(sigmaf = sigmaf,ell1=ell1,ell2,d)))
}

# -------------------------------------------------------------------------


sigmaf=20 ; ell1=1 ; ell2=10 ; d=0.5766923
X=c(1,3,4) ; Xstar=c(2,3,4) # input vectors


SquaredExpKernel(X,Xstar,sigmaf,ell1)
kernelMatrix(SEkernel(sigmaf,ell1), x = X, y = Xstar)


PerExpKernel(X,Xstar,sigmaf,ell1,ell2,d)
kernelMatrix(perKernel(sigmaf,ell1,ell2,d), x = X, y = Xstar)








