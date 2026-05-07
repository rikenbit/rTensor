###Tensor Decompositions

#'(Truncated-)Higher-order SVD
#'
#'Higher-order SVD of a K-Tensor. Write the K-Tensor as a (m-mode) product of a core Tensor (possibly smaller modes) and K orthogonal factor matrices. Truncations can be specified via \code{ranks} (making them smaller than the original modes of the K-Tensor will result in a truncation). For the mathematical details on HOSVD, consult Lathauwer et. al. (2000).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors.
#'@name hosvd
#'@rdname hosvd
#'@aliases hosvd
#'@param tnsr Tensor with K modes
#'@param ranks a vector of desired modes in the output core tensor, default is \code{tnsr@@modes}
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{core tensor with modes speficied by \code{ranks}}
#'\item{\code{U}}{a list of orthogonal matrices, one for each mode}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)} - if there was no truncation, then this is on the order of mach_eps * fnorm. }
#'}
#'@seealso \code{\link{tucker}}
#'@references L. Lathauwer, B.Moor, J. Vanderwalle "A multilinear singular value decomposition". Journal of Matrix Analysis and Applications 2000.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(6,7,8))
#'hosvdD <-hosvd(tnsr)
#'hosvdD$fnorm_resid
#'hosvdD2 <-hosvd(tnsr,ranks=c(3,3,4))
#'hosvdD2$fnorm_resid
hosvd <- function(tnsr,ranks=NULL){
	stopifnot(is(tnsr,"Tensor"))
	if (sum(ranks<=0)!=0) stop("ranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
	
	num_modes <- tnsr@num_modes
	#no truncation if ranks not provided
	if(is.null(ranks)){
		ranks <- tnsr@modes
	}else{
		if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=num_modes,style=3)
	#loops through and performs SVD on mode-m matricization of tnsr
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		temp_mat <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
		setTxtProgressBar(pb,m)
	}
	close(pb)
	#computes the core tensor
	Z <- ttl(tnsr,lapply(U_list,t),ms=1:num_modes)
	est <- ttl(Z,U_list,ms=1:num_modes)
	resid <- fnorm(est-tnsr)
	#put together the return list, and returns
	list(Z=Z,U=U_list,est=est,fnorm_resid=resid)	
}

#'Canonical Polyadic Decomposition
#'
#'Canonical Polyadic (CP) decomposition of a tensor, aka CANDECOMP/PARAFRAC. Approximate a K-Tensor using a sum of \code{num_components} rank-1 K-Tensors. A rank-1 K-Tensor can be written as an outer product of K vectors. There are a total of \code{num_compoents *tnsr@@num_modes} vectors in the output, stored in \code{tnsr@@num_modes} matrices, each with \code{num_components} columns. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on CP decomposition, consult Kolda and Bader (2009).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name cp
#'@rdname cp
#'@aliases cp
#'@param tnsr Tensor with K modes
#'@param num_components the number of rank-1 K-Tensors to use in approximation
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following \describe{
#'\item{\code{lambdas}}{a vector of normalizing constants, one for each component}
#'\item{\code{U}}{a list of matrices - one for each mode - each matrix with \code{num_components} columns}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{tucker}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@examples
#'subject <- faces_tnsr[,,14,]
#'cpD <- cp(subject,num_components=10) 
#'cpD$conv 
#'cpD$norm_percent 
#'plot(cpD$all_resids) 
cp <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	tnsr_norm <- fnorm(tnsr)
	for(m in 1:num_modes){
		unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
	}
	est <- tnsr
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(est){
		curr_resid <- fnorm(est - tnsr)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{ return(FALSE)}
	}	
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	norm_vec <- function(vec){
	norm(as.matrix(vec))
	}
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)
		for(m in 1:num_modes){
			V <- hadamard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
			V_inv <- solve(V)			
			tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
			lambdas <- apply(tmp,2,norm_vec)
			U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
			Z <- .superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
			est <- ttl(Z,U_list,ms=1:num_modes)
		}
		#checks convergence
		if(CHECK_CONV(est)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			 }
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
}

#'Tucker Decomposition
#'
#'The Tucker decomposition of a tensor. Approximates a K-Tensor using a n-mode product of a core tensor (with modes specified by \code{ranks}) with orthogonal factor matrices. If there is no truncation in one of the modes, then this is the same as the MPCA, \code{\link{mpca}}. If there is no truncation in all the modes (i.e. \code{ranks = tnsr@@modes}), then this is the same as the HOSVD, \code{\link{hosvd}}. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the Tucker decomposition, consult Kolda and Bader (2009).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure also known as Higher-Order Orthogonal Iteration (HOOI). Intialized using a (Truncated-)HOSVD. A progress bar is included to help monitor operations on large tensors.
#'@name tucker
#'@rdname tucker
#'@aliases tucker
#'@param tnsr Tensor with K modes
#'@param ranks a vector of the modes of the output core Tensor
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{Z}}{the core tensor, with modes specified by \code{ranks}}
#'\item{\code{U}}{a list of orthgonal factor matrices - one for each mode, with the number of columns of the matrices given by \code{ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{hosvd}}, \code{\link{mpca}}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes}.
#'@examples
#'tnsr <- rand_tensor(c(4,4,4,4))
#'tuckerD <- tucker(tnsr,ranks=c(2,2,2,2))
#'tuckerD$conv 
#'tuckerD$norm_percent
#'plot(tuckerD$all_resids)
tucker <- function(tnsr,ranks=NULL,max_iter=25,tol=1e-5){
	stopifnot(is(tnsr,"Tensor"))
	if(is.null(ranks)) stop("ranks must be specified")
	if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
	if (sum(ranks<=0)!=0) stop("ranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	#initialization via truncated hosvd
	num_modes <- tnsr@num_modes
	U_list <- vector("list",num_modes)
	for(m in 1:num_modes){
		temp_mat <- rs_unfold(tnsr,m=m)@data
		U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
	}
	tnsr_norm <- fnorm(tnsr)
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(Z,U_list){
		est <- ttl(Z,U_list,ms=1:num_modes)
		curr_resid <- fnorm(tnsr - est)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)	
	modes <- tnsr@modes
	modes_seq <- 1:num_modes
		for(m in modes_seq){
			#core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-m],t),ms=modes_seq[-m])
			#truncated SVD of X
			#U_list[[m]] <- (svd(rs_unfold(X,m=m)@data,nu=ranks[m],nv=prod(modes[-m]))$u)[,1:ranks[m]]
			U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
		}
		#compute core tensor Z
		Z <- ttm(X,mat=t(U_list[[num_modes]]),m=num_modes)

		#checks convergence
		if(CHECK_CONV(Z, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)	
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
	est <- ttl(Z,U_list,ms=1:num_modes)
	invisible(list(Z=Z, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}

#'Multilinear Principal Components Analysis
#'
#'This is basically the Tucker decomposition of a K-Tensor, \code{\link{tucker}}, with one of the modes uncompressed. If K = 3, then this is also known as the Generalized Low Rank Approximation of Matrices (GLRAM). This implementation assumes that the last mode is the measurement mode and hence uncompressed. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on the MPCA of tensors, consult Lu et al. (2008).
#'@export
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name mpca
#'@rdname mpca
#'@aliases mpca
#'@param tnsr Tensor with K modes
#'@param ranks a vector of the compressed modes of the output core Tensor, this has length K-1
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{Z_ext}}{the extended core tensor, with the first K-1 modes given by \code{ranks}}
#'\item{\code{U}}{a list of K-1 orthgonal factor matrices - one for each compressed mode, with the number of columns of the matrices given by \code{ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{tucker}}, \code{\link{hosvd}}
#'@references H. Lu, K. Plataniotis, A. Venetsanopoulos, "Mpca: Multilinear principal component analysis of tensor objects". IEEE Trans. Neural networks, 2008.
#'@note The length of \code{ranks} must match \code{tnsr@@num_modes-1}.
#'@examples
#'subject <- faces_tnsr[,,21,]
#'mpcaD <- mpca(subject,ranks=c(10,10))
#'mpcaD$conv
#'mpcaD$norm_percent
#'plot(mpcaD$all_resids)
mpca <- function(tnsr, ranks = NULL, max_iter = 25, tol=1e-5){
	if(is.null(ranks)) stop("ranks must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if (sum(ranks>tnsr@modes)!=0) stop("ranks must be smaller than the corresponding mode")
	if (sum(ranks<=0)!=0) stop("ranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	#initialization via hosvd of M-1 modes
	num_modes <- tnsr@num_modes
	stopifnot(length(ranks)==(num_modes-1))
	ranks <- c(ranks,1)
	modes <- tnsr@modes
	U_list <- vector("list",num_modes)
	unfolded_mat <- vector("list",num_modes)
	for(m in 1:(num_modes-1)){
		unfolded_mat <- rs_unfold(tnsr,m=m)@data
		mode_m_cov <- unfolded_mat%*%t(unfolded_mat)
		U_list[[m]] <- svd(mode_m_cov, nu=ranks[m])$u
	}
	Z_ext <- ttl(tnsr,lapply(U_list[-num_modes],t),ms=1:(num_modes-1))
	tnsr_norm <- fnorm(tnsr)
	curr_iter <- 1
	converged <- FALSE
	#set up convergence check
	fnorm_resid <- rep(0, max_iter)
	CHECK_CONV <- function(Z_ext,U_list){
		est <- ttl(Z_ext,U_list[-num_modes],ms=1:(num_modes-1))
		curr_resid <- fnorm(tnsr - est)
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}
	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (until convergence or max_iter)
	while((curr_iter < max_iter) && (!converged)){
	setTxtProgressBar(pb,curr_iter)
	modes <-tnsr@modes
	modes_seq <- 1:(num_modes-1)
		for(m in modes_seq){
			#extended core Z minus mode m
			X <- ttl(tnsr,lapply(U_list[-c(m,num_modes)],t),ms=modes_seq[-m])
			#truncated SVD of X
			U_list[[m]] <- svd(rs_unfold(X,m=m)@data,nu=ranks[m])$u
		}
		#compute core tensor Z_ext
		Z_ext <- ttm(X,mat=t(U_list[[num_modes-1]]),m=num_modes-1)
		#checks convergence
		if(CHECK_CONV(Z_ext, U_list)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
			}
	}
	close(pb)
	#end of main loop
	#put together return list, and returns
	est <- ttl(Z_ext,U_list[-num_modes],ms=1:(num_modes-1))
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent<-(1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(Z_ext=Z_ext, U=U_list, conv=converged, est=est, norm_percent = norm_percent, fnorm_resid=tail(fnorm_resid,1), all_resids=fnorm_resid))
}

#'Population Value Decomposition
#'
#'The default Population Value Decomposition (PVD) of a series of 2D images. Constructs population-level matrices P, V, and D to account for variances within as well as across the images. Structurally similar to Tucker (\code{\link{tucker}}) and GLRAM (\code{\link{mpca}}), but retains crucial differences. Requires \code{2*n3 + 2} parameters to specified the final ranks of P, V, and D, where n3 is the third mode (how many images are in the set). Consult Crainiceanu et al. (2013) for the construction and rationale behind the PVD model.
#'@export
#'@details The PVD is not an iterative method, but instead relies on \code{n3 + 2}separate PCA decompositions. The third mode is for how many images are in the set.
#'@name pvd
#'@rdname pvd
#'@aliases pvd
#'@param tnsr 3-Tensor with the third mode being the measurement mode
#'@param uranks ranks of the U matrices
#'@param wranks ranks of the W matrices
#'@param a rank of \code{P = U\%*\%t(U)}
#'@param b rank of \code{D = W\%*\%t(W)}
#'@return a list containing the following:\describe{
#'\item{\code{P}}{population-level matrix \code{P = U\%*\%t(U)}, where U is constructed by stacking the truncated left eigenvectors of slicewise PCA along the third mode}
#'\item{\code{V}}{a list of image-level core matrices}
#'\item{\code{D}}{population-leve matrix \code{D = W\%*\%t(W)}, where W is constructed by stacking the truncated right eigenvectors of slicewise PCA along the third mode}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'}
#'@references C. Crainiceanu, B. Caffo, S. Luo, V. Zipunnikov, N. Punjabi, "Population value decomposition: a framework for the analysis of image populations". Journal of the American Statistical Association, 2013.
#'@examples
#'subject <- faces_tnsr[,,8,]
#'pvdD<-pvd(subject,uranks=rep(46,10),wranks=rep(56,10),a=46,b=56)
pvd <- function(tnsr,uranks=NULL,wranks=NULL,a=NULL,b=NULL){
	if(tnsr@num_modes!=3) stop("PVD only for 3D")
	if (sum(uranks<=0)!=0) stop("uranks must be positive")
	if (sum(wranks<=0)!=0) stop("wranks must be positive")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	if(is.null(uranks)||is.null(wranks)) stop("U and V ranks must be specified")
	if(is.null(a)||is.null(b)) stop("a and b must be specified")
	modes <- tnsr@modes
	n <- modes[3]
	if(length(uranks)!=n||length(wranks)!=n) stop("ranks must be of length n3")
	pb <- txtProgressBar(min=0,max=(n+3),style=3)
	x <- tnsr@data
	Us <- vector('list',n)
	Vs <- vector('list',n)
	S <- vector('list',n)
	for(i in 1:n){
		svdz <- svd(x[,,i],nu=uranks[i],nv=wranks[i])
		Us[[i]] <- svdz$u
		Vs[[i]] <- svdz$v
		S[[i]] <- svdz$d[1:min(uranks[i],wranks[i])]
		setTxtProgressBar(pb,i)
	}
	U <- matrix(unlist(Us),nrow=modes[1],ncol=sum(uranks)*n)
	#eigenU <- eigen(U%*%t(U))
	P <- eigen(U%*%t(U))$vectors[,1:a] #E-vecs of UU^T
	setTxtProgressBar(pb,n+1)
	V <- matrix(unlist(Vs),nrow=modes[2],ncol=sum(wranks)*n)
	#eigenV <- eigen(V%*%t(V))
	Dt <- eigen(V%*%t(V))$vectors[,1:b] #E-vecs of VV^T
	D <- t(Dt)
	setTxtProgressBar(pb,n+2)
	V2 <- vector('list',n)
	est <- array(0,dim=modes)
	for(i in 1:n){
		V2[[i]] <- (t(P)%*%Us[[i]])%*%diag(S[[i]],nrow=uranks[i],ncol=wranks[i])%*%(t(Vs[[i]])%*%Dt)
		est[,,i] <- P%*%V2[[i]]%*%D
	}
	est <- as.tensor(est)
	fnorm_resid <- fnorm(est-tnsr)	
	setTxtProgressBar(pb,n+3)
	norm_percent<-(1-(fnorm_resid/fnorm(tnsr)))*100
	invisible(list(P=P,D=D,V=V2,est=est,norm_percent=norm_percent,fnorm_resid=fnorm_resid))
}

#'Tensor Singular Value Decomposition
#'
#'TSVD for a 3-Tensor. Constructs 3-Tensors \code{U, S, V} such that \code{tnsr = t_mult(t_mult(U,S),t(V))}. \code{U} and \code{V} are orthgonal 3-Tensors with orthogonality defined in Kilmer et al. (2013), and \code{S} is a 3-Tensor consists of facewise diagonal matrices. For more details on the TSVD, consult Kilmer et al. (2013).
#'@export
#'@name t_svd
#'@rdname t_svd
#'@aliases t_svd
#'@param tnsr 3-Tensor to decompose via TSVD
#'@return a list containing the following:\describe{
#'\item{\code{U}}{the left orthgonal 3-Tensor}
#'\item{\code{V}}{the right orthgonal 3-Tensor}
#'\item{\code{S}}{the middle 3-Tensor consisting of face-wise diagonal matrices}
#'}
#'@seealso \code{\link{t_mult}}, \code{\link{t_svd_reconstruct}}
#'@references M. Kilmer, K. Braman, N. Hao, and R. Hoover, "Third-order tensors as operators on matrices: a theoretical and computational framework with applications in imaging". SIAM Journal on Matrix Analysis and Applications 2013.
#'@note Computation involves complex values, but if the inputs are real, then the outputs are also real. Some loss of precision occurs in the truncation of the imaginary components during the FFT and inverse FFT.
#'@examples
#'tnsr <- rand_tensor()
#'tsvdD <- t_svd(tnsr)
t_svd<-function(tnsr){
	if(tnsr@num_modes!=3) stop("T-SVD only implemented for 3d so far")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	modes <- tnsr@modes
	n1 <- modes[1]
	n2 <- modes[2]
	n3 <- modes[3]
	#progress bar
	pb <- txtProgressBar(min=0,max=n3,style=3)
	#define ifft
	#.ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
	#fft for each of the n1n2 vectors (of length n3) along mode 3
	fftz <- aperm(apply(tnsr@data,MARGIN=1:2,fft),c(2,3,1))
	#svd for each face (svdz is a list of the results)
	U_arr <- array(0,dim=c(n1,n1,n3))
	V_arr <- array(0,dim=c(n2,n2,n3))
	m <- min(n1,n2)		
	S_arr <- array(0,dim=c(n1,n2,n3))
	#Think of a way to avoid a loop in the beginning
	#Problem is that svd returns a list but ideally we want 3 arrays
	#Even with unlist this doesn't seem possible
	for (j in 1:n3){
		setTxtProgressBar(pb,j)
		decomp <- svd(fftz[,,j],nu=n1,nv=n2)
		U_arr[,,j] <- decomp$u
		V_arr[,,j] <- decomp$v
		S_arr[,,j] <- diag(decomp$d,nrow=n1,ncol=n2) #length is min(n1,n2)
	}	
	close(pb)
	#for each svd result, we want to apply ifft
	U <- as.tensor(aperm(apply(U_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	V <- as.tensor(aperm(apply(V_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	S <- as.tensor(aperm(apply(S_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	invisible(list(U=U,V=V,S=S))
}

#'Reconstruct Tensor From TSVD
#'
#'Reconstruct the original 3-Tensor after it has been decomposed into \code{U, S, V} via \code{\link{t_svd}}.
#'@export
#'@name t_svd_reconstruct
#'@rdname t_svd_reconstruct
#'@aliases t_svd_reconstruct
#'@param L list that is an output from \code{\link{t_svd}}
#'@return a 3-Tensor 
#'@seealso \code{\link{t_svd}}
#'@examples
#'tnsr <- rand_tensor(c(10,10,10))
#'tsvdD <- t_svd(tnsr)
#'1 - fnorm(t_svd_reconstruct(tsvdD)-tnsr)/fnorm(tnsr)
t_svd_reconstruct <- function(L){
	t_mult(t_mult(L$U,L$S),t(L$V))
}

#####
.is_zero_tensor <- function(tnsr){
	if (sum(tnsr@data==0)==prod(tnsr@modes)) return(TRUE)
	return(FALSE)
}

#'INDSCAL Decomposition
#'
#'Individual Differences Scaling (INDSCAL) decomposition of a 3-Tensor. Decomposes a symmetric 3-Tensor into a shared factor matrix \code{A} and slice-specific diagonal weight matrices \code{D_k} such that each frontal slice \code{X_k = A \%*\% D_k \%*\% t(A)}. Uses the Alternating Least Squares (ALS) estimation procedure. For more details on INDSCAL, consult Carroll and Chang (1970).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors. The input tensor must be 3-dimensional with square frontal slices (i.e. the first two modes must be equal).
#'@name indscal
#'@rdname indscal
#'@aliases indscal
#'@param tnsr 3-Tensor with square frontal slices
#'@param num_components the number of components for the decomposition
#'@param max_iter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{A}}{the shared factor matrix with \code{num_components} columns}
#'\item{\code{D}}{a list of diagonal weight matrices, one for each frontal slice}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after decomposition}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{cp}}, \code{\link{rescal}}
#'@references J. Carroll, J. Chang, "Analysis of individual differences in multidimensional scaling via an N-way generalization of Eckart-Young decomposition". Psychometrika 1970.
#'@note The first two modes of \code{tnsr} must be equal (square frontal slices).
#'@examples
#'tnsr <- rand_tensor(c(4,4,3))
#'# make symmetric frontal slices
#'for(k in 1:3) tnsr[,,k] <- (tnsr[,,k] + t(tnsr[,,k]@data)) / 2
#'indscalD <- indscal(tnsr, num_components=2)
#'indscalD$conv
#'indscalD$norm_percent
#'plot(indscalD$all_resids)
indscal <- function(tnsr, num_components=NULL, max_iter=25, tol=1e-5){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if(tnsr@num_modes!=3) stop("INDSCAL only for 3D tensors")
	if(tnsr@modes[1]!=tnsr@modes[2]) stop("INDSCAL requires square frontal slices")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	modes <- tnsr@modes
	I <- modes[1]
	K <- modes[3]
	R <- num_components
	x <- tnsr@data
	tnsr_norm <- fnorm(tnsr)

	#initialization
	A <- qr.Q(qr(matrix(rnorm(I * R), I, R)))
	D_list <- vector("list", K)

	curr_iter <- 1
	converged <- FALSE
	fnorm_resid <- rep(0, max_iter)

	CHECK_CONV <- function(est_arr){
		curr_resid <- sqrt(sum((est_arr - x)^2))
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}

	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop
	while((curr_iter < max_iter) && (!converged)){
		setTxtProgressBar(pb,curr_iter)

		#update D_k for each slice
		AtA <- t(A) %*% A
		for(k in 1:K){
			d_k <- diag(t(A) %*% x[,,k] %*% A) / diag(AtA * AtA)
			D_list[[k]] <- diag(d_k, nrow=R, ncol=R)
		}

		#update A: stack all slices into a big least-squares problem
		# X_(1) = [X_1; X_2; ...;X_K] = [A D_1; A D_2; ...; A D_K] A^T
		# => solve for A using stacked normal equations
		lhs <- matrix(0, I, R)
		rhs <- matrix(0, R, R)
		for(k in 1:K){
			lhs <- lhs + x[,,k] %*% A %*% D_list[[k]]
			rhs <- rhs + D_list[[k]] %*% AtA %*% D_list[[k]]
		}
		A <- lhs %*% solve(rhs)

		#reconstruct estimate
		est_arr <- array(0, dim=modes)
		for(k in 1:K){
			est_arr[,,k] <- A %*% D_list[[k]] %*% t(A)
		}

		if(CHECK_CONV(est_arr)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
		}
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)

	#put together return list
	est <- as.tensor(est_arr)
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent <- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(A=A, D=D_list, conv=converged, est=est,
		norm_percent=norm_percent, fnorm_resid=tail(fnorm_resid,1),
		all_resids=fnorm_resid))
}

#'RESCAL Decomposition
#'
#'RESCAL decomposition of a 3-Tensor for relational data. Decomposes a 3-Tensor into a shared entity factor matrix \code{A} and slice-specific core matrices \code{R_k} such that each frontal slice \code{X_k = A \%*\% R_k \%*\% t(A)}. Uses the Alternating Least Squares (ALS) estimation procedure. For more details on RESCAL, consult Nickel et al. (2011).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors. The input tensor must be 3-dimensional with square frontal slices (i.e. the first two modes must be equal). Unlike \code{\link{indscal}}, the core matrices \code{R_k} are dense (not restricted to be diagonal), which allows modeling of asymmetric relations.
#'@name rescal
#'@rdname rescal
#'@aliases rescal
#'@param tnsr 3-Tensor with square frontal slices
#'@param num_components the number of components for the decomposition
#'@param max_iter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{A}}{the shared entity factor matrix with \code{num_components} columns}
#'\item{\code{R}}{a list of core matrices (one per frontal slice), each of size \code{num_components} by \code{num_components}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after decomposition}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{indscal}}, \code{\link{dedicom}}
#'@references M. Nickel, V. Tresp, H. Kriegel, "A Three-Way Model for Collective Learning on Multi-Relational Data". Proceedings of the 28th International Conference on Machine Learning 2011.
#'@note The first two modes of \code{tnsr} must be equal (square frontal slices).
#'@examples
#'tnsr <- rand_tensor(c(5,5,3))
#'rescalD <- rescal(tnsr, num_components=2)
#'rescalD$conv
#'rescalD$norm_percent
#'plot(rescalD$all_resids)
rescal <- function(tnsr, num_components=NULL, max_iter=25, tol=1e-5){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if(tnsr@num_modes!=3) stop("RESCAL only for 3D tensors")
	if(tnsr@modes[1]!=tnsr@modes[2]) stop("RESCAL requires square frontal slices")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	modes <- tnsr@modes
	N <- modes[1]
	K <- modes[3]
	R <- num_components
	x <- tnsr@data
	tnsr_norm <- fnorm(tnsr)

	#initialization
	A <- qr.Q(qr(matrix(rnorm(N * R), N, R)))
	R_list <- vector("list", K)

	curr_iter <- 1
	converged <- FALSE
	fnorm_resid <- rep(0, max_iter)

	CHECK_CONV <- function(est_arr){
		curr_resid <- sqrt(sum((est_arr - x)^2))
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}

	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop
	while((curr_iter < max_iter) && (!converged)){
		setTxtProgressBar(pb,curr_iter)

		#update R_k for each slice: R_k = (A^T A)^{-1} A^T X_k A (A^T A)^{-1}
		AtA <- t(A) %*% A
		AtA_inv <- solve(AtA)
		for(k in 1:K){
			R_list[[k]] <- AtA_inv %*% (t(A) %*% x[,,k] %*% A) %*% AtA_inv
		}

		#update A: solve stacked normal equations
		# sum_k X_k A R_k^T + X_k^T A R_k = sum_k A R_k AtA R_k^T + A R_k^T AtA R_k
		# Simplified: vec(A) via Kronecker, but easier to do direct gradient step
		# Actually use direct least squares: stack [X_1; ...; X_K] = stack [A R_1; ...; A R_K] A^T
		# => A = (sum_k X_k A R_k^T + X_k^T A R_k) (sum_k R_k AtA R_k^T + R_k^T AtA R_k)^{-1} / 2
		# Simpler: solve sum_k X_k %*% A %*% R_k^T = A %*% sum_k R_k %*% AtA %*% R_k^T
		lhs <- matrix(0, N, R)
		rhs <- matrix(0, R, R)
		for(k in 1:K){
			lhs <- lhs + x[,,k] %*% A %*% t(R_list[[k]])
			lhs <- lhs + t(x[,,k]) %*% A %*% R_list[[k]]
			rhs <- rhs + R_list[[k]] %*% AtA %*% t(R_list[[k]])
			rhs <- rhs + t(R_list[[k]]) %*% AtA %*% R_list[[k]]
		}
		A <- lhs %*% solve(rhs)

		#reconstruct estimate
		est_arr <- array(0, dim=modes)
		for(k in 1:K){
			est_arr[,,k] <- A %*% R_list[[k]] %*% t(A)
		}

		if(CHECK_CONV(est_arr)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
		}
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)

	#put together return list
	est <- as.tensor(est_arr)
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent <- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(A=A, R=R_list, conv=converged, est=est,
		norm_percent=norm_percent, fnorm_resid=tail(fnorm_resid,1),
		all_resids=fnorm_resid))
}

#'DEDICOM Decomposition
#'
#'Decomposition into Directional Components (DEDICOM) of a 3-Tensor. Decomposes a 3-Tensor into a shared factor matrix \code{A}, a shared asymmetric relation matrix \code{R}, and slice-specific diagonal weight matrices \code{D_k} such that each frontal slice \code{X_k = A \%*\% D_k \%*\% R \%*\% D_k \%*\% t(A)}. Uses the Alternating Least Squares (ALS) estimation procedure. For more details on DEDICOM, consult Bader et al. (2007).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors. The input tensor must be 3-dimensional with square frontal slices (i.e. the first two modes must be equal).
#'@name dedicom
#'@rdname dedicom
#'@aliases dedicom
#'@param tnsr 3-Tensor with square frontal slices
#'@param num_components the number of components for the decomposition
#'@param max_iter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{A}}{the shared factor matrix with \code{num_components} columns}
#'\item{\code{R}}{the asymmetric relation matrix of size \code{num_components} by \code{num_components}}
#'\item{\code{D}}{a list of diagonal weight matrices, one for each frontal slice}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after decomposition}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{rescal}}, \code{\link{indscal}}
#'@references B. Bader, R. Harshman, T. Kolda, "Temporal analysis of semantic graphs using ASALSAN". Proceedings of the 7th IEEE International Conference on Data Mining 2007.
#'@note The first two modes of \code{tnsr} must be equal (square frontal slices).
#'@examples
#'tnsr <- rand_tensor(c(5,5,3))
#'dedicomD <- dedicom(tnsr, num_components=2)
#'dedicomD$conv
#'dedicomD$norm_percent
#'plot(dedicomD$all_resids)
dedicom <- function(tnsr, num_components=NULL, max_iter=25, tol=1e-5){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if(tnsr@num_modes!=3) stop("DEDICOM only for 3D tensors")
	if(tnsr@modes[1]!=tnsr@modes[2]) stop("DEDICOM requires square frontal slices")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	modes <- tnsr@modes
	N <- modes[1]
	K <- modes[3]
	R <- num_components
	x <- tnsr@data
	tnsr_norm <- fnorm(tnsr)

	#initialization via reparametrization: let M_k = A D_k, so X_k ~ M_k R M_k^T
	A <- qr.Q(qr(matrix(rnorm(N * R), N, R)))
	R_mat <- diag(1, R)
	D_list <- vector("list", K)
	M_list <- vector("list", K)
	for(k in 1:K){
		D_list[[k]] <- diag(rep(1, R), nrow=R, ncol=R)
		M_list[[k]] <- A
	}

	curr_iter <- 1
	converged <- FALSE
	fnorm_resid <- rep(0, max_iter)

	CHECK_CONV <- function(est_arr){
		curr_resid <- sqrt(sum((est_arr - x)^2))
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}

	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop
	while((curr_iter < max_iter) && (!converged)){
		setTxtProgressBar(pb,curr_iter)

		#update R: given M_k, solve for R
		# min sum_k || X_k - M_k R M_k^T ||^2
		# => R = (sum_k M_k^T M_k kron M_k^T M_k)^{-1} vec(sum_k M_k^T X_k M_k)
		gram_sum <- matrix(0, R*R, R*R)
		rhs_vec <- rep(0, R*R)
		for(k in 1:K){
			MtM <- t(M_list[[k]]) %*% M_list[[k]]
			gram_sum <- gram_sum + kronecker(MtM, MtM)
			rhs_vec <- rhs_vec + as.vector(t(M_list[[k]]) %*% x[,,k] %*% M_list[[k]])
		}
		gram_sum_reg <- gram_sum + diag(1e-8, R*R)
		R_vec <- solve(gram_sum_reg, rhs_vec)
		R_mat <- matrix(R_vec, R, R)

		#update M_k: for each slice, solve M_k given R
		# X_k ~ M_k R M_k^T => vectorize and solve via Kronecker
		for(k in 1:K){
			# Use gradient-based update: dL/dM_k = -2 X_k M_k R^T - 2 X_k^T M_k R + 2 M_k R M_k^T M_k R^T + 2 M_k R^T M_k^T M_k R
			# Simplified ALS: fix M_k on RHS, solve for M_k on LHS
			# X_k M_k R^T + X_k^T M_k R = M_k (R M_k^T M_k R^T + R^T M_k^T M_k R)
			MtM <- t(M_list[[k]]) %*% M_list[[k]]
			lhs_M <- x[,,k] %*% M_list[[k]] %*% t(R_mat) + t(x[,,k]) %*% M_list[[k]] %*% R_mat
			rhs_M <- R_mat %*% MtM %*% t(R_mat) + t(R_mat) %*% MtM %*% R_mat + diag(1e-8, R)
			M_list[[k]] <- lhs_M %*% solve(rhs_M)
		}

		#extract A and D_k from M_k: A is shared, D_k is per-slice scaling
		# Compute weighted average of M_k column spaces
		M_avg <- Reduce("+", M_list) / K
		qr_avg <- qr(M_avg)
		A_new <- qr.Q(qr_avg)
		# D_k: project M_k onto A to get diagonal scaling
		AtA_inv <- solve(t(A_new) %*% A_new + diag(1e-10, R))
		for(k in 1:K){
			# M_k ~ A D_k => D_k = diag of (A^T A)^{-1} A^T M_k
			proj <- AtA_inv %*% t(A_new) %*% M_list[[k]]
			D_list[[k]] <- diag(diag(proj), nrow=R, ncol=R)
			# reconstruct M_k with DEDICOM structure
			M_list[[k]] <- A_new %*% D_list[[k]]
		}
		A <- A_new

		#reconstruct estimate
		est_arr <- array(0, dim=modes)
		for(k in 1:K){
			est_arr[,,k] <- A %*% D_list[[k]] %*% R_mat %*% D_list[[k]] %*% t(A)
		}

		if(CHECK_CONV(est_arr)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
		}
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)

	#put together return list
	est <- as.tensor(est_arr)
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent <- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(A=A, R=R_mat, D=D_list, conv=converged, est=est,
		norm_percent=norm_percent, fnorm_resid=tail(fnorm_resid,1),
		all_resids=fnorm_resid))
}


#'PARAFAC2 Decomposition
#'
#'PARAFAC2 decomposition of a 3-Tensor. Decomposes a 3-Tensor into slice-specific orthogonal matrices \code{H_k}, a shared profile matrix \code{B}, a mode-3 factor matrix \code{C}, and slice-specific diagonal weight matrices \code{D_k} such that each frontal slice \code{X_k = H_k \%*\% B \%*\% D_k \%*\% t(C)}. The key PARAFAC2 constraint is that \code{t(H_k) \%*\% H_k} is constant across all slices. Uses the ALS estimation procedure of Kiers et al. (1999). For more details on PARAFAC2, consult Harshman (1972).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors. The input tensor must be 3-dimensional. Unlike \code{\link{cp}}, PARAFAC2 allows different row spaces per slice while constraining their cross-products to be equal.
#'@name parafac2
#'@rdname parafac2
#'@aliases parafac2
#'@param tnsr 3-Tensor to decompose
#'@param num_components the number of components for the decomposition
#'@param max_iter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{H}}{a list of orthogonal matrices (one per frontal slice), each of size \code{mode1} by \code{num_components}}
#'\item{\code{B}}{the shared profile matrix of size \code{num_components} by \code{num_components}}
#'\item{\code{C}}{the mode-3 factor matrix of size \code{mode2} by \code{num_components}}
#'\item{\code{D}}{a list of diagonal weight matrices, one for each frontal slice}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after decomposition}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{cp}}, \code{\link{tucker}}
#'@references R. Harshman, "PARAFAC2: Mathematical and technical notes". UCLA Working Papers in Phonetics 1972.
#'@references H. Kiers, J. ten Berge, R. Bro, "PARAFAC2 - Part I. A direct fitting algorithm for the PARAFAC2 model". Journal of Chemometrics 1999.
#'@note The input tensor must be 3-dimensional.
#'@examples
#'tnsr <- rand_tensor(c(6,5,4))
#'pf2D <- parafac2(tnsr, num_components=2)
#'pf2D$conv
#'pf2D$norm_percent
#'plot(pf2D$all_resids)
parafac2 <- function(tnsr, num_components=NULL, max_iter=25, tol=1e-5){
	if(is.null(num_components)) stop("num_components must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if(tnsr@num_modes!=3) stop("PARAFAC2 only for 3D tensors")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	modes <- tnsr@modes
	I <- modes[1]
	J <- modes[2]
	K <- modes[3]
	R <- num_components
	x <- tnsr@data
	tnsr_norm <- fnorm(tnsr)

	#initialization
	C <- qr.Q(qr(matrix(rnorm(J * R), J, R)))
	H_list <- vector("list", K)
	D_list <- vector("list", K)
	B <- diag(R)
	for(k in 1:K){
		H_list[[k]] <- qr.Q(qr(matrix(rnorm(I * R), I, R)))
		D_list[[k]] <- diag(rep(1, R), nrow=R, ncol=R)
	}

	curr_iter <- 1
	converged <- FALSE
	fnorm_resid <- rep(0, max_iter)

	CHECK_CONV <- function(est_arr){
		curr_resid <- sqrt(sum((est_arr - x)^2))
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}

	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	#main loop (Kiers et al. 1999 ALS)
	while((curr_iter < max_iter) && (!converged)){
		setTxtProgressBar(pb,curr_iter)

		#step 1: update H_k via Procrustes rotation
		# For each k, H_k = U_k V_k^T where SVD(X_k C D_k B^T) = U_k S_k V_k^T
		for(k in 1:K){
			target <- x[,,k] %*% C %*% D_list[[k]] %*% t(B)
			svd_res <- svd(target, nu=I, nv=R)
			H_list[[k]] <- svd_res$u[,1:R,drop=FALSE] %*% t(svd_res$v[,1:R,drop=FALSE])
		}

		#step 2: update B and D_k
		# Define Y_k = H_k^T X_k C, then Y_k ~ B D_k (for each k)
		# Stack: [Y_1; Y_2;...; Y_K] = [B D_1; B D_2;...; B D_K]
		# This is a CP-like subproblem: Y_k = B * diag(d_k)
		Y_list <- vector("list", K)
		for(k in 1:K){
			Y_list[[k]] <- t(H_list[[k]]) %*% x[,,k] %*% C
		}

		# Update B: given D_k, solve B = (sum_k Y_k D_k) (sum_k D_k^2)^{-1}
		lhs_B <- matrix(0, R, R)
		rhs_B <- matrix(0, R, R)
		for(k in 1:K){
			lhs_B <- lhs_B + Y_list[[k]] %*% D_list[[k]]
			rhs_B <- rhs_B + D_list[[k]] %*% D_list[[k]]
		}
		B <- lhs_B %*% solve(rhs_B + diag(1e-10, R))

		# Update D_k: D_k = diag(diag(B^{-1} Y_k)) approximately
		# More precisely: d_k = diag(solve(B^T B) B^T Y_k)
		BtB_inv_Bt <- solve(t(B) %*% B + diag(1e-10, R)) %*% t(B)
		for(k in 1:K){
			d_k <- diag(BtB_inv_Bt %*% Y_list[[k]])
			D_list[[k]] <- diag(d_k, nrow=R, ncol=R)
		}

		#step 3: update C
		# min sum_k || X_k - H_k B D_k C^T ||^2 w.r.t. C
		# => C = (sum_k X_k^T H_k B D_k) (sum_k D_k B^T B D_k)^{-1}
		lhs_C <- matrix(0, J, R)
		rhs_C <- matrix(0, R, R)
		for(k in 1:K){
			lhs_C <- lhs_C + t(x[,,k]) %*% H_list[[k]] %*% B %*% D_list[[k]]
			rhs_C <- rhs_C + D_list[[k]] %*% t(B) %*% B %*% D_list[[k]]
		}
		C <- lhs_C %*% solve(rhs_C + diag(1e-10, R))

		#reconstruct estimate
		est_arr <- array(0, dim=modes)
		for(k in 1:K){
			est_arr[,,k] <- H_list[[k]] %*% B %*% D_list[[k]] %*% t(C)
		}

		if(CHECK_CONV(est_arr)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
		}
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)

	#put together return list
	est <- as.tensor(est_arr)
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	norm_percent <- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(H=H_list, B=B, C=C, D=D_list, conv=converged, est=est,
		norm_percent=norm_percent, fnorm_resid=tail(fnorm_resid,1),
		all_resids=fnorm_resid))
}

#'Two-Dimensional Linear Discriminant Analysis
#'
#'Two-Dimensional Linear Discriminant Analysis (2DLDA) for a 3-Tensor of matrix observations with class labels. Finds left and right projection matrices \code{L} and \code{R} that maximize class separability in the projected space. Each frontal slice is treated as a matrix observation. This is an iterative algorithm. For more details on 2DLDA, consult Ye et al. (2005).
#'@export
#'@details A progress bar is included to help monitor operations on large tensors. The input tensor must be 3-dimensional, where the third mode indexes the observations. Unlike other decompositions in this package, 2DLDA is a supervised method that requires class labels.
#'@name twodlda
#'@rdname twodlda
#'@aliases twodlda
#'@param tnsr 3-Tensor where the third mode indexes observations
#'@param labels a vector of class labels (length must equal the third mode)
#'@param r_ranks number of left projection vectors (columns of \code{L})
#'@param c_ranks number of right projection vectors (columns of \code{R})
#'@param max_iter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following:\describe{
#'\item{\code{L}}{the left projection matrix of size \code{mode1} by \code{r_ranks}}
#'\item{\code{R}}{the right projection matrix of size \code{mode2} by \code{c_ranks}}
#'\item{\code{Z}}{a list of projected feature matrices (one per observation), each of size \code{r_ranks} by \code{c_ranks}}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{est}}{estimate of \code{tnsr} after projection and back-projection}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@seealso \code{\link{mpca}}, \code{\link{tucker}}
#'@references J. Ye, R. Janardan, Q. Li, "Two-Dimensional Linear Discriminant Analysis". Advances in Neural Information Processing Systems 2005.
#'@note The length of \code{labels} must match the third mode of \code{tnsr}.
#'@examples
#'tnsr <- rand_tensor(c(5,4,10))
#'labels <- rep(c(1,2), each=5)
#'twodldaD <- twodlda(tnsr, labels=labels, r_ranks=2, c_ranks=2)
#'twodldaD$conv
#'twodldaD$norm_percent
twodlda <- function(tnsr, labels=NULL, r_ranks=NULL, c_ranks=NULL, max_iter=25, tol=1e-5){
	if(is.null(labels)) stop("labels must be specified")
	stopifnot(is(tnsr,"Tensor"))
	if(tnsr@num_modes!=3) stop("2DLDA only for 3D tensors")
	if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")

	modes <- tnsr@modes
	I <- modes[1]
	J <- modes[2]
	K <- modes[3]
	if(length(labels)!=K) stop("length of labels must match the third mode")
	classes <- unique(labels)
	num_classes <- length(classes)
	if(num_classes < 2) stop("at least 2 classes are required")

	x <- tnsr@data
	tnsr_norm <- fnorm(tnsr)

	if(is.null(r_ranks)) r_ranks <- min(I, num_classes)
	if(is.null(c_ranks)) c_ranks <- min(J, num_classes)

	#compute grand mean
	grand_mean <- apply(x, 1:2, mean)

	#compute class means
	class_means <- vector("list", num_classes)
	class_sizes <- integer(num_classes)
	for(g in 1:num_classes){
		idx <- which(labels == classes[g])
		class_sizes[g] <- length(idx)
		if(length(idx)==1){
			class_means[[g]] <- x[,,idx]
		}else{
			class_means[[g]] <- apply(x[,,idx,drop=FALSE], 1:2, mean)
		}
	}

	#initialize R with random orthonormal columns
	R_mat <- qr.Q(qr(matrix(rnorm(J * c_ranks), J, c_ranks)))

	curr_iter <- 1
	converged <- FALSE
	fnorm_resid <- rep(0, max_iter)

	CHECK_CONV <- function(est_arr){
		curr_resid <- sqrt(sum((est_arr - x)^2))
		fnorm_resid[curr_iter] <<- curr_resid
		if (curr_iter==1) return(FALSE)
		if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
		else{return(FALSE)}
	}

	#progress bar
	pb <- txtProgressBar(min=0,max=max_iter,style=3)
	L_mat <- NULL
	#main loop
	while((curr_iter < max_iter) && (!converged)){
		setTxtProgressBar(pb,curr_iter)

		#step 1: compute row scatter matrices given R_mat
		# Project columns: X_k R -> I x c_ranks, then compute scatter in I-space
		Sb_row <- matrix(0, I, I)
		Sw_row <- matrix(0, I, I)
		grand_proj_row <- grand_mean %*% R_mat
		for(g in 1:num_classes){
			mean_proj <- class_means[[g]] %*% R_mat
			diff_row <- mean_proj - grand_proj_row
			Sb_row <- Sb_row + class_sizes[g] * (diff_row %*% t(diff_row))
			idx <- which(labels == classes[g])
			for(i in idx){
				diff_i <- x[,,i] %*% R_mat - mean_proj
				Sw_row <- Sw_row + diff_i %*% t(diff_i)
			}
		}
		# Solve generalized eigenvalue problem: Sb_row v = lambda Sw_row v
		Sw_row_reg <- Sw_row + diag(1e-8, I)
		eig_row <- eigen(solve(Sw_row_reg) %*% Sb_row)
		L_mat <- Re(eig_row$vectors[, 1:r_ranks, drop=FALSE])

		#step 2: compute column scatter matrices given L_mat
		Sb_col <- matrix(0, J, J)
		Sw_col <- matrix(0, J, J)
		grand_proj_col <- t(L_mat) %*% grand_mean
		for(g in 1:num_classes){
			mean_proj <- t(L_mat) %*% class_means[[g]]
			diff_col <- mean_proj - grand_proj_col
			Sb_col <- Sb_col + class_sizes[g] * (t(diff_col) %*% diff_col)
			idx <- which(labels == classes[g])
			for(i in idx){
				diff_i <- t(L_mat) %*% x[,,i] - mean_proj
				Sw_col <- Sw_col + t(diff_i) %*% diff_i
			}
		}
		Sw_col_reg <- Sw_col + diag(1e-8, J)
		eig_col <- eigen(solve(Sw_col_reg) %*% Sb_col)
		R_mat <- Re(eig_col$vectors[, 1:c_ranks, drop=FALSE])

		#reconstruct estimate: est_k = L Z_k R^T where Z_k = L^T X_k R
		est_arr <- array(0, dim=modes)
		for(k in 1:K){
			z_k <- t(L_mat) %*% x[,,k] %*% R_mat
			est_arr[,,k] <- L_mat %*% z_k %*% t(R_mat)
		}

		if(CHECK_CONV(est_arr)){
			converged <- TRUE
			setTxtProgressBar(pb,max_iter)
		}else{
			curr_iter <- curr_iter + 1
		}
	}
	if(!converged){setTxtProgressBar(pb,max_iter)}
	close(pb)

	#compute projected features
	Z_list <- vector("list", K)
	est_arr <- array(0, dim=modes)
	for(k in 1:K){
		Z_list[[k]] <- t(L_mat) %*% x[,,k] %*% R_mat
		est_arr[,,k] <- L_mat %*% Z_list[[k]] %*% t(R_mat)
	}

	est <- as.tensor(est_arr)
	fnorm_resid <- fnorm_resid[fnorm_resid!=0]
	if(length(fnorm_resid)==0) fnorm_resid <- fnorm(est - tnsr)
	norm_percent <- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
	invisible(list(L=L_mat, R=R_mat, Z=Z_list, conv=converged, est=est,
		norm_percent=norm_percent, fnorm_resid=tail(fnorm_resid,1),
		all_resids=fnorm_resid))
}


###t-compress (Not Supported)
.t_compress <- function(tnsr,k){
	modes <- tnsr@modes
	n1 <- modes[1]
	n2 <- modes[2]
	n3 <- modes[3]
	#progress bar
	pb <- txtProgressBar(min=0,max=n3,style=3)
	#define ifft
	#.ifft <- function(x){suppressWarnings(as.numeric(fft(x,inverse=TRUE))/length(x))}
	#fft for each of the n1n2 vectors (of length n3) along mode 3
	fftz <- aperm(apply(tnsr@data,MARGIN=1:2,fft),c(2,3,1))
	#svd for each face (svdz is a list of the results)
	U_arr <- array(0,dim=c(n1,n1,n3))
	V_arr <- array(0,dim=c(n2,n2,n3))
	m <- min(n1,n2)		
	S_arr <- array(0,dim=c(n1,n2,n3))
	#Think of a way to avoid a loop in the beginning
	#Problem is that svd returns a list but ideally we want 3 arrays
	#Even with unlist this doesn't seem possible
	for (j in 1:n3){
		setTxtProgressBar(pb,j)
		decomp <- svd(fftz[,,j],nu=n1,nv=n2)
		U_arr[,,j] <- decomp$u
		V_arr[,,j] <- decomp$v
		S_arr[,,j] <- diag(decomp$d,nrow=n1,ncol=n2) #length is min(n1,n2)
	}	
	close(pb)
	#for each svd result, we want to apply ifft
	U <- as.tensor(aperm(apply(U_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	V <- as.tensor(aperm(apply(V_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	S <- as.tensor(aperm(apply(S_arr,MARGIN=1:2, .ifft),c(2,3,1)))
	
	est <- as.tensor(array(0,dim=modes))
	for (i in 1:k){
		est <- est + t_mult(t_mult(U[,i,,drop=FALSE],S[i,i,,drop=FALSE]),t(V[,i,,drop=FALSE]))
	}
	resid <- fnorm(est-tnsr)
	invisible(list(est=est, fnorm_resid = resid, norm_percent = (1-resid/fnorm(tnsr))*100))
}

###t-compress2 (Not Supported)
.t_compress2 <- function(tnsr,k1,k2){
	A = modeSum(tnsr,m=3,drop=TRUE)
	svdz <- svd(A@data,nu=k1,nv=k2)
	Util <- svdz$u
	Vtil <- svdz$v
	modes <- tnsr@modes
	n3 <- modes[3]
	core <- array(0,dim=c(k1,k2,n3))
	for(i in 1:n3){
	core[,,i]<-t(Util)%*%tnsr[,,i]@data%*%Vtil
	}
	est <- array(0,dim=modes)
	for(i in 1:k1){
		for (j in 1:k2){
			est = est + Util[,i] %o% Vtil[,j] %o% core[i,j,]
		}	
	}
	resid <- fnorm(tnsr - est)
	invisible(list(core = as.tensor(core), est=est, fnorm_resid = resid, norm_percent = (1-resid/fnorm(tnsr))*100))
}
