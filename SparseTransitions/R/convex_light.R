jl_pkg_setup <- function (...){
  julia <- JuliaCall::julia_setup(...)
  JuliaCall::julia_library("SparseTransitions")
}


alternating_minim_light <- function(rcount, max_iter=7L, lambda=0.1, weighted = TRUE,
                                     post_fit=TRUE, t_cutoff=10^(-4)){

  max_iter <- as.integer(max_iter)

  JuliaCall::julia_assign("rcount",rcount)
  JuliaCall::julia_assign("max_iter", max_iter)
  JuliaCall::julia_assign("lambda", lambda)
  JuliaCall::julia_assign("weighted", weighted)
  JuliaCall::julia_assign("post_fit", post_fit)
  JuliaCall::julia_assign("t_cutoff", t_cutoff)


  res_init <-  JuliaCall::julia_eval("fast_minim(rcount; max_iter=max_iter, λ=lambda, weighted=weighted, threshold=t_cutoff, refit=post_fit)")



  res <- with(res_init$res, list( obj_vals=obj_vals, transition_P = list(P1,P2), theta_m = theta_m))

  colnames(res$transition_P[[1]]) <- rownames(rcount)
  colnames(res$transition_P[[2]]) <- rownames(rcount)
  rownames(res$transition_P[[1]]) <- rownames(rcount)
  rownames(res$transition_P[[2]]) <- rownames(rcount)


  if (post_fit){

    post_res <- with(res_init$post_res, list( obj_vals=obj_vals, transition_P = list(P1,P2), theta_m = theta_m))

    colnames(post_res$transition_P[[1]]) <- rownames(rcount)
    colnames(post_res$transition_P[[2]]) <- rownames(rcount)
    rownames(post_res$transition_P[[1]]) <- rownames(rcount)
    rownames(post_res$transition_P[[2]]) <- rownames(rcount)


    res <- list(res=res, post_res=post_res)

  }

  class(res) <- "ForwardForwardFit"

  return (res)
}




predict_one_step <- function(fit, r2_count, desparsified=TRUE, weighted=TRUE, wts=NULL){

  r2 <- scale(r2_count,F,colSums(r2_count))

  if (is.null(wts)){
    wts <- 1/pmax(r2_count,2)
    wts <- wts/sum(wts)*length(wts)
  }

  if (desparsified) {
    p1 <- fit$post_res$transition_P[[1]]
    p2 <- fit$post_res$transition_P[[2]]
  } else {
    p1 <- fit$res$transition_P[[1]]
    p2 <- fit$res$transition_P[[2]]
  }


  # pass with p1
  r2_fit_1 <- r2[,1]
  r2_fit_2 <- p1 %*% r2[,1]
  r2_fit_imputed_1 <- p1 %*% r2[,2]
  r2_fit_3 <- p1 %*% r2_fit_imputed_1
  r2_fit_imputed_2 <- p1 %*% r2[,3]
  r2_fit_4 <- p1 %*% r2_fit_imputed_2

  # pass with p2
  r2_fit_5 <- p2 %*% r2[,4]
  r2_fit_imputed_3 <- p2 %*% r2[,5]
  r2_fit_6 <- p2 %*% r2_fit_imputed_3
  r2_fit_imputed_4 <- p2 %*% r2[,6]
  r2_fit_7 <- p2 %*% r2_fit_imputed_4

  preds <- cbind(r2_fit_1, r2_fit_2, r2_fit_3, r2_fit_4, r2_fit_5,r2_fit_6, r2_fit_7)
  if (weighted){
    mse <- sum( wts*((preds - r2)^2))
  } else {
    mse <- sum((preds-r2)^2)
  }
  mse
}


full_forward_pass <- function(fit, r2_count, desparsified=FALSE, weighted=TRUE, wts=NULL){
  r2 <- scale(r2_count,F,colSums(r2_count))

  if (is.null(wts)){
    wts <- 1/pmax(r2_count,2)
    wts <- wts/sum(wts)*length(wts)
  }


  if (desparsified) {
    p1 <- fit$post_res$transition_P[[1]]
    p2 <- fit$post_res$transition_P[[2]]
  } else {
    p1 <- fit$res$transition_P[[1]]
    p2 <- fit$res$transition_P[[2]]
  }

  # pass with p1
  r2_fit_1 <- r2[,1]
  r2_fit_2 <- p1 %*% r2_fit_1
  r2_fit_imputed_1 <- p1 %*% r2_fit_2
  r2_fit_3 <- p1 %*% r2_fit_imputed_1
  r2_fit_imputed_2 <- p1 %*% r2_fit_3
  r2_fit_4 <- p1 %*% r2_fit_imputed_2

  # pass with p2
  r2_fit_5 <- p2 %*% r2_fit_4
  r2_fit_imputed_3 <- p2 %*% r2_fit_5
  r2_fit_6 <- p2 %*% r2_fit_imputed_3
  r2_fit_imputed_4 <- p2 %*% r2_fit_6
  r2_fit_7 <- p2 %*% r2_fit_imputed_4

  preds <- cbind(r2_fit_1, r2_fit_2, r2_fit_3, r2_fit_4, r2_fit_5,r2_fit_6, r2_fit_7)

  if (weighted){
    mse <- sum( wts*((preds - r2)^2))
  } else {
    mse <- sum((preds-r2)^2)
  }
  mse
}
