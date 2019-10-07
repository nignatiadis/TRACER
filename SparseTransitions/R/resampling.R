split_into_two <- function(rcount, prop=0.5){
  rcount_left <- matrix(0, nrow(rcount), ncol(rcount))
  rcount_right <- matrix(0, nrow(rcount), ncol(rcount))

  for (i in 1:nrow(rcount)){
    for (j in 1:ncol(rcount)){
      rcount_left[i,j] <- rbinom(1, size=rcount[i,j], prob=prop)
      rcount_right[i,j] <- rcount[i,j] - rcount_left[i,j]
    }
  }

  rownames(rcount_left) <- rownames(rcount)
  rownames(rcount_right) <- rownames(rcount)

  colnames(rcount_left) <- colnames(rcount)
  colnames(rcount_right) <- colnames(rcount)

  rprop_left =scale(rcount_left,F,colSums(rcount_left))  #make colsums eq to 1
  rprop_right =scale(rcount_right,F,colSums(rcount_right))  #make colsums eq to 1

  list( rcount_left = rcount_left,
        rprop_left = rprop_left,
        rcount_right = rcount_right,
        rprop_right = rprop_right)
}

poisson_boot <- function(rcount){
  rcount_pois <- matrix(0, nrow(rcount), ncol(rcount))

  for (i in 1:nrow(rcount_pois)){
    for (j in 1:ncol(rcount_pois)){
      rcount_pois[i,j] <- rpois(1, rcount[i,j])
    }
  }

  rownames(rcount_pois) <- rownames(rcount)
  colnames(rcount_pois) <- colnames(rcount)

  rcount_pois
}


column_multinomial_boot <- function(rcount){
  rcount_mult <- matrix(0, nrow(rcount), ncol(rcount))

  n_time <-  ncol(rcount)
  n_clusters <- nrow(rcount)

  for (i in 1:n_time) {
    Ns <- sum(rcount[,i])
    pis <- rcount[,i]/Ns
    rcount_mult[,i] = rmultinom(1, Ns, pis)
  }
  rownames(rcount_mult) <- rownames(rcount)
  colnames(rcount_mult) <- colnames(rcount)

  rcount_mult
}

simulate_cell_proportions <- function(theta_matrix, Ns){
  n_time <-  ncol(theta_matrix)
  n_clusters <- nrow(theta_matrix)
  rs <- matrix(NA, n_clusters, n_time)

  # let us be lazy and use a loop..
  for (i in 1:n_time) {
    rs[,i] = rmultinom(1, Ns[i], theta_matrix[,i])/Ns[i]
  }

  rs
}



boot_run <- function(rcount, seed_idx, sampling="column multinomial", nseeds_cv=5, ...){
  set.seed(seed_idx)
  #bootstrap
  if (sampling=="column multinomial"){
    count_pois <- column_multinomial_boot(rcount)
  } else if (sampling=="poisson"){
    count_pois <- poisson_boot(rcount)
  }

  cv_test_split <- split_into_two(count_pois, p=0.5)

  lambda_grid <- c(0,10^seq(-6,-2, length=40))
  cv_res_list <- run_lambda_split_grid(cv_test_split$rcount_left,
                                       nseeds=nseeds_cv,
                                       lambdas=lambda_grid)
  cv_res_df <- outer_res_list_to_df(cv_res_list)


  lambda_cv_df <- filter(cv_res_df,  predict_mode=="one step", desparsified==TRUE) %>%
    group_by(lambda) %>%
    summarize(mean_err = mean(error), sd_error=sd(error)/sqrt(2)) %>%
    arrange(mean_err)

  se_rule <- lambda_cv_df$mean_err[1] + lambda_cv_df$sd_error[1]
  lambda_cv_df_filt <- filter(lambda_cv_df, mean_err <= se_rule)

  lambda_cv <- lambda_cv_df$lambda[1]
  lambda_1se <- max(lambda_cv_df_filt$lambda)


  train_res <- alternating_minim_light(cv_test_split$rcount_right, lambda=lambda_1se, post_fit = TRUE, max_iter=7)

  list(count_pois=count_pois, cv_test_split=cv_test_split,
       train_res=train_res, cv_res_list=cv_res_list,
       cv_res_df=cv_res_df, lambda_cv=lambda_cv, lambda_1se=lambda_1se)
}
