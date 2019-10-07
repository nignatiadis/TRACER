run_lambda_split_grid <- function(rcount, nseeds =2, lambdas= c( 0, round(10^seq(-5,-3.5, length=9), 5))){
  res_list_outer <- list()

  for (seed_idx in 1:nseeds){
    print(paste0("seed:",seed_idx))
    #set.seed(seeds[seed_idx])
    spl <- split_into_two(rcount, p=0.5)

    # Remember to update

    res_list <- list()

    for (i in 1:length(lambdas)){
      print(paste0("lambda iter:", i, " with lambda:", lambdas[i]))
      total_cts <- spl$rcount_left + spl$rcount_right
      wts <- 1/pmax(total_cts,2)
      wts <- wts/sum(wts)*length(wts)

      res_l <- alternating_minim_light(spl$rcount_left, lambda=lambdas[i], post_fit = TRUE, max_iter=7)
      res_r <- alternating_minim_light(spl$rcount_right,lambda=lambdas[i], post_fit = TRUE, max_iter=7)

      res_df <-  rbind(data.frame(split=1, desparsified=FALSE, predict_mode="one step",  error= predict_one_step(res_l, spl$rcount_right, wts= wts,desparsified=FALSE)),
                       data.frame(split=1, desparsified=TRUE, predict_mode="one step", error= predict_one_step(res_l, spl$rcount_right, wts= wts, desparsified=TRUE)),
                       data.frame(split=2, desparsified=FALSE, predict_mode="one step", error= predict_one_step(res_r, spl$rcount_left, wts= wts,desparsified=FALSE)),
                       data.frame(split=2, desparsified=TRUE, predict_mode="one step", error= predict_one_step(res_r, spl$rcount_left, wts= wts, desparsified=TRUE)),
                       data.frame(split=1, desparsified=FALSE, predict_mode="full forward",  error=full_forward_pass(res_l, spl$rcount_right, wts= wts,desparsified=FALSE)),
                       data.frame(split=1, desparsified=TRUE, predict_mode="full forward",  error=full_forward_pass(res_l, spl$rcount_right, wts= wts, desparsified=TRUE)),
                       data.frame(split=2, desparsified=FALSE, predict_mode="full forward", error=full_forward_pass(res_r, spl$rcount_left, wts= wts,desparsified=FALSE)),
                       data.frame(split=2, desparsified=TRUE, predict_mode="full forward", error=full_forward_pass(res_r, spl$rcount_left, wts= wts, desparsified=TRUE)))

      res_df <- dplyr::mutate(res_df, seed_idx=seed_idx, lambda=lambdas[i])
      res_res <- list(res_l=res_l,res_r=res_r,res_df=res_df, spl=spl) # just to have it available
      res_list[[i]] <- res_res
    }
    res_list_outer[[seed_idx]] <- res_list
  }
  res_list_outer
}

outer_res_list_to_df <- function(res_list_outer){
  all_df_list <- list()
  for (seed_idx in 1:length(res_list_outer)){
    all_df_list[[seed_idx]] <- list()
    for (i in 1:length(res_list_outer[[seed_idx]])){
      all_df_list[[seed_idx]][[i]] <- res_list_outer[[seed_idx]][[i]]$res_df
    }
  }
  all_df <- bind_rows(unlist(all_df_list, recursive=FALSE)) %>% mutate(split=as.factor(split)) %>%
    mutate(split_seed = paste0("seed:",seed_idx," & split:", split))
  all_df
}
