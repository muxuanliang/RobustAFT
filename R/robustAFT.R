robustAFT <- function(surv_obj, covariate, loss = 'square', penalty = TRUE, group_index = NULL, lambda_lasso = 60, lambda_group = NULL, standardize = TRUE,...){
  # check
  penalty_type <- 'lasso'
  if (!is.null(group_index)){
    stopifnot(NROW(covariate)==length(group_index))
    stopifnot(!is.null(lambda_group))
    penalty_type <- 'glasso'
  }
  if (penalty == FALSE){
    penalty_type <- 'refit'
  }
  stopifnot(loss %in% c("square", "absolute", "huber", "tukey"))

  fit <- switch(penalty_type,
                lasso=bje_ly(covariate, surv_obj, method=loss, lambda=lambda_lasso, standardize=standardize, ...),
                glasso=bje_sgl(covariate, surv_obj, group_index, lambda_lasso , lambda_group, method=loss, standardize=standardize, ...),
                refit=bje_refit(covariate, surv_obj, method=loss))
  fit
}
