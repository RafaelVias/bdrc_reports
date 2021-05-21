tournament.get_report <- function(t_obj,directory=NULL,type=1){
    h_dat <- t_obj$winner$data[[all.vars(t_obj$winner$formula)[2]]]
    q_dat <- t_obj$winner$data[[all.vars(t_obj$winner$formula)[1]]]
    posterior_list <- lapply(t_obj$contestants,function(m){
        post_dat <- spread_draws(m,'rating_curve','f','sigma_eps','c')
        post_dat[post_dat$h>=min(h_dat) & post_dat$h<=max(h_dat),]
    })
    max_res <- lapply(names(posterior_list), function(x) {
                  mod_obj <- t_obj$contestants[[x]]
                  c_param <- if(is.null(mod_obj$run_info$c_param)) median(mod_obj$c_posterior) else mod_obj$run_info$c_param
                  resid_dat <- merge(mod_obj$rating_curve[,c('h','median')],mod_obj$data,by.x='h',by.y=all.vars(mod_obj$formula)[2])
                  resid_dat[,'log(h-c_hat)'] <- log(resid_dat$h-c_param)
                  max(log(resid_dat$Q)-log(resid_dat$median))
               })
    max_res <- max(abs(do.call('rbind',max_res)))
    lim_list <- lapply(posterior_list,function(df){
                    sigma_eps_median <- quantile(df$sigma_eps,0.5)
                    data.frame(rating_curve_x_min=quantile(df$rating_curve,0.025),rating_curve_x_max=1.01*max(quantile(df[df$h==max(df$h),]$rating_curve,0.975),max(q_dat)),
                               rating_curve_y_min=min(df$h),rating_curve_y_max=1.01*max(df$h)-0.01*min(df$h),
                               residuals_y_min=1.1*min((-1.96*sigma_eps_median),-max_res),residuals_y_max=1.1*max((1.96*sigma_eps_median),max_res),
                               sigma_eps_x_min=min(df$h),sigma_eps_x_max=max(df$h),
                               sigma_eps_y_min=0,sigma_eps_y_max=max(df$sigma_eps),
                               f_x_min=min(df$h),f_x_max=max(df$h),
                               f_y_min=min(df$f,1),f_y_max=max(df$f,3.5))
                 })
    lim_dat <- do.call('rbind',lim_list)
    main_plot_types <- c('rating_curve','residuals','sigma_eps','f')
    main_plot_list <- lapply(t_obj$contestants,function(m){
        pt_plot_list <- lapply(main_plot_types,function(pt){
            if(pt=="residuals") {
              autoplot(m,type=pt) + 
                scale_y_continuous(limits = c(min(lim_dat[[paste0(pt,'_y_min')]]),max(lim_dat[[paste0(pt,'_y_max')]])))
            }else{
              autoplot(m,type=pt) +
                scale_x_continuous(limits = c(min(lim_dat[[paste0(pt,'_x_min')]]),max(lim_dat[[paste0(pt,'_x_max')]]))) +
                scale_y_continuous(limits = c(min(lim_dat[[paste0(pt,'_y_min')]]),max(lim_dat[[paste0(pt,'_y_max')]])))
            }
        })
        do.call('grid.arrange',pt_plot_list)
    })
    mcmc_hist_list <- lapply(t_obj$contestants,function(m){
        autoplot(m,type='histogram',param=c('latent_parameters','hyperparameters'),transformed = T)
    })
    MCMC_table <- lapply(t_obj$contestants,function(m){
                    data.frame(nr_iter=m$run_info$nr_iter,
                               nr_chains=m$run_info$num_chains,
                               burnin=m$run_info$burnin,
                               thin=m$run_info$thin,
                               nr_eff_param=format(m$num_effective_param,digits=2),
                               acceptance_rate=format(m$acceptance_rate,digits=2))
                  })
    MCMC_table <- t(do.call('rbind',MCMC_table))


}
