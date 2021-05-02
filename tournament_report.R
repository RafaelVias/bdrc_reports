tournament.get_report <- function(t_obj,directory=NULL,type=1){
    h_dat <- t_obj$winner$data[[all.vars(t_obj$winner$formula)[2]]]
    q_dat <- t_obj$winner$data[[all.vars(t_obj$winner$formula)[1]]]
    posterior_list <- lapply(t_obj$contestants,function(m){
        post_dat <- spread_draws(m,'rating_curve','f','sigma_eps','c')
        post_dat[post_dat$h>=min(h_dat) & post_dat$h<=max(h_dat),]
    })
    main_plot_types <- c('rating_curve','residuals','sigma_eps','f')
    lim_list <- lapply(posterior_list,function(df){
                    log_h_minus_c <- if(is.null(t_obj$winner$run_info$c_param)) log(df$h-quantile(df$c,0.5)) else log(df$h-t_obj$winner$run_info$c_param)
                    sigma_eps_median <- quantile(df$sigma_eps,0.5)
                    data.frame(rating_curve_x_min=quantile(df$rating_curve,0.025),rating_curve_x_max=1.01*quantile(df$rating_curve,0.975),
                               rating_curve_y_min=min(df$h),rating_curve_y_max=1.01*max(df$h) - 0.01*min(df$h),
                               residuals_x_min=min(log_h_minus_c),residuals_x_max=max(log_h_minus_c),
                               residuals_y_min=1.1*(-1.96*sigma_eps_median),residuals_y_max=1.1*(1.96*sigma_eps_median),
                               sigma_eps_x_min=min(df$h),sigma_eps_x_max=max(df$h),
                               sigma_eps_y_min=0,sigma_eps_y_max=max(df$sigma_eps),
                               f_x_min=min(df$h),f_x_max=max(df$h),
                               f_y_min=min(df$f,1),f_y_max=max(df$f,3.5))
                })
    lim_dat <- do.call('rbind',lim_list)
    plot_list <- lapply(t_obj$contestants,function(m){
        pt_plot_list <- lapply(plot_types,function(pt){
            autoplot(m,type=pt) +
            scale_x_continuous(limits = c(min(lim_dat[[paste0(pt,'_x_min')]]),max(lim_dat[[paste0(pt,'_x_max')]]))) +
            scale_y_continuous(limits = c(min(lim_dat[[paste0(pt,'_y_min')]]),max(lim_dat[[paste0(pt,'_y_max')]])))
        })
        do.call('grid.arrange',model_plot_list)
    })
    mcmc_hist_list <- lapply(t_obj$contestants,function(m){
        plot(m,type='histogram',param=c('latent_parameters','hyperparameters'))
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
