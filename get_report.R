



rhat_ready_draws <- function(draws){
  rhat_ready <- split(draws$value,draws$chain,drop=TRUE)
  rhat_ready <- do.call('cbind',rhat_ready)
  return(rhat_ready)
}




make_rhat_plots <- function(m_obj){
  draw_list <- lapply(m_obj,function(m){
    lapply(get_hist_params(class(m)),function(x){
      draws  <- gather_draws(m,x,transformed = T)
    })
  })
  rhat_plot_list <- lapply(names(m_obj),function(m){
    model_draws <- lapply(1:length(draw_list[[m]]),function(x){
      draws  <- rhat_ready_draws(draw_list[[m]][[x]])
      param <- get_hist_params(m,c_param = c_param)
      data.frame('iterations'=4:nrow(draws),'parameters'=rep(param[[x]],length(4:nrow(draws))),'Rhat'=sapply(4:nrow(draws), function(r) R_hat(draws[1:r,])))
    })
    model_draws <- do.call('rbind',model_draws)
    ggplot(data=model_draws, aes(x=iterations,y=Rhat,color=parameters))+
      geom_hline(yintercept = 1.1, linetype='dashed')+
      geom_line()+
      scale_y_continuous(expand = c(0,0),limits = c(1,2))+
      scale_x_continuous(expand = c(0,0))+
      xlab('')+
      ggtitle(m)+
      theme_bdrc() 
  })
  p <- grid.arrange(grobs=rhat_plot_list,nrow=3,ncol=2)
  return(p)
}




R_hat <- function(four_chains){ 
  if(nrow(four_chains)%%2!=0) four_chains <- four_chains[-1,]
  psi <- lapply(1:4,function(x){
    chain <- data.frame('chain'=rep(1:2,each=nrow(four_chains)/2),'value'=four_chains[,x])
    chain <- split(chain$value, chain$chain , drop=TRUE)
    chain <- do.call('cbind',chain)
  })
  psi <- do.call('cbind',psi)
  m <- dim(psi)[[2]]
  n <- dim(psi)[[1]]
  psi_bar_dj <- sapply(1:m, function(x) mean(psi[,x]))
  psi_bar_dd <- mean(psi_bar_dj)
  B <- ((n)/(m-1))*sum((psi_bar_dj-psi_bar_dd)^2)
  s_sqrd <- sapply(1:8, function(x) ((1)/(n-1))*sum((psi[,x]-psi_bar_dj[x])^2))
  W <- mean(s_sqrd)
  var_hat_p <- ((n-1)/(n))*W+((1)/(n))*B
  rhat <- sqrt(var_hat_p/W)
  return(rhat)
}




rhat_vector <- function(m){
  c_param <- if(is.null(m$run_info$c_param)) F else T
  params <- get_hist_params(class(m),c_param = c_param)
  rhat_vec <- sapply(params,function(x){
    long_chain <- gather_draws(m,x,transformed = T)
    four_chains <- rhat_ready_draws(long_chain)
    R_hat(four_chains)
  })
  return(rhat_vec)
}




predict_matrix <- function(x){
  c_param <- if(is.null(x$run_info$c_param)) median(x$c_posterior) else x$run_info$c_param
  c_param <- ceiling(c_param*100)/100
  grid_max <- x$run_info$h_max
  p_dat <- predict(x,newdata=seq(c_param,grid_max,by=0.01))[c('h','median')]
  p_dat$decimal <- floor(p_dat$h*10)/10
  first_decimal <- length(p_dat$decimal[p_dat$decimal==p_dat$decimal[1]])
  if(first_decimal!=10) {
    n <- 10-first_decimal
    top_rows <- data.frame(h=sapply(n:1,function(x) p_dat$h[1]-0.01*x),median=rep(0,n),decimal=rep(p_dat$decimal[1],n))
    p_dat <- rbind(top_rows,p_dat)
  }
  last_decimal <- length(p_dat$decimal[p_dat$decimal==p_dat$decimal[length(p_dat$decimal)]])
  if(last_decimal!=10){
    m <- 10-last_decimal
    bot_rows <- data.frame(h=sapply(1:m,function(x) p_dat$h[length(p_dat$h)]+0.01*x),median=rep(NA,m),decimal=rep(p_dat$decimal[length(p_dat$decimal)],m))
    p_dat <- rbind(p_dat,bot_rows)
  }
  p_mat <- lapply(unique(p_dat$decimal),function(d) p_dat$median[p_dat$decimal==d]) %>% do.call('rbind',.)
  rownames(p_mat) <- unique(p_dat$decimal)
  colnames(p_mat) <- seq(0,0.09,by=0.01)
  p_mat <- round(p_mat,digits=3)
  p_mat <- tableGrob(p_mat,theme=ttheme_minimal(base_family = "Times",
                                                core=list(bg_params = list(fill = blues9[1:2], col=NA),fg_params=list(fontface=3)),
                                                colhead=list(fg_params=list(col="black",fontface=2L)),
                                                rowhead=list(fg_params=list(col="black",fontface=2L))))
  return(p_mat)
}




get_hist_params <- function(m_class,c_param=FALSE){
  param_sets_T <- list('gplm'=c(1,2,3,5,6,7,8,9,10,11,12,13),
                       'gplm0'=c(1,2,3,4,5,6),
                       'plm'=c(1,2,3,7,8,9,10,11,12,13),
                       'plm0'=c(1,2,3,4))
  parameter_vec <- c('a','b','c','sigma_eps','sigma_beta',
                     'phi_beta','sigma_eta','eta_1','eta_2',
                     'eta_3','eta_4','eta_5','eta_6')
  if(c_param==T){
    param_sets_T <- lapply(param_sets_T,function(x) x[-3])
  } 
  params <- parameter_vec[ param_sets_T[[m_class]] ]
  return(params)
}




get_report.gplm <- function(x,directory=NULL,type=1,...){
  get_report(x,directory=directory,type=type,...)
}




get_report.gplm0 <- function(x,directory=NULL,type=1,...){
  get_report(x,directory=directory,type=type,...)
}




get_report.plm <- function(x,directory=NULL,type=1,...){
  get_report(x,directory=directory,type=type,...)
}




get_report.plm0 <- function(x,directory=NULL,type=1,...){
  get_report(x,directory=directory,type=type,...)
}


get_report.tournament <- function(x,directory=NULL,type=1){
  get_report(x,directory=directory,type=type,...)
}




get_report <- function(...,directory=NULL,type=1,rhat_plots=FALSE){
  args <- list(...)
  error_msg1 <- 'Please provide either a single tournament object or a single model object of types gplm, gplm0, plm and plm0.'
  legal_types <- c('gplm','gplm0','plm','plm0','tournament')
  if(!(type %in% c(1,2))){
    stop('Please input an integer value of 1 or 2 to indicate which type of report is to be produce.')
  }else{
    args_class <- unlist(lapply(args,class))
    if(!(args_class%in%legal_types)){
      stop(error_msg1)
    }
  }
  if(!is.null(directory)){
    if(!dir.exists(directory)){
      stop('Please provide a valid directory to store the .pdf report.')
    }
  }else{
    directory <- getwd() 
  }
  if(args_class=='tournament'){
    if(type==1){
      m_obj <- list(args[[1]]$winner)
      names(m_obj) <- class(args[[1]]$winner)
    }else{
      m_obj <- args[[1]]$contestants
      t_obj <- args[[1]]
    }
  }else{
    if(type==2){
      stop('It is only possible to produce a type 1 report for a single model object of type gplm, gplm0, plm or plm0.')
    }
    names(args) <- args_class
    m_obj <- args
  }
  h_dat <- m_obj[[1]]$data[[all.vars(m_obj[[1]]$formula)[2]]]
  q_dat <- m_obj[[1]]$data[[all.vars(m_obj[[1]]$formula)[1]]]
  posterior_list <- lapply(m_obj,function(m){
    post_dat <- spread_draws(m,'rating_curve','f','sigma_eps','c')          # ma nota c her ??
    post_dat[post_dat$h>=min(h_dat) & post_dat$h<=max(h_dat),]
  })
  max_res <- lapply(names(posterior_list), function(x) {
    c_param <- if(is.null(m_obj[[x]]$run_info$c_param)) median(m_obj[[x]]$c_posterior) else m_obj[[x]]$run_info$c_param
    resid_dat <- merge(m_obj[[x]]$rating_curve[,c('h','median')],m_obj[[x]]$data,by.x='h',by.y=all.vars(m_obj[[x]]$formula)[2])
    resid_dat[,'log(h-c_hat)'] <- log(resid_dat$h-c_param)
    max(log(resid_dat$Q)-log(resid_dat$median))
  })
  max_res <- max(abs(do.call('rbind',max_res)))
  lim_list <- lapply(posterior_list,function(df){
    sigma_eps_median <- quantile(df$sigma_eps,0.5)
    data.frame(rating_curve_x_min=quantile(df$rating_curve,0.025),rating_curve_x_max=1.01*max(quantile(df[df$h==max(df$h),]$rating_curve,0.975),max(q_dat)),   
               rating_curve_y_min=min(df$h),rating_curve_y_max=1.01*max(df$h)-0.01*min(df$h),
               residuals_y_min=1.1*min((-1.96*sigma_eps_median),-max_res),residuals_y_max=1.1*max((1.96*sigma_eps_median),max_res),  
               residuals_x_min=NA,residuals_x_max=NA,
               sigma_eps_x_min=min(df$h),sigma_eps_x_max=max(df$h),
               sigma_eps_y_min=0,sigma_eps_y_max=max(df$sigma_eps),
               f_x_min=min(df$h),f_x_max=max(df$h),
               f_y_min=min(df$f,1),f_y_max=max(df$f,3.5))
  })
  lim_dat <- do.call('rbind',lim_list)
  main_plot_types <- c('rating_curve','residuals','sigma_eps','f')  
  main_plot_list <- lapply(m_obj,function(m){
    pt_plot_list <- lapply(main_plot_types,function(pt){
      autoplot(m,type=pt) +
        scale_x_continuous(limits = c(min(lim_dat[[paste0(pt,'_x_min')]]),max(lim_dat[[paste0(pt,'_x_max')]]))) +
        scale_y_continuous(limits = c(min(lim_dat[[paste0(pt,'_y_min')]]),max(lim_dat[[paste0(pt,'_y_max')]])))
    })
    do.call('grid.arrange',pt_plot_list)
  })
  main_table_list <- lapply(m_obj,function(m){
    table <- rbind(m$param_summary,m$Deviance_summary)
    row.names(table) <- c(sapply(1:(nrow(table)-1), function(x) get_param_expression(row.names(table)[[x]])[[1]]),"Deviance")
    rhat_col <- c(rhat_vector(m),NA)   
    table <- cbind(table,'R_hat'=rhat_col)
    table <- format(round(table,digits=3),nsmall=3)
    table[nrow(table),4] <- ''
    table <- tableGrob(table,theme=ttheme_minimal(base_family="Times",rowhead=list(fg_params=list(parse=TRUE))))
  })
  if(type==1){ 
    p_mat <- predict_matrix(m_obj[[1]])
  }else{  
    mcmc_hist_list <- lapply(m_obj,function(m){
      params <- get_hist_params(class(m))
      hist_plot_list <- lapply(1:length(params), function(j){ 
        autoplot(m,type='histogram',param=params[j],transformed = T)
      })
    })
    MCMC_table <- lapply(m_obj,function(m){
      data.frame(nr_iter=m$run_info$nr_iter,
                 nr_chains=m$run_info$num_chains,
                 burnin=m$run_info$burnin,
                 thin=m$run_info$thin,
                 nr_eff_param=format(m$num_effective_param,digits=2),
                 acceptance_rate=format(m$acceptance_rate,digits=2))
    })
    MCMC_table <- t(do.call('rbind',MCMC_table))
    MCMC_table <- tableGrob(MCMC_table,theme=ttheme_minimal(base_family="Times"))
    tour_table <- t_obj$summary[c("round","game","model","DIC","P","winner")]
    tour_table[c('DIC','P')] <- round(tour_table[c('DIC','P')],digits=2)
    tour_table <- tableGrob(tour_table,theme=ttheme_minimal(base_family="Times"),rows=NULL)
    dev_boxplot <- autoplot(t_obj,type='deviance')
    if(rhat_plots) rhat_plot <- make_rhat_plots(m_obj)           ### <- MUNA !!! ###
  }
  pdf(file=paste0(directory,'/report.pdf'),paper='a4',width=8,height=11)
  if(type==1){
      grid.arrange(main_plot_list[[1]],main_table_list[[1]],nrow=2,as.table=TRUE,heights=c(5,3),top=textGrob(class(m_obj[[1]]),gp=gpar(fontfamily="Times",fontsize=22,facetype='bold')))
      grid.arrange(p_mat,nrow=1,as.table=TRUE,heights=c(1),top=textGrob(paste0('Rating curve predictions for ',class(m_obj[[1]])),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
  }else{
    lapply(names(posterior_list),function(x){ 
      grid.arrange(main_plot_list[[x]],main_table_list[[x]],nrow=2,as.table=TRUE,heights=c(5,3),top=textGrob(x,gp=gpar(fontfamily="Times",fontsize=22,facetype='bold')))
    })
    ### MUNA !!! ###
    if(rhat_plots){
      grid.arrange(arrangeGrob(MCMC_table,tour_table,nrow=1),arrangeGrob(dev_boxplot,ncol=2),nrow=2,as.table=TRUE,heights=c(1,1),top=textGrob('Model comparison and MCMC diagnostics',gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
      grid.arrange(rhat_plot,as.table=TRUE,top=textGrob('',gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
    }else{
      grid.arrange(tour_table,arrangeGrob(dev_boxplot,MCMC_table,ncol=2),nrow=2,as.table=TRUE,heights=c(1,2),top=textGrob('Model comparison',gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
    }
    lapply(names(posterior_list), function(x) {
      grid.arrange(arrangeGrob(grobs=mcmc_hist_list[[x]],nrow=4,ncol=3),top=textGrob(paste0('Estimated parameters of ',x),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
    })
  }
  dev.off()
}



