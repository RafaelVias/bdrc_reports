



extract_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




smaller_legend <- function(p, pointSize = 1, textSize = 11, spaceLegend = 0.8) {
  p +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize+5), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}



disjoint_four_chains <- function(draws){
  if(nrow(draws)%%2!=0){
    draws <- draws[-1,]
  } 
  draws <- split(draws$value,draws$chain,drop=TRUE)
  do.call('cbind',draws)
}



split_chains <- function(four_chains){
  if(nrow(four_chains)%%2!=0){
    four_chains <- four_chains[-1,]
  } 
  psi <- lapply(1:ncol(four_chains),function(x){
    chain <- data.frame('chain'=rep(1:2,each=nrow(four_chains)/2),'value'=four_chains[,x])
    chain <- split(chain$value,chain$chain,drop=TRUE)
    chain <- do.call('cbind',chain)
  })
  psi <- do.call('cbind',psi)
  return(psi)
}



param_draws_list <- function(m,param,split_chains=TRUE){
  draws_list <- lapply(param,function(x){
    draws <- gather_draws(m,x,transformed=T)
    draws <- disjoint_four_chains(draws)
    if(split_chains==TRUE){
      split_chains(draws)
    }else{
      draws
    }
  })
  names(draws_list) <- param
  return(draws_list)
}




chain_statistics <- function(psi){
  m <- ncol(psi)
  n <- nrow(psi)
  psi_bar_dj <- sapply(1:m, function(x) mean(psi[,x]))
  psi_bar_dd <- mean(psi_bar_dj)
  B <- ((n)/(m-1))*sum((psi_bar_dj-psi_bar_dd)^2)
  s_sqrd <- sapply(1:m, function(x) (n-1)^(-1)*sum((psi[,x]-psi_bar_dj[x])^2))
  W <- mean(s_sqrd)
  var_hat <- ((n-1)/(n))*W+((1)/(n))*B
  return(list('W'=W,'var_hat'=var_hat))
}




variogram <- function(t,psi){
  m <- ncol(psi)
  n <- nrow(psi)
  double_sum <- sum(sapply(1:m,function(j){ sum((psi[(t+1):(n),j]-psi[(1):(n-t),j])^2) }))    
  V_t <- (m*(n-t))^(-1)*double_sum
  return(V_t)
}




autocorrelation <- function(t,psi){
  V_t <- variogram(t,psi)
  var_hat <- chain_statistics(psi)$var_hat
  rho_hat_t <- 1-(V_t/(2*var_hat))
  return(rho_hat_t)
}




n_eff_samples <- function(psi){
  m <- ncol(psi)
  n <- nrow(psi)
  rho_vec <- autocorrelation(t=1,psi)
  stop <- FALSE
  t <- 2
  while(stop==FALSE & t<nrow(psi)){
    rho_t <- autocorrelation(t,psi)
    rho_vec <- append(rho_vec,rho_t)
    if((t-2)%%2!=0){
      if( rho_vec[t-1]+rho_vec[t]<0 ){
        stop <- TRUE
      }
    }
    t <- t+1
  }
  rho_vec <- rho_vec[1:(length(rho_vec)-2)]
  n_eff <- (n*m)/(1+2*sum(rho_vec))
  n_eff <- round(n_eff)
  return(n_eff)
}




R_hat <- function(four_chains){ 
  psi <- split_chains(four_chains)
  numbers <- chain_statistics(psi)
  W <- numbers$W
  var_hat <- numbers$var_hat
  rhat <- sqrt(var_hat/W)
  return(rhat)
}




rhat_vector <- function(m,param){
  draws_list <- param_draws_list(m,param,split_chains=FALSE)
  rhat_vec <- sapply(param,function(x){
    R_hat(draws_list[[x]])
  })
  return(rhat_vec)
}




neff_vector <- function(m,param){
  draws_list <- param_draws_list(m,param)
  neff_vec <- sapply(param,function(x){
    n_eff_samples(draws_list[[x]])
  })
  return(neff_vec)
}




get_rhat_dat <- function(m,param,smoothness=40){
  thin <- m$run_info$thin
  burnin <- m$run_info$burnin
  draws_list <- param_draws_list(m,param,split_chains=FALSE)
  rhat_dat <- lapply(param,function(x){    
    draws  <- draws_list[[x]]
    real_iter <- seq(4*thin+burnin,((nrow(draws))*thin)+burnin,by=smoothness*thin)    
    by_n <- seq(4,nrow(draws),by=smoothness)
    data.frame('iterations'=real_iter,'parameters'=rep(x,length(by_n)),'Rhat'=sapply(by_n, function(r) R_hat(draws[1:r,]))) 
  })
  rhat_dat <- do.call('rbind',rhat_dat)
  return(rhat_dat)
}




get_autocorrelation_dat <- function(m,param,lags=30){
  draws_list <- param_draws_list(m,param)
  auto_dat <- lapply(param,function(x){
    draws  <- draws_list[[x]]
    data.frame('lags'=0:lags,'parameters'=rep(x,lags+1),'autocorrelation'=sapply(0:lags, function(t) autocorrelation(t,draws)))
  })
  auto_dat <- do.call('rbind',auto_dat)
  return(auto_dat)
}




get_conv_diagnostics_plots <- function(m_obj){
  m_class <- class(m_obj)
  model_type <- c('gplm','gplm0','plm','plm0')
  if(m_class%in%model_type){
    m_obj <- list(m_obj)
    names(m_obj) <- m_class
  }
  plot_list <- lapply(m_obj,function(m){
    param <- get_param_names(class(m),m$run_info$c_param)
    p1 <- autoplot(m,type='rhat',param=param)
    p2 <- autoplot(m,type='autocorrelation',param=param)
    p1 <- smaller_legend(p1)
    p2 <- smaller_legend(p2)
    my_legend <- extract_legend(p1)
    grid.arrange(arrangeGrob(p1+theme(legend.position="none"),p2+theme(legend.position="none"),nrow=1),
                 my_legend,ncol=2,widths=c(4,1),top=textGrob(class(m),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
  })
  names(plot_list) <- names(m_obj)
  if(m_class=='list'){
    p <- grid.arrange(grobs=plot_list,nrow=4)
  }else{
    p <- plot_list[[m_class]]
  }
  return(p)
}




get_residuals_dat <- function(m){
  resid_dat <- merge(m$rating_curve[,c('h','median')],m$data,by.x='h',by.y=all.vars(m$formula)[2],)
  if('sigma_eps_summary' %in% names(m)){
    resid_dat <- merge(resid_dat,m$sigma_eps_summary[,c('h','median')],by = 'h')
    names(resid_dat) <- c('h','median','Q','sigma_eps')
  }else{
    resid_dat$sigma_eps <- m$param_summary['sigma_eps','median']
  }
  c_hat <- if(is.null(m$run_info$c_param)) median(m$c_posterior) else m$run_info$c_param
  resid_dat[,'log(h-c_hat)'] <- log(resid_dat$h-c_hat)
  resid_dat$r_median <- log(resid_dat$Q)-log(resid_dat$median)
  resid_dat$r_lower <- -1.96*resid_dat$sigma_eps
  resid_dat$r_upper <- 1.96*resid_dat$sigma_eps
  return(resid_dat)
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
  p_mat <- lapply(unique(p_dat$decimal),function(d) p_dat$median[p_dat$decimal==d])
  p_mat <- do.call('rbind',p_mat)
  rownames(p_mat) <- unique(p_dat$decimal)
  colnames(p_mat) <- seq(0,0.09,by=0.01)
  p_mat <- round(p_mat,digits=3)
  p_mat <- tableGrob(p_mat,theme=ttheme_minimal(base_family = "Times",
                                                core=list(bg_params = list(fill = blues9[1:2], col=NA),fg_params=list(fontface=3)),
                                                colhead=list(fg_params=list(col="black",fontface=2L)),
                                                rowhead=list(fg_params=list(col="black",fontface=2L))))
  return(p_mat)
}




get_report.gplm <- function(x,directory=NULL,report_title=NULL,type=1,...){
  get_report(x,directory=directory,report_title=report_title,type=type,...)
}




get_report.gplm0 <- function(x,directory=NULL,report_title=NULL,type=1,...){
  get_report(x,directory=directory,report_title=report_title,type=type,...)
}




get_report.plm <- function(x,directory=NULL,report_title=NULL,type=1,...){
  get_report(x,directory=directory,report_title=report_title,type=type,...)
}




get_report.plm0 <- function(x,directory=NULL,report_title=NULL,type=1,...){
  get_report(x,directory=directory,report_title=report_title,type=type,...)
}




get_report.tournament <- function(x,directory=NULL,report_title=NULL,type=1){
  get_report(x,directory=directory,report_title=report_title,type=type,...)
}




get_report <- function(...,directory=NULL,report_title=NULL,type=1){
  args <- list(...)
  error_msg1 <- 'Please provide either a single tournament object or a single model object of types gplm, gplm0, plm and plm0.'
  error_msg2 <- 'Please provide one report title that is type character.'
  legal_types <- c('gplm','gplm0','plm','plm0','tournament')
  if(!(type %in% c(1,2))){
    stop('Please input an integer value of 1 or 2 to indicate which type of report is to be produce.')
  }else{
    args_class <- unlist(lapply(args,class))
    if(!(args_class%in%legal_types)){
      stop(error_msg1)
    }
  }
  if(is.null(directory)){
    directory <- getwd() 
  }else{
    if(!dir.exists(directory)){
      stop('Please provide a valid directory to store the .pdf report.')
    }
  }
  if(is.null(report_title)){
    report_title <- '/report.pdf'
  }else{
    if(class(report_title)!='character'){
      stop(error_msg2)
    }else if(length(report_title)!=1){
      stop(error_msg2)
    }else{
      report_title <- paste0('/',report_title,'.pdf')
    }
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
    if(is.null(m$run_info$c_param)){
      post_dat <- spread_draws(m,'rating_curve','f','sigma_eps','c') 
    }else{
      post_dat <- spread_draws(m,'rating_curve','f','sigma_eps') 
      post_dat$c <- m$run_info$c_param
    }
    post_dat[post_dat$h>=min(h_dat) & post_dat$h<=max(h_dat),]
  })
  res_dat <- sapply(names(m_obj),function(x){
    max(abs(get_residuals_dat(m_obj[[x]])[,c('r_median','r_lower','r_upper')]))
  })
  max_res <- max(res_dat)
  lim_list <- lapply(posterior_list,function(df){
    data.frame(rating_curve_x_min=quantile(df$rating_curve,0.025),rating_curve_x_max=1.01*max(quantile(df[df$h==max(df$h),]$rating_curve,0.975),max(q_dat)),   
               rating_curve_y_min=min(df$h),rating_curve_y_max=1.01*max(df$h)-0.01*min(df$h),
               residuals_y_min=1.1*(-max_res),residuals_y_max=1.1*max_res,  
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
    param <- get_param_names(class(m),m$run_info$c_param)
    table <- rbind(m$param_summary,m$Deviance_summary)
    row.names(table) <- c(sapply(1:length(param), function(x) get_param_expression(x)[[1]]),"Deviance")
    rhat_col <- c(rhat_vector(m,param),NA)       
    table <- cbind(table,'R_hat'=rhat_col)
    table <- format(round(table,digits=3),nsmall=3)
    table[nrow(table),4] <- ''
    table <- tableGrob(table,theme=ttheme_minimal(base_family="Times",rowhead=list(fg_params=list(parse=TRUE))))
  })
  if(type==1){ 
    p_mat <- predict_matrix(m_obj[[1]])
  }else{  
    mcmc_hist_list <- lapply(m_obj,function(m){
      params <- get_param_names(class(m),m$run_info$c_param)
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
    conv_diag_plots <- get_conv_diagnostics_plots(m_obj)           
  }
  pdf(file=paste0(directory,report_title),paper='a4',width=8,height=11)
  if(type==1){
      grid.arrange(main_plot_list[[1]],main_table_list[[1]],nrow=2,as.table=TRUE,heights=c(5,3),top=textGrob(class(m_obj[[1]]),gp=gpar(fontfamily="Times",fontsize=22,facetype='bold')))
      grid.arrange(p_mat,nrow=1,as.table=TRUE,heights=c(1),top=textGrob(paste0('Rating curve predictions for ',class(m_obj[[1]])),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
  }else{
    lapply(names(posterior_list),function(x){ 
      grid.arrange(main_plot_list[[x]],main_table_list[[x]],nrow=2,as.table=TRUE,heights=c(5,3),top=textGrob(x,gp=gpar(fontfamily="Times",fontsize=22,facetype='bold')))
    })
    grid.arrange(arrangeGrob(MCMC_table,tour_table,nrow=1),arrangeGrob(dev_boxplot,ncol=2),nrow=2,as.table=TRUE,heights=c(1,1),top=textGrob('Model comparison and MCMC diagnostics',gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
    grid.arrange(conv_diag_plots,as.table=TRUE,top=textGrob('',gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
    lapply(names(posterior_list), function(x) {
      grid.arrange(arrangeGrob(grobs=mcmc_hist_list[[x]],nrow=4,ncol=3),top=textGrob(paste0('Estimated parameters of ',x),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
    })
  }
  dev.off()
}



