


R_hat <- function(posterior,num_ch=4){ 
  l <- length(posterior)/(num_ch)
  r <- 1
  s <- l
  #if(length(posterior)%%2!=0){posterior <- posterior[-1]}
  sims <- c()
  for (i in 1:num_ch) {
    sims <- cbind(sims,posterior[r:s])
    r <- r+l
    s <- s+l
  }
  if(length(sims[,1])%%2!=0){sims <- sims[-1,]}
  l <- length(sims[,i])/2
  r <- 1
  s <- l
  psi <- c()
  for (i in 1:num_ch) {
    psi <- cbind(psi,sims[(r:l),i],sims[((r+l):(s+l)),i])
  }
  m <- num_ch*2
  n <- l
  psi_bar_dj <- sapply(1:m, function(x) mean(psi[,x]))
  psi_bar_dd <- mean(psi_bar_dj)
  B <- ((n)/(m-1))*sum((psi_bar_dj-psi_bar_dd)^2)
  s_sqrd <- sapply(1:8, function(x) ((1)/(n-1))*sum((psi[,x]-psi_bar_dj[x])^2))
  W <- mean(s_sqrd)
  var_hat_p <- ((n-1)/(n))*W+((1)/(n))*B
  rhat <- sqrt(var_hat_p/W)
  return(rhat)
}


rhat_vector <- function(fit){
  if(is.null(fit$c_posterior)){
    if(is.null(fit$phi_beta_posterior)){
      if(is.null(fit$sigma_eta_posterior)){
        rhat_vec <- sapply(1:3,function(x) R_hat(fit[[x]]))
      }else{
        rhat_vec <- sapply(1:9,function(x) R_hat(fit[[x]]))
      }
    }else{
      if(is.null(fit$sigma_eta_posterior)){
        rhat_vec <- sapply(1:5,function(x) R_hat(fit[[x]]))
      }else{
        rhat_vec <- sapply(1:11,function(x) R_hat(fit[[x]]))
      }
    }
  }else{
    if(is.null(fit$phi_beta_posterior)){
      if(is.null(fit$sigma_eta_posterior)){
        rhat_vec <- sapply(1:4,function(x) R_hat(fit[[x]]))
      }else{
        rhat_vec <- sapply(1:10,function(x) R_hat(fit[[x]]))
      }
    }else{
      if(is.null(fit$sigma_eta_posterior)){
        rhat_vec <- sapply(1:6,function(x) R_hat(fit[[x]]))
      }else{
        rhat_vec <- sapply(1:12,function(x) R_hat(fit[[x]]))
      }
    }
  }
  return(rhat_vec)
}


main_list <- function(formula, data, min_obs=10){
  id_var <- colnames(data)[-which(colnames(data) %in% formula_args)]     
  if(class(data[[all_of(id_var)]])!='factor') data[[all_of(id_var)]] <- as.factor(data[[all_of(id_var)]])
  data <- data %>% arrange(data[[all_of(id_var)]]) 
  data <- select(data, id_var, formula_args[1],formula_args[2])
  snames <- unique(data[[all_of(id_var)]])
  main <- list()
  log_file <- list()
  n <- 1
  #create the mainlist
  for (i in 1:length(snames)){
    if(nrow(filter(data, data[[all_of(id_var)]] == snames[i]))<min_obs){
      log_file <- rbind(log_file,paste0(Sys.time(),'  [WARN]  Data set ',snames[i],' was excluded. Minimum requirement of ',min_obs,' observations was not fulfilled.' ))
    }else{
      main$sets[n] <- list(filter(data, data[[all_of(id_var)]] == snames[i])) 
      n <- n+1
    }
  }  
  error_msg7 <- 'No data sets have 10 or more observations. Minimum number of observations required per data set is 10.'
  if(n==1) stop(error_msg7) 
  #new sets vector
  snames <- snames[1:(n-1)]
  for (i in 1:(n-1))  snames[i] <- main$sets[[i]][[all_of(id_var)]][[1]] 
  return(list('snames'=snames,'main'=main,'log_file'=log_file))
}


c_list <- function(c_values,snames){
  if(!any(c_values[,1] %in% snames | c_values[,2] %in% snames)){
    c_values <- NULL
    log_file <- rbind(log_file,paste0(Sys.time(),'  [WARN]  All known c values from the "c_values" data frame were excluded. No names, corresponding to these c values, were found in the main data frame of paired stage and discharge observations.' )) 
  }else{
    n_0 <- nrow(c_values)
    if(any(c_values[,1] %in% snames)){name_col <- 1}else{name_col <- 2}
    c_values <- data.frame(name=c_values[[name_col]],c=c_values[[1+name_col%%2]])
    error_msg7 <- 'The c values need to be numerical.'
    if(!any(sapply(1:length(c_values$c),function(x) is.numeric(c_values$c[x])))){
      stop(error_msg7)
    }
    c_values <- c_values[which(c_values$name %in% snames),]
    c_values <- c_values[sapply(1:length(c_values$c),function(x) is.numeric(c_values$c[x])&!is.na(c_values$c[x])),]           
    n_1 <- nrow(c_values)
    if(n_0!=n_1){
      N <- n_0-n_1
      log_file <- rbind(log_file,paste0(Sys.time(),'  [WARN]  ',N,' known c values from the "c_values" data frame were excluded. Either the names, corresponding to these c values, were not found in the main data frame of paired stage and discharge observations or the actual values were not numerical.' )) 
    }
  }
  return(c_values)
}



plot_c_posterior <- function(fits,snames,c_values){
  c_plots <- list()
  for (i in 1:length(c_values$name)) {
    j <- which(snames %in% c_values$name[[i]])
    fit <- fits$bgplm[[j]]
    # if(all(c_tables$tour[[j]]$winner[c(1,5)] %in% T)){
    #   fit <- fits$bgplm[[j]]
    # }else if(all(c_tables$tour[[j]]$winner[c(2,5)] %in% T)){
    #   fit <- fits$bgplm0[[j]]
    # }else if(all(c_tables$tour[[j]]$winner[c(3,6)] %in% T)){
    #   fit <- fits$bplm[[j]]
    # }else{
    #   fit <- fits$bplm0[[j]]
    # }
    d <- data.frame(c=fit$c_posterior,name=rep('c',length(fit$c_posterior)))
    p <- ggplot(data=d,aes(x=c))+
      geom_density(fill='gray70',color='black',alpha=0.5) +
      geom_vline(xintercept =  c_values$c[[i]], color='red',size=1.5,linetype='solid')+
      theme_bw() +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 13),
            plot.title = element_text(size = 14,face = "bold",hjust = 0.45),
            text = element_text(family="Times", face="plain"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black")) +
      ggtitle(paste0('P(c|Data) and the known c (red line).'))+
      ylab('')
    c_plots[[i]] <- grid.arrange(p)
  }  
  return(list('c_plots'=c_plots))
}

predict_tables <- function(fits,snames,c_tables){
  p_tables <- list()
  for (i in 1:length(snames)) {
    if(all(c_tables$tour[[i]]$winner[c(1,5)] %in% T)){
      fit <- fits$bgplm[[i]]
    }else if(all(c_tables$tour[[i]]$winner[c(2,5)] %in% T)){
      fit <- fits$bgplm0[[i]]
    }else if(all(c_tables$tour[[i]]$winner[c(3,6)] %in% T)){
      fit <- fits$bplm[[i]]
    }else{
      fit <- fits$bplm0[[i]]
    }
    h <- fit$rating_curve[[1]]
    step <- (max(h)-min(h))/2^5
    h_grid <- seq(min(h),max(h),by=step)
    rating_curve_h_grid <- predict(fit,newdata=h_grid)
    p_tables[[i]] <- format(round(rating_curve_h_grid,digits = 3),nsmall=3)
  }
  return(p_tables)
}


comparison_table <- function(fits,snames){
  c_tables <- list()
  for (i in 1:length(snames)) {
    tour <- tournament(fits$bgplm[[i]],fits$bgplm0[[i]],fits$bplm[[i]],fits$bplm0[[i]])
    c_tables$tour[[i]] <- tour$summary
  } 
  return(c_tables)
}


mcmc_diag_table <- function(fits,snames){
  mcmc_tables <- list()
  for (i in 1:length(snames)) {
    mcmc_table <- data.frame()
    mcmc_table[1,1:4] <- c(fits$bgplm[[i]]$run_info$nr_iter,fits$bgplm0[[i]]$run_info$nr_iter,fits$bplm[[i]]$run_info$nr_iter,fits$bplm0[[i]]$run_info$nr_iter)
    mcmc_table[2,1:4] <- c(fits$bgplm[[i]]$run_info$num_chains,fits$bgplm0[[i]]$run_info$num_chains,fits$bplm[[i]]$run_info$num_chains,fits$bplm0[[i]]$run_info$num_chains)
    mcmc_table[3,1:4] <- c(fits$bgplm[[i]]$run_info$burnin,fits$bgplm0[[i]]$run_info$burnin,fits$bplm[[i]]$run_info$burnin,fits$bplm0[[i]]$run_info$burnin)
    mcmc_table[4,1:4] <- c(fits$bgplm[[i]]$run_info$thin,fits$bgplm0[[i]]$run_info$thin,fits$bplm[[i]]$run_info$thin,fits$bplm0[[i]]$run_info$thin)
    mcmc_table[5,1:4] <- c(fits$bgplm[[i]]$num_effective_param,fits$bgplm0[[i]]$num_effective_param,fits$bplm[[i]]$num_effective_param,fits$bplm0[[i]]$num_effective_param)
    mcmc_table[6,1:4] <- c(fits$bgplm[[i]]$acceptance_rate,fits$bgplm0[[i]]$acceptance_rate,fits$bplm[[i]]$acceptance_rate,fits$bplm0[[i]]$acceptance_rate)
    mcmc_table[5:6,] <- round(mcmc_table[5:6,],digits = 3)
    mcmc_table[1:4,] <- format(mcmc_table[1:4,],nsmall = 0)
    colnames(mcmc_table) <- c('gplm','gplm0','plm','plm0')
    rownames(mcmc_table) <- c('nr_iter','nr_chains','burnin','thin','nr_eff_param','accept_rate')
    mcmc_tables[[i]] <- mcmc_table
  }
  return(mcmc_tables)
}


dev_boxplot <- function(fits,snames){
  d_plots <- list()
  for (i in 1:length(snames)) {
    val_gplm <- fits$bgplm[[i]]$Deviance_posterior
    d1 <- cbind(rep("gplm", length(val_gplm)),val_gplm)
    val_gplm0 <- fits$bgplm0[[i]]$Deviance_posterior
    d2 <- cbind(rep("gplm0", length(val_gplm0)),val_gplm0)
    val_plm <- fits$bplm[[i]]$Deviance_posterior
    d3 <- cbind(rep("plm", length(val_plm)),val_plm)
    val_plm0 <- fits$bplm0[[i]]$Deviance_posterior
    d4 <- cbind(rep("plm0", length(val_plm0)),val_plm0)
    dev_df <- data.frame(rbind(d1,d2,d3,d4))
    colnames(dev_df) <- c("model","values")
    dev_df$values <- as.numeric(dev_df$values)
    DIC_p <- data.frame(model=c('gplm','gplm0','plm','plm0'), DIC=c(fits$bgplm[[i]]$DIC,fits$bgplm0[[i]]$DIC,fits$bplm[[i]]$DIC,fits$bplm0[[i]]$DIC))
    p <- ggplot(data = dev_df, aes(x=model, y=values)) +
      geom_boxplot(size=.4,color="black",outlier.size=0.1,outlier.shape=21,outlier.fill="gray90",fill="gray90") +
      stat_boxplot(geom ='errorbar') +
      geom_line(data = DIC_p, aes(x=model, y=DIC, group = 1),color='gray30') +
      geom_point(data = DIC_p, aes(x=model, y=DIC), size=3,shape=23,fill='red2',color='black') +
      theme_bw() +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 13),
            axis.text.x = element_text(size = 13),
            plot.title = element_text(size = 14,face = "bold",hjust = 0.45),
            text = element_text(family="Times", face="plain"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black")) +
      ggtitle('Deviance posterior and DIC ')
    d_plots[[i]] <- p
  }
  return(d_plots)
}


histogram_plots <- function(fit){
  if(!is.null(fit$c_posterior)){
    if(is.null(fit$sigma_eta_posterior)&is.null(fit$phi_beta_posterior)){
      p <- plot(fit,type='histogram','a','b','c','sigma_eps',transformed=F)
    }else if(!is.null(fit$sigma_eta_posterior)&is.null(fit$phi_beta_posterior)){
      p <- plot(fit,type='histogram','a','b','c','sigma_eta','eta_1','eta_2','eta_3','eta_4','eta_5','eta_6',transformed=F)
    }else if(is.null(fit$sigma_eta_posterior)&!is.null(fit$phi_beta_posterior)){
      p <- plot(fit,type='histogram','a','b','c','sigma_eps','sigma_beta','phi_beta',transformed=F)
    }else{
      p <- plot(fit,type='histogram','a','b','c','sigma_eta','eta_1','eta_2','eta_3','eta_4','eta_5','eta_6','sigma_beta','phi_beta',transformed=F)
    } 
  }else{
    if(is.null(fit$sigma_eta_posterior)&is.null(fit$phi_beta_posterior)){
      p <- plot(fit,type='histogram','a','b','sigma_eps',transformed=F)
    }else if(!is.null(fit$sigma_eta_posterior)&is.null(fit$phi_beta_posterior)){
      p <- plot(fit,type='histogram','a','b','sigma_eta','eta_1','eta_2','eta_3','eta_4','eta_5','eta_6',transformed=F)
    }else if(is.null(fit$sigma_eta_posterior)&!is.null(fit$phi_beta_posterior)){
      p <- plot(fit,type='histogram','a','b','sigma_eps','sigma_beta','phi_beta',transformed=F)
    }else{
      p <- plot(fit,type='histogram','a','b','sigma_eta','eta_1','eta_2','eta_3','eta_4','eta_5','eta_6','sigma_beta','phi_beta',transformed=F)
    } 
  }
  return(p)
}


mcmc_plots <- function(fits,snames){
  hist_plots <- list()
  for (i in 1:length(snames)) {
    hist_plots$bgplm[[i]] <- histogram_plots(fits$bgplm[[i]])
    hist_plots$bgplm0[[i]] <- histogram_plots(fits$bgplm0[[i]])
    hist_plots$bplm[[i]] <- histogram_plots(fits$bplm[[i]])
    hist_plots$bplm0[[i]] <- histogram_plots(fits$bplm0[[i]])
  }
  return(hist_plots)
}


multiplier <- function (Q_median, h_i) {
  Q_median_temp1 <- Q_median[0,]
  Q_median_temp2 <- Q_median
  account <- data.frame(table(h_i[which( h_i %in% h_i[duplicated(h_i)])]))
  if (nrow(account)!=0) {
    b <- 0
    for (i in 1:nrow(account)) {
      row_num <- as.numeric(rownames(Q_median[which(Q_median[[2]] %in% account[[1]][i] ),]))
      row_num <- row_num + b
      Q_median_temp1 <- slice(Q_median,rep(row_num-b,each=account[[2]][i])) 
      if (row_num==1) {
        Q_median_temp2  <- rbind(Q_median_temp1,Q_median[Q_median[[2]]>Q_median_temp1[1,2],])
      } else {
        Q_median_temp2  <- rbind(Q_median_temp2[(1:(row_num-1)),],Q_median_temp1,Q_median[Q_median[[2]]>Q_median_temp1[1,2],])
      }
      row.names(Q_median_temp2) <- NULL
      b <- b+(account$Freq[i]-1)
    }
  }
  return(Q_median_temp2)
}


logQmed <- function(fit){
  df1 <- as.data.frame(fit$data[[2]])
  df2 <- as.data.frame(fit$rating_curve$median)
  df3 <- as.data.frame(fit$rating_curve$h)
  df2and3 <- cbind(df2,df3)
  Q_median <- df2and3 %>% filter(df2and3[[2]] %in% df1[[1]]) 
  Q_median <- multiplier(Q_median,fit$data[[2]])
  logQ_median <- log(Q_median[[1]])
  return(logQ_median)
}



max_res <- function(fit) {
  logQ_i <- log(fit$data[[1]])
  logQ_median <- logQmed(fit)                    
  max_value <- max(abs(logQ_i - logQ_median))
  return(max_value)
}



page_prep <- function(fit1,fit2,fit3,fit4) {
  rate_upr_lim <- min(fit1$data[[2]]) + abs(min(fit1$data[[2]])-max(fit1$data[[2]]))*1.01
  rate_lwr_lim <- min(fit1$data[[2]])
  rate_lim <- c(rate_lwr_lim,rate_upr_lim)
  rate_upr_lim_x <- max(fit1$data[[1]])*1.2
  rate_lwr_lim_x <- min(min(fit1$rating_curve$lower),
                        min(fit2$rating_curve$lower),
                        min(fit3$rating_curve$lower),
                        min(fit4$rating_curve$lower))
  rate_lim_x <- c(rate_lwr_lim_x,rate_upr_lim_x)
  maxRes <- max(max_res(fit1),max_res(fit2),max_res(fit3),max_res(fit4))
  res_upr_lim <- max(max(fit1$sigma_eps_summary$median*1.96*1.1),
                     fit2$param_summary['sigma_eps',]$median*1.96*1.1,
                     max(fit3$sigma_eps_summary$median*1.96*1.1),
                     fit4$param_summary['sigma_eps',]$median*1.96*1.1,
                     maxRes*1.1)
  res_lwr_lim <- min(min(fit1$sigma_eps_summary$median*(-1.96)*1.1),
                     fit2$param_summary['sigma_eps',]$median*(-1.96)*1.1,
                     min(fit3$sigma_eps_summary$median*(-1.96)*1.1),
                     fit4$param_summary['sigma_eps',]$median*(-1.96)*1.1,
                     (-maxRes)*1.1)
  res_lim <- c(res_lwr_lim,res_upr_lim)
  beta_upr_lim <- max(max(fit1$f_summary$upper)*1.01,
                      max(fit2$f_summary$upper)*1.01,
                      fit3$param_summary['b',]$upper*1.01,
                      fit4$param_summary['b',]$upper*1.01,
                      3.5)
  beta_lwr_lim <- min(min(fit1$f_summary$lower)*0.99,
                      min(fit2$f_summary$lower)*0.99,
                      fit3$param_summary['b',]$lower*0.99,
                      fit4$param_summary['b',]$lower*0.99,
                      1)
  beta_lim <- c(beta_lwr_lim,beta_upr_lim)
  beta_upr_lim_x <- max(max(fit1$f_summary$h),
                        max(fit2$f_summary$h))
  beta_lwr_lim_x <- min(min(fit1$f_summary$h),
                        min(fit2$f_summary$h))
  beta_lim_x <- c(beta_lwr_lim_x,beta_upr_lim_x)
  sig_upr_lim <- max(max(fit1$sigma_eps_summary$upper),
                     fit2$param_summary['sigma_eps',]$upper,
                     max(fit3$sigma_eps_summary$upper),
                     fit4$param_summary['sigma_eps',]$upper)
  sig_lwr_lim <- 0
  sig_lim <- c(sig_lwr_lim,sig_upr_lim)
  sig_upr_lim_x <- max(fit1$data[[2]])
  sig_lwr_lim_x <- min(fit1$data[[2]])
  sig_lim_x <- c(sig_lwr_lim_x,sig_upr_lim_x)
  return(list('res_lim'=res_lim,'rate_lim'=rate_lim,'rate_lim_x'=rate_lim_x,'beta_lim'=beta_lim,'beta_lim_x'=beta_lim_x,'sig_lim'=sig_lim,'sig_lim_x'=sig_lim_x))
}



rc_plotter <- function(fit,rclim,rclimx){
  p <- ggplot(data=fit$rating_curve) +
    geom_point(data=fit$data,aes(fit$data[[1]],fit$data[[2]]), size=.9, shape=21, fill="gray60", color="black") +
    geom_line(aes(median,h)) +
    geom_path(aes(lower,h),linetype='dashed') +
    geom_path(aes(upper,h),linetype='dashed') +
    xlab(expression(paste('Q[',m^{3},'/s]'))) + 
    ylab("h[m]") +
    scale_y_continuous(expand = c(0,0), limits = c(rclim[1],rclim[2]) ) +
    scale_x_continuous(expand = c(0,0), limits = c(rclimx[1],rclimx[2]) ) +
    theme_bw() +
    theme( text = element_text(family="Times", face="plain"),
           axis.text = element_text(face="plain"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           panel.border = element_rect(colour = "black"))
  return(p)
}

res_plotter <- function(fit,rlim,c_value){
  if(is.null(c_value)){
    c <- fit$param_summary['c',]$median
    c_exp <- expression(paste('log(h - ',hat(c),')'))
  }else{
    c <- c_value
    c_exp <- expression(paste('log(h - c)'))
  }
  h <- fit$rating_curve$h
  h_i <- fit$data[[2]]
  logQ_i <- log(fit$data[[1]])
  logQ_median <- logQmed(fit)
  if(is.null(fit$sigma_eps_summary)){
    upr_r <- rep(fit$param_summary['sigma_eps',]$median*1.96, length(h))
    lwr_r <- rep(fit$param_summary['sigma_eps',]$median*(-1.96), length(h))
    res_2sd <- data.frame(upr_r,lwr_r,h)
  }else{
    upr_r <- fit$sigma_eps_summary$median*1.96
    lwr_r <- fit$sigma_eps_summary$median*(-1.96)
    res_2sd <- data.frame(upr_r,lwr_r,h)
  }
  p <- ggplot(fit$data, aes( log( h_i - c ) )) + 
    geom_point(aes(y=( logQ_i - logQ_median ) ), size=.9, shape=21, fill="gray60", color="black") +
    geom_hline(yintercept = 0, linetype="solid") +
    geom_line(data=filter(res_2sd,h>=min(fit$data[[2]]),h<=max(fit$data[[2]])),aes(log(h-c),upr_r),linetype='dashed') +
    geom_line(data=filter(res_2sd,h>=min(fit$data[[2]]),h<=max(fit$data[[2]])),aes(log(h-c),lwr_r),linetype='dashed') +
    ylim( rlim[1],rlim[2] ) +
    labs(x = c_exp, y = expression(paste('log(Q)-log(',hat(Q),')'))) +
    theme_bw() +
    theme( text = element_text(family="Times", face="plain"),
           axis.text = element_text(face="plain"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           panel.border = element_rect(colour = "black")) 
  return(p)
}

b_plus_beta_plotter <- function(fit,blim,blimx){
  p <- ggplot(data=fit$f_summary)  +
    geom_line(aes(h, median)) +
    geom_line(aes(h, lower),linetype='dashed') +
    geom_line(aes(h, upper),linetype='dashed') +
    ylim( blim[1],blim[2] ) +
    xlim( blimx[1],blimx[2] ) +
    theme_bw() +
    theme( text = element_text(family="Times", face="plain"),
           axis.text = element_text(face="plain"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           panel.border = element_rect(colour = "black")) +
    labs( x = "h[m]",
          y = expression(paste('b + ',beta,'(h)')) ) 
  return(p)
}


b_plotter <- function(fit,blim,blimx){
  h <- fit$rating_curve$h
  lwr_b <- rep(fit$param_summary['b',]$lower, length(h))
  med_b <- rep(fit$param_summary['b',]$median, length(h))
  upr_b <- rep(fit$param_summary['b',]$upper, length(h))
  b_df <- data.frame(h, upr_b, med_b, lwr_b)
  p <- ggplot(data=b_df)  +
    geom_line(aes(h, med_b)) +
    geom_line(aes(h, lwr_b),linetype='dashed') +
    geom_line(aes(h, upr_b),linetype='dashed') +
    ylim( blim[1],blim[2] ) +
    xlim( blimx[1],blimx[2] ) +
    theme_bw() +
    theme( text = element_text(family="Times", face="plain"),
           axis.text = element_text(face="plain"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           panel.border = element_rect(colour = "black")) +
    labs( x = 'h[m]', y = 'b' ) 
  return(p)
}

sigma_eps_h_plotter <- function(fit,slim,slimx){
  p <- ggplot(data=fit$sigma_eps_summary)  +
    geom_line(aes(h, median)) +
    geom_line(aes(h, lower),linetype='dashed') +
    geom_line(aes(h, upper),linetype='dashed') +
    labs(x = "h[m]", y = expression(paste(sigma[epsilon],'(h)'))) +
    scale_y_continuous(expand = c(0, 0), limits = c(slim[1],slim[2]*1.05)) +
    xlim( slimx[1],slimx[2] ) +
    theme_bw() +
    theme( text = element_text(family="Times", face="plain"),
           axis.text = element_text(face="plain"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           panel.border = element_rect(colour = "black")) 
  return(p)
}


sigma_eps_plotter <- function(fit,slim,slimx){
  h <- fit$rating_curve$h
  lwr_s <- rep(fit$param_summary['sigma_eps',]$lower, length(h))
  med_s <- rep(fit$param_summary['sigma_eps',]$median, length(h))
  upr_s <- rep(fit$param_summary['sigma_eps',]$upper, length(h))
  sigma_eps_df <- data.frame(h, upr_s, med_s, lwr_s)
  p <- ggplot(data=sigma_eps_df)  +
    geom_line(aes(h, med_s)) +
    geom_line(aes(h, lwr_s),linetype='dashed') +
    geom_line(aes(h, upr_s),linetype='dashed') +
    labs(x = 'h[m]', y = expression(sigma[epsilon])) +
    scale_y_continuous(expand = c(0, 0), limits = c(slim[1],slim[2]*1.05)) +
    xlim( slimx[1],slimx[2] ) +
    theme_bw() +
    theme( text = element_text(family="Times", face="plain"),
           axis.text = element_text(face="plain"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           panel.border = element_rect(colour = "black")) 
  return(p)
}




bgplm_page_maker <- function(fit,snames,rlim,rclim,rclimx,blim,blimx,slim,slimx,c_value) {
  p1 <- rc_plotter(fit,rclim,rclimx)
  p2 <- res_plotter(fit,rlim,c_value)
  p3 <- b_plus_beta_plotter(fit,blim,blimx)
  p4 <- sigma_eps_h_plotter(fit,slim,slimx)
  plots <- grid.arrange(p1,p2,p3,p4,top=textGrob(paste(as.character(snames),as.character("   -   "),as.character("gplm"),sep=""),gp=gpar(fontfamily="Times",fontsize=19,facetype='bold')))
  DEV <- as.data.frame(fit$Deviance_summary)
  row.names(DEV) <- "Deviance"
  table <- as.data.frame(fit$param_summary)
  rownames(table) <- sapply(1:nrow(table), function(x) get_param_expression( rownames(table)[[x]] )[[1]])
  table <- rbind(table, DEV)
  rhat_col <- data.frame(Rhat=c(rhat_vector(fit),NA))
  table <- cbind(table,rhat_col)
  table <- format(round(table,digits = 3),nsmall=3)
  table[length(table[,4]),4] <- ''
  table <- tableGrob(table, theme=ttheme_minimal(base_family = "Times",rowhead=list(fg_params = list(parse=TRUE))))
  return(list("plots" = plots, "table" = table))
}


bgplm0_page_maker <- function(fit,snames,rlim,rclim,rclimx,blim,blimx,slim,slimx,c_value) {
  p1 <- rc_plotter(fit,rclim,rclimx)
  p2 <- res_plotter(fit,rlim,c_value)
  p3 <- b_plus_beta_plotter(fit,blim,blimx)
  p4 <- sigma_eps_plotter(fit,slim,slimx)
  plots <- grid.arrange(p1,p2,p3,p4,top=textGrob(paste(as.character(snames),as.character("   -   "),as.character("gplm0"),sep=""),gp=gpar(fontfamily="Times",fontsize=19,facetype='bold')))
  DEV <- as.data.frame(fit$Deviance_summary)
  row.names(DEV) <- "Deviance"
  table <- as.data.frame(fit$param_summary)
  rownames(table) <- sapply(1:nrow(table), function(x) get_param_expression( rownames(table)[[x]] )[[1]])
  table <- rbind(table, DEV)
  rhat_col <- data.frame(Rhat=c(rhat_vector(fit),NA))
  table <- cbind(table,rhat_col)
  table <- format(round(table,digits = 3),nsmall=3)
  table[length(table[,4]),4] <- ''
  table <- tableGrob(table, theme=ttheme_minimal(base_family = "Times",rowhead=list(fg_params = list(parse=TRUE))))
  return(list("plots" = plots, "table" = table))
}


bplm_page_maker <- function(fit,snames,rlim,rclim,rclimx,blim,blimx,slim,slimx,c_value) {
  p1 <- rc_plotter(fit,rclim,rclimx)
  p2 <- res_plotter(fit,rlim,c_value)
  p3 <- b_plotter(fit,blim,blimx)
  p4 <- sigma_eps_h_plotter(fit,slim,slimx)
  plots <- grid.arrange(p1,p2,p3,p4,top=textGrob(paste(as.character(snames),as.character("   -   "),as.character("plm"),sep=""),gp=gpar(fontfamily="Times",fontsize=19,facetype='bold')))
  DEV <- as.data.frame(fit$Deviance_summary)
  row.names(DEV) <- "Deviance"
  table <- as.data.frame(fit$param_summary)
  rownames(table) <- sapply(1:nrow(table), function(x) get_param_expression( rownames(table)[[x]] )[[1]])
  table <- rbind(table, DEV)
  rhat_col <- data.frame(Rhat=c(rhat_vector(fit),NA))
  table <- cbind(table,rhat_col)
  table <- format(round(table,digits = 3),nsmall=3)
  table[length(table[,4]),4] <- ''
  table <- tableGrob(table, theme=ttheme_minimal(base_family = "Times",rowhead=list(fg_params = list(parse=TRUE))))
  return(list("plots" = plots, "table" = table))
}



bplm0_page_maker <- function(fit,snames,rlim,rclim,rclimx,blim,blimx,slim,slimx,DICt,c_value) {
  p1 <- rc_plotter(fit,rclim,rclimx)
  p2 <- res_plotter(fit,rlim,c_value)
  p3 <- b_plotter(fit,blim,blimx)
  p4 <- sigma_eps_plotter(fit,slim,slimx)
  plots <- grid.arrange(p1,p2,p3,p4,top=textGrob(paste( as.character(snames),as.character("   -   "),as.character("plm0"),sep=""),gp=gpar(fontfamily="Times",fontsize=19,facetype='bold')))
  DEV <- as.data.frame(fit$Deviance_summary)
  row.names(DEV) <- "Deviance"
  table <- as.data.frame(fit$param_summary)
  rownames(table) <- sapply(1:nrow(table), function(x) get_param_expression( rownames(table)[[x]] )[[1]])
  table <- rbind(table, DEV)
  rhat_col <- data.frame(Rhat=c(rhat_vector(fit),NA))
  table <- cbind(table,rhat_col)
  table <- format(round(table,digits = 3),nsmall=3)
  table[length(table[,4]),4] <- ''
  table <- tableGrob(table, theme=ttheme_minimal(base_family = "Times",rowhead=list(fg_params = list(parse=TRUE))))
  return(list("plots" = plots, "table" = table))
}

