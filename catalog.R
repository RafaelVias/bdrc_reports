
source('/Users/Vias/Desktop/R/HI/Bayes_Sumar/Rennsli_2/functions/catalog_functions.R')


catalog <- function(formula,data,c_values=NULL,directory,type=1){
  begin <- Sys.time()
  #argument checking
  error_msg1 <- 'Please provide a formula and data.'
  error_msg2 <- 'The data should have 3 variables: (1) Water stage (2) Water discharge (3) Name variable for seperating the different data sets.'
  error_msg3 <- 'The data should only have 3 variables: (1) Water stage (2) Water discharge (3) Name variable for seperating the different data sets.'
  error_msg4 <- 'Directory does not exist. Please enter a valid directory.'
  error_msg5 <- 'Invalid catalog type. Only three possible catalog types: 1, 2 and 3.'
  error_msg6 <- 'The known c values should be input as a data frame with 2 variables: (1) A name variable and (2) a corresponding c value.'
  if(is.null(formula) | is.null(data)){
    stop(error_msg1)  
  }else if(ncol(data)<3){
    stop(error_msg2)
  }else if(ncol(data)>3){
    stop(error_msg3)
  }else if(!dir.exists(directory)){
    stop(error_msg4)
  }else if(type!=1 & type!=2 & type!=3){
    stop(error_msg5)
  }else if(!is.null(c_values)&ncol(c_values)!=2){
    stop(error_msg6)
  }
  stopifnot('formula' %in% class(formula))
  stopifnot('data.frame' %in% class(data))
  stopifnot('data.frame' %in% class(c_values))
  stopifnot('character' %in% class(directory))
  formula_args <- all.vars(formula)
  stopifnot(length(formula_args)==2 & all(formula_args %in% names(data)))
  #start seting up the data
  main_list_output <- main_list(formula,data)
  snames <- main_list_output[[1]]
  main <- main_list_output[[2]]
  log_file <- main_list_output[[3]]
  if(!is.null(c_values)){
    c_values <- c_list(c_values,snames)
  }
  #run the models
  fits <- list()
  print_load <- T
  b <- 0
  for (i in 1:length(snames)) {
    rc_dat <- data.frame(main$sets[i])
    set.seed(3004)
    fits$bgplm[[i]] <- bgplm(formula,rc_dat)
    set.seed(3004)
    fits$bgplm0[[i]] <- bgplm0(formula,rc_dat)
    set.seed(3004)
    fits$bplm[[i]] <- bplm(formula,rc_dat)
    set.seed(3004)
    fits$bplm0[[i]] <- bplm0(formula,rc_dat)
    b <- b+1
    if(print_load==T) print(paste(100*round(b/(length(snames)+nrow(c_values)),digits = 3),"%",sep = ""))
  }
  if(!is.null(c_values)){
    fits_c <- list()
    for (i in 1:nrow(c_values)) {
      j <- which(snames %in% c_values$name[i])
      rc_dat <- data.frame(main$sets[j])
      set.seed(3004)
      fits_c$bgplm[[j]] <- bgplm(formula,rc_dat,c_param = c_values$c[i])
      set.seed(3004)
      fits_c$bgplm0[[j]] <- bgplm0(formula,rc_dat,c_param = c_values$c[i])
      set.seed(3004)
      fits_c$bplm[[j]] <- bplm(formula,rc_dat,c_param = c_values$c[i])
      set.seed(3004)
      fits_c$bplm0[[j]] <- bplm0(formula,rc_dat,c_param = c_values$c[i])
      b <- b+1
      if(print_load==T) print(paste(100*round(b/(length(snames)+nrow(c_values)),digits = 3),"%",sep = ""))
    }
    kc <- which(sapply(1:length(fits_c$bgplm),function(x) !is.null(fits_c$bgplm[[x]])))
    fits_com <- list()
    for (i in 1:length(snames)) {
      if(i %in% kc){
        fits_com$bgplm[[i]] <- fits_c[[1]][[i]]
        fits_com$bgplm0[[i]] <- fits_c[[2]][[i]]
        fits_com$bplm[[i]] <- fits_c[[3]][[i]]
        fits_com$bplm0[[i]] <- fits_c[[4]][[i]]
      }else{
        fits_com$bgplm[[i]] <- fits[[1]][[i]]
        fits_com$bgplm0[[i]] <- fits[[2]][[i]]
        fits_com$bplm[[i]] <- fits[[3]][[i]]
        fits_com$bplm0[[i]] <- fits[[4]][[i]]
      }
    }
  }else{ 
    fits_com <- fits
  }
  #start creating the objects to be presented in the pdf
  pmo <- list()
  for (i in 1:length(snames)) {
    fit1 <- fits_com$bgplm[[i]]
    fit2 <- fits_com$bgplm0[[i]]
    fit3 <- fits_com$bplm[[i]]
    fit4 <- fits_com$bplm0[[i]]
    if(snames[[i]] %in% c_values$name){c_value <- c_values$c[which(c_values$name %in% snames[[i]])]}else{c_value <- NULL}
    ppo <- page_prep(fit1,fit2,fit3,fit4)
    DICt <- data.frame(gplm=round(fit1$DIC,digits=2),gplm0=round(fit2$DIC,digits=2),plm=round(fit3$DIC,digits=2),plm0=round(fit4$DIC,digits=2))
    pmo$bgplm[[i]] <- bgplm_page_maker(fit1,snames[i],ppo[[1]],ppo[[2]],ppo[[3]],ppo[[4]],ppo[[5]],ppo[[6]],ppo[[7]],c_value)
    pmo$bgplm0[[i]] <- bgplm0_page_maker(fit2,snames[i],ppo[[1]],ppo[[2]],ppo[[3]],ppo[[4]],ppo[[5]],ppo[[6]],ppo[[7]],c_value)
    pmo$bplm[[i]] <- bplm_page_maker(fit3,snames[i],ppo[[1]],ppo[[2]],ppo[[3]],ppo[[4]],ppo[[5]],ppo[[6]],ppo[[7]],c_value)
    pmo$bplm0[[i]] <- bplm0_page_maker(fit4,snames[i],ppo[[1]],ppo[[2]],ppo[[3]],ppo[[4]],ppo[[5]],ppo[[6]],ppo[[7]],DICt,c_value)
  }
  if(type==1) {
    c_tables <- comparison_table(fits_com,snames)
    p_tables <- predict_tables(fits_com,snames,c_tables)
  }else if(type==2){
    c_tables <- comparison_table(fits_com,snames)           
    mcmc_tables <- mcmc_diag_table(fits_com,snames)
    d_plots <- dev_boxplot(fits_com,snames)
    c_plots <- plot_c_posterior(fits,snames,c_values)   
  }else if(type==3){
    h_plots <- mcmc_plots(fits_com,snames)
    c_tables <- comparison_table(fits_com,snames)
    mcmc_tables <- mcmc_diag_table(fits_com,snames)
    d_plots <- dev_boxplot(fits_com,snames)
    c_plots <- plot_c_posterior(fits,snames,c_values) 
  }
  #start printing the pdf
  pdf(file=paste(directory,'/catalog.pdf',sep = ''),paper = 'a4',width=8,height=11)
  for (i in 1:length(snames)) {
    if(type==1){
      winner <- as.character(c_tables$tour[[i]]$model[5:6][c_tables$tour[[i]]$winner[5:6]])
      grid.arrange(pmo[[winner]][[i]]$plots,pmo[[winner]][[i]]$table,nrow=2,as.table=TRUE,heights=c(5,3))
      p_table <- tableGrob(p_tables[[i]],theme=ttheme_minimal(base_family = "Times"),rows = NULL)
      grid.arrange(p_table,nrow=1,as.table=TRUE,heights=c(1),
                   top=textGrob(paste0('Rating curve predictions for ',as.character(snames[i])),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
    }else if(type==2){
      grid.arrange(pmo$bgplm[[i]]$plots, pmo$bgplm[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      grid.arrange(pmo$bgplm0[[i]]$plots, pmo$bgplm0[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      grid.arrange(pmo$bplm[[i]]$plots, pmo$bplm[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      grid.arrange(pmo$bplm0[[i]]$plots, pmo$bplm0[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      c_table <- tableGrob(c_tables$tour[[i]], theme=ttheme_minimal(base_family = "Times"),rows = NULL)
      mcmc_table <- tableGrob(mcmc_tables[[i]], theme=ttheme_minimal(base_family = "Times"))
      if(snames[[i]] %in% c_values$name){
        j <- which(c_values$name %in% snames[[i]])
        grid.arrange(c_table,arrangeGrob(d_plots[[i]],arrangeGrob( c_plots$c_plots[[j]], mcmc_table,nrow=2) ,ncol = 2),nrow=2,as.table=TRUE,heights=c(1,2),
                     top=textGrob(paste0('Model comparison for ',as.character(snames[i])),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
      }else{
        grid.arrange(c_table,arrangeGrob(d_plots[[i]],mcmc_table,ncol = 2),nrow=2,as.table=TRUE,heights=c(1,2),
                     top=textGrob(paste0('Model comparison for ',as.character(snames[i])),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
      }
    }else if(type==3){
      grid.arrange(pmo$bgplm[[i]]$plots, pmo$bgplm[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      grid.arrange(pmo$bgplm0[[i]]$plots, pmo$bgplm0[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      grid.arrange(pmo$bplm[[i]]$plots, pmo$bplm[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      grid.arrange(pmo$bplm0[[i]]$plots, pmo$bplm0[[i]]$table, nrow=2,as.table=TRUE,heights=c(5,3))
      c_table <- tableGrob(c_tables$tour[[i]], theme=ttheme_minimal(base_family = "Times"),rows = NULL)
      mcmc_table <- tableGrob(mcmc_tables[[i]], theme=ttheme_minimal(base_family = "Times"))
      if(snames[[i]] %in% c_values$name){
        j <- which(c_values$name %in% snames[[i]])
        grid.arrange(c_table,arrangeGrob(d_plots[[i]],arrangeGrob( c_plots$c_plots[[j]], mcmc_table,nrow=2) ,ncol = 2),nrow=2,as.table=TRUE,heights=c(1,2),
                     top=textGrob(paste0('Model comparison for ',as.character(snames[i])),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
      }else{
        grid.arrange(c_table,arrangeGrob(d_plots[[i]],mcmc_table,ncol = 2),nrow=2,as.table=TRUE,heights=c(1,2),
                     top=textGrob(paste0('Model comparison for ',as.character(snames[i])),gp=gpar(fontsize=22,facetype='bold',fontfamily="Times")))
      }
      grid.arrange(h_plots$bgplm[[i]],nrow=2,as.table=TRUE,heights=c(3,1),
                   top=textGrob(paste0('gplm MCMC sampling for ',as.character(snames[i])),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
      grid.arrange(h_plots$bgplm0[[i]],nrow=2,as.table=TRUE,heights=c(1,1),
                   top=textGrob(paste0('gplm0 MCMC sampling for ',as.character(snames[i])),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
      if(snames[[i]] %in% c_values$name){
        grid.arrange(h_plots$bplm[[i]],nrow=2,as.table=TRUE,heights=c(3,1),
                     top=textGrob(paste0('plm MCMC sampling for ',as.character(snames[i])),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
        grid.arrange(h_plots$bplm0[[i]],nrow=2,as.table=TRUE,heights=c(1,3),
                     top=textGrob(paste0('plm0 MCMC sampling for ',as.character(snames[i])),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
      }else{
        grid.arrange(h_plots$bplm[[i]],nrow=2,as.table=TRUE,heights=c(3,1),
                     top=textGrob(paste0('plm MCMC sampling for ',as.character(snames[i])),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
        grid.arrange(h_plots$bplm0[[i]],nrow=2,as.table=TRUE,heights=c(2,2),
                     top=textGrob(paste0('plm0 MCMC sampling for ',as.character(snames[i])),gp=gpar(fontsize=20,facetype='bold',fontfamily="Times")))
      }
    }
  }
  dev.off()
  end <- Sys.time()
  t_diff <- end-begin
  print(t_diff)
}




######################################################################################################################
#######################
# EXAMPLE 1 ###########
formula <-  Q~W
data <- rbind(mainlist$rivers[[2]],mainlist$rivers[[31]],mainlist$rivers[[26]],mainlist$rivers[[25]][1:9,])
data <- data.frame(data)
c_values <- data.frame(name=c('SUNDSTORP','GREDEBY','RANSTA'),c_val=c(65.7,7.15,23.4))
directory <- "/Users/Vias/Desktop/R/HI/Bayes_Sumar/Rennsli_2"
type <- 3
# EXAMPLE 2 ###########
formula <-  Q~W
data <- data.frame(mainlist$rivers[[1]])[0,]
for (i in 1:65) data <- rbind(data,mainlist$rivers[[i]])
c_values <- data.frame(raf=mainlist$value_c,dan=all_rnames)
directory <- "/Users/Vias/Desktop/R/HI/Bayes_Sumar/Rennsli_2"
type <- 3
# EXAMPLE 3 ###########
formula <-  Q~W
data <- rbind(mainlist$rivers[[2]])
data <- data.frame(data)
c_values <- data.frame(name=c('SUNDSTORP','GREDEBY','RANSTA'),c_val=c(65.7,7.15,23.4))
directory <- "/Users/Vias/Desktop/R/HI/Bayes_Sumar/Rennsli_2"
type <- 3
#######################
# catalog(formula = Q~W, data=d, c_values=md, directory="/Users/Vias/Desktop/R/HI/Bayes_Sumar/Rennsli_2",type=2)
# catalog(formula = Q~W, data=d, c_values=md, directory="/Users/Vias/Desktop/R/HI/Bayes_Sumar",type=1)
######################################################################################################################
