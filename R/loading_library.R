#### libraries   ----------------------------

library(tensorflow)
library(ggplot2)
library(patchwork)
library(keras)
library(ggthemes)
library(quantmod)
library(forecast)
library(KernSmooth)  
# library(locpol)  
library(locfit) 
require(stats); require(graphics)
library(splines)

library(tseries)
#library(rugarch)
library(prophet)
library(tsfknn)
library(neuralnet)
library(randomForest)
library(e1071)
library(xgboost)
library(Metrics)



days_before_date <- function(input_date, days_before) {
  
  date_obj <- as.Date(input_date)
  
  result_date <- date_obj - days_before
  
  return(format(result_date, "%Y-%m-%d"))
}
stretch_y <- function(y, y_mid_min, y_mid_max, y_new_min, y_new_max) {
  ifelse(y < y_mid_min, (y_new_min - y_mid_min) * y / y_mid_min + y_mid_min,
         ifelse(y > y_mid_max, (y_new_max - y_mid_max) * (y - y_mid_max) / (1 - y_mid_max) + y_mid_max,
                (y_new_max - y_new_min) * (y - y_mid_min) / (y_mid_max - y_mid_min) + y_new_min))
}
##################################################################################################################################################################
####   make the model 
model_two <-function(num_lags,laglag,pred_length,train_time,data_B,data_p1,data_p2,data_p3){
  # Create lagged predictors for each mutation site
  # set up the new_data set
  columns <- c("y",paste0("mutation_rate_lag_", 1:num_lags))
  new_data_B <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B) <- columns
  #    all mutations and some time as the training data set
  train_mu = ncol(data_B)
  for (i in 1:train_mu) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_B[j,i], data_B[(j-num_lags):(j-1),i])
      new_data_B[nrow(new_data_B) + 1,] <- newline
    }
  }
  #########################################################################################################   training data pattern1
  new_data_B1 <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B1) <- columns
  #    all mutations and some time as the training data set
  for (i in 1:min(train_mu,ncol(data_p1))) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_p1[j,i], data_p1[(j-num_lags):(j-1),i])
      new_data_B1[nrow(new_data_B1) + 1,] <- newline
    }
  }
  #########################################################################################################   training data pattern2
  new_data_B2 <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B2) <- columns
  #    all mutations and some time as the training data set
  for (i in 1:min(train_mu,ncol(data_p2))) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_p2[j,i], data_p2[(j-num_lags):(j-1),i])
      new_data_B2[nrow(new_data_B2) + 1,] <- newline
    }
  }
  #########################################################################################################   training data pattern3
  new_data_B3 <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B3) <- columns
  #    all mutations and some time as the training data set
  for (i in 1:min(train_mu,ncol(data_p3))) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_p3[j,i], data_p3[(j-num_lags):(j-1),i])
      new_data_B3[nrow(new_data_B3) + 1,] <- newline
    }
  }
  model_one_train = list()
  model_one_train[[1]] = new_data_B
  model_one_train[[2]] = new_data_B1
  model_one_train[[3]] = new_data_B2
  model_one_train[[4]] = new_data_B3
  names(model_one_train) = c('all','pattern1','pattern2','pattern3')
  return(model_one_train)
}

##################################################################################################################################################################
######    first,    use  the changed data  functions
## get mse,mae + plot together!!!!!!


test_0926_changed <- function(start, site,pred_length,num_lags,nn,test_pattern,test_pattern_minus,model,train_pattern,datac,date,plots=FALSE) { # function to run testing
  # start--> the position where test data begin to predict
  # site --> which mutation in test data
  # pred_length --> how may days u want to predict 
  # num_lags --> the number of days which we use to predict the next (usually 21 here, show be match the training data)
  # nn  --> mutation name
  # test_pattern --> processed data
  # test_pattern_minus --> just test part use true minus
  # model --> different model input here, default use LSTM
  # plots --> plot or not
  
  ############################################################################## loading time series using data first
  days_true = pred_length
  
  data=cbind(date,datac)
  data = as.data.frame(data)
  s2 <- data[, !duplicated(t(data))]
  data = s2
  j=1
  
  pred_true_matrix_list = list()
  
  
  
  ###############################################################################   only need to change here  , how many days u want to predict
  
  # result_date <- days_before_date(input_date, days_true)
  result_date <-data[start,1]
  print(result_date)
  
  ###Interval time
  I_time = abs(as.numeric(difftime(data[1,1], data[2,1], units = "days"))) ###Interval time
  days = ceiling(days_true/I_time)
  day_tag = paste(days_true,sep = " ",'days')
  DD = data[,c(1,which(colnames(data)%in%nn))]
  DD = as.data.frame(DD)
  name = colnames(DD)[2]
  colnames(DD)[2] = c('M')
  colnames(DD)[1] = c('Date')
  
  max_value <- max(DD$M)
  min_value <- min(DD$M)
  spread <- max_value - min_value
  
  #dataset <- (DD$M - min_value) / spread
  dataset = DD$M
  ############################################################################################################################
  ########################################################################################################################
  ##################################################################    replace place  ############################
  train_size = which(data[,1]==result_date)+pred_length  ###################################### change the time
  test_size <- length(dataset) - train_size
  
  train <- dataset[1:(train_size-1)]
  test <- dataset[(train_size):c(train_size+pred_length-1)]
  
  # cat(length(train), length(test))
  look_back <- length(test)
  ddf <- data.frame(ds = DD$Date[1:(train_size-1)],y = as.numeric(train))
  skirtsts<-ts(data[,which(colnames(data)%in%nn)],frequency=365)
  
  
  
  if(class(model)[1]=="character"){
    if(model=="LSTM"){
      ###################################   Lstm Mechine
      library(keras)
      model <- keras_model_sequential()
      model %>%
        layer_lstm(units = 5, input_shape = c(num_lags, 1)) %>%
        layer_dense(units = 1)
      model %>% compile(
        loss = 'mean_squared_error',
        optimizer = optimizer_adam()
      )
      train_data = as.matrix(train_pattern[,-1])
      train_labels = train_pattern[,1]
      # ??????????????????????????????LSTM???????????????
      train_data <- array(train_data, dim = c(dim(train_data)[1], dim(train_data)[2], 1))
      # ????????????
      history <- model %>% fit(
        train_data, train_labels,
        epochs = 10,
        batch_size = 128,
        verbose = 2
      )
      test_data <- test_pattern[start:(start+(num_lags-1)),site]
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]
      # test_data <- test_pattern[,site]
      test_data <- array(test_data, dim = c(1, num_lags, 1))
      test_data_minus <- array(test_data_minus, dim = c(1, num_lags, 1))
      
      # test_data =t(as.data.frame(test_data))
      pred_data <- c()
      pred_data_com<-c()
      for (j in c(1:pred_length)) {
        predictions <- model %>% predict(test_data)
        test_data[,1:num_lags,] = c(test_data[,2:num_lags,],predictions)
        # print(predictions)
        pred_data <- c(pred_data,predictions)
        pred_data_com <- c(pred_data_com,(test_data_minus[num_lags]+predictions))
        test_data_minus[,1:num_lags,] = c(test_data_minus[,2:num_lags,],(test_data_minus[num_lags]+predictions))
        
      }
    }
    if(model=="ARI"){
      
      
      # skirtsarima<-auto.arima(skirtsts[1:length(train)],trace=F)    ############-----------  better
      #skirtsarima <- arima(skirtsts[1:length(train)], order=c(1,0,0))  # AR(1) model
      # skirtsarima <- auto.arima(skirtsts[1:length(train)], seasonal=FALSE, trace=F)
      # skirtsarima <- auto.arima(skirtsts[1:length(train)], seasonal=FALSE, trace=F, stepwise=FALSE)
      skirtsarima <- auto.arima(skirtsts[1:length(train)], seasonal=FALSE, trace=F, d=0, D=0)
      
      
      
      skirtsarimaforecast<-forecast(skirtsarima,h=days,level=c(99.5))
      a1a1 = skirtsarima[["fitted"]]
      a2a2 = as.double(skirtsarimaforecast$mean)
      pred_data_com=a2a2
    }
    if(model=="PRO"){
      
      prophet_pred = prophet(ddf)
      future = make_future_dataframe(prophet_pred,periods=days)
      fcastprophet = predict(prophet_pred,future)
      #Creating train prediction dataset to compare real data
      dataprediction1 = data.frame(fcastprophet$ds,fcastprophet$yhat)
      trainlen = length(train)
      dataprediction = dataprediction1[c(1:trainlen),]
      
      c1c1 = dataprediction1[c(1:trainlen),2]
      c2c2 = dataprediction1[-c(1:trainlen),2]
      pred_data_com=c2c2
    }
    if(model=="KNN"){
      df_knn <- ddf
      # predknn <- knn_forecasting(df_knn$y, h = days, lags = 1:30, k = 50, msas = "MIMO")
      predknn <- knn_forecasting(df_knn$y, h = days, lags = 1:21, k = 21, msas = "MIMO")
      d1d1 = c(predknn[["model"]][["ts"]])
      d2d2 = c(predknn[["prediction"]])
      pred_data_com=d2d2
    }
    if(model=="FNN_1"){
      
      lambda = BoxCox.lambda(train)
      # dnn_fit = nnetar(train,lambda=lambda,size=10, repeats=10)
      dnn_fit = nnetar(train,lambda=lambda,size=5, repeats=2)
      fcast = forecast(dnn_fit,PI=T,h=days)
      #autoplot(fcast)
      f1f1 = c(fcast[["fitted"]])
      f2f2 = c(fcast[["mean"]])
      pred_data_com=f2f2
    }
    if(model=="LSTM_1"){
      minn  = -1
      while (minn<0) {
        ##################################################################    LSTM
        iteration  =50
        future_length  =days_true
        rr = my_LSTM(train,look_back,iteration,future_length)
        h1h1 = rr[[1]]
        h2h2 = rr[[2]]
        pred_data_com = h2h2
        minn = min(h2h2)
      }
      
    }
    if(model=="B_S"){
      x <- c(1:length(train))
      x2.new <- c((length(train)+1):(length(train)+pred_length))
      y2 <- skirtsts[1:length(train)]
      n0=c(5,10,15,20)
      fit.temp <- lm(y2 ~ bs(x, df=n0))
      j1j1 = fit.temp[["fitted.values"]]
      j2j2 = predict(fit.temp, newdata=data.frame(x=x2.new))
      pred_data_com=as.matrix(j2j2)
      
    }
    if(model=="C_S"){
      x <- c(1:length(train))
      x.new <- c((length(train)+1):(length(train)+pred_length))
      y <- skirtsts[1:length(train)]
      
      # My n.s fit:
      n0=c(5,10,15,20)
      fit.temp <- lm(y ~ ns(x, df=n0))
      i1i1 = fit.temp[["fitted.values"]]
      i2i2 = predict(fit.temp, newdata=data.frame(x=x.new))
      pred_data_com=as.matrix(i2i2)
      
    }
    
    
  }else{
    
    if(class(model)[1] =="randomForest"){
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # test_data <- test_pattern[start:(start+(num_lags-1)),site]
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        new_pred <- predict(model,test_data)      # that is changed ponits
        pred_data <- c(pred_data,new_pred)        #
        test_data <- c(test_data[-1], new_pred)
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        test_data_minus <- c(test_data_minus[-1], (test_data_minus[length(test_data_minus)]+new_pred))
        
      }
    }else{
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # plot(c(1:21),pred_data_com,type = 'line')
      test_data =t(as.data.frame(test_data))
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        if(class(model)[1]=='lm'){
          names_c = names(model[["coefficients"]])[-1]
          colnames(test_data)= names_c
          new_pred <- predict(model,as.data.frame(test_data))
        }else{
          new_pred <- predict(model,as.matrix(test_data))
        }
        
        pred_data <- c(pred_data,new_pred)
        test_data[1:(num_lags-1)] <- test_data[2:num_lags]
        test_data[num_lags] = new_pred
        
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        #test_pattern_minus <- c(test_pattern_minus[-1], (test_pattern_minus[length(test_pattern_minus)]+new_pred))
        
        test_data_minus[1:(num_lags-1)] <- test_data_minus[2:num_lags]
        test_data_minus[num_lags] = (test_data_minus[length(test_data_minus)]+new_pred)
        
      }
      
    }
  }
  score = matrix(data=NA,nrow = pred_length,ncol = 2)
  t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site]+test_pattern_minus[(start+num_lags):(start+(num_lags-1)+pred_length),site]
  
  # check accuracy
  
  rmse <- rmse(pred_data_com, t_true)
  # rse <- rse(pred_data, t_true)
  mse <- mean((pred_data_com - t_true)^2)
  # 
  mae <- mean(abs(pred_data_com - t_true))
  # 
  resi = sum(abs(pred_data_com - t_true))
  score_mse <- c(mse,mae)
  names(score_mse) = c('mse','mae')
  
  
  
  
  
  
  score[,1]=t_true
  # print(length(pred_data))
  score[,2]=pred_data_com
  if(plots==TRUE){
    # Plot actual vs. predicted mutation rate
    library(ggplot2)
    x <- seq(1:pred_length)+start+num_lags-1
    plot(seq(1:nrow(test_pattern)),test_pattern[,site]+test_pattern_minus[,site],type="l",col="#00A497", xlab = nn,ylab="mutation frequency",lwd=2)
    
    lines(x,t_true,col="#0095D9",lwd=2)
    lines(x,pred_data_com,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#00A497", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    ylim(c(min(c(t_true,pred_data_com)), max(c(t_true,pred_data_com))))
    plot(x,t_true,type="l",col="#0095D9", xlab = nn,ylab="mutation frequency",lwd=2)
    lines(x,pred_data_com,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#0095D9", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    
  }
  result = list()
  
  result[[1]]=score
  result[[2]]=score_mse
  names(result)[1]=c("score")
  names(result)[2]=c("score_mse")
  return(result)
}



test_0926_raw <- function(start, site,pred_length,num_lags,nn,test_pattern,model,train_pattern,datac,date,plots=FALSE) { # function to run testing
  # start--> the position where test data begin to predict
  # site --> which mutation in test data
  # pred_length --> how may days u want to predict 
  # num_lags --> the number of days which we use to predict the next (usually 21 here, show be match the training data)
  # nn  --> mutation name
  # test_pattern --> processed data
  # model --> different model input here, default use LSTM
  # plots --> plot or not
  
  ############################################################################## loading time series using data first
  
  if(class(model)[1]=="character"){
    if(model=="LSTM"){
      ###################################   Lstm Mechine
      library(keras)
      model <- keras_model_sequential()
      model %>%
        layer_lstm(units = 5, input_shape = c(num_lags, 1)) %>%
        layer_dense(units = 1)
      model %>% compile(
        loss = 'mean_squared_error',
        optimizer = optimizer_adam()
      )
      train_data = as.matrix(train_pattern[,-1])
      train_labels = train_pattern[,1]
      # ??????????????????????????????LSTM???????????????
      train_data <- array(train_data, dim = c(dim(train_data)[1], dim(train_data)[2], 1))
      # ????????????
      history <- model %>% fit(
        train_data, train_labels,
        epochs = 10,
        batch_size = 128,
        verbose = 2
      )
      test_data <- test_pattern[start:(start+(num_lags-1)),site]
      #test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]
      # test_data <- test_pattern[,site]
      test_data <- array(test_data, dim = c(1, num_lags, 1))
      #test_data_minus <- array(test_data_minus, dim = c(1, num_lags, 1))
      
      # test_data =t(as.data.frame(test_data))
      pred_data <- c()
      #pred_data_com<-c()
      for (j in c(1:pred_length)) {
        predictions <- model %>% predict(test_data)
        test_data[,1:num_lags,] = c(test_data[,2:num_lags,],predictions)
        # print(predictions)
        pred_data <- c(pred_data,predictions)
        #pred_data_com <- c(pred_data_com,(test_data_minus[num_lags]+predictions))
        #test_data_minus[,1:num_lags,] = c(test_data_minus[,2:num_lags,],(test_data_minus[num_lags]+predictions))
        
      }
    }
    
    
  }else{
    
    if(class(model)[1] =="randomForest"){
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      #test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # test_data <- test_pattern[start:(start+(num_lags-1)),site]
      pred_data <- c()
      #pred_data_com<-c()
      for (i in 1: pred_length) {
        new_pred <- predict(model,test_data)      # that is changed ponits
        pred_data <- c(pred_data,new_pred)        #
        test_data <- c(test_data[-1], new_pred)
        #pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        #test_data_minus <- c(test_data_minus[-1], (test_data_minus[length(test_data_minus)]+new_pred))
        
      }
    }else{
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      #test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # plot(c(1:21),pred_data_com,type = 'line')
      test_data =t(as.data.frame(test_data))
      pred_data <- c()
      #pred_data_com<-c()
      for (i in 1: pred_length) {
        if(class(model)[1]=='lm'){
          names_c = names(model[["coefficients"]])[-1]
          colnames(test_data)= names_c
          new_pred <- predict(model,as.data.frame(test_data))
        }else{
          new_pred <- predict(model,as.matrix(test_data))
        }
        
        pred_data <- c(pred_data,new_pred)
        test_data[1:(num_lags-1)] <- test_data[2:num_lags]
        test_data[num_lags] = new_pred
        
        #pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        #test_pattern_minus <- c(test_pattern_minus[-1], (test_pattern_minus[length(test_pattern_minus)]+new_pred))
        
        #test_data_minus[1:(num_lags-1)] <- test_data_minus[2:num_lags]
        #test_data_minus[num_lags] = (test_data_minus[length(test_data_minus)]+new_pred)
        
      }
      
    }
  }
  score = matrix(data=NA,nrow = pred_length,ncol = 2)
  t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site]#+test_pattern_minus[(start+num_lags):(start+(num_lags-1)+pred_length),site]
  
  # check accuracy
  
  rmse <- rmse(pred_data, t_true)
  # rse <- rse(pred_data, t_true)
  mse <- mean((pred_data - t_true)^2)
  # 
  mae <- mean(abs(pred_data - t_true))
  # 
  resi = sum(abs(pred_data - t_true))
  score_mse <- c(rmse,mse,mae,resi)
  names(score_mse) = c('rmse','mse','mae','resi')
  
  
  
  
  
  
  score[,1]=t_true
  # print(length(pred_data))
  score[,2]=pred_data
  if(plots==TRUE){
    # Plot actual vs. predicted mutation rate
    library(ggplot2)
    x <- seq(1:pred_length)+start+num_lags-1
    plot(seq(1:nrow(test_pattern)),test_pattern[,site]+test_pattern_minus[,site],type="l",col="#00A497", xlab = nn,ylab="mutation frequency",lwd=2)
    
    lines(x,t_true,col="#0095D9",lwd=2)
    lines(x,pred_data,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#00A497", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    ylim(c(min(c(t_true,pred_data)), max(c(t_true,pred_data))))
    plot(x,t_true,type="l",col="#0095D9", xlab = nn,ylab="mutation frequency",lwd=2)
    lines(x,pred_data,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#0095D9", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    
  }
  result = list()
  
  result[[1]]=score
  result[[2]]=score_mse
  names(result)[1]=c("score")
  names(result)[2]=c("score_mse")
  return(result)
}







##################################################################################################################################################################
########    second,  use the raw data   functions

model_two <-function(num_lags,laglag,pred_length,train_time,data_B,data_p1,data_p2,data_p3){
  # Create lagged predictors for each mutation site
  # set up the new_data set
  columns <- c("y",paste0("mutation_rate_lag_", 1:num_lags))
  new_data_B <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B) <- columns
  #    all mutations and some time as the training data set
  train_mu = ncol(data_B)
  for (i in 1:train_mu) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_B[j,i], data_B[(j-num_lags):(j-1),i])
      new_data_B[nrow(new_data_B) + 1,] <- newline
    }
  }
  #########################################################################################################   training data pattern1
  new_data_B1 <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B1) <- columns
  #    all mutations and some time as the training data set
  for (i in 1:min(train_mu,ncol(data_p1))) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_p1[j,i], data_p1[(j-num_lags):(j-1),i])
      new_data_B1[nrow(new_data_B1) + 1,] <- newline
    }
  }
  #########################################################################################################   training data pattern2
  new_data_B2 <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B2) <- columns
  #    all mutations and some time as the training data set
  for (i in 1:min(train_mu,ncol(data_p2))) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_p2[j,i], data_p2[(j-num_lags):(j-1),i])
      new_data_B2[nrow(new_data_B2) + 1,] <- newline
    }
  }
  #########################################################################################################   training data pattern3
  new_data_B3 <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(new_data_B3) <- columns
  #    all mutations and some time as the training data set
  for (i in 1:min(train_mu,ncol(data_p3))) {
    for (j in seq((num_lags+1), train_time, by = laglag)) {
      #for (j in seq((num_lags+1), 300, by = laglag)) {
      newline <- c(data_p3[j,i], data_p3[(j-num_lags):(j-1),i])
      new_data_B3[nrow(new_data_B3) + 1,] <- newline
    }
  }
  model_one_train = list()
  model_one_train[[1]] = new_data_B
  model_one_train[[2]] = new_data_B1
  model_one_train[[3]] = new_data_B2
  model_one_train[[4]] = new_data_B3
  names(model_one_train) = c('all','pattern1','pattern2','pattern3')
  return(model_one_train)
}





test_no_changed <- function(start, site,pred_length,num_lags,nn,test_pattern,model='LSTM',train_pattern,plots=FALSE) { # function to run testing
  # start--> the position where test data begin to predict
  # site --> which mutation in test data
  # pred_length --> how may days u want to predict 
  # num_lags --> the number of days which we use to predict the next (usually 21 here, show be match the training data)
  # nn  --> mutation name
  # test_pattern --> processed data
  # model --> different model input here, default use LSTM
  # plots --> plot or not
  if(class(model)=="character"){
    ###################################   Lstm
    library(keras)
    model <- keras_model_sequential()
    model %>%
      layer_lstm(units = 5, input_shape = c(num_lags, 1)) %>%
      layer_dense(units = 1)
    model %>% compile(
      loss = 'mean_squared_error',
      optimizer = optimizer_adam()
    )
    train_data = as.matrix(train_pattern[,-1])
    train_labels = train_pattern[,1]
    # 将训练数据重塑为适合LSTM的输入形状
    train_data <- array(train_data, dim = c(dim(train_data)[1], dim(train_data)[2], 1))
    # 训练模型
    history <- model %>% fit(
      train_data, train_labels,
      epochs = 10,
      batch_size = 128,
      verbose = 2
    )
    test_data <- test_pattern[start:(start+(num_lags-1)),site]
    # test_data <- test_pattern[,site]
    test_data <- array(test_data, dim = c(1, num_lags, 1))
    # test_data =t(as.data.frame(test_data))
    pred_data <- c()
    pred_data=c()
    for (j in c(1:pred_length)) {
      predictions <- model %>% predict(test_data)
      test_data[,1:num_lags,] = c(test_data[,2:num_lags,],predictions)
      # print(predictions)
      pred_data <- c(pred_data,predictions)
    }
    
  }else{
    test_data <- test_pattern[start:(start+(num_lags-1)),site]
    if(class(model) =="randomForest"){
      # test_data <- test_pattern[start:(start+(num_lags-1)),site]
      pred_data <- c()
      for (i in 1: pred_length) {
        new_pred <- predict(model,test_data)
        pred_data <- c(pred_data,new_pred)
        test_data <- c(test_data[-1], new_pred)
      }
    }else{
      
      # test_data <- test_pattern[,site]
      test_data =t(as.data.frame(test_data))
      pred_data <- c()
      for (i in 1: pred_length) {
        if(class(model)=='lm'){
          names_c = names(model[["coefficients"]])[-1]
          colnames(test_data)= names_c
          new_pred <- predict(model,as.data.frame(test_data))
        }else{
          new_pred <- predict(model,as.matrix(test_data))
        }
        
        pred_data <- c(pred_data,new_pred)
        test_data[1:(num_lags-1)] <- test_data[2:num_lags]
        test_data[num_lags] = new_pred
      }
      
    }
  }
  # check accuracy
  library(Metrics)
  score = matrix(data=NA,nrow = pred_length,ncol = 2)
  t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site]
  score[,1]=t_true
  print(length(pred_data))
  score[,2]=pred_data
  if(plots==TRUE){
    # Plot actual vs. predicted mutation rate
    library(ggplot2)
    x <- seq(1:pred_length)+start+num_lags-1
    plot(seq(1:nrow(test_pattern)),test_pattern[,site],type="l",col="#00A497", xlab = nn,ylab="mutation frequency",lwd=2)
    
    lines(x,t_true,col="#0095D9",lwd=2)
    lines(x,pred_data,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#00A497", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    ylim(c(min(c(t_true,pred_data)), max(c(t_true,pred_data))))
    plot(x,t_true,type="l",col="#0095D9", xlab = nn,ylab="mutation frequency",lwd=2)
    lines(x,pred_data,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#0095D9", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    
  }
  return(score)
}

#############################################################


test_new_changed <- function(start, site,pred_length,num_lags,nn,test_pattern,test_pattern_minus,model='LSTM',train_pattern,plots=FALSE) { # function to run testing
  # start--> the position where test data begin to predict
  # site --> which mutation in test data
  # pred_length --> how may days u want to predict 
  # num_lags --> the number of days which we use to predict the next (usually 21 here, show be match the training data)
  # nn  --> mutation name
  # test_pattern --> processed data
  # test_pattern_minus --> just test part use true minus
  # model --> different model input here, default use LSTM
  # plots --> plot or not
  if(class(model)=="character"){
    ###################################   Lstm
    library(keras)
    model <- keras_model_sequential()
    model %>%
      layer_lstm(units = 5, input_shape = c(num_lags, 1)) %>%
      layer_dense(units = 1)
    model %>% compile(
      loss = 'mean_squared_error',
      optimizer = optimizer_adam()
    )
    train_data = as.matrix(train_pattern[,-1])
    train_labels = train_pattern[,1]
    # 将训练数据重塑为适合LSTM的输入形状
    train_data <- array(train_data, dim = c(dim(train_data)[1], dim(train_data)[2], 1))
    # 训练模型
    history <- model %>% fit(
      train_data, train_labels,
      epochs = 10,
      batch_size = 128,
      verbose = 2
    )
    test_data <- test_pattern[start:(start+(num_lags-1)),site]
    test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]
    # test_data <- test_pattern[,site]
    test_data <- array(test_data, dim = c(1, num_lags, 1))
    test_data_minus <- array(test_data_minus, dim = c(1, num_lags, 1))
    
    # test_data =t(as.data.frame(test_data))
    pred_data <- c()
    pred_data_com<-c()
    for (j in c(1:pred_length)) {
      predictions <- model %>% predict(test_data)
      test_data[,1:num_lags,] = c(test_data[,2:num_lags,],predictions)
      # print(predictions)
      pred_data <- c(pred_data,predictions)
      pred_data_com <- c(pred_data_com,(test_data_minus[num_lags]+predictions))
      test_data_minus[,1:num_lags,] = c(test_data_minus[,2:num_lags,],(test_data_minus[num_lags]+predictions))
      
    }
    
  }else{
    
    if(class(model) =="randomForest"){
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # test_data <- test_pattern[start:(start+(num_lags-1)),site]
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        new_pred <- predict(model,test_data)      # that is changed ponits
        pred_data <- c(pred_data,new_pred)        #
        test_data <- c(test_data[-1], new_pred)
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        test_data_minus <- c(test_data_minus[-1], (test_data_minus[length(test_data_minus)]+new_pred))
        
      }
    }else{
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # plot(c(1:21),pred_data_com,type = 'line')
      test_data =t(as.data.frame(test_data))
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        if(class(model)=='lm'){
          names_c = names(model[["coefficients"]])[-1]
          colnames(test_data)= names_c
          new_pred <- predict(model,as.data.frame(test_data))
        }else{
          new_pred <- predict(model,as.matrix(test_data))
        }
        
        pred_data <- c(pred_data,new_pred)
        test_data[1:(num_lags-1)] <- test_data[2:num_lags]
        test_data[num_lags] = new_pred
        
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        #test_pattern_minus <- c(test_pattern_minus[-1], (test_pattern_minus[length(test_pattern_minus)]+new_pred))
        
        test_data_minus[1:(num_lags-1)] <- test_data_minus[2:num_lags]
        test_data_minus[num_lags] = (test_data_minus[length(test_data_minus)]+new_pred)
        
      }
      
    }
  }
  # check accuracy
  library(Metrics)
  score = matrix(data=NA,nrow = pred_length,ncol = 2)
  t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site]+test_pattern_minus[(start+num_lags):(start+(num_lags-1)+pred_length),site]
  score[,1]=t_true
  print(length(pred_data))
  score[,2]=pred_data_com
  if(plots==TRUE){
    # Plot actual vs. predicted mutation rate
    library(ggplot2)
    x <- seq(1:pred_length)+start+num_lags-1
    plot(seq(1:nrow(test_pattern)),test_pattern[,site]+test_pattern_minus[,site],type="l",col="#00A497", xlab = nn,ylab="mutation frequency",lwd=2)
    
    lines(x,t_true,col="#0095D9",lwd=2)
    lines(x,pred_data_com,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#00A497", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    ylim(c(min(c(t_true,pred_data_com)), max(c(t_true,pred_data_com))))
    plot(x,t_true,type="l",col="#0095D9", xlab = nn,ylab="mutation frequency",lwd=2)
    lines(x,pred_data_com,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#0095D9", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    
  }
  return(score)
}

##########################################################################################################    with   mse

test_new_changed_mse <- function(start, site,pred_length,num_lags,nn,test_pattern,test_pattern_minus,model='LSTM',train_pattern,plots=FALSE) { # function to run testing
  # start--> the position where test data begin to predict
  # site --> which mutation in test data
  # pred_length --> how may days u want to predict 
  # num_lags --> the number of days which we use to predict the next (usually 21 here, show be match the training data)
  # nn  --> mutation name
  # test_pattern --> processed data
  # test_pattern_minus --> just test part use true minus
  # model --> different model input here, default use LSTM
  # plots --> plot or not
  if(class(model)=="character"){
    ###################################   Lstm
    library(keras)
    model <- keras_model_sequential()
    model %>%
      layer_lstm(units = 5, input_shape = c(num_lags, 1)) %>%
      layer_dense(units = 1)
    model %>% compile(
      loss = 'mean_squared_error',
      optimizer = optimizer_adam()
    )
    train_data = as.matrix(train_pattern[,-1])
    train_labels = train_pattern[,1]
    # 将训练数据重塑为适合LSTM的输入形状
    train_data <- array(train_data, dim = c(dim(train_data)[1], dim(train_data)[2], 1))
    # 训练模型
    history <- model %>% fit(
      train_data, train_labels,
      epochs = 10,
      batch_size = 128,
      verbose = 2
    )
    test_data <- test_pattern[start:(start+(num_lags-1)),site]
    test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]
    # test_data <- test_pattern[,site]
    test_data <- array(test_data, dim = c(1, num_lags, 1))
    test_data_minus <- array(test_data_minus, dim = c(1, num_lags, 1))
    
    # test_data =t(as.data.frame(test_data))
    pred_data <- c()
    pred_data_com<-c()
    for (j in c(1:pred_length)) {
      predictions <- model %>% predict(test_data)
      test_data[,1:num_lags,] = c(test_data[,2:num_lags,],predictions)
      # print(predictions)
      pred_data <- c(pred_data,predictions)
      pred_data_com <- c(pred_data_com,(test_data_minus[num_lags]+predictions))
      test_data_minus[,1:num_lags,] = c(test_data_minus[,2:num_lags,],(test_data_minus[num_lags]+predictions))
      
    }
    
  }else{
    
    if(class(model) =="randomForest"){
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # test_data <- test_pattern[start:(start+(num_lags-1)),site]
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        new_pred <- predict(model,test_data)      # that is changed ponits
        pred_data <- c(pred_data,new_pred)        #
        test_data <- c(test_data[-1], new_pred)
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        test_data_minus <- c(test_data_minus[-1], (test_data_minus[length(test_data_minus)]+new_pred))
        
      }
    }else{
      test_data <- test_pattern[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),site]    ### e.g.   1~21
      # plot(c(1:21),pred_data_com,type = 'line')
      test_data =t(as.data.frame(test_data))
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        if(class(model)=='lm'){
          names_c = names(model[["coefficients"]])[-1]
          colnames(test_data)= names_c
          new_pred <- predict(model,as.data.frame(test_data))
        }else{
          new_pred <- predict(model,as.matrix(test_data))
        }
        
        pred_data <- c(pred_data,new_pred)
        test_data[1:(num_lags-1)] <- test_data[2:num_lags]
        test_data[num_lags] = new_pred
        
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        #test_pattern_minus <- c(test_pattern_minus[-1], (test_pattern_minus[length(test_pattern_minus)]+new_pred))
        
        test_data_minus[1:(num_lags-1)] <- test_data_minus[2:num_lags]
        test_data_minus[num_lags] = (test_data_minus[length(test_data_minus)]+new_pred)
        
      }
      
    }
  }
  score = matrix(data=NA,nrow = pred_length,ncol = 2)
  t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site]+test_pattern_minus[(start+num_lags):(start+(num_lags-1)+pred_length),site]
  
  # check accuracy
  library(Metrics)
  rmse <- rmse(pred_data_com, t_true)
  # rse <- rse(pred_data, t_true)
  mse <- mean((pred_data_com - t_true)^2)
  # 
  mae <- mean(abs(pred_data_com - t_true))
  # 
  resi = sum(abs(pred_data_com - t_true))
  score_mse <- c(rmse,mse,mae,resi)
  names(score_mse) = c('rmse','mse','mae','resi')
  
  
  
  
  
  
  score[,1]=t_true
  print(length(pred_data))
  score[,2]=pred_data_com
  if(plots==TRUE){
    # Plot actual vs. predicted mutation rate
    library(ggplot2)
    x <- seq(1:pred_length)+start+num_lags-1
    plot(seq(1:nrow(test_pattern)),test_pattern[,site]+test_pattern_minus[,site],type="l",col="#00A497", xlab = nn,ylab="mutation frequency",lwd=2)
    
    lines(x,t_true,col="#0095D9",lwd=2)
    lines(x,pred_data_com,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#00A497", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    ylim(c(min(c(t_true,pred_data_com)), max(c(t_true,pred_data_com))))
    plot(x,t_true,type="l",col="#0095D9", xlab = nn,ylab="mutation frequency",lwd=2)
    lines(x,pred_data_com,col="#EA5506" ,lwd=2)
    # legend('topright', legend=c("Actual", "Predicted"), col=c("#0095D9", "#EA5506"), pch=19,cex=1.7,lwd = 1.5)
    
    
  }
  result = list()
  
  result[[1]]=score
  result[[2]]=score_mse
  names(result)[1]=c("score")
  names(result)[2]=c("score_mse")
  return(result)
}



break_down_test <- function(start,pred_length,num_lags,nn,test_pattern,test_pattern_minus,model='LSTM',train_pattern) { # function to run testing
  # start--> the position where test data begin to predict
  # site --> which mutation in test data
  # pred_length --> how may days u want to predict 
  # num_lags --> the number of days which we use to predict the next (usually 21 here, show be match the training data)
  # nn  --> mutation name
  # test_pattern --> processed data
  # test_pattern_minus -->data add back
  # model --> different model input here, default use LSTM
  LL=0
  LM=0
  LH=0
  ML=0
  MM=0
  MH=0
  HL=0
  HM=0
  HH=0
  if(class(model)=="character"){
    ###################################   Lstm
    library(keras)
    model <- keras_model_sequential()
    model %>%
      layer_lstm(units = 5, input_shape = c(num_lags, 1)) %>%
      layer_dense(units = 1)
    model %>% compile(
      loss = 'mean_squared_error',
      optimizer = optimizer_adam()
    )
    train_data = as.matrix(train_pattern[,-1])
    train_labels = train_pattern[,1]
    # 将训练数据重塑为适合LSTM的输入形状
    train_data <- array(train_data, dim = c(dim(train_data)[1], dim(train_data)[2], 1))
    # 训练模型
    history <- model %>% fit(
      train_data, train_labels,
      epochs = 5,
      batch_size = 128,
      verbose = 2
    )
    
    # test_data <- test_pattern[(start-20):start]#[start:(start+(num_lags-1)),site]
    # test_data <- array(test_data, dim = c(1, num_lags, 1))
    # 
    test_data <- test_pattern[start:(start+(num_lags-1))]
    test_data_minus<-test_pattern_minus[start:(start+(num_lags-1))]
    # test_data <- test_pattern[,site]
    test_data <- array(test_data, dim = c(1, num_lags, 1))
    test_data_minus <- array(test_data_minus, dim = c(1, num_lags, 1))
    
    # test_data =t(as.data.frame(test_data))
    pred_data <- c()
    pred_data_com<-c()
    for (j in c(1:num_lags)) {
      predictions <- model %>% predict(test_data)
      test_data[,1:num_lags,] = c(test_data[,2:num_lags,],predictions)
      print(predictions)
      pred_data <- c(pred_data,predictions)
      
      pred_data_com <- c(pred_data_com,(test_data_minus[num_lags]+predictions))
      test_data_minus[,1:num_lags,] = c(test_data_minus[,2:num_lags,],(test_data_minus[num_lags]+predictions))
      
    }
    
  }else{
    # par(mfrow=c(1,2))
    if(class(model) =="randomForest"){
      test_data <- test_pattern[start:(start+(num_lags-1))]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1))]    ### e.g.   1~21
      # test_data <- test_pattern[start:(start+(num_lags-1)),site]
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        new_pred <- predict(model,test_data)      # that is changed ponits
        pred_data <- c(pred_data,new_pred)        #
        test_data <- c(test_data[-1], new_pred)
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        test_data_minus <- c(test_data_minus[-1], (test_data_minus[length(test_data_minus)]+new_pred))
        
      }
    }else{
      test_data <- test_pattern[start:(start+(num_lags-1))]    ### e.g.   1~21
      test_data_minus<-test_pattern_minus[start:(start+(num_lags-1))]    ### e.g.   1~21
      # plot(c(1:21),pred_data_com,type = 'line')
      test_data =t(as.data.frame(test_data))
      pred_data <- c()
      pred_data_com<-c()
      for (i in 1: pred_length) {
        if(class(model)=='lm'){
          names_c = names(model[["coefficients"]])[-1]
          colnames(test_data)= names_c
          new_pred <- predict(model,as.data.frame(test_data))
        }else{
          new_pred <- predict(model,as.matrix(test_data))
        }
        
        pred_data <- c(pred_data,new_pred)
        test_data[1:(num_lags-1)] <- test_data[2:num_lags]
        test_data[num_lags] = new_pred
        
        pred_data_com<-c(pred_data_com,(test_data_minus[length(test_data_minus)]+new_pred))   #  that is the true what we predict data
        #test_pattern_minus <- c(test_pattern_minus[-1], (test_pattern_minus[length(test_pattern_minus)]+new_pred))
        
        test_data_minus[1:(num_lags-1)] <- test_data_minus[2:num_lags]
        test_data_minus[num_lags] = (test_data_minus[length(test_data_minus)]+new_pred)
        
      }
      
      
      
    }
  }
  
  p_result = pred_data_com
  t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site]+test_pattern_minus[(start+num_lags):(start+(num_lags-1)+pred_length),site]
  #########################
  #########################
  if(t_true[length(t_true)]<=0.01){
    if(p_result[length(p_result)]<=0.01){
      LL=LL+1
    }else if(p_result[length(p_result)]>0.1){
      LH=LH+1
    }else{
      LM=LM+1
    }
  }else if(t_true[length(t_true)]>0.1){
    if(p_result[length(p_result)]<=0.01){
      HL=HL+1
    }else if(p_result[length(p_result)]>0.1){
      HH=HH+1
    }else{
      HM=HM+1
    }
  }else{
    if(p_result[length(p_result)]<=0.01){
      ML=ML+1
    }else if(p_result[length(p_result)]>0.1){
      MH=MH+1
    }else{
      MM=MM+1
    }
  }
  # # check accuracy
  # library(Metrics)
  # rmse <- rmse(pred_data, test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site])
  # rse <- rse(pred_data, test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site])
  # mse <- mean((pred_data - test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site])^2)
  # 
  # mae <- mean(abs(pred_data - test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),site]))
  # 
  # score <- c(rmse, mse,mae)
  
  
  conf_matrix = matrix(0,nrow = 3,ncol=3)
  rownames(conf_matrix)=c('true_low','true_mid','true_high')
  colnames(conf_matrix)=c('pred_low','pred_mid','pred_high')
  conf_matrix[1,]=c(LL,LM,LH)
  conf_matrix[2,]=c(ML,MM,MH)
  conf_matrix[3,]=c(HL,HM,HH)
  
  return(conf_matrix)
}









###   new   changed   compare

library(ggplot2)
library(patchwork)
library(keras)
library(ggthemes)
library(quantmod)
library(forecast)
library(KernSmooth)  
# library(locpol)  
library(locfit) 
require(stats); require(graphics)
library(splines)
library(keras)
library(quantmod)
library(tseries)
#library(rugarch)
library(prophet)
library(tsfknn)
##########################################################    plot all [(  mechine learning +  time series )] together, for per mutation 