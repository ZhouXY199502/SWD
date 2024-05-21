#### some small function

# R/days_before_date.R

#' days_before_date
#'
#' This function sets a pattern.
#'
#' @param x The input parameter.
#' @return The pattern.
#' @export
days_before_date <- function(input_date, days_before) {

  date_obj <- as.Date(input_date)

  result_date <- date_obj - days_before

  return(format(result_date, "%Y-%m-%d"))
}


####  loading and making data
# R/loading.R

#' loading
#'
#' This function sets a pattern.
#'
#' @param x The input parameter.
#' @return The pattern.
#' @export
loading<-function(path,name){
  # path = '/Users/xinyuzhou/Library/CloudStorage/OneDrive-IndianaUniversity/shared/covid/github_packages/'
  data <- read_tsv(paste(path,sep = '',name))
  colnames(data)[1]=c('date')
  # data <- data[, !duplicated(t(data))]
  date = data[,1]
  ### datac without non-zero
  datac <- data.frame(matrix(nrow = nrow(data), ncol = 0))
  max <- c()
  index <- c() # index of non-zero columns
  t=1
  for (i in 2:ncol(data)) {
    if (sum(data[,i] >= 0.0001)) { # only keep those columns with non-zeros
      datac <- cbind(datac,data[,i])
      max <- c(max,max(data[,i]))
      index <- c(index,i)
      colnames(datac)[t]=colnames(data)[i]
      t=t+1
    }
  }
  datac2_m = rbind(rep(0,ncol(datac)),datac[-nrow(datac),])
  data_B=datac-datac2_m            #    changed data
  results = list(data_B,datac2_m,datac,data,date)
  names(results)=c("Changed","Minus","Raw","Original","date")
  return(results)
}

#####threshold，default=0.05
# R/set_pattern.R

#' set_pattern
#'
#' This function sets a pattern.
#'
#' @param x The input parameter.
#' @return The pattern.
#' @export
set_pattern<-function(data_raw,data_changed,data_ori,data_minus,threshold=0.05,cut_time="2021-01-01"){
  index_0.05 <- as.numeric(which(apply(data_raw, 2, max)<threshold))            # pattern one, maximum less than 0.05
  data_p1=data_changed[,index_0.05]
  index_big <- as.numeric(which(apply(data_raw, 2, max)>=threshold))
  data_pre=data_raw[,index_big]
  index_date <- max(which(as.POSIXct(data_ori$date) < as.POSIXct(cut_time)))
  idx <- as.numeric(apply(data_pre >= threshold, 2, which.max))
  index_before = index_big[which(idx<index_date)]
  index_after = index_big[which(idx>=index_date)]
  data_p2 = data_changed[,index_before]                                    # pattern two, peak happend at a earlier time
  data_p3 = data_changed[,index_after]                                     # pattern three, peak happend at a later time
  rm(data_pre)

  # be minused data
  data_m1 = data_minus[,index_0.05]         #pattern 1
  data_m2 = data_minus[,index_before]       #pattern 2
  data_m3 = data_minus[,index_after]        #pattern 3
  data_p1_raw=data_raw[,index_0.05]                                        # pattern one, maximum less than 0.05
  data_p2_raw = data_raw[,index_before]                                    # pattern two, peak happend at a earlier time
  data_p3_raw = data_raw[,index_after]
  results = list(data_p1,data_m1,data_p2,data_m2,data_p3,data_m3,data_p1_raw,data_p2_raw,data_p3_raw)
  names(results)=c("pattern1_changed","pattern1_minus","pattern2_changed","pattern2_minus","pattern3_changed","pattern3_minus",
                   "data_p1_raw","data_p2_raw","data_p3_raw")
  return(results)
}


#######################
# num_lags=21        # slide
# laglag=1
# pred_length=30     # dates what to predict
# timee = c("2022-09-01")       #  before that day, regurd as train data
# meth_name = c("RF","S_L","S_R","XG","LR","FNN","ARI","PRO","KNN","C_S","B_S")
# ## data should be from  "loading"
# ############----------------  which mutation you want to test   ********
# removed_names=c('C_28311_T',"AG_28878_TC")        # removed mutation name


# R/get_train_matrix.R

#' get_train_matrix
#'
#' This function sets a pattern.
#'
#' @param x The input parameter.
#' @return The pattern.
#' @export
get_train_matrix<-function(num_lags,laglag,pred_length,timee,data,removed_names,data_pattern){
  train_time = max(which(as.POSIXct(data$Original$date) < as.POSIXct(timee)))
  removed_n = which(colnames(data$Raw)%in% removed_names)
  data_B_removed = data$Changed[,-removed_n]         # not inclued changed data
  train_matrix <-model_two(num_lags,laglag,pred_length,train_time,data_B_removed,data_pattern$pattern1_changed,
                                     data_pattern$pattern2_changed,data_pattern$pattern3_changed)
  return(train_matrix)
}









######  用来训练使用的 不同pattern and all的数据

#   使用时，查看raw里的data_B_removed应该使用什么？？
#model_two_data_changed_raw <-model_two(num_lags,laglag,pred_length,train_time,data_B_removed,data_p1_raw,data_p2_raw,data_p3_raw)



# save(model_two_data_changed,file = '/Users/xinyuzhou/Library/CloudStorage/OneDrive-IndianaUniversity/shared/covid/github_packages/model_two_data_all.RData')



###################   Going mechine learning training function    -------------------------------------
# MLE_meth_name_list = c("RF","SL","SR","XG","LR")   #,"FNN","ARI","PRO","KNN","C_S","B_S")

# R/Train_ML.R

#' Train_ML
#'
#' This function sets a pattern.
#'
#' @param x The input parameter.
#' @return The pattern.
#' @export
Train_ML<-function(train_patterns_list,path,MLE_meth_name,svm_type="eps-regression",XG_eta = 0.3,XG_max_depth = 6,
                   XG_eval_metric = "rmse",XG_num_rounds=300){
  all_pattern_trained=list()
  for (pattern_name in names(train_patterns_list)) {
    train_pattern <- train_patterns_list[[pattern_name]]
    ###########################    get train model

    # ######################################    Random Forest
    # model_list=list(rf_model,svm_L_model,svm_model_rbf,xgmodel,LL_model)
    model_list=list()
    for (mm in c(1:length(MLE_meth_name))) {
      if(MLE_meth_name[mm]=="RF"){
        system.time({
          model_list[[mm]] <- randomForest(x=train_pattern[,-1], y=train_pattern[,1])
        })
      }else if(MLE_meth_name[mm]=="SL"){
        system.time({
          model_list[[mm]] <- svm(y ~ ., data = train_pattern, type = svm_type,kernel = "linear")
        })
      }else if(MLE_meth_name[mm]=="SR"){
        system.time({
          model_list[[mm]] <- svm(y ~ ., data = train_pattern,type = svm_type, kernel = "radial")
        })
      }else if(MLE_meth_name[mm]=="XG"){
        X <- as.matrix(train_pattern[,-1])
        y <- as.matrix(train_pattern[,1])
        params <- list(
          objective = "reg:squarederror",  # 回归任务
          eta = XG_eta,
          max_depth = XG_max_depth,
          eval_metric = XG_eval_metric
        )
        num_rounds <- XG_num_rounds
        model_list[[mm]] <- xgboost(data = X, label = y, params = params, nrounds = num_rounds)

      }else if(MLE_meth_name[mm]=="LR"){
        model_list[[mm]] <- lm(y ~ ., data = train_pattern)
      }else if(MLE_meth_name[mm]=="FNN"){
        system.time({
          model_list[[mm]] <- neuralnet(y ~ ., data = train_pattern, hidden = c(7, 5), linear.output = TRUE)
        })
      }



    }
    names(model_list) = MLE_meth_name#,'nn1_model')
    all_pattern_trained[[pattern_name]]=model_list
    #########################################################################################################    get a model list
    # save_filename <- paste0(path, "/model_list_", pattern_name, ".RData")
    # save(model_list, file = save_filename)
  }
  return(all_pattern_trained)
}
######################-----------------------------------   Going to predict
# train_pattern =  train_matrix$all
#nn = "AG_28878_TC"
#start =880
# site = mu_num
# pred_length = 30
# num_lags =21
# # data_input
# date = data[,1]
# test_pattern = data_B
# test_pattern_minus = datac2_m
# position  =c(270,410,500,600,710,880,920,980,1000)
# datasum=data
# model_list kinds of mechine learning trained model inside

# R/Predicting.R

#' Predicting
#'
#' This function sets a pattern.
#'
#' @param x The input parameter.
#' @return The pattern.
#' @export
Predicting<-function(train_pattern,nn,pred_length,num_lags,datasum,position,model_list){
  date = datasum[["date"]]$date
  test_pattern = datasum$Changed
  test_pattern_minus = datasum$Minus
  dddd = position
  result_list = vector("list",length(dddd))
  names(result_list) <-as.character(dddd)
  ##### going to predict from different position in dddd
  for (i in c(1:length(dddd))) {
    start=dddd[i]
    # if want to plot
    # pdf_namesss <- paste0("F:/OneDrive - Indiana University/shared/covid/one_code_clean/a_final_1204/diff/2//", paste0(nn,start,"--222222-"), ".pdf")
    {
      ############################################################################################################################################################################
      ###########################################    this part for time serises      ###########################################
      ############################################################################################################################################################################
      #   time serises using the raw data
      days_true = pred_length
      data=cbind(date,datasum$Raw)
      data = as.data.frame(data)
      s2 <- data[, !duplicated(t(data))]
      data = s2
      j=1
      # MSE_list = list()
      # MSE_matrix1 = matrix(NA,nrow = 7,ncol=(length(data[1,])-1))
      #
      # rownames(MSE_matrix1)=c('Arima','Prophet','KNN Regression','Feed-Forward Neural Networks',"LSTM","B_spline","Cubic")
      # colnames(MSE_matrix1) = colnames(data)[2:length(data[1,])]
      #
      # MSE_matrix2 = matrix(NA,nrow = 7,ncol=(length(data[1,])-1))
      #
      # rownames(MSE_matrix2)=c('Arima','Prophet','KNN Regression','Feed-Forward Neural Networks',"LSTM","B_spline","Cubic")
      # colnames(MSE_matrix2) = colnames(data)[2:length(data[1,])]
      # MSE_matrix3 = matrix(NA,nrow = 7,ncol=(length(data[1,])-1))
      #
      # rownames(MSE_matrix3)=c('Arima','Prophet','KNN Regression','Feed-Forward Neural Networks',"LSTM","B_spline","Cubic")
      # colnames(MSE_matrix3) = colnames(data)[2:length(data[1,])]
      # MSE_matrix4 = matrix(NA,nrow = 7,ncol=(length(data[1,])-1))
      #
      # rownames(MSE_matrix4)=c('Arima','Prophet','KNN Regression','Feed-Forward Neural Networks',"LSTM","B_spline","Cubic")
      # colnames(MSE_matrix4) = colnames(data)[2:length(data[1,])]
      ####
      pred_true_matrix_list = list()


      ###############################################################################   only need to change here  , how many days u want to predict

      result_date <-data[start,1]
      print(result_date)

      ###Interval time
      I_time = abs(as.numeric(difftime(data[1,1], data[2,1], units = "days"))) ###Interval time
      days = ceiling(days_true/I_time)
      day_tag = paste(days_true,sep = " ",'days')
      #   make the mutation which we want to test predict
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

      train <- dataset[1:(train_size-1)]   ###选定的日期起，+之后x天. 一起作为training data
      test <- dataset[(train_size):c(train_size+pred_length-1)]   ### train size  之后再+x 天作为test

      # cat(length(train), length(test))
      look_back <- length(test)

      pred_true_matrix = matrix(NA,nrow = 6,ncol=look_back)

      # rownames(pred_true_matrix)=c('True','Arima','Prophet','KNN Regression','Feed-Forward Neural Networks',"LSTM","B_spline","Cubic")
      rownames(pred_true_matrix)=c('True','ARI','PRO','KNN',"B_S","C_S")

      pred_true_matrix[1,] = c(test)
      # ##################################################################    LSTM
      # # data=datas[1:1109,96]
      #
      # # look_back = 30
      # iteration  =50
      # future_length  =days_true
      # rr = my_LSTM(train,look_back,iteration,future_length)
      # h1h1 = rr[[1]]
      # h2h2 = rr[[2]]
      # test_na = which(h2h2%in%NA)
      # train_na = which(h1h1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[5,j] = mean(abs(test-h2h2))/mean(abs(test))
      #   MSE_matrix3[5,j] = mean((test-h2h2)^2)
      # }else{
      #   MSE_matrix1[5,j] = mean(abs(test[-test_na]-h2h2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[5,j] = mean((test[-test_na]-h2h2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[5,j] = mean(abs(train-h1h1))/mean(abs(train))      # train
      #   MSE_matrix4[5,j] = mean((train-h1h1)^2)
      # }else{
      #   MSE_matrix2[5,j] = mean(abs(train[-train_na]-h1h1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[5,j] = mean((train[-train_na]-h1h1[-train_na])^2)
      # }
      ##############################################################   Arima

      skirtsts<-ts(data[,which(colnames(data)%in%nn)],frequency=365)
      skirtsarima<-auto.arima(skirtsts[1:length(train)],trace=F)
      skirtsarimaforecast<-forecast(skirtsarima,h=days,level=c(99.5))
      a1a1 = skirtsarima[["fitted"]]
      a2a2 = as.double(skirtsarimaforecast$mean)

      # test_na = which(a2a2%in%NA)
      # train_na = which(a1a1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[1,j] = mean(abs(test-a2a2))/mean(abs(test))
      #   MSE_matrix3[1,j] = mean((test-a2a2)^2)
      # }else{
      #   MSE_matrix1[1,j] = mean(abs(test[-test_na]-a2a2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[1,j] = mean((test[-test_na]-a2a2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[1,j] = mean(abs(train-a1a1))/mean(abs(train))      # train
      #   MSE_matrix4[1,j] = mean((train-a1a1)^2)
      # }else{
      #   MSE_matrix2[1,j] = mean(abs(train[-train_na]-aa1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[1,j] = mean((train[-train_na]-a1a1[-train_na])^2)
      # }
      # MSE_matrix[1,j] = mean(abs(train-a1a1))/mean(abs(train))      # train
      pred_true_matrix[2,] = c(a2a2)

      ##############################################################   Prophet
      ddf <- data.frame(ds = DD$Date[1:(train_size-1)],
                        y = as.numeric(train))
      prophet_pred = prophet(ddf)
      future = make_future_dataframe(prophet_pred,periods=days)
      fcastprophet = predict(prophet_pred,future)
      #Creating train prediction dataset to compare real data
      dataprediction1 = data.frame(fcastprophet$ds,fcastprophet$yhat)
      trainlen = length(train)
      dataprediction = dataprediction1[c(1:trainlen),]

      c1c1 = dataprediction1[c(1:trainlen),2]
      c2c2 = dataprediction1[-c(1:trainlen),2]
      # test_na = which(c2c2%in%NA)
      # train_na = which(c1c1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[2,j] = mean(abs(test-c2c2))/mean(abs(test))
      #   MSE_matrix3[2,j] = mean((test-c2c2)^2)
      # }else{
      #   MSE_matrix1[2,j] = mean(abs(test[-test_na]-c2c2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[2,j] = mean((test[-test_na]-c2c2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[2,j] = mean(abs(train-c1c1))/mean(abs(train))      # train
      #   MSE_matrix4[2,j] = mean((train-c1c1)^2)
      # }else{
      #   MSE_matrix2[2,j] = mean(abs(train[-train_na]-c1c1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[2,j] = mean((train[-train_na]-c1c1[-train_na])^2)
      # }
      # MSE_matrix[3,j] = mean(abs(train-c1c1))/mean(abs(train))      # train
      pred_true_matrix[3,] = c(c2c2)

      ##############################################################   KNN Regression
      df_knn <- ddf

      predknn <- knn_forecasting(df_knn$y, h = days, lags = 1:30, k = 50, msas = "MIMO")
      d1d1 = c(predknn[["model"]][["ts"]])
      d2d2 = c(predknn[["prediction"]])
      #
      # test_na = which(d2d2%in%NA)
      # train_na = which(d1d1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[3,j] = mean(abs(test-d2d2))/mean(abs(test))
      #   MSE_matrix3[3,j] = mean((test-d2d2)^2)
      # }else{
      #   MSE_matrix1[3,j] = mean(abs(test[-test_na]-d2d2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[3,j] = mean((test[-test_na]-d2d2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[3,j] = mean(abs(train-d1d1))/mean(abs(train))      # train
      #   MSE_matrix4[3,j] = mean((train-d1d1)^2)
      # }else{
      #   MSE_matrix2[3,j] = mean(abs(train[-train_na]-d1d1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[3,j] = mean((train[-train_na]-d1d1[-train_na])^2)
      # }
      # MSE_matrix[4,j] = mean(abs(train-d1d1))/mean(abs(train))      # train
      pred_true_matrix[4,] = c(d2d2)

      ##############################################################   Feed-Forward Neural Networks
      #
      #
      # lambda = BoxCox.lambda(train)
      # dnn_fit = nnetar(train,lambda=lambda)
      # dnn_fit
      # fcast = forecast(dnn_fit,PI=T,h=days)
      # #autoplot(fcast)
      # f1f1 = c(fcast[["fitted"]])
      # f2f2 = c(fcast[["mean"]])
      #
      # test_na = which(f2f2%in%NA)
      # train_na = which(f1f1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[4,j] = mean(abs(test-f2f2))/mean(abs(test))
      #   MSE_matrix3[4,j] = mean((test-f2f2)^2)
      # }else{
      #   MSE_matrix1[4,j] = mean(abs(test[-test_na]-f2f2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[4,j] = mean((test[-test_na]-f2f2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[4,j] = mean(abs(train-f1f1))/mean(abs(train))      # train
      #   MSE_matrix4[4,j] = mean((train-f1f1)^2)
      # }else{
      #   MSE_matrix2[4,j] = mean(abs(train[-train_na]-f1f1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[4,j] = mean((train[-train_na]-f1f1[-train_na])^2)
      # }
      # pred_true_matrix[5,] = c(f2f2)

      #########################   Cubic Spline Forecast    ##################
      # x <- c(1:length(train))
      # x.new <- c((length(train)+1):(length(train)+pred_length))
      # y <- skirtsts[1:length(train)]
      #
      # # My n.s fit:
      # n0=c(5,10,15,20)
      # fit.temp <- lm(y ~ ns(x, df=n0))
      # i1i1 = fit.temp[["fitted.values"]]
      # i2i2 = predict(fit.temp, newdata=data.frame(x=x.new))

      # x <- c(1:length(train))
      # x.new <- c((length(train)+1):(length(train)+pred_length))
      # y <- skirtsts[1:length(train)]
      #
      # # Adjusting the degrees of freedom for a smoother curve
      # # Higher degrees of freedom can lead to a smoother curve
      # n0 <- c(30, 40, 50, 60) # Increased df values for smoother fit
      #
      # fit.temp <- lm(y ~ ns(x, df=n0))
      # i1i1 = fit.temp[["fitted.values"]]
      # i2i2 = predict(fit.temp, newdata=data.frame(x=x.new))
      #

      x <- c(1:length(train))
      x.new <- c((length(train)+1):(length(train)+pred_length))
      y <- skirtsts[1:length(train)]

      # Start with a low number of df for a smoother fit
      df_value <- 4 # Adjust this value as needed

      fit.temp <- lm(y ~ ns(x, df=df_value))
      i1i1 = fit.temp[["fitted.values"]]
      i2i2 = predict(fit.temp, newdata=data.frame(x=x.new))
      pred_true_matrix[6,] = c(i2i2)

      # # Optionally, you can plot the original and smoothed predictions to compare
      # plot(x, y, type = "l", col = "red", lwd = 2) # Original data
      # lines(x, i1i1, col = "blue", lwd = 2) # Fitted values
      # lines(x.new, i2i2, col = "green", lwd = 2) # Predicted values
      #
      # MSE_matrix4[7,j] = mean((train-i1i1)^2)
      # test_na = which(i2i2%in%NA)
      # train_na = which(i1i1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[7,j] = mean(abs(test-i2i2))/mean(abs(test))
      #   MSE_matrix3[7,j] = mean((test-i2i2)^2)
      # }else{
      #   MSE_matrix1[7,j] = mean(abs(test[-test_na]-i2i2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[7,j] = mean((test[-test_na]-i2i2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[7,j] = mean(abs(train-i1i1))/mean(abs(train))      # train
      #   MSE_matrix4[7,j] = mean((train-i1i1)^2)
      # }else{
      #   MSE_matrix2[7,j] = mean(abs(train[-train_na]-i1i1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[7,j] = mean((train[-train_na]-i1i1[-train_na])^2)
      # }
      #########################   B Spline Forecast    ######################
      x2 <- c(1:length(train))
      x2.new <- c((length(train)+1):(length(train)+pred_length))
      y2 <- skirtsts[1:length(train)]
      n0=c(5,10,15,20)
      fit.temp <- lm(y2 ~ bs(x2, df=n0))
      j1j1 = fit.temp[["fitted.values"]]
      j2j2 = predict(fit.temp, newdata=data.frame(x2=x2.new))
      pred_true_matrix[5,] = c(j2j2)
      # test_na = which(j2j2%in%NA)
      # train_na = which(j1j1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[6,j] = mean(abs(test-j2j2))/mean(abs(test))
      #   MSE_matrix3[6,j] = mean((test-j2j2)^2)
      # }else{
      #   MSE_matrix1[6,j] = mean(abs(test[-test_na]-j2j2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[6,j] = mean((test[-test_na]-j2j2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[6,j] = mean(abs(train-j1j1))/mean(abs(train))      # train
      #   MSE_matrix4[6,j] = mean((train-j1j1)^2)
      # }else{
      #   MSE_matrix2[6,j] = mean(abs(train[-train_na]-j1j1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[6,j] = mean((train[-train_na]-j1j1[-train_na])^2)
      # }


      ########################################################################################################################################################################
      ########################################################################################################################################################################
      ##########  mechine learning predicting



      mm = length(model_list)+1   #   want to plot all model together
      # par(mfrow=c(1,1))
      score_list = list()


      output_data = matrix(data=NA,nrow = 11,ncol=pred_length)

      for (m in c(1:mm)) {
        #

        if(m<mm){
          model=model_list[[m]]

          if(class(model) =="randomForest"){
            test_data <- c(test_pattern[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)] )   ### e.g.   1~21
            test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)]    ### e.g.   1~21
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
            test_data <- test_pattern[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)]    ### e.g.   1~21
            test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)]    ### e.g.   1~21
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



        # score = matrix(data=NA,nrow = pred_length,ncol = 2)
        # t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),which(colnames(test_pattern)%in%nn)]+test_pattern_minus[(start+num_lags):(start+(num_lags-1)+pred_length),which(colnames(test_pattern)%in%nn)]
        # score[,1]=t_true
        # print(length(pred_data))
        # score[,2]=pred_data_com
        #
        #
        #
        # x <- seq(1:pred_length)+start+num_lags-1

        output_data[m,]=pred_data_com
      }
      # print("$$$$")
      # print(dim(output_data))
      # print(a2a2)
      # print(c2c2)
      # print(d2d2)
      # print(i2i2)
      # print(j2j2)
      # print(paste(names(model_list),c("ARI","PRO","KNN","C_S","B_S")))
      #

      output_data[mm,] = a2a2
      output_data[mm+1,] = c2c2
      output_data[mm+2,] = d2d2
      output_data[mm+3,] = i2i2
      output_data[mm+4,] = j2j2#f2f2
      rownames(output_data)=c(names(model_list),c("ARI","PRO","KNN","C_S","B_S"))
      # print(c(names(model_list),c("ARI","PRO","KNN","C_S","B_S")))
      print(output_data)
    }
    result_list[[as.character(start)]] = output_data
  }
  return(result_list)
}

######################-----------------------------------   Going to predict with parameter
# train_pattern =  train_matrix$all
#nn = "AG_28878_TC"
#start =880
# site = mu_num
# pred_length = 30
# num_lags =21
# # data_input
# date = data[,1]
# test_pattern = data_B
# test_pattern_minus = datac2_m
# position  =c(270,410,500,600,710,880,920,980,1000)
# datasum=data
# model_list kinds of mechine learning trained model inside


# R/Predicting2.R

#' Predicting2
#'
#' This function sets a pattern.
#'
#' @param x The input parameter.
#' @return The pattern.
#' @export
Predicting2<-function(train_pattern,nn,pred_length,num_lags,datasum,position,model_list){
  date = datasum[["date"]]$date
  test_pattern = datasum$Changed
  test_pattern_minus = datasum$Minus
  dddd = position
  result_list = vector("list",length(dddd))
  names(result_list) <-as.character(dddd)
  ##### going to predict from different position in dddd
  for (i in c(1:length(dddd))) {
    start=dddd[i]
    # if want to plot
    # pdf_namesss <- paste0("F:/OneDrive - Indiana University/shared/covid/one_code_clean/a_final_1204/diff/2//", paste0(nn,start,"--222222-"), ".pdf")
    {
      ############################################################################################################################################################################
      ###########################################    this part for time serises      ###########################################
      ############################################################################################################################################################################
      #   time serises using the raw data
      days_true = pred_length
      data=cbind(date,datasum$Raw)
      data = as.data.frame(data)
      s2 <- data[, !duplicated(t(data))]
      data = s2
      j=1

      pred_true_matrix_list = list()


      ###############################################################################   only need to change here  , how many days u want to predict

      result_date <-data[start,1]
      print(result_date)

      ###Interval time
      I_time = abs(as.numeric(difftime(data[1,1], data[2,1], units = "days"))) ###Interval time
      days = ceiling(days_true/I_time)
      day_tag = paste(days_true,sep = " ",'days')
      #   make the mutation which we want to test predict
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

      train <- dataset[1:(train_size-1)]   ###选定的日期起，+之后x天. 一起作为training data
      test <- dataset[(train_size):c(train_size+pred_length-1)]   ### train size  之后再+x 天作为test

      # cat(length(train), length(test))
      look_back <- length(test)

      pred_true_matrix = matrix(NA,nrow = 6,ncol=look_back)

      # rownames(pred_true_matrix)=c('True','Arima','Prophet','KNN Regression','Feed-Forward Neural Networks',"LSTM","B_spline","Cubic")
      rownames(pred_true_matrix)=c('True','ARI','PRO','KNN',"B_S","C_S")

      pred_true_matrix[1,] = c(test)
      # ##################################################################    LSTM
      # # data=datas[1:1109,96]
      #
      # # look_back = 30
      # iteration  =50
      # future_length  =days_true
      # rr = my_LSTM(train,look_back,iteration,future_length)
      # h1h1 = rr[[1]]
      # h2h2 = rr[[2]]
      # test_na = which(h2h2%in%NA)
      # train_na = which(h1h1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[5,j] = mean(abs(test-h2h2))/mean(abs(test))
      #   MSE_matrix3[5,j] = mean((test-h2h2)^2)
      # }else{
      #   MSE_matrix1[5,j] = mean(abs(test[-test_na]-h2h2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[5,j] = mean((test[-test_na]-h2h2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[5,j] = mean(abs(train-h1h1))/mean(abs(train))      # train
      #   MSE_matrix4[5,j] = mean((train-h1h1)^2)
      # }else{
      #   MSE_matrix2[5,j] = mean(abs(train[-train_na]-h1h1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[5,j] = mean((train[-train_na]-h1h1[-train_na])^2)
      # }
      ##############################################################   Arima

      skirtsts<-ts(data[,which(colnames(data)%in%nn)],frequency=365)
      skirtsarima<-auto.arima(skirtsts[1:length(train)],trace=F)
      skirtsarimaforecast<-forecast(skirtsarima,h=days,level=c(99.5))
      a1a1 = skirtsarima[["fitted"]]
      a2a2 = as.double(skirtsarimaforecast$mean)

      # test_na = which(a2a2%in%NA)
      # train_na = which(a1a1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[1,j] = mean(abs(test-a2a2))/mean(abs(test))
      #   MSE_matrix3[1,j] = mean((test-a2a2)^2)
      # }else{
      #   MSE_matrix1[1,j] = mean(abs(test[-test_na]-a2a2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[1,j] = mean((test[-test_na]-a2a2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[1,j] = mean(abs(train-a1a1))/mean(abs(train))      # train
      #   MSE_matrix4[1,j] = mean((train-a1a1)^2)
      # }else{
      #   MSE_matrix2[1,j] = mean(abs(train[-train_na]-aa1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[1,j] = mean((train[-train_na]-a1a1[-train_na])^2)
      # }
      # MSE_matrix[1,j] = mean(abs(train-a1a1))/mean(abs(train))      # train
      pred_true_matrix[2,] = c(a2a2)

      ##############################################################   Prophet
      ddf <- data.frame(ds = DD$Date[1:(train_size-1)],
                        y = as.numeric(train))
      prophet_pred = prophet(ddf)
      future = make_future_dataframe(prophet_pred,periods=days)
      fcastprophet = predict(prophet_pred,future)
      #Creating train prediction dataset to compare real data
      dataprediction1 = data.frame(fcastprophet$ds,fcastprophet$yhat)
      trainlen = length(train)
      dataprediction = dataprediction1[c(1:trainlen),]

      c1c1 = dataprediction1[c(1:trainlen),2]
      c2c2 = dataprediction1[-c(1:trainlen),2]
      # test_na = which(c2c2%in%NA)
      # train_na = which(c1c1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[2,j] = mean(abs(test-c2c2))/mean(abs(test))
      #   MSE_matrix3[2,j] = mean((test-c2c2)^2)
      # }else{
      #   MSE_matrix1[2,j] = mean(abs(test[-test_na]-c2c2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[2,j] = mean((test[-test_na]-c2c2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[2,j] = mean(abs(train-c1c1))/mean(abs(train))      # train
      #   MSE_matrix4[2,j] = mean((train-c1c1)^2)
      # }else{
      #   MSE_matrix2[2,j] = mean(abs(train[-train_na]-c1c1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[2,j] = mean((train[-train_na]-c1c1[-train_na])^2)
      # }
      # MSE_matrix[3,j] = mean(abs(train-c1c1))/mean(abs(train))      # train
      pred_true_matrix[3,] = c(c2c2)

      ##############################################################   KNN Regression
      df_knn <- ddf

      predknn <- knn_forecasting(df_knn$y, h = days, lags = 1:30, k = 50, msas = "MIMO")
      d1d1 = c(predknn[["model"]][["ts"]])
      d2d2 = c(predknn[["prediction"]])
      #
      # test_na = which(d2d2%in%NA)
      # train_na = which(d1d1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[3,j] = mean(abs(test-d2d2))/mean(abs(test))
      #   MSE_matrix3[3,j] = mean((test-d2d2)^2)
      # }else{
      #   MSE_matrix1[3,j] = mean(abs(test[-test_na]-d2d2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[3,j] = mean((test[-test_na]-d2d2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[3,j] = mean(abs(train-d1d1))/mean(abs(train))      # train
      #   MSE_matrix4[3,j] = mean((train-d1d1)^2)
      # }else{
      #   MSE_matrix2[3,j] = mean(abs(train[-train_na]-d1d1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[3,j] = mean((train[-train_na]-d1d1[-train_na])^2)
      # }
      # MSE_matrix[4,j] = mean(abs(train-d1d1))/mean(abs(train))      # train
      pred_true_matrix[4,] = c(d2d2)

      ##############################################################   Feed-Forward Neural Networks
      #
      #
      # lambda = BoxCox.lambda(train)
      # dnn_fit = nnetar(train,lambda=lambda)
      # dnn_fit
      # fcast = forecast(dnn_fit,PI=T,h=days)
      # #autoplot(fcast)
      # f1f1 = c(fcast[["fitted"]])
      # f2f2 = c(fcast[["mean"]])
      #
      # test_na = which(f2f2%in%NA)
      # train_na = which(f1f1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[4,j] = mean(abs(test-f2f2))/mean(abs(test))
      #   MSE_matrix3[4,j] = mean((test-f2f2)^2)
      # }else{
      #   MSE_matrix1[4,j] = mean(abs(test[-test_na]-f2f2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[4,j] = mean((test[-test_na]-f2f2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[4,j] = mean(abs(train-f1f1))/mean(abs(train))      # train
      #   MSE_matrix4[4,j] = mean((train-f1f1)^2)
      # }else{
      #   MSE_matrix2[4,j] = mean(abs(train[-train_na]-f1f1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[4,j] = mean((train[-train_na]-f1f1[-train_na])^2)
      # }
      # pred_true_matrix[5,] = c(f2f2)

      #########################   Cubic Spline Forecast    ##################
      # x <- c(1:length(train))
      # x.new <- c((length(train)+1):(length(train)+pred_length))
      # y <- skirtsts[1:length(train)]
      #
      # # My n.s fit:
      # n0=c(5,10,15,20)
      # fit.temp <- lm(y ~ ns(x, df=n0))
      # i1i1 = fit.temp[["fitted.values"]]
      # i2i2 = predict(fit.temp, newdata=data.frame(x=x.new))

      # x <- c(1:length(train))
      # x.new <- c((length(train)+1):(length(train)+pred_length))
      # y <- skirtsts[1:length(train)]
      #
      # # Adjusting the degrees of freedom for a smoother curve
      # # Higher degrees of freedom can lead to a smoother curve
      # n0 <- c(30, 40, 50, 60) # Increased df values for smoother fit
      #
      # fit.temp <- lm(y ~ ns(x, df=n0))
      # i1i1 = fit.temp[["fitted.values"]]
      # i2i2 = predict(fit.temp, newdata=data.frame(x=x.new))
      #

      x <- c(1:length(train))
      x.new <- c((length(train)+1):(length(train)+pred_length))
      y <- skirtsts[1:length(train)]

      # Start with a low number of df for a smoother fit
      df_value <- 4 # Adjust this value as needed

      fit.temp <- lm(y ~ ns(x, df=df_value))
      i1i1 = fit.temp[["fitted.values"]]
      i2i2 = predict(fit.temp, newdata=data.frame(x=x.new))
      pred_true_matrix[6,] = c(i2i2)

      # # Optionally, you can plot the original and smoothed predictions to compare
      # plot(x, y, type = "l", col = "red", lwd = 2) # Original data
      # lines(x, i1i1, col = "blue", lwd = 2) # Fitted values
      # lines(x.new, i2i2, col = "green", lwd = 2) # Predicted values
      #
      # MSE_matrix4[7,j] = mean((train-i1i1)^2)
      # test_na = which(i2i2%in%NA)
      # train_na = which(i1i1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[7,j] = mean(abs(test-i2i2))/mean(abs(test))
      #   MSE_matrix3[7,j] = mean((test-i2i2)^2)
      # }else{
      #   MSE_matrix1[7,j] = mean(abs(test[-test_na]-i2i2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[7,j] = mean((test[-test_na]-i2i2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[7,j] = mean(abs(train-i1i1))/mean(abs(train))      # train
      #   MSE_matrix4[7,j] = mean((train-i1i1)^2)
      # }else{
      #   MSE_matrix2[7,j] = mean(abs(train[-train_na]-i1i1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[7,j] = mean((train[-train_na]-i1i1[-train_na])^2)
      # }
      #########################   B Spline Forecast    ######################
      x2 <- c(1:length(train))
      x2.new <- c((length(train)+1):(length(train)+pred_length))
      y2 <- skirtsts[1:length(train)]
      n0=c(5,10,15,20)
      fit.temp <- lm(y2 ~ bs(x2, df=n0))
      j1j1 = fit.temp[["fitted.values"]]
      j2j2 = predict(fit.temp, newdata=data.frame(x2=x2.new))
      pred_true_matrix[5,] = c(j2j2)
      # test_na = which(j2j2%in%NA)
      # train_na = which(j1j1%in%NA)
      # if(length(test_na)==0){
      #   MSE_matrix1[6,j] = mean(abs(test-j2j2))/mean(abs(test))
      #   MSE_matrix3[6,j] = mean((test-j2j2)^2)
      # }else{
      #   MSE_matrix1[6,j] = mean(abs(test[-test_na]-j2j2[-test_na]))/mean(abs(test[-test_na]))
      #   MSE_matrix3[6,j] = mean((test[-test_na]-j2j2[-test_na])^2)
      # }
      # if(length(train_na)==0){
      #   MSE_matrix2[6,j] = mean(abs(train-j1j1))/mean(abs(train))      # train
      #   MSE_matrix4[6,j] = mean((train-j1j1)^2)
      # }else{
      #   MSE_matrix2[6,j] = mean(abs(train[-train_na]-j1j1[-train_na]))/mean(abs(train[-train_na]))      # train
      #   MSE_matrix4[6,j] = mean((train[-train_na]-j1j1[-train_na])^2)
      # }


      ########################################################################################################################################################################
      ########################################################################################################################################################################
      ##########  mechine learning predicting



      mm = length(model_list)+1   #   want to plot all model together
      # par(mfrow=c(1,1))
      score_list = list()


      output_data = matrix(data=NA,nrow = 11,ncol=pred_length)

      for (m in c(1:mm)) {
        #

        if(m<mm){
          model=model_list[[m]]

          if(class(model) =="randomForest"){
            test_data <- c(test_pattern[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)] )   ### e.g.   1~21
            test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)]    ### e.g.   1~21
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
            test_data <- test_pattern[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)]    ### e.g.   1~21
            test_data_minus<-test_pattern_minus[start:(start+(num_lags-1)),which(colnames(test_pattern)%in%nn)]    ### e.g.   1~21
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



        # score = matrix(data=NA,nrow = pred_length,ncol = 2)
        # t_true = test_pattern[(start+num_lags):(start+(num_lags-1)+pred_length),which(colnames(test_pattern)%in%nn)]+test_pattern_minus[(start+num_lags):(start+(num_lags-1)+pred_length),which(colnames(test_pattern)%in%nn)]
        # score[,1]=t_true
        # print(length(pred_data))
        # score[,2]=pred_data_com
        #
        #
        #
        # x <- seq(1:pred_length)+start+num_lags-1

        output_data[m,]=pred_data_com
      }
      # print("$$$$")
      # print(dim(output_data))
      # print(a2a2)
      # print(c2c2)
      # print(d2d2)
      # print(i2i2)
      # print(j2j2)
      # print(paste(names(model_list),c("ARI","PRO","KNN","C_S","B_S")))
      #

      output_data[mm,] = a2a2
      output_data[mm+1,] = c2c2
      output_data[mm+2,] = d2d2
      output_data[mm+3,] = i2i2
      output_data[mm+4,] = j2j2#f2f2
      rownames(output_data)=c(names(model_list),c("ARI","PRO","KNN","C_S","B_S"))
      # print(c(names(model_list),c("ARI","PRO","KNN","C_S","B_S")))
      print(output_data)
    }
    result_list[[as.character(start)]] = output_data
  }
  return(result_list)
}



