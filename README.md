# Steps for usage & Function explanation

## loading function

path = '/Users/..../'

name= "result_matrix_1D_202306.tsv"

data = loading(path,name)

### output
![image](https://github.com/ZhouXY199502/SWD/assets/70504823/8f612ff4-6942-4fcc-80d6-edc268b190e8)


## set_pattern function

data_raw = data$Raw

data_changed = data$Changed

data_ori = data$Original

data_minus = data$Minus

threshold=0.05

cut_time="2021-01-01"

data_pattern = set_pattern(data_raw,data_changed,data_ori,data_minus,threshold,cut_time)

### output
![image](https://github.com/ZhouXY199502/SWD/assets/70504823/a3c4038f-7038-4de4-8315-bd843c56c038)


## train_matrix function

num_lags=21        

laglag=1

pred_length=30     

timee = c("2022-09-01")       

meth_name = c("RF","S_L","S_R","XG","LR","FNN","ARI","PRO","KNN","C_S","B_S")

removed_names=c('C_28311_T',"AG_28878_TC")        

train_matrix = get_train_matrix(num_lags,laglag,pred_length,timee,data,removed_names,data_pattern)
### output
<img width="1035" alt="image" src="https://github.com/ZhouXY199502/SWD/assets/70504823/0a373b39-d26f-4822-bf44-08f27bdffbaf">


## Train_ML function

MLE_meth_name_list = c("RF","SL","SR","XG","LR","FNN")

ML_trained = Train_ML(train_patterns_list,path,MLE_meth_name=MLE_meth_name_list,svm_type="eps-regression",XG_eta = 0.3,XG_max_depth = 6,
                      XG_eval_metric = "rmse",XG_num_rounds=300)



<img width="709" alt="image" src="https://github.com/ZhouXY199502/SWD/assets/70504823/91f53d21-8ae1-4d7f-9b0a-d54a2ca1e860">

## Predicting function

train_pattern =  train_matrix$all
nn = "AG_28878_TC"    # the mutation what we want to test

pred_length = 30   #  predict future length
num_lags =21   # use to input and predict the future
datasum = data
position  =c(270,410,500,600,710,880,920,980,1000)
model_list = ML_trained$train_all
Pred_results<- Predicting(train_pattern,nn,pred_length,num_lags,datasum,position,model_list)


<img width="1182" alt="image" src="https://github.com/ZhouXY199502/SWD/assets/70504823/ddeb8923-5aa3-4978-8b72-4825ee5bae79">















