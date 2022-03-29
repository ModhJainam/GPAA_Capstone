#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import io
from sklearn.utils import shuffle
from sklearn.naive_bayes import BernoulliNB
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedShuffleSplit
from collections import defaultdict 
from sklearn.model_selection import StratifiedKFold
import random
from matplotlib import pyplot as plt
import seaborn as sns
from collections import defaultdict
import itertools
from subprocess import call
from sklearn import tree
from sklearn.metrics import fbeta_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
#from sklearn.linear_model import RandomizedLasso
from sklearn.linear_model import Lasso
from sklearn.feature_selection import SelectFromModel
from scipy.stats import rankdata

with open('/project/gpaa/machine_learning/jainam_capstone/preprocessing/LAD_labeled_batches0-15.csv', 'r') as f:
  df_labeled = pd.read_csv(f)
  

df_binary_labeled = df_labeled.drop(df_labeled.columns[0:28], axis=1)
df_binary_labeled = df_binary_labeled.drop(df_labeled.columns[-6:], axis=1)

  
X = df_binary_labeled.drop('binary_pathology', axis=1)
target = df_binary_labeled['binary_pathology']

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

def randomForest(column, dataDummy,PrecolNum,rfTrees,threads):
    x_train = dataDummy[list(dataDummy.columns[0:PrecolNum-1])]
    names = dataDummy.columns[0:PrecolNum-1]
    y_train =dataDummy[column]
    rf = RandomForestClassifier(n_estimators = rfTrees, n_jobs = threads, random_state=123456)
    rf.fit(x_train, y_train)
    Ranked_Features=sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), names),reverse=True)
    return Ranked_Features

def rankInformative(Ranked_Features):
    RankedList = []
    midcounter = 0
    for x in Ranked_Features:
        midcounter +=1
        RankedList.append(x[1])
        rankedDict[column] = RankedList
        if midcounter==30:
            break
    return RankedList  

def rankInformativeAll(Ranked_Features):
    RankedList = []
    midcounter = 0
    for x in Ranked_Features:
        midcounter +=1
        RankedList.append(x[1])
        rankedDict[column] = RankedList
        if midcounter==10000:
            break
    return RankedList             

def negativeOut(x, column, medianValues,Median_Expression_Level):
    Positive_RankedList_Complete = []
    for i in x:
        if len(Positive_RankedList_Complete) == 240:
          break
        if medianValues.loc[column, i] > Median_Expression_Level:
            #print(i)
            #print(medianValues.loc[column, i])
            Positive_RankedList_Complete.append(i)

    return Positive_RankedList_Complete

def binaryScore(Positive_RankedList_Complete, informativeGenes, medianValues, column):
    Positive_RankedList=list(Positive_RankedList_Complete[0:InformativeGenes])
    Median_RF_Subset=medianValues.loc[:, Positive_RankedList]
    Rescaled_Matrix=pd.DataFrame()
    
    for i in Positive_RankedList:
        Target_value=medianValues.loc[column, i]
        Rescaled_values=Median_RF_Subset[[i]].divide(Target_value)
        Rescaled_Matrix=pd.concat([Rescaled_Matrix,Rescaled_values],axis=1)
    difference_matrix=Rescaled_Matrix.apply(lambda x: 1-x, axis=1)
    difference_matrix_clean = difference_matrix.where(difference_matrix > 0, 0)
    ColumnSums=difference_matrix_clean.sum(0)
    rescaled = ColumnSums/clusters2Loop

    # Double sort so that for ties, the RF ranking prevails!     
    Ranked_Features_df=pd.DataFrame(Ranked_Features)
    Ranked_Features_df.rename(columns={1: 'Symbol'}, inplace=True)
    Ranked_Features_df_indexed=Ranked_Features_df.set_index("Symbol")
    rescaled_df=pd.DataFrame(rescaled)
    binaryAndinformation_Ranks=rescaled_df.join(Ranked_Features_df_indexed,lsuffix='_scaled', rsuffix='_informationGain')
    binaryAndinformation_Ranks.sort_values(by=['0_scaled','0_informationGain'],ascending= [False, False], inplace = True)
    Binary_ranked_Genes=binaryAndinformation_Ranks.index.tolist()
    Binary_RankedList=list(Binary_ranked_Genes[0:Genes_to_testing])
    Binary_scores=rescaled.to_dict()
    global Binary_store_DF 
    Binary_store_DF = Binary_store_DF.append(binaryAndinformation_Ranks)
    return Binary_RankedList

def DT_cutOffs(x, column):
    cut_dict = {}
    for i in x:
        filename=str(i)
        y_train =dataDummy[column]
        x_train = dataDummy[i]
        X = x_train[:, None]
        clf = tree.DecisionTreeClassifier(max_leaf_nodes=2)
        clf = clf.fit(X, y_train)
        threshold = clf.tree_.threshold
        cut_dict[i]=threshold[0]
    return cut_dict

def queryGenerator(x, cut_dict):
    queryList = []
    for i in x:
        str1 = i
        current_value = cut_dict.get(str1)
        queryString1 = str(str1)+'>='+ str(current_value)
        queryList.append(queryString1)
    return queryList

def permutor(x):
    binarylist2 = x
    combs = []
    for i in range(1, len(x)+1):
        els = [list(x) for x in itertools.permutations(binarylist2, i)]
        combs.extend(els)
    return combs

def fbetaTest(x, column,testArray, betaValue):
    fbeta_dict = {}
    for list in x:
        testArray['y_pred']= 0
        betaQuery = '&'.join(list)
        Ineq1=dataFull.query(betaQuery)
        testList=Ineq1.index.tolist()
        testArray.loc[testList, 'y_pred'] = 1
        f1 = fbeta_score(testArray['y_true'], testArray['y_pred'], average= 'binary', beta=betaValue)        
        dictName = column+"&"+betaQuery
        fbeta_dict[dictName] = f1 
    return fbeta_dict

'''
############################ RLASSO ##########################################
split =1
for train_index, test_index in skf.split(X, target):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = target[train_index], target[test_index]
    
    randomized_lasso = RandomizedLasso(random_state=101, normalize=False, verbose=True)
    randomized_lasso.fit(X_train, y_train)
    
    randomized_lasso_rank = randomized_lasso.scores_.tolist()
    
    ranks2 = rankdata(randomized_lasso_rank)
    RLASSO_to_write = ranks2
    np.savetxt("/project/gpaa/machine_learning/jainam_capstone/feature_selection/feature_selection_results/top_240_indices_LAD_RLASSO_split_"+str(split)+".csv", ranks2.argsort()[-240:][::-1], delimiter=",")
    split += 1
'''
############################ LASSO ##########################################
split =1
for train_index, test_index in skf.split(X, target):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = target[train_index], target[test_index]
    
    lasso = Lasso(random_state=101, selection = 'random', normalize=False)
    lasso.fit(X_train, y_train)
    
    lasso_rank = lasso.coef_.tolist()
    
    ranks2 = rankdata(lasso_rank)
    LASSO_to_write = ranks2
    np.savetxt("/project/gpaa/machine_learning/jainam_capstone/feature_selection/feature_selection_results/top_240_indices_LAD_LASSO_split_"+str(split)+".csv", ranks2.argsort()[-240:][::-1], delimiter=",")
    split += 1
########################### Random Forest ####################################
split = 1    
for train_index, test_index in skf.split(X, target):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = target[train_index], target[test_index]
    
    rfc = RandomForestClassifier(random_state=101, n_jobs = -1)
    rfc.fit(X_train, y_train)
    
    RF = SelectFromModel(estimator=rfc)
    RF_rank = RF.estimator.feature_importances_.tolist()
    
    ranks3 = rankdata(RF_rank)
    RF_to_write = ranks3
    np.savetxt("/project/gpaa/machine_learning/jainam_capstone/feature_selection/feature_selection_results/top_240_indices_LAD_random_forest_split_"+str(split)+".csv", ranks3.argsort()[-240:][::-1], delimiter=",")
    split += 1

'''    
############################## NSF ###########################################

for train_index, test_index in skf.split(X, target):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = target[train_index], target[test_index]
    
    dataFull = df_labeled
    dataFull = dataFull.set_index('sample_id')
    dataFull = dataFull.rename(columns={"binary_pathology" : "Clusters"})
    dataFull = dataFull.drop(['age', 'white', 'black', 'hispanic', 'male', 'female'], axis=1)
    dataFull['Clusters'] = dataFull.Clusters.astype(str)
    dataFull['Clusters'] = dataFull.Clusters.str.replace('1','A') #CAD
    dataFull['Clusters'] = dataFull.Clusters.str.replace('0','B') #control
    dataFull = dataFull.iloc[train_indices]
    #dataFull.head()
    
    #Creates dummy columns for one vs all Random Forest modeling
    dataDummy = pd.get_dummies(dataFull, columns=["Clusters"], prefix = "", prefix_sep = "")
    
    #Creates matrix of cluster median expression values
    medianValues = dataFull.groupby(by="Clusters").median()
    medianValues.to_csv('Function_medianValues.csv')     
    
    #Finding the number of clusters and printing that to screen (sanity check)
    PrecolNum = len(dataFull.columns)
    PostcolNum = len(dataDummy.columns)
    adjustedColumns = PrecolNum-1
    clusters2Loop=PostcolNum-PrecolNum
    #print clusters2Loop
    
    
    ####Random Forest parameters
    rfTrees=1000 #Number of trees
    threads=1     #Number of threads to use, -1 is the greedy option where it will take all available CPUs/RAM
    
    ####Filtering and ranking of genes from random forest parameters
    
    Median_Expression_Level = 0
    InformativeGenes = 240 #How many top genes from the Random Forest ranked features will be evaluated for binariness 
    Genes_to_testing = 240    #How many top genes ranked by binary score will be evaluated in permutations by fbeta-score (as the number increases the number of permutation rises exponentially!)
    
    #### fbeta-score parameters                   
    
    betaValue = 0.5 #Set values for fbeta weighting. 1 is default f-measure. close to zero is Precision, greater than 1 weights toward Recall
    
    
    #Core analysis - pathology
    rankedDict =  {}  ###gives us the top ten features from RF
    f1_store_1D = {}
    Binary_score_store_DF=pd.DataFrame()
    DT_cutoffs_store={}
    
    NSF_pathology_1 = []
    NSF_pathology_0 = []
    
    random.seed(42)
    for column in dataDummy.columns[PrecolNum-1:PostcolNum]:
            print(column)
            ## Run Random Forest and get a ranked list 
            Ranked_Features= randomForest(column, dataDummy, PrecolNum, rfTrees, threads)
            RankedList = rankInformative(Ranked_Features)
            RankedListAll = rankInformativeAll(Ranked_Features)
            print(len(RankedListAll))
            
            ## Setup testArray for f-beta evaluation
            testArray = dataDummy[[column]]
            testArray.columns = ['y_true']
            
            #Rerank according to expression level and binary score
            Positive_RankedList_Complete = negativeOut(RankedList, column, medianValues, Median_Expression_Level) 
            Positive_RankedList_Complete240 = negativeOut(RankedListAll, column, medianValues, Median_Expression_Level)
            print(len(Positive_RankedList_Complete240))
            Binary_store_DF = pd.DataFrame()
            Binary_RankedList = binaryScore(Positive_RankedList_Complete240, InformativeGenes, medianValues, column)
            print(len(Binary_RankedList))
            if column == 'A':
              NSF_pathology_1 = Binary_RankedList
            if column == 'B':
              NSF_pathology_0 = Binary_RankedList
              
    np.savetxt("/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/feature_selection/holdout_feature_selections_results/top_120_indices_AA_NSF1_b13_split_"+str(split)+".csv", NSF_pathology_1, delimiter=",", fmt="%s")
    np.savetxt("/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/feature_selection/holdout_feature_selections_results/top_120_indices_AA_NSF0_b13_split_"+str(split)+".csv", NSF_pathology_0, delimiter=",", fmt="%s")
    
'''    
############################## RFE ###########################################
split = 1
for train_index, test_index in skf.split(X, target):
    print("TRAIN:", train_index, "TEST:", test_index)
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = target[train_index], target[test_index]  
    
    rfc = RandomForestClassifier(random_state=101, n_jobs = -1)
    rfe = RFE(estimator=rfc, step=1, verbose=0, n_features_to_select=1)
    rfe.fit(X_train, y_train)
    
    rfe_rank = rfe.ranking_.tolist()
    ranks1 = rankdata(rfe_rank)
    
    rfe_to_write = ranks1
    np.savetxt("/project/gpaa/machine_learning/jainam_capstone/feature_selection/feature_selection_results/top_240_indices_LAD_rfe_split_"+str(split)+".csv", ranks1.argsort()[-240:][::-1], delimiter=",")
    split += 1
    
    