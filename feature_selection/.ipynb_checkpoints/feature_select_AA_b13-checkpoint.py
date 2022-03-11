# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np


with open('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/feature_selection/AA_labeled_b13_fixed.csv', 'r') as f:
  df_labled = pd.read_csv(f)
  

df_binary_labeled = df_labled.drop(['sample_id', 'age', 'white', 'black', 'hispanic', 'male', 'female'], axis=1)


from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE

X = df_binary_labeled.drop('binary_pathology', axis=1)
target = df_binary_labeled['binary_pathology']

rfc = RandomForestClassifier(random_state=101, n_jobs = -1)
rfe = RFE(estimator=rfc, step=1, verbose=0, n_features_to_select=1)
rfe.fit(X, target)


rfe_rank = rfe.ranking_.tolist()


from scipy.stats import rankdata
print("RFE")
print(len(rfe_rank))
print(rfe_rank)
print("Ranks:")
print(rankdata(rfe_rank).tolist())
print("----------")

# Random LASSO (RLASSO)

from sklearn.linear_model import RandomizedLasso
randomized_lasso = RandomizedLasso(random_state=101, normalize=False, verbose=True)
randomized_lasso.fit(X, target)

randomized_lasso_rank = randomized_lasso.scores_.tolist()

print("Randomized Lasso")
print(len(randomized_lasso_rank))
print(randomized_lasso_rank)
print("Ranks:")
print(rankdata(randomized_lasso_rank).tolist())
print("----------")

# Random Forrest (RF)

from sklearn.feature_selection import SelectFromModel

rfc.fit(X, target)
RF = SelectFromModel(estimator=rfc)

RF_rank = RF.estimator.feature_importances_.tolist()

print("Random Forrest")
print(len(RF_rank))
print(RF_rank)
print("Ranks:")
print(rankdata(RF_rank).tolist())
print("----------")

ranks1 = rankdata(rfe_rank)
ranks2 = rankdata(randomized_lasso_rank)
ranks3 = rankdata(RF_rank)
#ranks4 = rankdata(regr_rank)
added = np.add(np.add(ranks1, ranks2), ranks3)

ranks_avg = np.true_divide(added, 3)

ranks_ranked = rankdata(ranks_avg)
list_to_write = ranks_ranked.tolist()
print("average ranks")
print(list_to_write)

np.savetxt("average_ranks_AA_pathology_b13_fixed.csv", ranks_ranked, delimiter=",")

np.savetxt("top_240_indices_AA_averaged_b13_fixed.csv", ranks_ranked.argsort()[-240:][::-1], delimiter=",")

rfe_to_write = ranks1
np.savetxt("top_240_indices_AA_rfe_b13_fixed.csv", ranks1.argsort()[-240:][::-1], delimiter=",")

RF_to_write = ranks3
np.savetxt("top_240_indices_AA_random_forest_b13_fixed.csv", ranks3.argsort()[-240:][::-1], delimiter=",")

RLASSO_to_write = ranks2
np.savetxt("top_240_indices_AA_RLASSO_b13_fixed.csv", ranks2.argsort()[-240:][::-1], delimiter=",")

############################## NSF Function ##################################

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


############################## NSF ###########################################
    
dataFull = df_labled
dataFull = dataFull.set_index('sample_id')
dataFull = dataFull.rename(columns={"binary_pathology" : "Clusters"})
dataFull = dataFull.drop(['age', 'white', 'black', 'hispanic', 'male', 'female'], axis=1)
dataFull['Clusters'] = dataFull.Clusters.astype(str)
dataFull['Clusters'] = dataFull.Clusters.str.replace('1','A') #CAD
dataFull['Clusters'] = dataFull.Clusters.str.replace('0','B') #control
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
          
np.savetxt("/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/feature_selection/holdout_feature_selections_results/top_120_indices_AA_NSF1_b13.csv", NSF_pathology_1, delimiter=",", fmt="%s")
np.savetxt("/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/feature_selection/holdout_feature_selections_results/top_120_indices_AA_NSF0_b13.csv", NSF_pathology_0, delimiter=",", fmt="%s")
 