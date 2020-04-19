

# In[1]:


import pandas as pd
import sklearn
from sklearn.pipeline import Pipeline
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_score
import numpy as np
from scipy import stats
from sklearn.metrics import roc_curve, auc,f1_score,precision_score,recall_score,roc_auc_score,accuracy_score
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn import svm
from sklearn.svm import LinearSVC
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from __future__ import division
from sklearn.model_selection import train_test_split
from sklearn.cross_validation import train_test_split
from sklearn.feature_selection import RFE
from sklearn.svm import SVR
from sklearn import preprocessing
from sklearn.preprocessing import Imputer
from sklearn_pandas import DataFrameMapper
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, chi2

#load the data 
data = pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step37data/Cyto_mutations_rare.csv',sep=',')
data = pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step37data/ProteinDomain_mutations_all.csv',sep=',')
data = pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step37data/somatic_mutations_rare.csv',sep=',')
data = pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step37data/Gene_mutations_rare.csv',sep=',')

data=data.rename(columns={ data.columns[0]: 'Sample' })
Selected_Genetics = data.drop(['Sample','group'], axis=1)

RiskFactor=pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/RiskFactor.csv',sep=',')
RiskFactor = RiskFactor.drop(['Study_ID','TissueStatus','BirthDate','Nonsense_mutation'], axis=1)
Selected_RiskFactor=RiskFactor[['Benign_Age','Class','AgeofMerche','MotherBC','Preg1Age_cate']]
Selected_RiskFactor=RiskFactor[['Benign_Age','Class','AgeofMerche','MotherBC','Preg1Age_cate','scar','adhe']]

Selected_RiskFactor=RiskFactor[['MotherBC','Preg1Age_cate','Signature_1','Signature_2','Signature_3','Nonsense_mutation']]

#X=Selected_RiskFactor
#X=Selected_Genetics
X=pd.concat([Selected_RiskFactor.reset_index(drop=True),Selected_Genetics], axis=1)
X=Selected_Genetics
X=Selected_RiskFactor
imputer = Imputer(missing_values='NaN',strategy='median',axis=0)
X["Benign_Age"] = imputer.fit_transform(X[["Benign_Age"]]).ravel()
imputer = Imputer(missing_values='NaN',strategy='median',axis=0)
X["AgeofMerche"] = imputer.fit_transform(X[["AgeofMerche"]]).ravel()


mapper3 = DataFrameMapper([
        ('Class',sklearn.preprocessing.LabelEncoder()),
        ('MotherBC',sklearn.preprocessing.LabelEncoder()),
        ('Preg1Age_cate',sklearn.preprocessing.LabelEncoder()),
        ('scar',sklearn.preprocessing.LabelEncoder()),
        ('adhe',sklearn.preprocessing.LabelEncoder()),       
      ], default=None)


mapper3 = DataFrameMapper([
        ('MotherBC',sklearn.preprocessing.LabelEncoder()),
        ('Preg1Age_cate',sklearn.preprocessing.LabelEncoder())
      ], default=None)

X_x=X
X_x=pd.DataFrame(mapper3.fit_transform(X.copy()))

y = data['group']
le=preprocessing.LabelEncoder()
y=le.fit_transform(y)

print('train data size: ', X_x.shape)


real_results=[]
pre_results=[]

result=[]
for i in range(1,10,1):
    X_train,X_test,y_train,y_test= train_test_split(X_x, y, test_size=0.3, random_state=i*33,stratify=y)
    #use chi saure to select 20% of data
    imputer = Imputer(missing_values='NaN',strategy='median',axis=0)
    imputer = imputer.fit(X_train)
    X_train= imputer.transform(X_train)
    imputer = Imputer(missing_values='NaN',strategy='mean',axis=0)
    imputer = imputer.fit(X_test)
    X_test= imputer.transform(X_test)
    ch2 = SelectKBest(chi2, k=int(X_x.shape[1]*0.2))
    X_train = ch2.fit_transform(X_train, y_train)
    X_test=ch2.transform(X_test)
    clf = LogisticRegression(penalty='l1')
    clf.fit(X_train, y_train)
    predictions = clf.predict_proba(X_test)
    score=(roc_auc_score(y_test, predictions[:,'1'])*100)
    print(score)
    result.append(score)
    pre_results.append((predictions[:,1]))
    real_results.append(y_test)


print("AUC: %0.2f (+/- %0.2f)" % (np.mean(result, axis=0), np.std(result, axis=0) ))


draw_multi_roc(real_results,pre_results,np.mean(result, axis=0)/100,'Genetic information')
draw_multi_roc(real_results,pre_results,np.mean(result, axis=0)/100,'Clinical information')



from sklearn import metrics

import matplotlib.pyplot as plt
plt.switch_backend('agg')
colors = ['aqua', 'darkorange', 'cornflowerblue','navy','deeppink','black','red','green','orange','yellow']

def draw_multi_roc(real_results,pre_results, auc_score_input,title):
  direct='/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step38data/'
  plt.title(title)
  i=0
  for real_result, pre_result in zip(real_results, pre_results):
      fpr, tpr, threshold = metrics.roc_curve(real_result, pre_result)
      roc_auc = metrics.auc(fpr, tpr)
      if i==0:
          plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % auc_score_input)
      else:
          plt.plot(fpr, tpr, 'b',linewidth=1,color=colors[i])
      i=i+1
  plt.legend(loc = 'lower right')
  plt.plot([0, 1], [0, 1],'r--')
  plt.xlim([0, 1])
  plt.ylim([0, 1])
  plt.ylabel('True Positive Rate')
  plt.xlabel('False Positive Rate')
  plt.savefig(direct+title+'.png',dpi=1000)
  #plt.show()
#draw_roc(real_result, pre_result)