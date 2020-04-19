

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
from sklearn.feature_selection import RFE
from sklearn.svm import SVR
from sklearn import preprocessing
from sklearn.preprocessing import Imputer
from sklearn_pandas import DataFrameMapper
from sklearn.linear_model import LogisticRegression

#load the data 
data = pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step35data/somatic_mutations_all.csv',sep=',')
data = pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step35data/somatic_mutations_all_removeSynonymous.csv',sep=',')
data = pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step35data/somatic_mutations_nonSense.csv',sep=',',encoding="utf-8")
data=data.rename(columns={ data.columns[0]: 'Sample' })
RiskFactor=pd.read_csv('/projects/p30007/Zexian/Alignment/BBCAR_NEW/administrative/Step10data/RiskFactor.csv',sep=',')

#select data as predictors
Selected_Genetics = data.drop(['Sample','group','age','menaposal'], axis=1)
Selected_Genetics=pd.concat([Selected_Genetics,RiskFactor[['Nonsense_mutation']]])
RiskFactor = RiskFactor.drop(['Study_ID','TissueStatus','BirthDate','Nonsense_mutation'], axis=1)
#data.columns.values[0]='Sample'

Selected_RiskFactor=RiskFactor
Selected_RiskFactor=RiskFactor[['Benign_Age','Class','AgeofMerche','FirstDegreeRelDiagBreastCancer','Preg1Age_cate']]
Selected_RiskFactor=RiskFactor[['Benign_Age','Class','AgeofMerche','FirstDegreeRelDiagBreastCancer','Preg1Age_cate','Nonsense_mutation']]

#combine data for prediction
X=Selected_RiskFactor
X=Selected_Genetics
X=pd.concat([Selected_RiskFactor.reset_index(drop=True),Selected_Genetics], axis=1)

imputer = Imputer(missing_values='NaN',strategy='median',axis=0)
X["Benign_Age"] = imputer.fit_transform(X[["Benign_Age"]]).ravel()
imputer = Imputer(missing_values='NaN',strategy='median',axis=0)
X["AgeofMerche"] = imputer.fit_transform(X[["AgeofMerche"]]).ravel()


mapper3 = DataFrameMapper([
	        (['Benign_Age'],sklearn.preprocessing.StandardScaler()),
        (['AgeofMerche'],sklearn.preprocessing.StandardScaler()),
        ('Class',sklearn.preprocessing.LabelEncoder()),
        ('FirstDegreeRelDiagBreastCancer',sklearn.preprocessing.LabelEncoder()),
        ('Preg1Age_cate',sklearn.preprocessing.LabelEncoder())
      ], default=None)

X_x=pd.DataFrame(mapper3.fit_transform(X.copy()))

y = data['group']
le=preprocessing.LabelEncoder()
y=le.fit_transform(y)

#the dimension of predictors
print('train data size: ', X_x.shape)
# summarize the number of rows and columns in the dataset
print(X.shape)

X_train,X_test,y_train,y_test= train_test_split(X_x, y, test_size=0.3, random_state=1,stratify=y)

imputer = Imputer(missing_values='NaN',strategy='median',axis=0)
imputer = imputer.fit(X_train)
X_train= imputer.transform(X_train)

imputer = Imputer(missing_values='NaN',strategy='mean',axis=0)
imputer = imputer.fit(X_test)
X_test= imputer.transform(X_test)


#five fold cross validation
clf = svm.SVC(kernel='linear', C=1)
clf=linear_model.Lasso(alpha=0.1)
clf=LogisticRegression(random_state=0, solver='lbfgs')
scores = cross_val_score(clf, X_train, y_train, cv=3,scoring='roc_auc')
print scores
print("AUC: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
#held-out test
clf.fit(X_train, y_train)
predictions = clf.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)


predictions = clf.predict(X_test)
print(roc_auc_score(y_test, predictions)*100)



print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)
print(accuracy_score(y_test, predictions)*100)
predictions = best_linear_model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)



#search for models 
#try logistic regression 
linear_pipeline = Pipeline([('std_scaler', StandardScaler()),
                            ("model", SGDClassifier( class_weight='balanced', n_iter=3000, random_state=131))
                           ])
loss_options = [ 'log']
penalty_options = ['l1','l2', 'elasticnet']
linear_params = {
        'model__loss': loss_options,
        'model__penalty': penalty_options
    }
grid = GridSearchCV(linear_pipeline, cv=5, n_jobs=-1, verbose=1, param_grid=linear_params, scoring='roc_auc')

grid.fit(X_train, y_train)
print('crossvalidation mean is: '+str(grid.best_score_*100))
print('crossvalidation sd is: '+str(grid.cv_results_['std_test_score'][grid.best_index_]*100))

best_linear_model = grid.best_estimator_
print('Best Linear Model: \n')
print(best_linear_model)
predictions = best_linear_model.predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)
print(accuracy_score(y_test, predictions)*100)
predictions = best_linear_model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)


#train SVM
linear_pipeline = Pipeline([('std_scaler', StandardScaler()),
                            ("model", SGDClassifier( class_weight='balanced', n_iter=3000, random_state=1314))
                           ])
loss_options = ['hinge', 'squared_hinge', 'perceptron']
penalty_options = ['l2', 'elasticnet']
linear_params = {
        'model__loss': loss_options,
        'model__penalty': penalty_options
    }
grid = GridSearchCV(linear_pipeline, cv=5, n_jobs=-1, verbose=1, param_grid=linear_params, scoring='roc_auc')


grid.fit(X_train, y_train)
print('crossvalidation mean is: '+str(grid.best_score_*100))
print('crossvalidation sd is: '+str(grid.cv_results_['std_test_score'][grid.best_index_]*100))

best_linear_model = grid.best_estimator_
print('Best Linear Model: \n')
print(best_linear_model)
predictions = best_linear_model.predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)
print(accuracy_score(y_test, predictions)*100)

predictions = best_linear_model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)

#start the nonlinear 


# Returns the best paramsuration for a model using crosvalidation and grid search
def best_params(model, name, parameters, X, y):
    print('Start search: ', name)
    grid = GridSearchCV(model, parameters, cv=5,n_jobs=-1,
                       scoring="f1_macro", verbose=1)
    grid.fit(X, y)
    best_estimator = grid.best_estimator_
    print('Best model: ', best_estimator)
    return [name, grid.best_score_, grid.cv_results_['std_test_score'][grid.best_index_], best_estimator]

# Returns the best model from a set of model families given training data using cross-validation.
def best_model(classifier_group, X, y):
    best_score = 0.0
    best_classifier = None
    classifiers = []
    for name, model, parameters in classifier_group:
        classifiers.append(best_params(model, name, parameters, X, y))
    for name, score, std, classifier in classifiers:
        print('Considering classifier... ' + name)
        if (score > best_score):
            best_score = score
            best_classifier = [name, classifier,score,std]
    return best_classifier
 
#random forest
def nonlinear_models():
    models = []
    rf_tuned_parameters = [{'criterion': ['gini', 'entropy'],
                            'max_depth' : [4, 6,8, 10,12,14,16],
                             'max_features': ['auto', 'sqrt', 'log2'],
                            'min_samples_leaf' : [1, 5, 10, 20]}]
    models.append(['RandomForest', RandomForestClassifier(n_jobs=-1), rf_tuned_parameters])
    return models

model=best_model(nonlinear_models(), X_train, y_train)
print('Best  Model: \n')
print(model[1])
print('crossvalidation mean is: '+str(model[2]*100))
print('crossvalidation sd is: '+str(model[3]*100))

predictions = model[1].predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)
predictions = best_linear_model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)

#GradientBoostingTree
def nonlinear_models():
    models = []
    gbt_tuned_parameters = [{       "loss":["deviance"],
                                    #"learning_rate": [0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2],
                                    "min_samples_split": np.linspace(0.1, 0.5, 12),
                                    #"min_samples_leaf": np.linspace(0.1, 0.5, 12),
                                    "max_depth":[3,5,8], 
                                    'min_samples_split':range(200,1001,200),
                                    #"max_features":["log2","sqrt"],
                                    #"criterion": ["friedman_mse",  "mae"],
                                    #"subsample":[0.5, 0.618, 0.8, 0.85, 0.9, 0.95, 1.0],
                                    "n_estimators":[10]
                                    }]
    models.append(['GradientBoostingTree', GradientBoostingClassifier(), gbt_tuned_parameters])
    return models

model=best_model(nonlinear_models(), X_train, y_train)
print('Best  Model: \n')
print(model[1])
print('crossvalidation mean is: '+str(model[2]*100))
print('crossvalidation sd is: '+str(model[3]*100))

predictions = model[1].predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)
predictions = best_linear_model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)

#kNN
def nonlinear_models():
    models = []
    knn_tuned_parameters = [{"n_neighbors": [5, 10, 20, 50]}]
    models.append(['kNN', KNeighborsClassifier(n_jobs=-1), knn_tuned_parameters])   
    return models

model=best_model(nonlinear_models(), X_train, y_train)
print('Best  Model: \n')
print(model[1])
print('crossvalidation mean is: '+str(model[2]*100))
print('crossvalidation sd is: '+str(model[3]*100))

predictions = model[1].predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)
predictions = best_linear_model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)

#SVM
def nonlinear_models():
    models = []
    SVM_tuned_params = [{'C': [ 1, 10, 100, 1000], 
                        'gamma':[ 1e-3, 1e-4],
                        'kernel': ['rbf']}]
    models.append(['SVM', svm.SVC(),  SVM_tuned_params]) 
    return models

model=best_model(nonlinear_models(), X_train, y_train)
print('Best  Model: \n')
print(model[1])
print('crossvalidation mean is: '+str(model[2]*100))
print('crossvalidation sd is: '+str(model[3]*100))

predictions = model[1].predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)
predictions = best_linear_model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)

#MLP
def nonlinear_models():
    models = []
    mlp_tuned_parameters = [{'solver': ['lbfgs'], 
                            # 'max_iter': [500,1000,1500],
                             'learning_rate': ['invscaling', 'adaptive'],
                             'activation': ['relu', 'logistic'],
                             'random_state':[0,1,2,3,4,5,6,7,8,9],
                             'alpha': 10.0 ** -np.arange(1, 7),
                             "hidden_layer_sizes": [(10,), (15,), (20,), (10,5), (15,5), (20,5)]}]
    models.append(['MLP', MLPClassifier(random_state=1314),  mlp_tuned_parameters])
    return models


def nonlinear_models():
    models = []
    mlp_tuned_parameters = [{'solver': ['lbfgs'], 
                            # 'max_iter': [500,1000,1500],
                             'learning_rate': ['invscaling'],
                             'activation': [ 'logistic'],
                             'random_state':[8],
                             'alpha': 10.0 ** -np.arange(3, 7),
                             "hidden_layer_sizes": [(10,5)]}]
    models.append(['MLP', MLPClassifier(random_state=1314),  mlp_tuned_parameters])
    return models


model=best_model(nonlinear_models(), X_train, y_train)
print('Best  Model: \n')
print(model[1])
print('crossvalidation mean is: '+str(model[2]*100))
print('crossvalidation sd is: '+str(model[3]*100))

predictions = model[1].predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)
print(roc_auc_score(y_test, predictions)*100)

predictions = model[1].predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)





#selected model
model=MLPClassifier(activation='logistic', alpha=9.9999999999999995e-07,
       batch_size='auto', beta_1=0.9, beta_2=0.999, early_stopping=False,
       epsilon=1e-08, hidden_layer_sizes=(10, 5),
       learning_rate='invscaling', learning_rate_init=0.001, max_iter=1500,
       momentum=0.9, nesterovs_momentum=True, power_t=0.5, random_state=8,
       shuffle=True, solver='lbfgs', tol=0.0001, validation_fraction=0.1,
       verbose=False, warm_start=False)
model.fit(X_train, y_train)
predictions = model.predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions)*100)
print(recall_score(y_test, predictions)*100)
print(f1_score(y_test, predictions)*100)

model.fit(X_train, y_train)
scores=cross_val_score(model, X_train, y_train, cv=5, scoring='precision')
scores.mean(),stats.sem(scores)

predictions = model.predict(X_test)
print('----------------\nTest data prediction:\n')
print(precision_score(y_test, predictions))
print(recall_score(y_test, predictions))
print(f1_score(y_test, predictions))

predictions = model.predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)




def plot_roc_auc(actual, prediction):
    fpr, tpr, thresholds = roc_curve(actual, prediction[:,1])
    plt.plot(fpr, tpr,'r')
    plt.plot([0,1],[0,1],'b')
    plt.title('AUC: {}'.format(auc(fpr,tpr)))
    plt.show()  
    
plot_roc_auc(y_test, best_linear_model_rm_vif.predict_proba(X_test_rm[cols_rm_vif]))









#old code



    gbt_tuned_parameters = [{"n_estimators": [50,200,400,1000],
                             'learning_rate': [1,0.5,0.1, 0.05, 0.01],
                              'max_depth': [4, 8, 16],
                              'min_samples_leaf': [1, 5, 10, 20]}]
    models.append(['GradientBoostingTree', GradientBoostingClassifier(), gbt_tuned_parameters])
    #kNN
    
    #SVM
    # (no of inputs + no of outputs) ^ 0.5 + range(1,10)
    SVM_tuned_params = [{'C': [0.001, 0.01, 0.1, 1, 10], 
                        'gamma':[0.001, 0.01, 0.1, 1],
                        'kernel': ['rbf']}]
    models.append(['SVM', svm.SVC(),  SVM_tuned_params])
    #LinearSVM
    SVM_tuned_params = [{'C': [0.001, 0.01, 0.1, 1, 10], 
                        'gamma':[0.001, 0.01, 0.1, 1],
                        'kernel': ['rbf']}]
    models.append(['SVM', svm.SVC(),  SVM_tuned_params])
    #MLP
    mlp_tuned_parameters = [{'solver': ['lbfgs'], 
                            # 'max_iter': [500,1000,1500],
                             'learning_rate': ['invscaling', 'adaptive'],
                             'activation': ['relu', 'logistic'],
                             'random_state':[0,1,2,3,4,5,6,7,8,9],
                             'alpha': 10.0 ** -np.arange(1, 7),
                             "hidden_layer_sizes": [(10,), (15,), (20,), (10,5), (15,5), (20,5)]}]
    models.append(['MLP', MLPClassifier(random_state=1314),  mlp_tuned_parameters])




def nonlinear_models():
    models = []
    mlp_tuned_parameters = [{'solver': ['lbfgs'], 
                            # 'max_iter': [500,1000,1500],
                             'learning_rate': ['invscaling'],
                             'activation': [ 'logistic'],
                             'random_state':[8],
                             'alpha': [1.00000000e-06],
                             "hidden_layer_sizes": [ (10,5)]}]
    models.append(['MLP', MLPClassifier(random_state=1314),  mlp_tuned_parameters])
    return models

# In[ ]:


#predictions = best_linear_model.predict(X_train)
#print('----------------\nTrain data prediction:\n')
#print(classification_report(y_train, predictions))
#print(roc_auc_score(y_train, predictions))
scores=cross_val_score(model, X_train, y_train, cv=5, scoring='f1')
scores.mean(),stats.sem(scores)

scores=cross_val_score(model, X_train, y_train, cv=5, scoring='precision')
scores.mean(),stats.sem(scores)

predictions = model.predict(X_test)
print('----------------\nTest data prediction:\n')
#print(classification_report(y_test, predictions))
print(precision_score(y_test, predictions))
print(recall_score(y_test, predictions))
print(f1_score(y_test, predictions))
print(roc_auc_score(y_test, predictions))



model=MLPClassifier(activation='logistic', alpha=9.9999999999999995e-07,
       batch_size='auto', beta_1=0.9, beta_2=0.999, early_stopping=False,
       epsilon=1e-08, hidden_layer_sizes=(10, 5),
       learning_rate='invscaling', learning_rate_init=0.001, max_iter=1500,
       momentum=0.9, nesterovs_momentum=True, power_t=0.5, random_state=8,
       shuffle=True, solver='lbfgs', tol=0.0001, validation_fraction=0.1,
       verbose=False, warm_start=False)
model.fit(X_train, y_train)
scores=cross_val_score(model, X_train, y_train, cv=5, scoring='precision')
scores.mean(),stats.sem(scores)

predictions = model.predict(X_test)
print('----------------\nTest data prediction:\n')
print(precision_score(y_test, predictions))
print(recall_score(y_test, predictions))
print(f1_score(y_test, predictions))

predictions = model[1].predict_proba(X_test)
print(roc_auc_score(y_test, predictions[:,'1'])*100)



cm1 = confusion_matrix(y_test, predictions)
total1=sum(sum(cm1))
accuracy1=(cm1[0,0]+cm1[1,1])/total1
print ('Accuracy : ', accuracy1)

sensitivity1 = cm1[0,0]/(cm1[0,0]+cm1[0,1])
print('Sensitivity : ', sensitivity1 )

specificity1 = cm1[1,1]/(cm1[1,0]+cm1[1,1])
print('Specificity : ', specificity1)

TP = np.diag(cm1)
FP = cm1.sum(axis=0) - np.diag(cm1)  
PPV = TP/(TP+FP)
print('PPV : ', PPV)






def plot_roc_auc(actual, prediction):
    fpr, tpr, thresholds = roc_curve(actual, prediction[:,1])
    plt.plot(fpr, tpr,'r')
    plt.plot([0,1],[0,1],'b')
    plt.title('AUC: {}'.format(auc(fpr,tpr)))
    plt.show()  
    
plot_roc_auc(y_test, best_linear_model_rm_vif.predict_proba(X_test_rm[cols_rm_vif]))
