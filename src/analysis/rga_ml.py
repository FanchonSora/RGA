from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

trainset = pd.read_csv('dataset_v5.csv')
testset = pd.read_csv('testset_v4.csv')

X_train = trainset.iloc[:, 2:].values
y_train = trainset.iloc[:, 1].values

X_test = testset.iloc[:, 3:].values
y_test = testset.iloc[:, 1].values

for i in range(X_train.shape[1]):
    curCol = X_train[:, i]
    diffs = []
    for j in range(len(curCol)):
        if curCol[j] != 0:
            diffs.append(abs(curCol[j] - y_train[j]))

    print("Col", i, ":", np.mean(diffs))

# sc = StandardScaler()
# X_train = sc.fit_transform(X_train)
# X_test = sc.transform(X_test)

# sc_mean = np.mean(X_train, axis=0)
# sc_std = np.std(X_train, axis=0)
#
# X_train = (X_train - sc_mean) / (sc_std * 1.0)
# X_test = (X_test - sc_mean) / (sc_std * 1.0)

print("Mean y_train:", np.mean(y_train))

ridge = Ridge(alpha=0.5, fit_intercept=True)
# ridge = KernelRidge(alpha=0.0005, kernel='rbf')
ridge.fit(X_train, y_train)

print(ridge.coef_, ridge.intercept_)

y_pred = ridge.predict(X_train)

print(r2_score(y_true=y_train, y_pred=y_pred))
print(mean_absolute_error(y_true=y_train, y_pred=y_pred))
print(mean_squared_error(y_true=y_train, y_pred=y_pred))

ytest_pred = []

for row in X_test:
    curW = 0
    curVal = 0
    for i in range(len(row)):
        col = row[i]
        if col != 0:
            curW += ridge.coef_[i]
            curVal += ridge.coef_[i] * col
            # curW += 1
            # curVal += col
    ytest_pred.append(curVal / (curW * 1.0))

ytest_pred = np.asarray(ytest_pred)
print(ytest_pred)
print(y_test)

# ytest_pred = ridge.predict(X_test)
#
# print(r2_score(y_true=y_test, y_pred=ytest_pred))
# print(mean_absolute_error(y_true=y_test, y_pred=ytest_pred))
# print(mean_squared_error(y_true=y_test, y_pred=ytest_pred))

print("Mean cross-validation absolute error:", abs(cross_val_score(ridge, X_train, y_train, cv=10, scoring='neg_mean_absolute_error').mean()))
# print abs(cross_val_score(ridge, X_train, y_train, cv=10, scoring='neg_mean_squared_error').mean())

print(np.mean(abs(ytest_pred - y_test)))
print(np.mean(abs(testset.iloc[:, 2].values - y_test)))
