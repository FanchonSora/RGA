from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

trainset = pd.read_csv('dataset_same_CBS.csv')
testset = pd.read_csv('testset_same_CBS.csv')

# trainset = pd.read_csv('dataset_same_AM1.csv')
# testset = pd.read_csv('testset_same_AM1.csv')

X_train = trainset.iloc[:, 3:].values
y_train = trainset.iloc[:, 1].values

X_test = testset.iloc[:, 3:].values
y_test = testset.iloc[:, 1].values

# sc = StandardScaler()
# X_train = sc.fit_transform(X_train)
# X_test = sc.transform(X_test)

print("Mean y_train:", np.mean(y_train))

models = []

for index in range(X_train.shape[1]):
    nonZeroIndices = np.where(X_train[:, index] != 0)[0]
    X_train_cur = X_train[nonZeroIndices, :index + 1]
    y_train_cur = y_train[nonZeroIndices]

    # ridge = Ridge(alpha=0.5, fit_intercept=True)
    ridge = KernelRidge()
    # ridge = KernelRidge(alpha=0.5, kernel='rbf')
    ridge.fit(X_train_cur, y_train_cur)
    print("Cur index:", (index + 1))
    # print ridge.coef_, ridge.intercept_
    print("-----")
    models.append(ridge)

ytrain_pred = []

for row in X_train:
    for index in range(X_train.shape[1] - 1, -1, -1):
        if row[index] != 0:
            ytrain_pred.append(models[index].predict(row[:index + 1].reshape(1, -1))[0])
            # ytrain_pred.append(np.mean([models[i].predict(row[:i + 1].reshape(1, -1))[0] for i in range(index + 1)]))

            # totalVar = 0
            # for i in range(index + 1):
            #     totalVar += np.std(row[:i + 1])
            #
            # pred = 0
            # for i in range(index + 1):
            #     pred += models[i].predict(row[:i + 1].reshape(1, -1))[0] * np.std(row[:i + 1]) / (totalVar * 1.0)
            #
            # ytrain_pred.append(pred)
            break

ytrain_pred = np.asarray(ytrain_pred)

print(np.mean(abs(ytrain_pred - y_train)))

ytest_pred = []

for row in X_test:
    for index in range(X_train.shape[1] - 1, -1, -1):
        if row[index] != 0:
            ytest_pred.append(models[index].predict(row[:index + 1].reshape(1, -1))[0])
            # ytest_pred.append(np.mean([models[i].predict(row[:i + 1].reshape(1, -1))[0] for i in range(index + 1)]))

            # totalVar = 0
            # for i in range(index + 1):
            #     totalVar += np.std(row[:i + 1])
            #
            # pred = 0
            # for i in range(index + 1):
            #     pred += models[i].predict(row[:i + 1].reshape(1, -1))[0] * np.std(row[:i + 1]) / (totalVar * 1.0)
            #
            # ytest_pred.append(pred)
            break

from sklearn.model_selection import KFold

kf = KFold(n_splits=10, random_state=None, shuffle=False)

cv_results = []

for train_index, test_index in kf.split(X_train):
    # print("TRAIN:", train_index, "TEST:", test_index)
    X_train_cur, X_test_cur = X_train[train_index], X_train[test_index]
    y_train_cur, y_test_cur = y_train[train_index], y_train[test_index]

    models = []

    for index in range(X_train_cur.shape[1]):
        nonZeroIndices = np.where(X_train_cur[:, index] != 0)[0]
        # X_train_cur2 = X_train_cur[nonZeroIndices, index - 1].reshape(-1, 1)
        X_train_cur2 = X_train_cur[nonZeroIndices, :index + 1]
        y_train_cur2 = y_train_cur[nonZeroIndices]

        # ridge = Ridge(alpha=0.5, fit_intercept=True)
        ridge = KernelRidge()

        ridge.fit(X_train_cur2, y_train_cur2)
        # print "Cur index:", (index + 1)
        # print ridge.coef_, ridge.intercept_
        # print "-----"
        models.append(ridge)

    curCV = []

    for row in X_test_cur:
        for index in range(X_train_cur.shape[1] - 1, -1, -1):
            if row[index] != 0:
                curCV.append(models[index].predict(row[:index + 1].reshape(1, -1))[0])
                # curCV.append(np.mean([models[i].predict(row[:i + 1].reshape(1, -1))[0] for i in range(index + 1)]))

                # totalVar = 0
                # for i in range(index + 1):
                #     totalVar += np.std(row[:i + 1])
                #
                # pred = 0
                # for i in range(index + 1):
                #     pred += models[i].predict(row[:i + 1].reshape(1, -1))[0] * np.std(row[:i + 1]) / (totalVar * 1.0)
                #
                # curCV.append(pred)

                break

    curCV = np.asarray(curCV)

    cv_results.append(np.mean(abs(curCV - y_test_cur)))

print(np.mean(np.asarray(cv_results)))

# print "Mean cross-validation absolute error:", abs(cross_val_score(ridge, X_train, y_train, cv=10, scoring='neg_mean_absolute_error').mean())
# print abs(cross_val_score(ridge, X_train, y_train, cv=10, scoring='neg_mean_squared_error').mean())

ytest_pred = np.asarray(ytest_pred)

# print ytest_pred
# print y_test

print(np.mean(abs(ytest_pred - y_test)))
print(np.mean(abs(testset.iloc[:, 2].values - y_test)))
