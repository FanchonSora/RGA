import pandas as pd
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

trainset = pd.read_csv('dataset_same_CBS.csv')
testset = pd.read_csv('testset_same_CBS.csv')

# trainset = pd.read_csv('dataset_same_AM1.csv')
# testset = pd.read_csv('testset_same_AM1.csv')

# trainset = pd.read_csv('dataset_same_CBS_no_filter.csv')
# trainset = pd.read_csv('dataset_same_AM1_no_filter.csv')
# testset = pd.read_csv('testset_same_CBS.csv')

X_train = trainset.iloc[:, 3:].values
y_train = trainset.iloc[:, 1].values

print("Error of rxn class in training set")
for cIndex in range(X_train.shape[1]):
    nonZeroIndices = np.where(X_train[:, cIndex] != 0)[0]
    print("Column", cIndex, ":", np.mean(abs(X_train[nonZeroIndices, cIndex] - y_train[nonZeroIndices])))
    print("Empty rows:", X_train.shape[0] - len(nonZeroIndices))

X_test = testset.iloc[:, 3:].values
y_test = testset.iloc[:, 1].values

print("Quantum method error training:", np.mean(np.abs(trainset.iloc[:, 2].values - y_train)))
print("Quantum method error testing:", np.mean(np.abs(testset.iloc[:, 2].values - y_test)))

print("\nError of rxn class in testing set")
for cIndex in range(X_test.shape[1]):
    nonZeroIndices = np.where(X_test[:, cIndex] != 0)[0]
    print("Column", cIndex, ":", np.mean(abs(X_test[nonZeroIndices, cIndex] - y_test[nonZeroIndices])))
    print("Empty rows:", X_train.shape[0] - len(nonZeroIndices))


# sc = StandardScaler()
# X_train = sc.fit_transform(X_train)
# X_test = sc.transform(X_test)

print("Mean y_train:", np.mean(y_train))

print(X_train.shape)

# Validation

from sklearn.model_selection import KFold

kf = KFold(n_splits=10, random_state=None, shuffle=False)

cv_results = []
cv_results2 = []

for train_index, test_index in kf.split(X_train):
    # print("TRAIN:", train_index, "TEST:", test_index)
    X_train_cur, X_test_cur = X_train[train_index], X_train[test_index]
    y_train_cur, y_test_cur = y_train[train_index], y_train[test_index]

    models = [None]

    X_train_cur_copy = X_train_cur[:]
    # print X_train_cur_copy

    for index in range(1, X_train_cur.shape[1]):
        nonZeroIndices = np.where(X_train_cur_copy[:, index] != 0)[0]

        # X_train_cur2 = X_train_cur_copy[nonZeroIndices, index - 1].reshape(-1, 1)
        X_train_cur2 = X_train_cur_copy[nonZeroIndices, :index]
        y_train_cur2 = X_train_cur_copy[nonZeroIndices, index]

        # ridge = Ridge(alpha=0.5, fit_intercept=True)
        ridge = KernelRidge()
        # ridge = KernelRidge(kernel='rbf')
        # ridge = GaussianProcessRegressor()
        ridge.fit(X_train_cur2, y_train_cur2)
        # print "Cur index:", (index + 1)
        # print ridge.coef_, ridge.intercept_
        # print "-----"
        models.append(ridge)

        zeroIndices = np.where(X_train_cur_copy[:, index] == 0)[0]

        # print "Zero indices:", zeroIndices

        if zeroIndices.size == 0:
            continue

        # X_test_cur2 = X_train_cur_copy[zeroIndices, index - 1].reshape(-1, 1)
        X_test_cur2 = X_train_cur_copy[zeroIndices, :index]
        y_test_cur2 = X_train_cur_copy[zeroIndices, index]

        y_mean = ridge.predict(X_test_cur2)

        X_train_cur_copy[zeroIndices, index] = y_mean.copy()

        # print X_test_cur2
        # print y_mean_cur2
        # print y_mean

    X_test_cur_copy = X_test_cur[:]

    for cIndex in range(1, X_test_cur_copy.shape[1]):

        zeroIndices = np.where(X_test_cur_copy[:, cIndex] == 0)[0]

        if zeroIndices.size == 0:
            continue

        # print cIndex

        # X_test_cur_copy[zeroIndices, cIndex] = models[cIndex].predict(X_test_cur_copy[zeroIndices, cIndex - 1].reshape(-1, 1))
        X_test_cur_copy[zeroIndices, cIndex] = models[cIndex].predict(X_test_cur_copy[zeroIndices, :cIndex])
        # print X_test_cur_copy[zeroIndices, cIndex - 1]
        # print X_test_cur_copy[zeroIndices, cIndex]

    # ridge = Ridge(alpha=0.5)
    ridge = KernelRidge()
    # ridge = KernelRidge(kernel='rbf')
    ridge.fit(X_train_cur_copy, y_train_cur)
    curCV = np.asarray(ridge.predict(X_test_cur_copy))

    cIndex = 5
    # ridge = Ridge(alpha=0.5)
    ridge = KernelRidge()
    # ridge = KernelRidge(kernel='rbf')
    ridge.fit(X_train_cur_copy[:, cIndex].reshape(-1, 1), y_train_cur)

    # print "y_train_cur", y_train_cur
    # print "X_train_cur_copy:", X_train_cur_copy

    curCV2 = np.asarray(ridge.predict(X_test_cur_copy[:, cIndex].reshape(-1, 1)))
    # print X_test_cur_copy[:, cIndex]
    # print curCV2
    # print y_test_cur

    cv_results.append(np.mean(np.abs(curCV - y_test_cur)))
    cv_results2.append(np.mean(np.abs(curCV2 - y_test_cur)))

print("Validation results all:", np.mean(np.asarray(cv_results)))
print("Validation results highest:", np.mean(np.asarray(cv_results2)))

# Training and testing

models = [None]

for index in range(1, X_train.shape[1]):

    nonZeroIndices = np.where(X_train[:, index] != 0)[0]

    # X_train_cur = X_train[nonZeroIndices, index - 1].reshape(-1, 1)
    X_train_cur = X_train[nonZeroIndices, :index]
    y_train_cur = X_train[nonZeroIndices, index]

    # print "Cur index:", (index + 1)

    # ridge = Ridge(alpha=0.5, fit_intercept=True)
    ridge = KernelRidge()
    # ridge = KernelRidge(kernel='linear')
    # ridge = GaussianProcessRegressor()
    ridge.fit(X_train_cur, y_train_cur)

    # print(ridge.coef_, ridge.intercept_)

    models.append(ridge)

    # print "-----"

    zeroIndices = np.where(X_train[:, index] == 0)[0]

    if zeroIndices.size == 0:
        continue

    # X_test_cur = X_train[zeroIndices, index - 1].reshape(-1, 1)
    X_test_cur = X_train[zeroIndices, :index]
    y_test_cur = X_train[zeroIndices, index]

    y_mean = ridge.predict(X_test_cur)
    X_train[zeroIndices, index] = y_mean.copy()

for cIndex in range(1, X_test.shape[1]):

    zeroIndices = np.where(X_test[:, cIndex] == 0)[0]

    if zeroIndices.size == 0:
        continue

    # X_test[zeroIndices, cIndex] = models[cIndex].predict(X_test[zeroIndices, cIndex - 1].reshape(-1, 1))
    X_test[zeroIndices, cIndex] = models[cIndex].predict(X_test[zeroIndices, :cIndex])

#  Predict HoF using all features
# ridge = Ridge(alpha=0.5, fit_intercept=True)
ridge = KernelRidge()
ridge.fit(X_train, y_train)
y_pred = ridge.predict(X_test)

# print(ridge.coef_, ridge.intercept_)

print("Testing error all rxns", mean_absolute_error(y_test, y_pred))

# Predict HoF using only last feature (hyperhomodesmotic)
# ridge = Ridge(alpha=0.5, fit_intercept=True)
ridge = KernelRidge()

cIndex = 5

ridge.fit(X_train[:, cIndex].reshape(-1, 1), y_train)
y_pred2 = ridge.predict(X_train[:, cIndex].reshape(-1, 1))
print("Training error:", mean_absolute_error(y_train, y_pred2))
y_pred2 = ridge.predict(X_test[:, cIndex].reshape(-1, 1))
# print(ridge.coef_, ridge.intercept_)

# print(y_pred2)
# print(y_test)

print("Testing error highest:", mean_absolute_error(y_test, y_pred2))

# Predict HoF directly from all rxns
# ridge = Ridge(alpha=10, fit_intercept=True)
# ridge = KernelRidge(kernel='rbf')
#
# ridge.fit(X_train[:, 0].reshape(-1, 1), y_train)
# y_pred3 = ridge.predict(X_test[:, 0].reshape(-1, 1))
# # print(ridge.coef_, ridge.intercept_)
#
# print(mean_absolute_error(y_test, y_pred3))
