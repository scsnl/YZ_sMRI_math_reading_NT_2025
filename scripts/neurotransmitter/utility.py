"""
Utility functions for neurotransmitter analysis
"""

import threadpoolctl
from brainsmash import mapgen
from brainspace.null_models import moran
import numpy as np
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy.stats import pearsonr
from scipy.spatial.distance import squareform, pdist

project_dir = '/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/'
# dist_fname = project_dir + 'scripts/spin/BNA_distmat_218.txt'
# n_perm = 5000  # number of permutations for null models
seed = 10000  # reproducibility
n_proc = 4  # number of parallel workers for surrogate generation


def make_surrogates(data, dist_fname, spatnull, n_perm):
    """
    Generates surrogates for `data` using `spatnull` method

    Parameters
    ----------
    data : (N,) pd.DataFrame
    dist_fname: distance matrix filename
    spatnull : {'burt2020', 'moran'}

    Returns
    -------
    surrogates : (N, `N_PERM`) np.ndarray
    """

    if spatnull not in ('burt2020', 'moran'):
        raise ValueError(f'Cannot make surrogates for null method {spatnull}')

    darr = np.asarray(data)
    dist = np.loadtxt(dist_fname)
    # print(dist.shape)
    # print(dist[0,:])

    surrogates = np.zeros((len(data), n_perm))

    if spatnull == 'burt2020':
        surrogates = mapgen.Base(darr, dist, seed=seed, n_jobs=n_proc)(n_perm).T
    elif spatnull == 'moran':
        mrs = moran.MoranRandomization(joint=True, n_rep=n_perm, tol=1e-6, random_state=seed)
        with threadpoolctl.threadpool_limits(limits=n_proc):
            np.fill_diagonal(dist, 1)
            dist **= -1
            surrogates = mrs.fit(dist).randomize(darr).T

    outputf = project_dir + 'scripts/spin/surrogatebrainmap_' + spatnull + '_nperm' + str(n_perm) + '_seed' + str(seed) + '.txt'
    np.savetxt(outputf, surrogates)

    return surrogates


def get_reg_r_sq(X, y):

    lin_reg = LinearRegression()
    lin_reg.fit(X, y)

    yhat = lin_reg.predict(X)
    SS_Residual = sum((y - yhat) ** 2)
    SS_Total = sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (float(SS_Residual)) / SS_Total
    adjusted_r_squared = 1 - (1 - r_squared) * \
        (len(y) - 1) / (len(y) - X.shape[1] - 1)

    return adjusted_r_squared, lin_reg.coef_


def get_reg_r_pval(X, y, surrogates):
    # print(surrogates.shape)
    emp = get_reg_r_sq(X, y)
    nspins = surrogates.shape[1]
    null = np.zeros((nspins, ))
    for s in range(nspins):
        null[s] = get_reg_r_sq(X, surrogates[:,s])
    return (1 + sum(null > emp))/(nspins + 1)


def get_cook_distance(X, y):

    X_with_const = sm.add_constant(X)
    # Fit the OLS model using statsmodels
    ols_model = sm.OLS(y, X_with_const).fit()
    # Get influence object to calculate Cook's distance
    influence = ols_model.get_influence()
    # Get Cook's distance
    cooks_d, p_values = influence.cooks_distance

    return cooks_d, p_values

def get_cook_distance_spin_pval(X, y, surrogates):
    # print(surrogates.shape)
    emp, _ = get_cook_distance(X, y)
    nspins = surrogates.shape[1]
    nroi = surrogates.shape[0]
    null = np.zeros((nroi, nspins))

    for s in range(nspins):
        cookd_null, _ = get_cook_distance(X, surrogates[:,s])
        null[:, s] = cookd_null

    comparison = null > emp[:, np.newaxis]

    return (1 + np.sum(comparison, axis=1))/(nspins + 1)


# def get_perm_p(emp, null):
#     return (1 + sum(abs(null - np.mean(null))
#                     > abs(emp - np.mean(null)))) / (len(null) + 1)

def cv_slr_distance_dependent(X, y, coords, train_pct=.75, metric='rsq'):
    '''
    cross validates linear regression model using distance-dependent method.
    X = n x p matrix of input variables
    y = n x 1 matrix of output variable
    coords = n x 3 coordinates of each observation
    train_pct (between 0 and 1), percent of observations in training set
    metric = {'rsq', 'corr'}
    '''

    P = squareform(pdist(coords, metric="euclidean"))
    train_metric = []
    test_metric = []

    for i in range(len(y)):
        distances = P[i, :]  # for every node
        idx = np.argsort(distances)

        train_idx = idx[:int(np.floor(train_pct * len(coords)))]
        test_idx = idx[int(np.floor(train_pct * len(coords))):]

        mdl = LinearRegression()
        mdl.fit(X[train_idx, :], y[train_idx])
        if metric == 'rsq':
            # get r^2 of train set
            train_metric.append(get_reg_r_sq(X[train_idx, :], y[train_idx]))

        elif metric == 'corr':
            rho, _ = pearsonr(mdl.predict(X[train_idx, :]), y[train_idx])
            train_metric.append(rho)

        yhat = mdl.predict(X[test_idx, :])
        if metric == 'rsq':
            # get r^2 of test set
            SS_Residual = sum((y[test_idx] - yhat) ** 2)
            SS_Total = sum((y[test_idx] - np.mean(y[test_idx])) ** 2)
            r_squared = 1 - (float(SS_Residual)) / SS_Total
            adjusted_r_squared = 1-(1-r_squared)*((len(y[test_idx]) - 1) /
                                                  (len(y[test_idx]) -
                                                   X.shape[1]-1))
            test_metric.append(adjusted_r_squared)

        elif metric == 'corr':
            rho, _ = pearsonr(yhat, y[test_idx])
            test_metric.append(rho)

    return train_metric, test_metric


# YZ added
# get p values of empirical CV test mean correlation against spin test CV test man correlation
def get_p_cv_slr_distance_dependent(X, y, coords, surrogates, train_pct=.75, metric='corr'):
    emp_train_metric, emp_test_metric = cv_slr_distance_dependent(X, y, coords, train_pct=train_pct, metric=metric)
    emp_test_mean = np.mean(emp_test_metric, axis=0)
    nspins = surrogates.shape[1]
    null = np.zeros((nspins, ))
    for s in range(nspins):
        null_train_metric, null_test_metric = cv_slr_distance_dependent(X, surrogates[:,s], coords, train_pct=train_pct, metric=metric)
        null[s] = np.mean(null_test_metric, axis=0)

    spin_p = (1 + sum(null > emp_test_mean)) / (nspins + 1)

    return emp_test_mean, null, spin_p





