#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import cvxpy as cp
import numpy as np
import gurobipy
import time

def loss(X, Y, beta1, beta2, intercept_value):
    beta = cp.hstack([beta1, beta2])
    intercept = np.repeat(intercept_value.value, Y.shape[0])
    return cp.norm2(Y - cp.matmul(X, beta) - intercept)**2

def modalityRegularizer(beta):
    return cp.norm1(beta)

def groupRegularizer(groups, *args):
    exprs = list()
    if len(args) != len(groups):                               # if we do not want/cannot provide a lambda for each group
        rest = np.ones(len(groups) - len(args))
        args = list(args)
        args.extend(rest)
    for gr, group in enumerate(groups):
        exprs.append(args[gr] * cp.norm2(cp.hstack(group[0])))    # we assume every arg is a lambda
    return sum(exprs)

def mapper(vec, indices):
    n_groups = len(set(indices))
    groups = list()
    for i in range(1, n_groups + 1):                             # we are assuming each group is given a number starting from 1
        idxs = (np.asarray(indices) == i).nonzero()
        groups.append(([vec[idx] for idx in idxs[0]], i))
    
    return groups

def objective_fn(X, Y, beta1, beta2, lambd1, lambd2, intercept_value, indices1, indices2):
    return loss(X, Y, beta1, beta2, intercept_value) + lambd1 * modalityRegularizer(beta1) + lambd2 * modalityRegularizer(beta2) + groupRegularizer(mapper(beta1, indices1), 1) + groupRegularizer(mapper(beta2, indices2), 1)

def mse(X, Y, beta1, beta2, intercept_value):
    return (1.0 / X.shape[0]) * loss(X, Y, beta1, beta2, intercept_value).value

def optimise(worker_id, lambd1, lambd2, intercept, X_train, Y_train, X_test, 
             Y_test, flux_dimensions_index, gene_dimensions_index, indices1, indices2):
    beta1 = cp.Variable(flux_dimensions_index + 1)
    beta2 = cp.Variable(gene_dimensions_index - flux_dimensions_index)
    lambd1_par = cp.Parameter(nonneg=True)
    lambd2_par = cp.Parameter(nonneg=True)
    intercept_par = cp.Parameter()
    lambd1_par.value = lambd1
    lambd2_par.value = lambd2
    intercept_par.value = intercept
    
    print("Worker: {}, Lambda1: {}, Lambda2: {}, Intercept: {}".format(worker_id, lambd1_par.value, 
                                                                       lambd2_par.value, intercept_par.value))
    substart_time = time.time()
    problem = cp.Problem(cp.Minimize(objective_fn(X_train, Y_train, beta1, beta2, lambd1_par, lambd2_par, 
                                                  intercept_par, indices1, indices2))) 
    problem.solve(solver=cp.GUROBI)
    train_error = mse(X_train, Y_train, beta1, beta2, intercept_par)
    test_error =  mse(X_test, Y_test, beta1, beta2, intercept_par)
    total_time = (time.time() - substart_time)/60
    print("Worker: {}, Total time in minutes: {}, Train error:{}, Test error: {}".format(
        worker_id, total_time, train_error, test_error))
    return (worker_id, lambd1_par.value, lambd2_par.value, intercept_par.value, list(beta1.value), list(beta2.value), total_time, train_error, test_error)
