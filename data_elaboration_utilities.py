#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# utilities for elaborating the data as per the code

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA as PCA
from matplotlib import pyplot as plt

# creating test and training set from scratch (one source)
def splitTrainingTestSet(path, percent_train):
    partition = dict() 
    labels = dict()
    dataset = pd.read_csv(path, encoding="utf-8")
    dataset.drop(["Row"], axis=1, inplace=True)                  # modify according to the structure of your data
    numpy_cleaned = dataset.values                               # converting to numpy array
    
    # creating the dictionary of labels
    for i in numpy_cleaned:
        labels[i[0]] = i[1]                                      

    split_idx = int(len(numpy_cleaned)*percent_train)
    partition["training"] = numpy_cleaned[0:split_idx]
    partition["test"] = numpy_cleaned[split_idx:]
    return partition, labels

# creating test, validation and training set from scratch (one source)
def splitTrainingValidationTestSet(path, percent_train, percent_validation):
    partition = dict() 
    labels = dict()
    dataset = pd.read_csv(path, encoding="utf-8")
    dataset.drop(["Row"], axis=1, inplace=True)                  # modify according to the structure of your data
    numpy_cleaned = dataset.values                               # converting to numpy array
    
    # creating the dictionary of labels
    for i in numpy_cleaned:
        labels[i[0]] = i[1]                                      
        
    split_idx = int(len(numpy_cleaned)*percent_train)
    split_idx2 = int(len(numpy_cleaned)*percent_validation) + split_idx
    partition["training"] = numpy_cleaned[0:split_idx]
    partition["validation"] = numpy_cleaned[split_idx:split_idx2] 
    partition["test"] = numpy_cleaned[split_idx2:]
    return partition, labels

# creating test and training set from two different sources
def loadTrainingTestSet(training_path, test_path):
    partition = dict() 
    labels = dict()
    for path in (training_path, test_path):
        dataset = pd.read_csv(path, encoding="utf-8")
        dataset.drop("Unnamed: 0", axis=1, inplace=True)        # modify according to the structure of your data
        dataset = dataset.values                                # converting to numpy array
    
        # creating the dictionary of labels
        for i in dataset:
            labels[i[0]] = i[1] 
        
        #creating partition
        if path == training_path:
            partition["training"] = dataset
        else:
            partition["test"] = dataset
    return partition, labels

# creating test, validation and training set from two different sources (validation is obtained from the test set)
def loadTrainingTestSetCreateValidation(training_path, test_path, percent_train, percent_validation):
    partition = dict() 
    labels = dict()
    for path in (training_path, test_path):
        dataset = pd.read_csv(path, encoding="utf-8")
        dataset.drop("Unnamed: 0", axis=1, inplace=True)        # modify according to the structure of your data
        dataset = dataset.values                                # converting to numpy array
    
        # creating the dictionary of labels
        for i in dataset:
            labels[i[0]] = i[1] 
        
        #creating partition
        if path == training_path:
            partition["training"] = dataset
        else:
            split_idx = int(len(dataset)*percent_validation/(1 - percent_train))
            partition["validation"] = dataset[:split_idx]
            partition["test"] = dataset[split_idx:]
    return partition, labels

def loadTrainingValidationTestSet(training_path, validation_path, test_path):
    partition = dict() 
    labels = dict()
    for path in (training_path, validation_path, test_path):
        dataset = pd.read_csv(path, encoding="utf-8")
        dataset.drop("Unnamed: 0", axis=1, inplace=True)        # modify according to the structure of your data
        dataset = dataset.values                                # converting to numpy array
    
        # creating the dictionary of labels
        for i in dataset:
            labels[i[0]] = i[1] 
        
        #creating partition
        if path == training_path:
            partition["training"] = dataset
        elif path == validation_path:
            partition["validation"] = dataset
        else:
            partition["test"] = dataset
    return partition, labels

# visualising how data (latent representation) evolves throughout the epochs
def visualizeLatentRepresentation(info, height, width, features, columns, rows, vis, plot): 
    
    plt.figure(figsize=(width, height))
    legend = list()
    for n in range(len(info)):
        plt.subplot(rows, columns, n+1)                     
        z, epoch = info[n]
        z = z.numpy()                                       # to use plt.plot 
        pca = PCA()                                         # PCA
        pca.fit(z);
        z_pca = pca.transform(z)
        for i in range(features):
            plt.plot(z_pca[:, i], vis)                     # vis defines clearly the rendering
            legend.append('z_pca[{}]'.format(i))
        plt.title("Epoch: {}".format(epoch))
        plt.legend(legend) 
    if plot:
        plt.show()
    
    
# visualising how losses evolve throughout the epochs
def visualizeLossesOverEpochs(first, second, third, first_name, second_name, third_name, height, width, vis, plot): 
    plt.figure(figsize=(width, height))
    plt.plot(first, vis)                     
    plt.plot(second, vis)                    
    plt.plot(third, vis)      
           
    plt.legend([first_name, second_name, third_name]) 
    plt.title("Losses evolution")
    if plot:
        plt.show()

        
# visualising the tensors saved in txt files    
def readLoss(path):                       
    previous_loss = open(path, "r")
    prev = previous_loss.read()
    prev = prev.replace("array", "").replace("(", "").replace(")", "").replace("dtype=float32", "") \
    .replace("[", "").replace("]", "").replace(" ", "").replace("tensor", "")
    prev = prev.split(",")
    compar = list()
    for j in prev:
        compar.append(j)
    first = [x for x in compar if compar.index(x)%2 == 0]
    second = [float(x) for x in first]                          # loss to compare to the current one
    return second   
    
    
    
    

