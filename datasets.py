#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# datasets of data used in the publication, change it according to the structure of your file

import torch
from torch.utils import data
import numpy as np

class YeastDataset(data.Dataset):                    
    def __init__(self, examples, labels):
        self.labels = labels
        self.examples = examples
        
    def __len__(self):
        return len(self.examples)
    
    def __getitem__(self, index):
        return self.examples[index][2:], self.labels[self.examples[index][0]] 
        # return only the features and the labels, separately

