# HybridGroupIPFLasso_pc2Lasso
The presented code is the one used in the publication "Multimodal regularised linear models and flux balance analysis for mechanistic integration of omic data".
The provided files can be used to train and test the methods **pc2Lasso** and **hybrid Group-IPF Lasso**.

### Prerequisites

It is required to install the following dependencies in order to be able to run the files:
For the **hybrid Group-IPF Lasso** method
```
python >= 3.5
jupyter notebook
gurobi 
cvxpy
seaborn
numpy >= 1.15
scipy >= 1.1.0
pandas >= 0.22
matplotlib >= 2.1.2
multiprocess
```
For the **pc2Lasso** method 
```
pclasso
```
### Running the code
File worker_function.py must be in the same directory of Hybrid Group Lasso - IPF Lasso.ipynb.
In Hybrid Group Lasso - IPF Lasso.ipynb and pc2Lasso, the data to used have to be loaded where indicated, as indicated. 
