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
In Hybrid Group Lasso - IPF Lasso.ipynb and pc2Lasso, the data to use have to be loaded where indicated, as indicated. 

File neural_network.ipynb contains the entire training pipeline of randomly-generated neural networks for regression, while final_neural_regressoripynb contains the code for the automatic selection, training and testing of a neural net defined by the best parameter combination for a neural-network regression given a file containing the results from previous simulations. The data to use have to be loaded where indicated in the files, as indicated. 

File Autoencoder.ipynb provides the code for the entire training pipeline of randomly-generated autoencoders, while final_autoencoder.ipynb allows the user to adopt warm start and test_set_evaluation.ipynb contains the code to evaluate the goodness of a previously-trained model. The data to use have to be loaded where indicated in the files, as indicated. 

### License
This is free software for academic use: you can redistribute it and/or modify it under the terms of the GNU Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Public License for more details.

Giuseppe Magazzù - September 2020
