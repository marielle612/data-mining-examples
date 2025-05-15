A program was developed to solve a multi-class classification problem using Bagging + QDA, based on user-specified data and variable input, with automated performance evaluation via cross-validation or a test dataset.

Consider only a classification problem. That is, there is a variable which indicates classes. The location of the class variable is not fixed. Make your program to handle more than two classes. You can assume that values of the class variable are integers starting with 1. Assume your data has both numerical and categorical variables. Further assume that the categorical variables are coded as integers starting with 1. 

1.	Prompt the user to type in the filename of the training data. 
2.	Prompt the user to enter the locations of the categorical variables and the class variables.
3.	Prompt the user whether a cross-validation is to be used or a test dataset is to be used as the model evaluation tool.
4.	If a cross-validation is chosen, perform v-fold cross-validation.  You need to ask the used to enter the value ‘v’.
5.	If a test dataset is chosen, ask the user the filename of the test dataset. (Assume the column location of the class variable is as same as that of the training dataset.)
6.	Perform Bagging depending on the choice by the user. For categorical variables, the program needs to create dummy variables to replace them automatically.
7.	For Bagging Ensemble method, use QDA as the classifier and 51 bootstraps as the number of re-sampled dat
