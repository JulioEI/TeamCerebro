# TeamCerebro
### Julio Esparza, Angelica Parra, Maria Royo, Zihao Wang <br />Dept. of Bioengineering, Imperial College London, London, UK
This repository is part of BMI class at Imperial College of London. 
The aim of its functions are to evaluate different methods to decode a monkey's hand position from neural activity.
The repository is subdivided into 3 main folders:
  1. Select_parameters: This folder contains the functions used to tune the different parameters the algorithms contain as to select the most appropiate ones ad hoc.<br />
          a. **CNB_kfold.m**: evaluate CNB accuracy through 10-kfold. <br />
          b. **SVM_del_kfold.m**: evaluate SVM (with delete cells as the preprocessing step) accuracy through 10-kfold. <br />
          c. **SVM_gauss_kfold.m**: evaluate SVM (with gaussian as the preprocessing step) accuracy through 10-kfold. <br />
          d. **parameters_performance.m**: general function that calls for the different classifiers parameter testing. This function may take several hours to run due to the                  iterative nature of the parameter optimization. <br />
          e. **fillformat.m** & **plotFill1.m**: functions used for plotting. Developed by Sara Mederos and Julio Esparza at Cajal Institute (Spanish Research Council), Madrid, Spain. <br />
  2. Method_evaluation: This folder contains the functions used to test the different algorithm combinations (direction classifier + trajectory predictor). <br />
          a. **evaluate_methods.m**: evaluation across different combinations of direction (classification) and trajectory (regression) prediction methods. This requires the                    **positionEstimatorTraining.m** & **positionEstimator.m** functions located in "Train&Predict" folder. <br />
  3. Train&Predict: This folder contains the functions used to train and predict the data for each algorithm (i.e. this folder contains the functions in the competition format).<br />
          a. **positionEstimatorTraining.m**: train different models for neural decoding.<br />
          b.  **positionEstimator.m**: predicts reaching direction and movement trajectory using the trained model.<br />
          c. **testFunction_for_students_MTb.m**: splits dataset, call training and testing functions.<br />
          
Last Updated: 17/05/2021 (DD/MM/YYYY)

Happy training!
