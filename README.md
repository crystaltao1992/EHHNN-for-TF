# EHHNN-for-TF
# Read me

This is the demo for efficient hinging hyperplanes neural networks (EHHNN) for traffic flow (TF) prediction.

(qtao@esat.kuleuven.be and zueslee.hitsz@foxmail.com)

## Main files

Three categories are mainly included: data, EHHNN as the TF Predictor, and the analysis of variance (ANOVA) decomposition for the  predictor.

### Data

The traffic data can be accessed online via https://pems.dot.ca.gov/ (PeMS), where the traffic factors of TF, average speed of vehicles (AVS) and road occupancy (OR) are adopted in our work for prediction and also analysis. 

The data pre-processing codes to process the raw data (.xlsx files) downloaded from (PeMS) are given together with 'process_data.py', where 'README.txt' gives the explanation. In addition, the processed datasets are uploaded for an illustration (see Section IV), where the datasets are named based on the rules in 'README.txt'.

### EHHNN for TF prediction

The EHHNN predictor is trained with the traffic inputs formulated in our work, which involve temporal information (with lags of time-series data) and spatial information (e.g., TF, AVS, and RO collected from different detectors locating in different road segments).

Run 'main_demo_single_layer.m' to do the prediction with single-layered EHHNN, where the ANOVA decomposition is done for variable selection and interpretation analysis is also given for an illustration.

Run 'main_demo.m' to do the prediction with the selected variables.

### ANOVA decomposition

The EHHNN predictor possesses sparse neuron connections, which make the variables and their interaction decomposable with ANOVA decomposition. With  the ANOVA decomposition and its corresponding $\sigma$​​ values, the relative importance of different so-called ANOVA function (e.g.,  traffic variables and their interactions).

See 'anova_ehh.m' applying ANOVA decomposition to the EHHNN predictor, and it is run in 'main_demo_single_layer.m' for varied variable selections and also in 'main_demo.m' after obtaining the predictor. With the ANOVA results, varied analysis can be done, such as spatial-temporal analysis (up to the users),  etc. The last block in  'main_demo_single_layer.m' gives an illustration of analysis.



 

 



 
