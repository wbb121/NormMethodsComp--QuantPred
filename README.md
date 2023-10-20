# NormMethodsComp--QuantPred
Evaluation of normalization methods for predicting quantitative phenotypes in metagenomic data analysis

### [helper.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/helper.R):
Including functions to merge count tables, normalize the data, predict quantitative phenotypes in simulated data, and predict quantitative phenotypes in real data.

### [sim_scenario1.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/sim_scenario1.R) / [sim_scenario2.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/sim_scenario2.R) / [sim_scenario3.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/sim_scenario3.R):
Validate the performance of different normalization methods on quantitative phenotype prediction in simulation scenarios 1 / 2 / 3.
Use the following command to run the script:
``` shell
Rscript sim_scenario1.R norm_method pred_method
Rscript sim_scenario2.R norm_method pred_method
Rscript sim_scenario3.R norm_method pred_method
```
- `norm_method`: normalization method
- `pred_method`: prediction method

### [curatedMetadgenomicData.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/curatedMetadgenomicData.R):
Obtain metadata and count tables of datasets from curatedMetagenomicData according to the inclusion criteria.

### [curatedMetagenomicData_analysis.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/curatedMetagenomicData_analysis.R):
Validate the performance of different normalization methods on quantitative phenotype prediction in datasets from curatedMetagenomicData.
Use the following command to run the script:
``` shell
Rscript curatedMetagenomicData_analysis.R meta count norm_method pred_method pred_cluster
```
- `meta`: metadata for included datasets
- `count`: count table for included datasets
- `norm_method`: normalization method
- `pred_method`: prediction method
- `pred_cluster`: number of clusters used for predictions

### [summ_pred_res.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/summ_pred_res.R):
Summarize the prediction results from simulation and real data.

### [draft_figures.R](https://github.com/wbb121/NormMethodsComp--QuantPred/blob/main/draft_figures.R)
Draw the figures in the manuscript.

