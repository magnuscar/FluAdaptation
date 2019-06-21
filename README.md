# Profiling the Host ANP32A Splicing Landscape to Predict Influenza A Virus Polymerase Adaptation

**Authors:**
Patricia Domingues, Davide Eletto, Carsten Magnus, Hannah L. Turkington, Stefan Schmutz, Osvaldo Zagordi, Matthias Lenk, Michael Huber, Martin Beer, Silke Stertz,  Benjamin G. Hale


-----

## Description

This repository contains code and data to repeat and extend the analyses of "Profiling the Host ANP32A Splicing Landscape to Predict Influenza A Virus Polymerase Adaptation" by Domingues, Eletto, Magnus et al (2019). In addition, it contains all necessary code to generate the [shiny app](https://magnuscar.shinyapps.io/FluAdaptation/) on your machine. 
The analysis was implemented in the freely available [`R`](https://cran.r-project.org/web/checks/check_results_drtmle.html) programming language.


-----

## Analyses

We stored the code into different folders according to the respective analyses:

1. Phylogenetic analysis of ANP splice regions (/phylogenetic_analysis/)
2. Estimation of virus production rates and bootstrapping. (/virus_production_rates/)
3. Passage prediction and risk scores (/passage_predictions/)
4. Heatmaps and shiny app (/shiny_app/)
5. Statistical analysis of surveillance data (/surveillance_analysis/)
6. Sensitvity analyses (/sensitivity_analyses/)

 Folder 1 contains the .xml file to repeat the phylogenetic analysis, this file also contains the sequence alignment. Folder 4 contains more detailed instructions on how to set up the shiny app for interactive heatmaps and on how to generate customized heatmaps. Additional data is stored in the [data folder](data/) which consists of one folder with input data and one folder of results of the modelling procedures that are time intensive to produce.

-----

## Figure plotting  

Besides being able to use our modeling framework for analyses of newly collected data on ANP32A splice variants or Influenza sequences found in mammalian/avian species, this repository also contains information on how to reproduce the figures of the paper "Profiling the Host ANP32A Splicing Landscape to Predict Influenza A Virus Polymerase Adaptation". Graphs in Figures 1, 2 and Supplement Figure 1 are produced in Prism based on the data accompying the original paper. The remaining figures include more complicated modeling and can be obtained in the following way:

| Figure | How to |
|---|---|
| 3e | Phylogenetic analysis with [BEAST2](https://www.beast2.org/) with the .xml file [Flu_ANP32A.xml](phylogenetic_analysis/Flu_ANP32A.xml). The maximum clade credibility tree can be obtained by TreeAnnotatar and [FigTree](http://tree.bio.ed.ac.uk/software/figtree/). |
| 4b-d, f, g | function calls explained in [display_bestestimates_CI.R](virus_production_rates/display_bestestimates_CI.R) |
| 5, SuppFig 2 | function calls explained in [passage_predictions_risk_scores.R](passage_predictions/passage_predictions_risk_scores.R) |
| 6a | function calls explained in [passage_predictions_risk_scores.R](passage_predictions/passage_predictions_risk_scores.R) |
| 6b | function calls explained in [heatmap.R](shiny_app/additional_heatmaps/heatmap.R) |
| 7 | function calls explained in [host_adaptation.R](surveillance_analysis/host_adaptation.R)|
| SuppFig 3| function calls explained in [sensitivityanalyses.R](sensitivity_analyses/sensitivityanalyses.R) |
-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/magnuscar/FluAdaptation).
