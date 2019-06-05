# Profiling the Host ANP32A Splicing Landscape to Predict Influenza A Virus Polymerase Adaptation

**Authors:**
Patricia Domingues, Davide Eletto, Carsten Magnus, Hannah L. Turkington, Stefan Schmutz, Osvaldo Zagordi, Matthias Lenk, Michael Huber, Martin Beer, Silke Stertz,  Benjamin G. Hale


-----

## Description

This repository contains code and data to reproduce the analyses of "Profiling the Host ANP32A Splicing Landscape to Predict Influenza A Virus Polymerase Adaptation" by Domingues, Eletto, Magnus et al (2019). In addition, it contains all necessary code to generate the [shiny app](https://magnuscar.shinyapps.io/FluAdaptation/) on your machine. The analysis was implemented in the freely available [`R`](https://cran.r-project.org/web/checks/check_results_drtmle.html) programming language.


-----

## Analyses

We stored the code into different folders according to the respective analyses:

1. Phylogenetic analysis of ANP splice regions (/phylogenetic_analysis/)
2. Estimation of virus production rates and bootstrapping. (/virus_production_rates/)
3. Passage prediction and risk scores (/passage_predictions/)
4. Heatmaps and shiny app (/shiny_app/)
5. Statistical analysis of surveillance data (/surveillance_analysis/)

 Folder 1 contains the .xml file to repeat the phylogenetic analysis, this file also contains the sequence alignment. Folder 4 contains more detailed instructions on how to set up the shiny app for interactive heatmaps and on how to generate customized heatmaps.

-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/magnuscar/FluAdaptation).
