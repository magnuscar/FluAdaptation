# Profiling the Host ANP32A Splicing Landscape to Predict Influenza A Virus Polymerase Adaptation

**Authors:**
Patricia Domingues, Davide Eletto, Carsten Magnus, Hannah L. Turkington, Stefan Schmutz, Osvaldo Zagordi, Matthias Lenk, Michael Huber, Martin Beer, Silke Stertz,  Benjamin G. Hale


-----

## Description

We deployed the web application [online](https://magnuscar.shinyapps.io/FluAdaptation/). However, if you want to re-build it locally you can also download the necessary files from this repository. The calculation of a heatmap for given parameter values needs some time and thus, we pre-calculated and stored heatmaps for the parameter values that one can select in the web app. We added the R-script [heatmap.R](additional_heatmaps/heatmap.R) to generate heatmap data for other parameter combinations, that can then either be displayed and stored directly or built into the shiny app.


-----

## How to use the app

With the [app](https://magnuscar.shinyapps.io/FluAdaptation/), one can test whether a newly found species with a certain ratio of ANP32A-splice variants would select for mammalian adapted viral variants (red in the heatmap) or avian adapted viral variants (blue in the heatmap). In addition, one can check whether this result is consistent for other parameters that are included in the model. For more details on the model, please check the original publication. 

-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/magnuscar/FluAdaptation).
