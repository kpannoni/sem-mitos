# sem-mitos
Quantifying morphometrics of dendritic mitochondria from segmented SEM images.

## MCU-enriched dendritic mitochondria regulate plasticity in distinct hippocampal circuits

This repository is for analysis codes related to the BioRxiv preprint by ***Pannoni et al., 2024***.<br>
[Read the preprint here.](https://doi.org/10.1101/2023.11.10.566606)

## Dataset

Dataset was generated using a custom [Biodock AI](https://biodock.ai) model to selectively segment dendritic mitochondria in scanning electron microscopy (SEM) images from mouse hippocampus, which were taken at 2nm resolution with a ThermoFisher Aprea Volumescope.

<!-- add an example image of the ROIs-->

### Training dataset:
The Biodock AI was trained on five 100 &micro;m<sup>2</sup> images from each dendritic layer of CA2 (basal, proximal and distal) from two control (CTL) and two MCU KO mice (cKO).  Training image ROIs were manually selected to represent a variety of features in the dataset. 

After running the AI on a test dataset, images that had >5 errors were edited and added to the training to improve AI performance.

**The final AI was trained on 80 images.**<br>
_Images used for training were not used in the analysis._

### Analysis dataset:
Dendritic mitochondria were segmented in basal, proximal and distal dendrites of CA2 in three hippocampal sections from three MCU cKO and three control mice.
Large 150 x 150 &micro;m ROIs from each layer were tiled into 100 &micro;m<sup>2</sup> tiles, and every 8th tile was sampled for the analysis.

#### Parameters of interest from Biodock:
- Area (pixels)
- Ferets Diameter (pixels)
- Length of major axis (pixels)
- Length of minor axis (pixels)
- (X, Y) coordinates

Get the exported data from Biodock [here](Biodock_AI_V6_output_den_mitos_dendrites_cKO_CTL_CA2.csv).

## Validation of AI performance

A manual spot check of about 1% of the dataset was performed to validate the performance of the trained AI. 

#### Three different types of errors were manually counted for each image in the spot check: 
1) Number of missed dendritic mitochondria  *(false negatives)*
2) Number of incorrectly segmented objects  *(false positives)*
3) Border errors
   
The spot check found an accuracy of about 97% across both genotypes.

|   AI Performance   | CTL             | cKO            |   
| ------------------ |:---------------:|:--------------:|
| total objects      | 1081            |      1337      |
| % correct          | 97.1 %          |     97.5 %      |
| missed objects / 100 um2   | 1.1             |      1.4       |

## Setup to run the analysis

*Ran using Spyder 5.4.3 with Python 3.11 packaged by Anaconda.*

#### Required packages:

+ [easygui 0.98.3](https://pypi.org/project/easygui/)
+ [pingouin 0.5.4](https://pingouin-stats.org/build/html/index.html)

#### Packages included with Anaconda:

    pandas 1.5.3
    matplotlib 3.7.1
    seaborn 0.12.2
    numpy 1.24.3
    scipy 1.10.1

### Main analysis code: 
[EM_image_mito_analysis_Biodock_6_20_2024_cleaned.py](EM_image_mito_analysis_Biodock_6_20_2024_cleaned.py)<br>
**Input:**  `Biodock_AI_V6_output_den_mitos_dendrites_cKO_CTL_CA2.csv`, the object level data exported from Biodock<br>
*User will be prompted to select the CSV file.*
#### Before you run:
- Download the [data file](Biodock_AI_V6_output_den_mitos_dendrites_cKO_CTL_CA2.csv) and codes
- Set the name of your save directory under `dirName` (this directory will be created)
- Make sure the meta data about the analysis in the first code block is correct (AI version, dataset, etc.)
- Add any relevant notes in the `notes` variable, which will be printed in `analysis_summary.txt` along with the other meta data
- Make sure the `mito_functions.py` code is in your working directory along with this analysis code

*If you are running the same analysis, you shouldn't have to change the meta data.*

---

### Supplemental functions: 
[mito_functions.py](mito_functions.py)<br>
Custom functions that will be called by the main analysis code. <br>
*You will need this code in your working directory for the analysis code to run properly.*
#### Custom functions used:
    WT_KO_violin()
    ANOVA_posthoc()
    get_stats_summary()
    prism_format()

---
  
### Bootstrap analysis code: 
[mito_bootstrap_median_for_Biodock_7_1_2024_cleaned.py](mito_bootstrap_median_for_Biodock_7_1_2024_cleaned.py)<br>
**Input:** CSV files `SEM_indiv_mito_data.csv` and `SEM_mito_tile_avgs.csv` saved by the main analysis code
#### Before you run:
- Run the `EM_image_mito_analysis_Biodock_6_20_2024_cleaned.py` code first to process and export the data
- Set the location of the data directory with your two CSV files under the `loc` variable
- Set the number of repetitions with `n_boot` (10000 for the analysis)
- If desired, the number of image tiles and mitochondria to be sampled can be changed in the custom **get_boot_sample()** function by changing the `tile_samp` and `mito_samp` variables


## Analysis Methods Overview

1. Convert pixels to microns, remove any objects larger than 2 &micro;m<sup>2</sup> or smaller than 0.01 &micro;m<sup>2</sup>
2. Calculate aspect ratio as the major axis / minor axis for each segmented mitochondria
3. Use KDTree from Scipy to get distance to nearest neighbor between each segmented mitochondria to every other mitochondria in the same image tile
4. Flag and remove any tiles with less than 2 mitochondria
5. Normalize the metrics of interest to the average of the control (Cre -) mice
6. Calculate animal, section, tile and group averages and save the data as separate CSV files
7. Create violin plots and run two-way ANOVAs with sidak post hoc on each metric of interest
8. Create correlation plots of mitochondria count versus total mitochondria area per 100 &micro;m<sup>2</sup> for each layer and each genotype
9. Compare mitochondrial size and count across layers and genotype with a hierarchical statistical bootstrap

## Results Summary

***See Figure 4 and Supplemental Figure 2 in our [preprint](https://doi.org/10.1101/2023.11.10.566606) for the main results of this analysis.***

### Wild-type mitochondria across CA2 layers

- There is considerable morphological heterogeneity of mitochondria across dendritic layers of CA2.
- Mitochondria are larger, longer and more numerous in the distal dendrites (SLM) than the proximal dendrites (SR) of CA2.

### Comparing mitochondria in MCU cKO (cKO) to control (CTL)

- Mitochondria are smaller and more numerous across all dendritic layers of CA2 in the MCU cKO mice compared to CTL, suggesting that MCU cKO leads to mitochondrial fragmentation.

## Hierarchical statistical bootstrap

<!-- describe the bootstrap and maybe include schematic. Include description of sampling at each level. -->
Due to the hierarchical nature of the dataset, a hierarchical bootstrap was performed on mitochondria area and mitochondria count per 100 &micro;m<sup>2</sup> based on [Saravanan et al, 2020](https://nbdt.scholasticahq.com/article/13927-application-of-the-hierarchical-bootstrap-to-multi-level-data-in-neuroscience). Data was randomly sampled at the level of animal, stub (hippocampal section) and then tile. The median of the resampled data was calculated for each group (layer and genotype). For the bootstrap of individual mitochondria area, the data was additionally resampled at the level of individual mitochondria. This process was repeated a total of 10,000 times to generate a population of 10,000 medians. The bootstrap medians were compared across layers and genotypes by calculating the proportion of bootstrap repetitions where group 1 was larger than group 2.

<!-- #### Mitochondria Area

<!-- include summary bar plot for mitochondria area in the cKO and CTL
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_mito_area_CTL_layers.tif" alt="Proportion of bootstrap wins comparing mitochondria area across layers in CTL CA2" width="40"/>
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_CTL_cKO_mito_area.tif" alt="Proportion of bootstrap wins comparing mitochondria area across genotypes" width="40"/>

#### Mitochondria count per 100 &micro;m<sup>2</sup>

<!-- include summary bar plot for mitochondria count in the cKO and CTL
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_mito_count_CTL_layers.tif" alt="Proportion of bootstrap wins comparing mitochondria count across layers in CTL CA2" width="40"/>
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_CTL_cKO_mito_count.tif" alt="Proportion of bootstrap wins comparing mitochondria count across genotypes" width="40"/> -->

***Detailed results of the bootstrap are shown in Supplemental Figures 3 & 4 in the preprint.***

