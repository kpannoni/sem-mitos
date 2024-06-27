# sem-mitos
Quantifying morphometrics of dendritic mitochondria from segmented SEM images.

## MCU-enriched dendritic mitochondria regulate plasticity in distinct hippocampal circuits

This repository is for analysis codes related to the preprint below:<br>
[Read the Preprint](https://doi.org/10.1101/2023.11.10.566606)

## Dataset:

This dataset was generated using a custom [Biodock AI model](https://biodock.ai) to selectively segment dendritic mitochondria in scanning electron microscopy (SEM) images.

Dendritic mitochondria were segmented in Serial electron microscopy (SEM) images from basal, proximal and distal CA2 in 3 MCU cKO and 3 control mice.
Large 150 x 150 um ROIs from each layer were tiled into 100 um2 tiles, and every 8th tile was sampled for the analysis.

#### Parameters exported from Biodock:
- Area (pixels)
- Ferets Diameter (pixels)
- Length of major axis (pixels)
- Length of minor axis (pixels)
- (X, Y) coordinates

View the exported data from Biodock [here](Biodock_AI_V6_output_den_mitos_dendrites_cKO_CTL_CA2.csv).

## Analysis Methods:

1. Convert pixels to microns, remove any objects larger than 2 um2 or smaller than 0.01 um2 (likely not mitochondria)
2. Calculate aspect ratio as the major axis / minor axis for each segmented mitochondria
3. Use KDTree from Scipy to get distance to nearest neighbor between each segmented mitochondria to every other mitochondria in the same image tile
4. Flag and remove any tiles with less than 2 mitochondria
5. Normalize the metrics of interest to the average of the control (Cre -) mice
6. Calculate animal, section, tile and group averages and save the data as separate CSV files
7. Create violin plots and run two-way ANOVAs with sidak post hoc on each metric of interest
8. Create correlation plots of mitochondria count versus total mitochondria area per 100 um2 for each layer and each genotype

## Validation of AI performance with manual spot check

A manual spot check of about 1% of the dataset found an accuracy of about 97% across both genotypes.

|   AI Performance   | CTL             | cKO            |   
| ------------------ |:---------------:|:--------------:|
| total objects      | 1081            |      1337      |
| % correct          | 97.1 %          |     97.5 %      |
| missed objects / 100 um2   | 1.1             |      1.4       |


## Results Summary:

### Wild-type mitochondria across CA2 layers

- Mitochondria are bigger in area and Feret's diameter in the distal dendrites (SLM) than the proximal dendrites (SR) of CA2

<!-- Violin plots of area and ferets in the CTL-->
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Area_um_sq_by_layer_violin.png" alt="violin plot of mitochondrial area in CTL CA2" width="40"/> 
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Feret_diam_um_by_layer_violin.png" alt="Violin plot of mitochondrial diameter in CTL CA2", width="40"/>
  
- Mitochondria in the distal dendrites are also more numerous and closer together than the proximal dendrites (SR)

<!-- Plots of count per tile and distance to nearest neighbor-->
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Count_by_layer_violin.png" alt="box plot of mitochondrial area in CTL CA2" width="40"/> 
<img src="https://github.com/kpannoni/sem-mitos/plots_images/NN_Dist_um_by_layer_violin.png" alt="violin plot of mitochondria nearest neighbor distance in CTL CA2" width="40"/> 

- The basal dendrites (SO) have more rounded mitochondria than either proximal or distal dendrites.
  
<!-- violin plot of Aspect Ratio-->
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Aspect_Ratio_by_layer_violin.png" alt="violin plot of mitochondrial area in CTL CA2" width="40"/> 

### Comparing mitochondria in MCU cKO (cKO) to control (CTL)

- Mitochondria are smaller across all dendritic layers of CA2 in the MCU cKO mice compared to CTL

<!-- violin plot of area and feret's diameter in cKO + control-->
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Norm_Area_by_Genotype_violin.png" alt="violin plot of mitochondrial area in CTL vs cKO CA2" width="40"/>
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Norm_Diam_by_Genotype_violin.png" alt="violin plot of mitochondrial diameter in CTL vs cKO CA2" width="40"/>

- MCU deletion did not alter the aspect ratio of mitochondria

<!-- violin plot of aspect ratio in cKO + control-->
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Norm_Aspect_by_Genotype_violin.png" alt="violin plot of mitochondrial aspect ratio in CTL vs cKO CA2" width="40"/>

- Mitochondria are closer together in the MCU cKO mice compared to CTL

<!-- violin plot of NN distance in cKO + control-->
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Norm_Dist_by_Genotype_violin.png" alt="violin plot of nearest neighbor distance in CTL vs cKO CA2" width="40"/>
  
- Mitochondria count per 100um is increased in the cKO compared to CTL across all dendritic layers
  
<!-- Correlation plot comparing cKO and CTL across layers-->
<img src="https://github.com/kpannoni/sem-mitos/plots_images/Mito_count_total_area_CTL_KO_corr_by_layer.png" alt="correlation plot of mitochondria count and total mitochondria area per 100um2 in each layer of CA2" width="40"/>

## Hierarchical statistical bootstrap:

<img src="https://github.com/kpannoni/sem-mitos/bootstrap/bootstrap_schematic.png" alt="schematic of the bootstrap analysis" width="40"/>

<!-- describe the bootstrap and maybe include schematic. Include description of sampling at each level. -->
Due to the hierarchical nature of the dataset, a hierarchical bootstrap was performed on mitochondria area and mitochondria count per 100um2 based on [Saravanan et al, 2020](https://nbdt.scholasticahq.com/article/13927-application-of-the-hierarchical-bootstrap-to-multi-level-data-in-neuroscience). Data was randomly sampled at the level of animal, stub (hippocampal section) and then tile. The median of the resampled data was calculated for each group (layer and genotype). For the bootstrap of individual mitochondria area, the data was additionally resampled at the level of mitochondria. This process was repeated a total of 10,000 times to generate a population of 10,000 medians. We then compared across layers and genotypes by calculating the proportion of bootstrap repetitions where the first group was larger than the second group.

#### Mitochondria Area

<!-- include summary bar plot for mitochondria area in the cKO and CTL-->
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_mito_area_CTL_layers.tif" alt="Proportion of bootstrap wins comparing mitochondria area across layers in CTL CA2" width="40"/>
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_CTL_cKO_mito_area.tif" alt="Proportion of bootstrap wins comparing mitochondria area across genotypes" width="40"/>

#### Mitochondria count per 100um2

<!-- include summary bar plot for mitochondria count in the cKO and CTL-->
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_mito_count_CTL_layers.tif" alt="Proportion of bootstrap wins comparing mitochondria count across layers in CTL CA2" width="40"/>
<img src="https://github.com/kpannoni/sem-mitos/bootstrap/plots_images/Bootstrap_bar_CTL_cKO_mito_count.tif" alt="Proportion of bootstrap wins comparing mitochondria count across genotypes" width="40"/>

