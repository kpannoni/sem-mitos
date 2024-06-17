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

## Analysis Methods:

1. Convert pixels to microns, remove any objects larger than 2 um2 or smaller than 0.01 um2 (likely not mitochondria)
2. Calculate aspect ratio as the major axis / minor axis for each segmented mitochondria
3. Use KDTree from Scipy to get distance to nearest neighbor between each segmented mitochondria to every other mitochondria in the same image tile
4. Flag and remove any tiles with less than 2 mitochondria
5. Normalize the metrics of interest to the average of the control (Cre -) mice
6. Calculate animal, section, tile and group averages and save the data as separate CSV files
7. Create violin plots and run two-way ANOVAs with sidak post hoc on each metric of interest
8. Create correlation plots of mitochondria count versus total mitochondria area per 100 um2 for each layer and each genotype

## Results Summary:

<!--

| Genotype           | Dendritic Layer | Median Area    |   Median Feret Diameter | Median Aspect Ratio   | Median Distance to NN | Median Count per 100um2 | N mitochondria |
| ------------------ |:---------------:|:--------------:|:-----------------------:|:---------------------:|:---------------------:|:-----------------------:|:--------------:|
| Cre - (CTL)        | Basal           |      0.15      |        0.59             |                       |                       |                         |                |
| Cre -              | Proximal        |      0.13      |        0.58             |                       |                       |                         |                |
| Cre -              | Distal          |      0.15      |        0.62             |                       |                       |                         |                |
| Cre + (cKO)        | Basal           |      0.13      |        0.54             |                       |                       |                         |                |
| Cre +              | Proximal        |      0.11      |        0.52             |                       |                       |                         |                |
| Cre +              | Distal          |      0.15      |        0.62             |                       |                       |                         |                |
-->

### Wild-type mitochondria across CA2 layers

- Mitochondria are bigger in area and Feret's diameter in the distal dendrites (SLM) than the proximal dendrites (SR) of CA2

<!-- Violin plots of area and ferets in the CTL-->
  
- Mitochondria in the distal dendrites are also more numerous and closer together than the proximal dendrites (SR)

<!-- Plots of count per tile and distance to nearest neighbor-->

- The basal dendrites (SO) have more rounded mitochondria than either proximal or distal dendrites.

<!-- violin plot of Aspect Ratio-->

### Comparing mitochondria in MCU cKO (cKO) to control (CTL)

- Mitochondria are smaller across all dendritic layers of CA2 in the MCU cKO mice compared to CTL

<!-- violin plot of area and feret's diameter in cKO + control-->

- MCU deletion did not alter the aspect ratio of mitochondria

<!-- violin plot of aspect ratio in cKO + control-->

- Mitochondria are closer together in the MCU cKO mice compared to CTL

<!-- violin plot of NN distance in cKO + control-->
  
- Mitochondria count per 100um is increased in the cKO compared to CTL across all dendritic layers
  
<!-- Correlation plot comparing cKO and CTL across layers-->

## Hierarchical statistical bootstrap:

<!-- describe the bootstrap and maybe include schematic. Include description of sampling at each level. -->

#### Mitochondria Area

<!-- include summary bar plot for mitochondria area in the cKO and CTL-->

#### Mitochondria count per 100um2

<!-- include summary bar plot for mitochondria count in the cKO and CTL-->

