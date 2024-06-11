# sem-mitos
Quantifying morphometrics of dendritic mitochondria from segmented SEM images.

## MCU-enriched dendritic mitochondria regulate plasticity in distinct hippocampal circuits

This repository is for analysis codes related to the preprint below:
[Read the Preprint](https://doi.org/10.1101/2023.11.10.566606)

## Dataset:

This dataset was generated using a custom [Biodock AI model](https://biodock.ai) to selectively segment dendritic mitochondria in scanning electron microscopy (SEM) images.

Dendritic mitochondria were segmented in Serial electron microscopy (SEM) images from basal, proximal and distal CA2 in MCU cKO and control mice.
Large 150 x 150 um ROIs from each layer were tiled into 100 um2 tiles, and every 8th tile was sampled for the analysis.

### Parameters exported from Biodock:
- Area (pixels)
- Ferets Diameter (pixels)
- Length of major axis (pixels)
- Length of minor axis (pixels)
- (X, Y) coordinates

## Methods:

1. Convert all lengths and areas from pixels to microns
2. Remove any objects larger than 2 um2 or smaller than 0.01 um2 (likely not mitochondria)
3. Calculate aspect ratio as the major axis / minor axis for each segmented mitochondria
4. Use KDTree from Scipy to get distance to nearest neighbor between each segmented mitochondria to every other mitochondria within the same image tile
5. Normalize the data to the average of the control (Cre -) mice
6. Calculate the number of mitochondria and total mitochondria area per tile
7. Flag and remove any tiles with less than 2 mitochondria
8. Calculate animal, section, tile and group averages and save the data as CSV files
9. 

## Results:

- Comparing across layers in the control, the cluster with the largest spine heads is over-represented in the distal dendrites compared to proximal dendrites.

- Comparing across genotypes, the cluster with the largest spine heads is reduced in the MCU cKO relative to control, while the cluster with the smallest spine heads is increased. This suggests that loss of MCU causes a selective decrease in very large, round spines and an increase in spines with small, round spine heads.

<!-- Our model achieves the following performance on :
## Results: 

### [Image Classification on ImageNet](https://paperswithcode.com/sota/image-classification-on-imagenet)

| Model name         | Top 1 Accuracy  | Top 5 Accuracy |
| ------------------ |---------------- | -------------- |
| My awesome model   |     85%         |      95%       |

>📋  Include a table of results from your paper, and link back to the leaderboard for clarity and context. If your main result is a figure, include that figure and link to the command or notebook to reproduce it. -->
