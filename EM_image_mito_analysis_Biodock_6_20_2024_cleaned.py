# -*- coding: utf-8 -*-
"""
Created on 9.1.2023

Custom code to work with the object data exported from Biodock, an AI based image analysis platform. Mitochondria were segmented in 6144x4096 pixel SBFEM tiles from the basal, proximal and distal dendrits of CA2 and CA1 of WT and MCU KO mice. The biodock output includes a row for each segmented object in the SEM images along with metrics such as area, X and Y position, length of major and minor axes, perimiter, etc. The file structure in Biodock should be a folder for CTL and cKO, with a subfolder for each layer, then a subfolder for each animal analyzed and each section. This is important, because this code will be getting a lot of the metadata information from the folder names (ie. dendritic layer, animal name, genotype).

The code will convert any lengths or areas from pixels to microns, apply a size filter, then calculate aspect ratio and the distance to nearest neighbor for each segmented mitochondria object. The total number of mitochondria and the total mitochondria area per 100 um2 will be calculated for each image. Any image tiles containing less than 2 mitochondria will be flagged and removed. The data of interest (area, Feret's diameter, perimeter, distance to nearest neighbor and count) will be normalized to the average of the control. 
                                                                                                                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                        The custom function WT_KO_Violin() from "mito_functions.py" is used to plot each metric with violin plots (or box plots) and run a two way ANOVA comparing the effect of genotype and layer. Make sure you have mito_functions.py in your working directory along with this code. 

Input: The CSV file exported from Biodock's analysis. User will be prompted to select this file.
    
Output: CSV files with the individual data, tile data, section data and animal data. There will be an output folder within the main output directory for the violin plots and stats, with a separate folder specifically for figure plots. There will also be a folder with CSV files formated for easy import and plotting in Prism, if desired. Finally, an analysis summary CSV file will be saved in the main output directory with some general information about the analysis.

BEFORE RUNNING THIS CODE: Set the name of the output directory you want, make sure the metadata information in this code is accurate, particularly the Biodock AI version used and the dataset analyzed. If those things change, you will have to update them manually in the code. They will be marked with an #!!! (like this line here).

6.20.24  This version of the code was cleaned up and edited to be published on GitHub along side the preprint of the manuscript.

@author: pannoni
"""

# import needed packages
import pandas as pd
import easygui
import os
from datetime import date
import numpy as np
from mito_functions import WT_KO_violin, get_stats_summary, prism_format # custom function
import seaborn as sns
import matplotlib.pyplot as plt
import inspect
from scipy.spatial import KDTree

# Expand the printed text of pd dataframes
pd.set_option('display.max_column', None)
pd.set_option('display.expand_frame_repr', False)

# import the excel file as a dataframe
# If easygui fails or you don't want to install it, change the data_loc variable on the line below to the full file pathway to your Biodock output CSV file. The rest of the code should work fine.
data_loc = easygui.fileopenbox(msg='Please select the CSV file with the object data.', title='Select CSV File', default='*')
mito_data = pd.read_csv(data_loc)

#Create a new folder in the working directory to save all the output files for this code. 
#!!! Set the name of the output folder where the data and plots will be saved here.
dirName = "SEM_mitos_CA2_Biodock_V6_6_20_24"
 
try:
    # Create target Directory
    os.mkdir(dirName)
    print("\nDirectory" , dirName ,  "Created.\n") 
except FileExistsError:
    print("\nDirectory" , dirName ,  "already exists, no need to create a new directory.")
    
# Make sure the meta data about the analysis is correct below
    
#!!! Enter the name an version of the Biodock AI used to segment the images here. This is important information to have if we end up training multiple AI versions.
AI_name = "Mitochondria in MCU KO and WT mice (V6)"

#!!! Dataset that the images are from:
dataset = "CTL and KO CA2"

# Set the tile size for the metadata
tile_size = "6144 x 4969 pixels"

# Tile sampling
sampling = "Every 8th tile"

#!!! Space to put notes about the dataset or any changes made since a previous run.
notes = "Analysis of CTL and MCU KO data in CA2 using version 6 of the Biodock AI. Denderitic mitochondria and dendrites were segmented with a confidence threshold of 0.4. This dataset was used for Figure 4 and supplemental figures 2-4."
    
# Parse out the image metadata from the CSV file by breaking the "image origin" column into its separate folders

if len(mito_data["Image origin"][1].split("/")) == 5: # check that there are the right number of levels in the data
    # order is genotype, dendritic layer, animal ID, "stub" and then image file
    mito_data["Genotype"] = [x.split("/")[0] for x in mito_data["Image origin"]] # Genotype or group should be first in the hierarchy
    mito_data["Layer"] = [x.split("/")[1] for x in mito_data["Image origin"]] # dendritic layer should be second
    mito_data["Animal"] = [x.split("/")[2] for x in mito_data["Image origin"]] # Animal ID third
    mito_data["Stub"] = [x.split("/")[3] for x in mito_data["Image origin"]] # Then stub name
    mito_data["Tile"] = [x.split("/")[4] for x in mito_data["Image origin"]] # image file name
else:
    print("Error: Can't get metadata information from the Biodock output file. Please make sure this information is included in the 'Image origin' column of your file. The Image origin should include the genotype, layer, animal, stub or section, then the image tile associated each object. It should be in that order, separated by /. If the image file structure in biodock is set up as described, it should be exported that way. If not, you will have to fix it and run the code again.")

# Make the layer column into category so we can sort them in the order we want (Basal, Prox, Distal)
sorter = ["Basal", "Proximal", "Distal"] # sorting list
# Set the Layer column to be a "category" type that we can sort by
mito_data.Layer = pd.Categorical(mito_data.Layer, categories=sorter, ordered=True)
mito_data = mito_data.sort_values(['Genotype', 'Layer', 'Animal', "Stub"], ascending = [False, True, True, True]).reset_index(drop=True) # sort

# Change "object ID" to string so it won't be treated as numerical
mito_data["Object ID"] = mito_data["Object ID"].astype("int").astype("str")
    
print("\nCSV file uploaded.")

# Separate the dendrite objects from the mitochondria objects. We aren't going to do much with the dendrite data now, but likely will in the future.
dendrite_data = mito_data[mito_data["Class"] == "Dendrite Object"].reset_index(drop=True)
mito_data = mito_data[mito_data["Class"] != "Dendrite Object"].reset_index(drop=True)


#%%% Process the data by converting pixels to microns and setting a size filter to remove non-mitochondria objects. Then calculate aspect ratio as the ratio of the major and minor axes and give each section a unique ID.

print("\nConverting area and length from pixels to microns...")

# First we need to convert Area and Length into um or um squared
# image scale is 2 nm per pixel
scale_nm_per_px = 2
# Convert scale to um instead of nm
scale_um_per_px = scale_nm_per_px / 1000 # use for as conversion factor for length
# get the scale in microns squared
scale_um_sq = scale_um_per_px **2 # convert nm to um and square, use as conversion factor for area

# convert the mitochondrial area into microns sq
mito_data["Area_um_sq"] = mito_data["Area"] * scale_um_sq
# convert the dendrite area into microns sq
dendrite_data["Den_Area_um_sq"] = dendrite_data["Area"] * scale_um_sq

# convert all lengths from pixels to microns
mito_data["Major_Length_um"] = mito_data["Length of Major Axis"] * scale_um_per_px
mito_data["Minor_Length_um"] = mito_data["Length of Minor Axis"] * scale_um_per_px
mito_data["Perimeter_um"] = mito_data["Perimeter"] * scale_um_per_px
mito_data["Feret_diam_um"] = mito_data["Feret Diameter Maximum"] * scale_um_per_px

dendrite_data["Den_Width_um"] = dendrite_data["Length of Minor Axis"] * scale_um_per_px
dendrite_data["Den_Length_um"] = dendrite_data["Length of Major Axis"] * scale_um_per_px

print("\nRemoving mitochondria larger than 2 microns squared or smaller than 0.01 microns squared...")

total_mitos = len(mito_data)

# Let's remove any mitos below a size of 0.1 um2, or larger than 2um2
mito_data = mito_data.query('0.01 < Area_um_sq < 2.1').reset_index(drop=True)

mitos_removed = total_mitos - len(mito_data)

if mitos_removed > 0:
    print(" ", mitos_removed, " mitochondria removed from size exclusion.")
else:
    print("No mitochondria removed from size exclusion.")
    
# Let's remove any dendrite segments below a size of 0.02 um2 (~5000 px), these are likely spine heads and not actually dendrites
print("\nRemoving dendrites smaller than 0.4 microns squared since these are likely spine heads...")
# get the current number of dendrites
total_dendrites = len(dendrite_data)
dendrite_data = dendrite_data.query('Den_Area_um_sq > 0.4').reset_index(drop=True)
# Get the number of dendrites that were removed, if any
den_removed = total_dendrites - len(dendrite_data)

if den_removed > 0:
    print(" ", den_removed, " dendrites removed as likely spine heads.")
else:
    print("No dendrite segments were removed.")

print("\nCalculating the mitochondria aspect ratio...")
# Calculate the aspect ratio using the major and minor axes, although we will probably just use the length of the major axes instead
mito_data["Aspect_Ratio"] = round(mito_data["Length of Major Axis"] / mito_data["Length of Minor Axis"],2)

# We'll have to make a unique name by combining animal ID and stub since some animals have stubs with the same names
mito_data["Section_ID"] = mito_data["Animal"] + "_" + mito_data["Stub"]

#%% Now let's do a nearest neighbor analysis to get the distance from each mitochondria to the closest mitochondria within the same image tile.

print("\nCalculating the nearest neighbors...")

# Find the nearest neighbor for each mitochondria with a KDtree
# Note this will have to be done within the same tile image

# Add an empty column to mito_data to put the distances and the nearest neighbor IDs
tile_list = []

# First combine the X and Y coordinates into pairs for each mito to get the locaion
mito_data["Location"] = mito_data[["X Position", "Y Position"]].values.tolist()

for tile in list(set(mito_data["Tile"])): # loop across each image
    # Get the data for the current image
    tile_data = mito_data.where(mito_data["Tile"] == tile).dropna(how="all")
    # Get the mitochondria locations and IDs for the current image
    mito_coor = list(tile_data["Location"])
    mito_IDs = list(tile_data["Object ID"].astype("int"))
    
    if len(tile_data) > 1: # if there's more than one mitochondria in the tile
        
        # Create a KD tree to locate nearest neighbors
        tree = KDTree(mito_coor)
        
        mito_dist = [] # empty list to put the mito distances for the current image
        neighbor_IDs = [] # empty list to put the object ID of the nearest neighbor
        
        # For each mito in the image, query the tree to get the distance to the nearest neighbor
        for mito in mito_coor:
            dist, obj = tree.query(mito, k=2) # get the two nearest neighbors
            # Note the query counts the same mitochondria itself as the first neighbor (dist = 0), so we actually want the second neighbor
            
            # Save the distance to the nearest mito in a column in mito_data (drops the distance to self)
            mito_dist.append(dist[1])
            
            # Save the object ID of the nearest mito in another column in mito_data
            NN_ID = mito_IDs[obj[1]]
            neighbor_IDs.append(NN_ID)
            
        # add the distances and IDs to the tile_data dataframe
        tile_data["NN_Dist"] = mito_dist
        tile_data["NN_ID"] = neighbor_IDs
        
        # Add the nearest neighbor data to the original mito_data dataframe in the correct location
        tile_list.append(tile_data)
    else:
        print("\nTile '" + str(tile) + "' has only one mitochondria. Nearest neighbor analysis cannot be done on this tile.\nThis tile should likely be removed from the analysis.")
    
# combine the tile data back into the original dataframe
mito_data = pd.concat(tile_list).sort_index()

# Change the ID column to string so it's not a numerical data type
mito_data["NN_ID"] = mito_data["NN_ID"].astype(str)

# Convert the nearest neighbor distance into microns like all the other distances
mito_data["NN_Dist_um"] = mito_data["NN_Dist"] * scale_um_per_px

#%% Get mitochondrial counts and average metrics per tile and remove any tiles with less than 2 mitochondria.

print("\nGetting the number and total area of mitochondria and dendrites in each tile...")

# group the mito data by tile to get averages for the metrics
mito_avgs = mito_data.groupby("Tile").mean(numeric_only=True) # get the tile averages
meta_cols = mito_data.groupby("Tile").agg({"Animal": "first", "Genotype": "first", "Layer": "first", "Stub": "first", "Section_ID": "first", "Tile": "count"}) # keep the metadata and count the number of objects (mitochondria) for each tile
mito_avgs = meta_cols.join(mito_avgs) # join the averaged data with the metadata columns + count

# Rename the tile column to "count" 
mito_avgs = mito_avgs.rename(columns = {"Tile": "Count"})

# get total mito area per tile
mito_avgs["Total_mito_Area_um_sq"] = mito_data[["Tile", "Area_um_sq"]].groupby("Tile").sum()

# Round the data in mito_avgs to 3 decimil points
mito_avgs = round(mito_avgs, 2)
mito_avgs.Layer = pd.Categorical(mito_avgs.Layer, categories=sorter, ordered=True) # make sure Layer is categorical

# Let's flag any tiles where the count is under 2. This likely means the tile is not useable or something is wrong with the tile image.
flagged = mito_avgs["Count"].where(mito_avgs["Count"] <= 2).dropna(how="all") # get any tiles where count is 5 or less
num_flag = len(flagged)

# print a warning if there are any flagged tiles
if num_flag >= 1:
    print("\n" + str(num_flag) + " tiles were flagged because they have fewer than 2 mitochondria. Recommended to visually check these tiles to see if they should be analyzed. These tiles will be removed from the data now.\n" + str(list(flagged.index)))
    
# Remove any flagged tiles from the main dataframe and the tile average dataframe
# Comment out this line if you don't want to remove the flagged tiles from the data
for tile in flagged.index:
    mito_data = mito_data[~mito_data["Tile"].str.contains(tile)]
    mito_avgs = mito_avgs[~mito_avgs.index.str.contains(tile)]
    
#%%% Normalize the data of interest to the CTL average in a separate dataframe. This data will be used for figure plots in Figure 4 of the manuscript.

print("\nNormalizing the data to the control average (cre -)...")

# Note that as the code is written, the normalization includes layer SO

# Empty dataframe for the normalized data
norm_to_ctrl = pd.DataFrame(columns = ["Object ID", "Animal", "Stub", "Tile", "Layer", "Genotype", "Norm_Area", "Norm_Diam", "Norm_Aspect", "Norm_Dist"])

# First we need to get the WT mean of all the data for each metric. This is what we will used as the normalization factor. 

# Get the CTL data only
WT_data = mito_data[mito_data["Genotype"] == "MCC Cre -"]
# Get the means
WT_avg = WT_data[["Area_um_sq", "Feret_diam_um", "Aspect_Ratio", "NN_Dist_um"]].mean()

# Normalize individual mito area by the average of the WT control
mito_data["Norm_Area"] = round(mito_data["Area_um_sq"] / WT_avg["Area_um_sq"], 2)
# Normalize mito diameter by the WT control
mito_data["Norm_Diam"] = round(mito_data["Feret_diam_um"] / WT_avg["Feret_diam_um"], 2)
# Normalize the mito aspect ratio by the WT control
mito_data["Norm_Aspect"] = round(mito_data["Aspect_Ratio"] / WT_avg["Aspect_Ratio"], 2)
# Normalize the NN distance by the WT control
mito_data["Norm_Dist"] = round(mito_data["NN_Dist_um"] / WT_avg["NN_Dist_um"], 2)

print(" done.")
    
# Subset the mito_avg dataframe to get just the columns of interest
# Later in the code, we will add the normalized count and total mito area to this dataframe and save it as a CSV file
mito_avgs_sub = mito_avgs[['Animal', 'Genotype', 'Layer', 'Stub', 'Section_ID', 'Area_um_sq', 'Major_Length_um', 'Minor_Length_um', 'Feret_diam_um', 'Aspect_Ratio', 'Perimeter_um', 'Count', 'Total_mito_Area_um_sq']]
# Save as CSV file
mito_avgs_sub.to_csv(os.path.join(dirName, str("SEM_mito_tile_avgs.csv")), index=True)

# Also save the individual mito data
# This includes some columns we don't really use, but may want to look at later.
mito_data_sub = mito_data[['Object ID', 'Tile', 'Animal', 'Genotype', 'Layer', 'Stub', 'Area', 'Area_um_sq', 'Length of Major Axis', 'Length of Minor Axis', 'Major_Length_um', 'Minor_Length_um', 'Feret_diam_um', 'Aspect_Ratio', 'Perimeter_um', 'Dendrite Object parent ID', 'NN_Dist_um', 'NN_ID', 'Eccentricity', 'Average Intensity (channel 1)', 'Norm_Area', 'Norm_Diam', 'Norm_Aspect', 'Norm_Dist']]
# Save as CSV file
mito_data_sub.to_csv(os.path.join(dirName, str("SEM_indiv_mito_data.csv")), index=False)


#%%% Let's also take a look at the animal averages and the overall averages by genotype and layer

print("\nGetting animal, section and group averages for metrics of interest...")

animal_avgs = round(mito_data_sub.groupby(["Genotype", "Animal", "Layer"], observed=True).median(numeric_only=True),3) # get the animal averages by layer

# Get the average number of mitochondria per 100um2 (tile) for each animal
animal_avgs["Count_100um2"] = round(mito_avgs[["Genotype", "Animal", "Layer", "Count"]].groupby(["Genotype", "Animal", "Layer"], observed=True).median(),1)

# Also get the total count from the tile average dataframe "mito_avgs"
animal_avgs["Total_Count"] = mito_avgs[["Genotype", "Animal", "Layer", "Count"]].groupby(["Genotype", "Animal", "Layer"], observed=True).sum()

# We need the total area analyzed, which we will use to normalize total mito area

# Number of tiles per animal for each layer
animal_avgs["Tiles_per_Animal"] = mito_avgs.groupby(["Genotype", "Animal", "Layer"], observed = True).agg({"Animal": "count"}).rename(columns={"Animal":"Tiles_per_Animal"})

# Calculate the total area analyzed for each stub (100 um2 per tile)
animal_avgs["Total_Area_um_sq"] = animal_avgs["Tiles_per_Animal"] * 100

# Get the total mito area by layer
animal_avgs["Total_mito_mass"] = mito_data[["Genotype", "Animal", "Layer", "Area_um_sq"]].groupby(["Genotype", "Animal", "Layer"], observed=True).sum()

# Normalize the total mito area by the total area analyzed
animal_avgs["Total_mito_mass_per_area"] = animal_avgs["Total_mito_mass"] / animal_avgs["Total_Area_um_sq"]

# reset index so that genotype, animal and layer are columns instead of just index
animal_avgs = animal_avgs.reset_index().sort_values(['Genotype', 'Animal', 'Layer'], ascending = [False, True, True]).reset_index(drop=True)

# Subset the animal df to just save the columns of interst
animal_avgs_sub = animal_avgs[["Genotype", "Animal", "Layer", "Area_um_sq", "Major_Length_um", "Minor_Length_um", "Feret_diam_um", "Aspect_Ratio", "Perimeter_um", "NN_Dist_um", "Norm_Area", "Norm_Diam", "Norm_Aspect", "Norm_Dist", "Count_100um2", "Total_Count", "Tiles_per_Animal", "Total_Area_um_sq", "Total_mito_mass", "Total_mito_mass_per_area"]]

# Save the animal averages to a CSV file
animal_avgs_sub.to_csv(os.path.join(dirName, str("SEM_mito_animal_avgs.csv")), index=False)

# Now let's get the averages by stub 
# Note that although we looked at the data by stub, this isn't included in the analyses in the paper.
section_avgs = round(mito_data[["Genotype", "Animal", "Section_ID", "Layer", "Area_um_sq", "Major_Length_um", "Minor_Length_um", "Feret_diam_um", "Aspect_Ratio", "Perimeter_um", "Eccentricity", "NN_Dist_um", "Norm_Area", "Norm_Diam", "Norm_Aspect", "Norm_Dist"]].groupby(["Genotype", "Animal", "Section_ID", "Layer"]).median(),3).dropna(how="all")

# Add count from the tile average dataframe "mito_avgs"
section_avgs["Total_Count"] = mito_avgs[["Genotype", "Animal", "Section_ID", "Layer", "Count"]].groupby(["Genotype", "Animal", "Section_ID", "Layer"]).sum().dropna(how="all")

section_avgs["Total_mito_mass"] = mito_data[["Genotype", "Animal", "Section_ID", "Layer", "Area_um_sq"]].groupby(["Genotype", "Animal", "Section_ID", "Layer"]).sum().dropna(how="all")

section_avgs["Tiles_per_Stub"] = mito_avgs.groupby(["Genotype", "Animal", "Section_ID", "Layer"]).agg({"Section_ID": "count"})

# Calculate the total area analyzed for each stub
section_avgs["Total_Area_um_sq"] = section_avgs["Tiles_per_Stub"] * 100

# Use the total area to normalize both count and total mito mass
section_avgs["Total_Count_per_area"] = section_avgs["Total_Count"] / section_avgs["Total_Area_um_sq"]
section_avgs["Total_mito_mass_per_area"] = section_avgs["Total_mito_mass"] / section_avgs["Total_Area_um_sq"]

# Save the section data to a CSV file
section_avgs.to_csv(os.path.join(dirName, str("SEM_mito_section_avgs.csv")))

# Now let's compare the different groups of the data (WT vs KO) across layers. Print some summary statistics for the metrics of interest.

# Get the number of groups in the data
group_list = list(set(mito_data["Genotype"]))
num_groups = len(group_list)

# Go back and add count before here
if num_groups >1:
    
    # Get the group medians by layer
    genotype_med = round(mito_data.groupby(["Genotype", "Layer"]).median(numeric_only=True),2)
    
    # Also look at the group standard deviations
    genotype_std = round(mito_data.groupby(["Genotype", "Layer"]).std(numeric_only=True),2)
    
    # Subset to get only the metrics of interest
    genotype_sub = genotype_med.loc[:, ("Area_um_sq", "Aspect_Ratio", "Perimeter_um", "Feret_diam_um", "NN_Dist_um")]
    
    # Subset the standard dev table as well with the same columns
    genotype_sub_std = genotype_std.loc[:, ("Area_um_sq", "Aspect_Ratio", "Perimeter_um", "Feret_diam_um", "NN_Dist_um")]
    
    # Get count from the mito_avgs dataframe
    genotype_sub["Count_100um2"] = round(mito_avgs[["Genotype", "Layer", "Count"]].groupby(["Genotype", "Layer"]).median(),1)
    genotype_sub["N_mitos"] = round(mito_data[["Genotype", "Layer","Object ID"]].groupby(["Genotype", "Layer"]).count(),1)
    
    # print to the console
    print("\n Group Averages:\n", genotype_sub)
    # Save as a CSV file
    genotype_sub.to_csv(os.path.join(dirName, str("SEM_mito_group_summary_medians.csv")))
    genotype_sub_std.to_csv(os.path.join(dirName, str("SEM_mito_group_summary_std.csv")))
    
    # Get 75% percentile stats for area and diameter to use to pick representative images
    desc_stats = mito_data[["Area", "Feret Diameter Maximum", "Genotype", "Layer"]].groupby(["Genotype", "Layer"]).describe()
    
#%%% Let's take a look at the distributions for mitochondrial area in the CTL and the KO. Plot histograms and look at the summary stats.

# Separate CTL and cKO data
WT_data = mito_data[mito_data["Genotype"] == "MCC Cre -"]
KO_data = mito_data[mito_data["Genotype"] == "MCC Cre +"]

# plot histograms with the get_stats_summary() custom function

# plot the CTL histogram by layer
stats_CTL, norm_CTL = get_stats_summary(data = WT_data, data_col = "Area_um_sq", x_label = "Area (um2)", save_dir = dirName, group_name = "CTL", binsize = 0.01, hist_comp = "proportion", figsize = (8,5), xlim = [0,0.4], ylim = [0,0.1])

# similar histogram for the cKO
stats_KO, norm_KO = get_stats_summary(data = KO_data, data_col = "Area_um_sq", x_label = "Area (um2)", save_dir = dirName, group_name = "cKO", binsize = 0.01, hist_comp = "proportion", figsize = (8,5), xlim = [0,0.4], ylim = [0,0.1])
    

#%%% Generate the CTL and cKO violin plots for Figure 4 and supplemental Figure 2. These will be saved in a separate folder called "Figure plots" along with the results of an ANOVA test for each metric of interest.

# Create a folder in the output JUST for the KO paper figures
fig_path = os.path.join(dirName, "Figure Plots")
# Create a subfolder for just the prism files
try:
    os.mkdir(fig_path) # make directory 
except FileExistsError:
    pass # skip if the folder already exists

# If you need the hex codes for CTL and cKO:
# WT_colors = ["#2c6fbb", "#a2cffe", "#004577"]
# KO_colors = ["#ff9a8a", "#ffe8e1", "#cf524e"]

# Set the color scheme
WT_colors = sns.cubehelix_palette(n_colors = 3, start=2.6, light=0.8, dark=0.25, rot=0, hue=1.5)

KO_colors = sns.cubehelix_palette(n_colors = 3, start=1.05, light=0.85, dark=0.35, rot=0, hue=2.5)

print("\nAnalyzing denditic mitochondria in the CTL (not normalized):")

# violin plot for individual mito area in the CTL
# If you want to run an ANOVA with any of these plots, add "stats = 1" to the line below and an ANOVA and post hoc tests will be run between your defined groups and saved as CSV files and significance added to the plots
WT_KO_violin(data = WT_data, data_col = "Area_um_sq", group=0, title ="CTL CA2", save_path = fig_path, units = "Mito Area (\u00B5m\u00B2)", y_range = [0,0.5], colors = WT_colors)

# violin plot for mito Feret's diameter in CTL:
WT_KO_violin(data = WT_data, data_col = "Feret_diam_um", group=0, title ="CTL CA2", save_path = fig_path, units = "Feret's Diameter (\u00B5m)", y_range = [0,1.5], colors = WT_colors)

# violin plot for mito aspect ratio in CTL:
WT_KO_violin(data = WT_data, data_col = "Aspect_Ratio", group=0, title ="CTL CA2", save_path = fig_path, units = "Aspect Ratio", y_range = [0,4.5], colors = WT_colors)

# violin plot for nearest neighbor distance in CTL:
WT_KO_violin(data = mito_data, data_col = "NN_Dist_um", group=0, title ="CTL CA2", save_path = fig_path, units = "NN Distance (\u00B5m)", y_range = [0,4], colors = WT_colors)

# violin plot for mito count in CTL:
WT_KO_violin(data = mito_avgs, data_col = "Count", group=0, title ="CTL CA2", plot_type="box", save_path = fig_path, units = "Count per 100 \u00B5m\u00B2", y_range = [0,60], colors = WT_colors)

print("\nAnalyzing denditic mitochondria in the CTL and MCU KO (normalized to CTL):")

# For plotting, combine WT and KO colors into one list alternating WT and KO
violin_colors = [item for pair in zip(WT_colors, KO_colors + [0]) for item in pair]

# violin plot for individual mito area in cKO and CTL normalized to CTL
WT_KO_violin(data = mito_data, data_col = "Norm_Area", title = "CA2 CTL vs cKO", save_path = fig_path, units = "Mito Area \nNorm. to CTL mean", y_range = [0,2], colors = violin_colors)

# violin plot for mito Feret's diameter normalized to CTL:
WT_KO_violin(data = mito_data, data_col = "Norm_Diam", title = "CA2 CTL vs cKO", save_path = fig_path, units = "Feret's Diameter \nNorm. to CTL mean", y_range = [0,1.8], colors = violin_colors)

# violin plot for mito aspect ratio normalized to CTL:
WT_KO_violin(data = mito_data, data_col = "Norm_Aspect", title = "CA2 CTL vs cKO", save_path = fig_path, units = "Aspect Ratio \nNorm. to CTL mean", y_range = [0,2], colors = violin_colors)

# violin plot for nearest neighbor distance normalized to CTL:
WT_KO_violin(data = mito_data, data_col = "Norm_Dist", title = "CA2 CTL vs cKO", save_path = fig_path, units = "NN Distance \nNorm. to CTL mean", y_range = [0,3], colors = violin_colors)


#%%% Now lets make a few more plots to look at mitochondrial content across layers in the cKO and CTL. We will plot mitochondrial count and mitochondrial total area.

# Let's create a folder in the output to save the plots and stats
plot_path = os.path.join(dirName, "Plots and Stats")
# Create a subfolder for just the prism files
try:
    os.mkdir(plot_path) # make directory 
except FileExistsError:
    pass # skip if the folder already exists

print("\nMitochondrial content (per tile):")

# Get metrics for mitochondrial content like the number of mitochondria and total mito area. These will be at the level of image tile and will be plotted as box plots instead of violin plots.

# Box plot of total mitochondria area per tile (100 um2)
WT_KO_violin(data = mito_avgs, data_col = "Total_mito_Area_um_sq", plot_type ="box", title = "CA2 CTL vs cKO", save_path = plot_path, units = "Total mito Area \nper 100 \u00B5m\u00B2", y_range = [0,8])

# Box plot of Number of mitos per tile (100 um2)
WT_KO_violin(data = mito_avgs, data_col = "Count", plot_type ="box", title = "CA2 CTL vs cKO", save_path = plot_path, units = "Mito Count / 100 \u00B5m\u00B2", y_range = [0,60])

# Generate a scatter plot with total mito area on the X axis and mitochondria count on the Y. This is part of Figure 4 of the manuscript.

# For this plot, we are removing the Basal layer and just looking proximal and distal dendrites
plot_data = mito_avgs.where(mito_avgs["Layer"] != "Basal").dropna(how="all")
plot_data["Layer"] = plot_data["Layer"].cat.remove_unused_categories()

# format the overall figure
sns.set_style("ticks")
sns.set_context("talk")

plot = sns.relplot(x= "Total_mito_Area_um_sq", y= "Count", data=mito_avgs, col="Genotype", hue="Layer", palette= WT_colors, alpha=0.8, s=50, height=4, aspect=1)

plot.set_axis_labels("Total Mito Area per 100\u00B5m\u00B2", "Count per 100\u00B5m\u00B2", fontweight="bold")
axes = plot.axes
axes[0,0].set(xlim=[0,8], yticks=np.arange(0, 50.1, 10), ylim=[0,50])
axes[0,0].set_ylabel("Count per 100 \u00B5m\u00B2", fontsize=18, fontweight="bold", labelpad=8) # count
axes[0,0].set_xlabel('Total mito area per 100 \u00B5m\u00B2', fontsize=18, fontweight="bold")
# axes[0,0].xticks(fontsize=16)  

plt.savefig(os.path.join(fig_path, "Mito_count_total_area_corr_plot.png"), format='png', dpi=300, bbox_inches='tight')

# We'll also look correlation plot for area and aspect ratio. This figure is not in the manuscript.

# First, set the layer order so they'll be in order on the plot
mito_data_ordered = mito_data.copy()
mito_data_ordered["Layer"]= mito_data_ordered["Layer"].cat.reorder_categories(["Distal", "Proximal", "Basal"])

# Separate the data by Layer for plotting (WT data only)
SO_data = WT_data.where(WT_data["Layer"]== "Basal").dropna(how="all").reset_index(drop=True)
SR_data = WT_data.where(WT_data["Layer"]== "Proximal").dropna(how="all").reset_index(drop=True)
SLM_data = WT_data.where(WT_data["Layer"]== "Distal").dropna(how="all").reset_index(drop=True)

fig, axes = plt.subplots(figsize=(5,4))
sns.despine(fig=fig, top=True, right=True)

# Plot SLM data
sns.scatterplot(x= "Aspect_Ratio", y= "Area_um_sq", data=SLM_data, ax=axes, color= WT_colors[2], alpha= 0.9, s=20, label="SLM")
# Plot SR data
sns.scatterplot(x= "Aspect_Ratio", y= "Area_um_sq", data=SR_data, ax=axes, color= WT_colors[1], alpha= 0.9, s=20, label="SR")
# Plot SO data
sns.scatterplot(x= "Aspect_Ratio", y= "Area_um_sq", data=SO_data, ax=axes, color= WT_colors[0], alpha= 0.9, s=20, label="SO")

plt.legend(fontsize='12', title_fontsize='14', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, frameon=False)

axes.set(ylabel ="Mito Area (\u00B5m\u00B2)", xlabel = "Aspect Ratio", xlim=[0,9], ylim=[0,1.25])

# Save the area / aspect ratio figure
plt.savefig(os.path.join(fig_path, "Mito_area_aspect_ratio_corr_plot.png"), format='png', dpi=300, bbox_inches='tight')

#%%% Normalize mito count and total mito area to the CTL mean and create a scatter plot of the normalized data for each layer with total mito area on the X and Count on the Y axis.

# Now we need to normalize count and total area to the control, like we normalized the other metrics.

# Get the metadata for each tile
norm_to_ctl_tile = mito_avgs.loc[:,["Animal", "Stub", "Layer", "Genotype"]]

# Get the means for count and total area in CTL
WT_tile_avgs = mito_avgs.where(mito_avgs["Genotype"] == "MCC Cre -").dropna(how="all").reset_index(drop = True) # get only WT data
WT_tile_mean = WT_tile_avgs[["Layer","Total_mito_Area_um_sq", "Count"]].mean(numeric_only=True) # get the WT mean

# Normalize total mito area by the average of the CTL
norm_to_ctl_tile["Norm Total Area"] = round(mito_avgs["Total_mito_Area_um_sq"] / WT_tile_mean["Total_mito_Area_um_sq"], 2)
# Normalize mito count by the WT control
norm_to_ctl_tile["Norm Count"] = round(mito_avgs["Count"] / WT_tile_mean["Count"], 2)

# Empty dataframe for plotting the center of mass and error bars of the normalized data
norm_plot = pd.DataFrame()

norm_plot["X"] = norm_to_ctl_tile["Norm Total Area"]
norm_plot["Y"] = norm_to_ctl_tile["Norm Count"]
#  Combine genotype and layer to create a "Group" column
norm_plot["Group"] = norm_to_ctl_tile["Genotype"] + "_" + norm_to_ctl_tile["Layer"].astype("object")
# Get the medians and st. deviation
norm_means = norm_plot[["X", "Y", "Group"]].groupby("Group").agg([np.median, np.std])

norm_to_ctl_grouped = round(norm_to_ctl_tile.groupby(["Genotype", "Layer"], observed=True).mean(numeric_only=True), 2)

# Create a dataframe for the normalized data
# grab the normalized data for the other metrics from the genotype_avg dataframe
norm_indiv_data = genotype_med.loc[:, ("Norm_Area", "Norm_Aspect", "Norm_Diam", "Norm_Dist")]

# combine the normalized data into a dataframe
norm_to_ctl_grouped = pd.concat([norm_indiv_data, norm_to_ctl_grouped], axis=1)

# Save as a CSV, in case we want to look at it later
norm_to_ctl_grouped.to_csv(os.path.join(dirName, str("SEM_mito_group_summary_norm_to_ctrl.csv")), index=True)

# split the data into the 4 groups so we can plot each group separately
plot_data_KO_SO = norm_to_ctl_tile.where((norm_to_ctl_tile["Genotype"]== "MCC Cre +") & (norm_to_ctl_tile["Layer"] == "Basal")).dropna().reset_index(drop=True)
plot_data_KO_SR = norm_to_ctl_tile.where((norm_to_ctl_tile["Genotype"]== "MCC Cre +") & (norm_to_ctl_tile["Layer"] == "Proximal")).dropna().reset_index(drop=True)
plot_data_KO_SLM = norm_to_ctl_tile.where((norm_to_ctl_tile["Genotype"]== "MCC Cre +") & (norm_to_ctl_tile["Layer"] == "Distal")).dropna().reset_index(drop=True)

plot_data_WT_SO = norm_to_ctl_tile.where((norm_to_ctl_tile["Genotype"]== "MCC Cre -") & (norm_to_ctl_tile["Layer"] == "Basal")).dropna().reset_index(drop=True)
plot_data_WT_SR = norm_to_ctl_tile.where((norm_to_ctl_tile["Genotype"]== "MCC Cre -") & (norm_to_ctl_tile["Layer"] == "Proximal")).dropna().reset_index(drop=True)
plot_data_WT_SLM = norm_to_ctl_tile.where((norm_to_ctl_tile["Genotype"]== "MCC Cre -") & (norm_to_ctl_tile["Layer"] == "Distal")).dropna().reset_index(drop=True)

# Set color palettes
WT_linecolor = sns.dark_palette(n_colors = 5, color='#0d516f')[2]
KO_linecolor = sns.dark_palette(n_colors = 5, color='#c83005')[2]
KO_edge_color = sns.light_palette(n_colors = 8, color='#c83005')[5]

# Set up the figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
fig.tight_layout(pad = 1.2)
# Remove top and right axes lines
sns.despine(fig=fig, top=True, right=True)

# Plot the groups on separate plots. This will make things easier
# Line is the linear regression line with 95% confidence interval
a = sns.regplot(x="Norm Total Area",y="Norm Count", color=WT_colors[2], data=plot_data_WT_SO, ax=ax1, label="CTL SO", fit_reg=False, scatter_kws = {"s": 55, "alpha":0.8, 'edgecolor':"white"}, line_kws = {"color": WT_linecolor})
b = sns.regplot(x="Norm Total Area",y="Norm Count", color="white", data=plot_data_KO_SO, ax=ax1, label="cKO SO", fit_reg=False, scatter_kws = {"s": 30, "alpha":0.7, 'edgecolor':KO_edge_color, 'linewidth':2}, line_kws = {"color": KO_linecolor})
# plot trend line
 
c = sns.regplot(x="Norm Total Area",y='Norm Count', color=WT_colors[2], data=plot_data_WT_SR, ax=ax2, label="CTL SR", fit_reg=False, scatter_kws = {"s": 55, "alpha":0.8, 'edgecolor':"white"}, line_kws = {"color": WT_linecolor})
d = sns.regplot(x="Norm Total Area",y='Norm Count', color="white", data=plot_data_KO_SR, ax=ax2, label="cKO SR", fit_reg=False, scatter_kws = {"s": 30, "alpha":0.7, 'edgecolor':KO_edge_color, 'linewidth':2}, line_kws = {"color": KO_linecolor})

e = sns.regplot(x="Norm Total Area",y="Norm Count", color=WT_colors[2], data=plot_data_WT_SLM, ax=ax3, label="CTL SLM", fit_reg=False, scatter_kws = {"s": 55, "alpha":0.8, 'edgecolor':"white"}, line_kws = {"color": WT_linecolor})
f = sns.regplot(x="Norm Total Area",y='Norm Count', color="white", data=plot_data_KO_SLM, ax=ax3, label="cKO SLM", fit_reg=False, scatter_kws = {"s": 30, "alpha":0.7, 'edgecolor':KO_edge_color, 'linewidth':2}, line_kws = {"color": KO_linecolor})

# reset index of norm_means
norm_means = norm_means.reset_index()

# plot the center of mass with error bars in both dimensions for each plot
# Note, these colors should be correct if the KO is listed first in the norm_means dataframe. I recommend checking that this is the case and visually checking that the error bar colors match their group.
# Plot the center of mass for SO first
SO_means = norm_means[norm_means["Group"].str.contains('Basal')]
ax1.errorbar(x= SO_means["X", "median"], y= SO_means["Y", "median"], xerr= SO_means["X", "std"], yerr= SO_means["Y", "std"], linestyle="None", ecolor=[KO_linecolor,WT_linecolor], zorder=1000, elinewidth=3)
# Plot the center of mass for SR
SR_means = norm_means[norm_means["Group"].str.contains('Proximal')]
ax2.errorbar(x= SR_means["X", "median"], y= SR_means["Y", "median"], xerr= SR_means["X", "std"], yerr= SR_means["Y", "std"], linestyle="None", ecolor=[KO_linecolor,WT_linecolor], zorder=1000, elinewidth=3)
# plot the center of pass for SLM
SLM_means = norm_means[norm_means["Group"].str.contains('Distal')]
ax3.errorbar(x= SLM_means["X", "median"], y= SLM_means["Y", "median"], xerr= SLM_means["X", "std"], yerr= SLM_means["Y", "std"], linestyle="None", ecolor=[KO_linecolor,WT_linecolor], zorder=1000, elinewidth=3)

# Set the axes parameters
for ax in [ax1,ax2,ax3]:
    ax.set_ylabel("Count per 100 \u00B5m\u00B2\n Norm. to CTL mean", fontsize=18, fontweight="bold", labelpad=8) # count
    ax.set_xlabel("Total mito area per 100 \u00B5m\u00B2\n Norm. to CTL mean", fontsize=18, fontweight="bold") # total area
    
    ax.set(xticks=np.arange(0, 4.1, 1), xlim=[0,4], yticks=np.arange(0, 4.1, 1), ylim=[0,4])

# show axes and legend just for the first axis
ax1.legend(labels=["CTL", "cKO"], frameon=False, markerscale=1.4, fontsize=14, labelspacing=0.3, handletextpad=0.1)
    
ax1.set_title('SO', size=20, y=0.95, fontweight="bold")
ax2.set_title('SR', size=20, y=0.95, fontweight="bold")
ax3.set_title('SLM', size=20, y=0.95, fontweight="bold")

# Save the figure
plt.savefig(os.path.join(fig_path, "Mito_count_total_area_CTL_KO_corr_by_layer.png"), format='png', dpi=300, bbox_inches='tight')


#%% Get some basic metadata information about this code and the analysis to save as a txt file with the data.

# Total number of tiles analyzed in each region
basal_count = list(mito_avgs["Layer"]).count("Basal")
prox_count = list(mito_avgs["Layer"]).count("Proximal")
dist_count = list(mito_avgs["Layer"]).count("Distal")

total_analyzed = str(basal_count) + " Basal; " + str(prox_count) + " Proximal; " + str(dist_count) + " Distal (Total = " + str(len(mito_avgs)) + ")"

# For reference, we'll add the name of the code used to analyse the data (this code)
# Note: Keep this updated if you update the name or date on the code
code_name = inspect.getfile(inspect.currentframe()).split("\\")[-1]

# Let's get how many animals and how many stubs are in the dataset
animal_list = list(set(mito_avgs["Animal"])) # list of animal IDs in the data
num_animals = len(animal_list) # number of animals analyzed
animals = ", ".join(animal_list) # animal IDs formatted for text output

# get how many animals in each group

mice_per_group = {x:len(list(set(mito_avgs["Animal"].where(mito_avgs["Genotype"] == x).dropna()))) for x in group_list}

groups = [str(mice_per_group[x]) + " " + str(x) + " mice" for x in mice_per_group]
groups_str = ", ".join(groups)

# Get how many stubs per animal
num_stub_per_animal = [len(list(set(mito_avgs["Stub"].where(mito_avgs["Animal"] == x).dropna()))) for x in animal_list]

if len(list(set(num_stub_per_animal))) > 1: # if animals have different numbers of stubs
    num_stub = str(min(num_stub_per_animal)) + "-" + str(max(num_stub_per_animal)) + " stubs / mouse"
else: # if all animals have same number of stubs
    num_stub = str(num_stub_per_animal[1]) + " stubs / mouse"
    
# List of animals and stubs in the data
stub_list = mito_data[["Animal", "Stub"]].drop_duplicates().reset_index(drop=True)
    
# Once we get WT and KO data we may want to separate out how many WT and KO mice
    
analyzed = str(num_animals) + " mice analyzed (" + str(num_stub) + ")"

# explain any processing of the data by this code
processing = "\n" + str(mitos_removed) + " mitochondria larger than 2.1 um2 or smaller than 0.01 um2 and " + str(den_removed) + " dendrite segments less than 0.4 um2 were removed.\nArea and diameter converted from pixels to micron squared.\n Count and mitochondrial mass were normalized to the total area analyzed. MCU KO data was normalized to the overal mean of the CTL."

# Put together the metadata information for the text file
meta_summary = {"Date": date.today().isoformat(), "Dataset":dataset, "Animals": animals, "Analyzed": analyzed, "Groups": groups_str, "Sampling": sampling, "Tile Size": tile_size, "N Tiles": total_analyzed, "Processing": processing, "Data CSV": data_loc.split("\\")[-1], "Biodock AI": AI_name, "Analysis Code": code_name}

# String of how many tiles were excluded
excluded = str(num_flag) + " tiles\n" + ", ".join(flagged.index)

if num_flag >= 1: # If any tiles were excluded due to low counts
    meta_summary["Excluded"] = excluded # add excluded column to the metadata

if len(notes) > 1: # if there are notes
    meta_summary["Notes"] = notes # add notes column to the metadata

# write the metadata dictionary to a text file
meta_file = open(os.path.join(dirName, 'analysis_summary.txt'), 'wt')
for line in meta_summary:
    meta_file.write(str(line) + ":  " + str(meta_summary[line]) + "\n\n")
meta_file.close()

#%% Export several CSV files in the right format for Prism. Right now we want to look at Area, Count, Major Axis, perimeter

print("\nFormatting and saving the data for easy input into Prism...")

# path to save the prism files (within the output directory)
prism_path = os.path.join(dirName, "Prism Files")
# Create a subfolder for just the prism files
try:
    os.mkdir(prism_path) # make directory 
except FileExistsError:
    pass # skip if the folder already exists

# We need to format the data for each metric of interest into a table for easy import into Prism and plotting. For this, we want a column for each layer, with each row a tile average. There can be a different numbers of tiles for each column, but in general they will be the same length as we should have an equal number of tiles for each dendritic layer.

# Now we can use the function to save the Prism tables for each metric of interest!

# First, save a prism table for the tile level data, particularly for count and total mito area

# Mito count per tile
mito_count_tile = prism_format(data = mito_avgs, data_col = "Count", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Mito_number_per_tile.csv")))

# Total mito area per tile
mito_count_animal = prism_format(data = mito_avgs, data_col = "Total_mito_Area_um_sq", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Total_mito_area_per_tile.csv")))

# Now save some Prism tables to look at the whole population of mitochondria for the other metrics of interest

# We also want to export the individual mito areas and diameter grouped by genotype and layer
mito_indiv_area = prism_format(data = mito_data, data_col = "Area_um_sq", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Indiv_Mito_area_um2.csv")))

# save a Prism table for perimeter
mito_indiv_peri = prism_format(data = mito_data, data_col = "Perimeter_um", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Indiv_Mito_perimeter_px.csv")))

# save a Prism table for Aspect Ratio
mito_indiv_aspect = prism_format(data = mito_data, data_col = "Aspect_Ratio", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Indiv_Mito_aspect_um.csv")))

# save a Prism table for distance to nearest neighbor
mito_indiv_dist = prism_format(data = mito_data, data_col = "NN_Dist_um", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Indiv_Mito_dist_NN.csv")))

# save a Prism table for Feret's diameter
mito_indiv_Feret = prism_format(data = mito_data, data_col = "Feret_diam_um", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Indiv_Mito_Feret_diam.csv")))

# Spit the individual mito data by animal to look at each animal's distribution. With this we can check for animal to animal variability and any outliers. Right now we'll just look at area.
# This will save a separate file for each animal in the dataset
[prism_format(data = mito_data_sub.where(mito_data_sub["Animal"]==a).dropna(), data_col = "Area_um_sq", col=["Genotype", "Layer"], file_path = os.path.join(prism_path, str("Animal_indiv_mito_area_" + str(a) + ".csv"))) for a in animal_list]

# Now print the animal means organized to plot the groups on the same plot.
# We want the columns to be the animals, organized with cre - animals first, and the rows to be dendritic layer.

mito_area_animal_trans = prism_format(data = animal_avgs, data_col = "Area_um_sq", col=["Genotype","Animal"], row="Layer", sort="Genotype", file_path = os.path.join(prism_path, str("Trans_Med_Mito_Area_um2.csv")))

mito_diam_animal_trans = prism_format(data = animal_avgs, data_col = "Feret_diam_um", col=["Genotype","Animal"], row="Layer", sort="Genotype", file_path = os.path.join(prism_path, str("Trans_Med_Mito_major_Ferets_diam.csv")))

mito_peri_animal_trans = prism_format(data = animal_avgs, data_col = "Perimeter_um", col=["Genotype","Animal"], row="Layer", sort="Genotype", file_path = os.path.join(prism_path, str("Trans_Med_Mito_perimeter.csv")))

mito_aspect_animal_trans = prism_format(data = animal_avgs, data_col = "Aspect_Ratio", col=["Genotype","Animal"], row="Layer", sort="Genotype", file_path = os.path.join(prism_path, str("Trans_Med_Mito_aspect.csv")))

mito_count_animal_trans = prism_format(data = animal_avgs, data_col = "Count_100um2", col=["Genotype","Animal"], row="Layer", sort="Genotype", file_path = os.path.join(prism_path, str("Trans_Med_Mito_count_100um2.csv")))

mito_dist_animal_trans = prism_format(data = animal_avgs, data_col = "NN_Dist_um", col=["Genotype","Animal"], row="Layer", sort="Genotype", file_path = os.path.join(prism_path, str("Trans_Med_Mito_dist_NN.csv")))

print("\nCode finished! \n\nSeveral main CSV files have been saved with the individual and averaged data, along with a txt file summary of the analysis. A subfolder was created for the Prism data files and for the main figure images. Additional figures and stats will be in the 'Figures and Stats' subfolder. \n\n File location: " + str(os.path.join(os.getcwd(), dirName)))


