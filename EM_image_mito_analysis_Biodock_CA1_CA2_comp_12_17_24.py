# -*- coding: utf-8 -*-
"""
Created on 9.1.2023

Custom code to work with the object data exported from Biodock, an AI based image analysis platform. Mitochondria were segmented in 6144x4096 pixel SBFEM tiles from the basal, proximal and distal dendrits of CA2 and CA1 of WT and MCU KO mice. The biodock output includes a row for each segmented object in the SEM images along with metrics such as area, X and Y position, length of major and minor axes, perimiter, etc. The file structure in Biodock should be a folder for CTL and cKO, with a subfolder for each layer, then a subfolder for each animal analyzed and each section. This is important, because this code will be getting a lot of the metadata information from the folder names (ie. dendritic layer, animal name, genotype). If you want to compare CA1 and CA2, make sure you have a combined CSV file with both CA1 and CA2 together and add a column for "Subregion" that says either CA1 or CA2 for each row.

The code will convert any lengths or areas from pixels to microns, apply a size filter, then calculate aspect ratio and the distance to nearest neighbor for each segmented mitochondria object. The total number of mitochondria and the total mitochondria area per 100 um2 will be calculated for each image. Any image tiles containing less than 2 mitochondria will be flagged and removed. 

This code will filter the data to include only mitochondria that have an associated segmented dendrite and remove dendrites too small as likely spine heads (< 0.04 um2 in area), or too wide as primary dendrites (> 1.9 um in width). The mitochondria areas will then be normalized to the parent dendrite's area and mitochondria Feret's diameter will be normalized to the Feret's diameter of the parent dendrite.
                                                                                                                                                                                                                                                                                                                                                                                                                        The custom function WT_KO_Violin() from "mito_functions.py" is used to plot metrics of interest with violin plots (or box plots) comparing CA1 and CA2 across dendritic layers. Make sure you have mito_functions.py in your working directory along with this code. 

Input: The CSV file exported from Biodock's analysis. User will be prompted to select this file. This file should contain the CA1 and CA2 individual mitochondria and dendrite objects and a column added for "Subregion". Note that you may have to combine the data and add the Subregion column manually.
    
Output: CSV files with the individual data, tile data, section data and animal data. There will be an output folder within the main output directory for the violin plots and stats, with a separate folder specifically for figure plots. There will also be a folder with CSV files formated for easy import and plotting in Prism, if desired. Finally, an analysis summary CSV file will be saved in the main output directory with some general information about the analysis.

BEFORE RUNNING THIS CODE: Set the name of the output directory you want, make sure the metadata information in this code is accurate, particularly the Biodock AI version used and the dataset analyzed. Make sure you have a CSV file that combines CA1 and CA2 together with a column for "Subregion". Also make sure you have the python file "mito_functions.py" in your working directory.

@author: pannoni
"""

# import needed packages
import pandas as pd
import easygui
import os
from datetime import date
from mito_functions import get_stats_summary, prism_format  # custom functions
import inspect
from scipy.spatial import KDTree

# Expand the printed text of pd dataframes
pd.set_option('display.max_column', None)
pd.set_option('display.expand_frame_repr', False)

# import the excel file as a dataframe
# If easygui fails or you don't want to install it, change the data_loc variable on the line below to the full file pathway to your Biodock output CSV file. The rest of the code should work fine.
data_loc = easygui.fileopenbox(msg='Please select the CSV file with the object data.', title='Select CSV File', default='*')
mito_data = pd.read_csv(data_loc)

# Create a new folder in the working directory to save all the output files for this code.
#!!! Set the name of the output folder where the data and plots will be saved here.
dirName = "SEM_mitos_CA1_CA2_comp_norm_dendrites_12_17_24"

try:
    # Create target Directory
    os.mkdir(dirName)
    print("\nDirectory", dirName,  "Created.\n")
except FileExistsError:
    print("\nDirectory", dirName, "already exists, no need to create a new directory.")

# Make sure the meta data about the analysis is correct below

#!!! Enter the name an version of the Biodock AI used to segment the images here. This is important information to have if we end up training multiple AI versions.
AI_name = "Mitochondria in MCU KO and WT mice (V6)"

#!!! Dataset that the images are from:
dataset = "CTL CA2 and CA1"

# Set the tile size for the metadata
tile_size = "6144 x 4969 pixels"

# Tile sampling
sampling = "Every 8th tile"

#!!! Set the main comparison you want to make below, either Subregion or Genotype. Use "Subregion" if you want to compare across CA1 and CA2 in your dataset. Set "Genotype" if you want to compare across genotypes. You can always run the code twice if you want to make both comparisons.
comparison = "Subregion"

#!!! Space to put notes about the dataset or any changes made since a previous run.
notes = "Analysis of CTL mitochondria and dendrites in CA2 and CA1 using version 6 of the Biodock AI, focusing on mitochondria within dendrites and controling for dendrite plane. Denderitic mitochondria and dendrites were segmented with a confidence threshold of 0.4. Dendrites were filtered to remove spine heads (< 0.4 area) and primary dendrites (> 1.9 width)."

# Parse out the image metadata from the CSV file by breaking the "image origin" column into its separate folders

# check that there are the right number of levels in the data
if len(mito_data["Image origin"][1].split("/")) == 5:
    # order is genotype, dendritic layer, animal ID, "stub" and then image file
    # Genotype or group should be first in the hierarchy
    mito_data["Genotype"] = [x.split("/")[0] for x in mito_data["Image origin"]]
    # dendritic layer should be second
    mito_data["Layer"] = [x.split("/")[1] for x in mito_data["Image origin"]]
    mito_data["Animal"] = [x.split("/")[2] for x in mito_data["Image origin"]]  # Animal ID third
    mito_data["Stub"] = [x.split("/")[3] for x in mito_data["Image origin"]]  # Then stub name
    mito_data["Tile"] = [x.split("/")[4] for x in mito_data["Image origin"]]  # image file name
else:
    print("Error: Can't get metadata information from the Biodock output file. Please make sure this information is included in the 'Image origin' column of your file. The Image origin should include the genotype, layer, animal, stub or section, then the image tile associated each object. It should be in that order, separated by /. If the image file structure in biodock is set up as described, it should be exported that way. If not, you will have to fix it and run the code again.")

# Check for a "Subregion" column and throw an error if there isn't one in the CSV file
if "Subregion" not in mito_data.columns:
    print("Error: Can't find the subregion information from the Biodock output file. Please make sure there is a column for 'Subregion' in your CSV data file and run the code again. You may have to manually add this column.")

# Make the layer column into category so we can sort them in the order we want (Basal, Prox, Distal)
sorter = ["Basal", "Proximal", "Distal"]  # sorting list
# Set the Layer column to be a "category" type that we can sort by
mito_data.Layer = pd.Categorical(mito_data.Layer, categories=sorter, ordered=True)
mito_data = mito_data.sort_values(['Subregion', 'Layer', 'Animal', "Stub"], ascending=[False, True, True, True]).reset_index(drop=True)  # sort

# Change "object ID" to string so it won't be treated as numerical
mito_data["Object ID"] = mito_data["Object ID"].astype("int").astype("str")
# remove any spaces from the animal numbers
mito_data["Animal"] = mito_data["Animal"].str.replace(' ', '')

print("\nCSV file uploaded.")

# Separate the dendrite objects from the mitochondria objects. We aren't going to do much with the dendrite data now, but likely will in the future.
dendrite_data = mito_data[mito_data["Class"] =="Dendrite Object"].reset_index(drop=True)
mito_data = mito_data[mito_data["Class"] != "Dendrite Object"].reset_index(drop=True)


# %%% Process the data by converting pixels to microns and setting a size filter to remove non-mitochondria objects. Then calculate aspect ratio as the ratio of the major and minor axes and give each section a unique ID.

print("\nConverting area and length from pixels to microns...")

# First we need to convert Area and Length into um or um squared
# image scale is 2 nm per pixel
scale_nm_per_px = 2
# Convert scale to um instead of nm
# use for as conversion factor for length
scale_um_per_px = scale_nm_per_px / 1000
# get the scale in microns squared
# convert nm to um and square, use as conversion factor for area
scale_um_sq = scale_um_per_px ** 2

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
dendrite_data["Den_Length_um"] = dendrite_data["Feret Diameter Maximum"] * scale_um_per_px

print("\nRemoving mitochondria larger than 2 microns squared or smaller than 0.01 microns squared...")

total_mitos = len(mito_data)  # total mitos before filtering

# Capture the mitos to be removed for size exclusion, if you want to look at them
mitos_removed_small = mito_data.query('0.01 > Area_um_sq').reset_index(drop=True)
mitos_removed_large = mito_data.query('Area_um_sq > 2.1').reset_index(drop=True)

# Now let's remove the mitos below a size of 0.1 um2, or larger than 2um2
mito_data = mito_data.query('0.01 < Area_um_sq < 2.1').reset_index(drop=True)

mitos_removed = total_mitos - len(mito_data)

if mitos_removed > 0:
    print(" ", mitos_removed, " mitochondria removed from size exclusion.")
else:
    print("No mitochondria removed from size exclusion.")

# Let's remove any dendrite segments below a size of 0.02 um2 (~5000 px), these are likely spine heads and not actually dendrites
print("\nRemoving dendrites smaller than 0.4 microns squared since these are likely spine heads and primary dendrites wider than 1.9 microns...")
# get the current number of dendrites
total_dendrites = len(dendrite_data)
# Captured the removed dendrites in case we need to look at them
removed_spine_heads = dendrite_data.query('Den_Area_um_sq < 0.4').reset_index(drop=True)
# Now we can remove them from the main dataframe
dendrite_data = dendrite_data.query('Den_Area_um_sq > 0.4').reset_index(drop=True)
# Get the number of dendrites that were removed, if any
den_removed_spines = total_dendrites - len(dendrite_data)

# Add another filter here to remove the primary dendrites
# Capture the primary dendrites being removed
dendrites_primary = dendrite_data.query("Den_Width_um > 1.9").reset_index(drop=True)
# Now filter out the primary dendrites from the original dataframe
dendrite_data = dendrite_data.query("Den_Width_um < 1.9").reset_index(drop=True)
# Get the number of dendrites that were removed as primary dendrites, if any
den_removed_primaries = den_removed_spines - len(dendrite_data)

if den_removed_spines > 0:
    print(" ", den_removed_spines, " dendrites removed as likely spine heads.")
else:
    print("No dendrite spine heads were removed.")

if den_removed_primaries > 0:
    print(" ", den_removed_primaries, " primary dendrites removed.")
else:
    print("No primary dendrites were removed.")

print("\nCalculating the mitochondria aspect ratio...")
# Calculate the aspect ratio using the major and minor axes, although we will probably just use the length of the major axes instead
mito_data["Aspect_Ratio"] = round(mito_data["Length of Major Axis"] / mito_data["Length of Minor Axis"], 2)

# We'll have to make a unique name by combining animal ID and stub since some animals have stubs with the same names
mito_data["Section_ID"] = mito_data["Animal"] + "_" + mito_data["Stub"]

# %% Now let's do a nearest neighbor analysis to get the distance from each mitochondria to the closest mitochondria within the same image tile.

print("\nCalculating the nearest neighbors...")

# Find the nearest neighbor for each mitochondria with a KDtree
# Note this will have to be done within the same tile image

# Add an empty column to mito_data to put the distances and the nearest neighbor IDs
tile_list = []

# First combine the X and Y coordinates into pairs for each mito to get the locaion
mito_data["Location"] = mito_data[["X Position", "Y Position"]].values.tolist()

for tile in list(set(mito_data["Tile"])):  # loop across each image
    # Get the data for the current image
    tile_data = mito_data.where(mito_data["Tile"] == tile).dropna(how="all")
    # Get the mitochondria locations and IDs for the current image
    mito_coor = list(tile_data["Location"])
    mito_IDs = list(tile_data["Object ID"].astype("int"))

    # get the current subregion
    curr_subregion = tile_data["Subregion"].unique()[0]

    if len(tile_data) > 1:  # if there's more than one mitochondria in the tile

        # Create a KD tree to locate nearest neighbors
        tree = KDTree(mito_coor)

        mito_dist = []  # empty list to put the mito distances for the current image
        neighbor_IDs = []  # empty list to put the object ID of the nearest neighbor

        # For each mito in the image, query the tree to get the distance to the nearest neighbor
        for mito in mito_coor:
            dist, obj = tree.query(mito, k=2)  # get the two nearest neighbors
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
        print("\nTile '" + str(tile) +
              f"' in {curr_subregion} has only one mitochondria. Nearest neighbor analysis cannot be done on this tile.\nThis tile should likely be removed from the analysis.")

# combine the tile data back into the original dataframe
mito_data = pd.concat(tile_list).sort_index()

# Change the ID column to string so it's not a numerical data type
mito_data["NN_ID"] = mito_data["NN_ID"].astype(str)

# Convert the nearest neighbor distance into microns like all the other distances
mito_data["NN_Dist_um"] = mito_data["NN_Dist"] * scale_um_per_px


# %% Get mitochondrial counts and average metrics per tile and remove any tiles with less than 2 mitochondria.

print("\nGetting the number and total area of mitochondria and dendrites in each tile...")

# group the mito data by tile to get averages for the metrics
mito_avgs = mito_data.groupby("Tile").mean(numeric_only=True)  # get the tile averages
meta_cols = mito_data.groupby("Tile").agg({"Animal": "first", "Subregion": "first", "Genotype": "first", "Layer": "first", "Stub": "first", "Section_ID": "first", "Tile": "count"})  # keep the metadata and count the number of objects (mitochondria) for each tile
# join the averaged data with the metadata columns + count
mito_avgs = meta_cols.join(mito_avgs)

# Rename the tile column to "count"
mito_avgs = mito_avgs.rename(columns={"Tile": "Count"})

# get total mito area per tile
mito_avgs["Total_mito_Area_um_sq"] = mito_data[["Tile", "Area_um_sq"]].groupby("Tile").sum()

# Round the data in mito_avgs to 3 decimal points
mito_avgs = round(mito_avgs, 2)
# make sure Layer is categorical
mito_avgs.Layer = pd.Categorical(mito_avgs.Layer, categories=sorter, ordered=True)

# Let's flag any tiles where the count is under 2. This likely means the tile is not useable or something is wrong with the tile image.
flagged = mito_avgs["Count"].where(mito_avgs["Count"] <= 2).dropna(how="all")  # get any tiles where count is 5 or less
num_flag = len(flagged)

# print a warning if there are any flagged tiles
if num_flag >= 1:
    print("\n" + str(num_flag) + " tiles were flagged because they have fewer than 2 mitochondria. Recommended to visually check these tiles to see if they should be analyzed. These tiles will be removed from the data now.\n" + str(list(flagged.index)))

# Remove any flagged tiles from the main dataframe and the tile average dataframe
# Comment out this line if you don't want to remove the flagged tiles from the data
for tile in flagged.index:
    mito_data = mito_data[~mito_data["Tile"].str.contains(tile)]
    mito_avgs = mito_avgs[~mito_avgs.index.str.contains(tile)]

#%%% From the dendrite data, we would like to get the number of mitochondria in each dendrite / length of the dendrite. This will be used to eventually compare to a different dataset.

# path to save the summary stats files (within the output directory)
stats_path = os.path.join(dirName, "Summary Statistics")
# Create a subfolder for just the prism files
try:
    os.mkdir(stats_path) # make directory 
except FileExistsError:
    pass # skip if the folder already exists

# Before subsetting, look at the frequency of mitochondria per dendrite counts in the data
# Get how many dendrites have 1, 2, 3 mitochondria for each group focusing on CTL SR and SLM

# Empty dataframe to store the counts
num_mito_per_den_counts = pd.DataFrame()

# Get the counts for the proximal dendrites and add to the dataframe
num_mito_per_den_counts["Proximal"] = dendrite_data[dendrite_data["Layer"] == "Proximal"]["Number of Dendritic mitochondria within Dendrite Object"].value_counts(dropna=False)

# Get the counts for the distal dendrites and add to the dataframe
num_mito_per_den_counts["Distal"] = dendrite_data[dendrite_data["Layer"] == "Distal"]["Number of Dendritic mitochondria within Dendrite Object"].value_counts(dropna=False)

# Sort them so that NaN comes first
num_mito_per_den_counts = num_mito_per_den_counts.sort_index(na_position='first')
num_mito_per_den_counts = num_mito_per_den_counts.rename_axis("Mitochondria per Dendrite")

num_mito_per_den_counts.to_csv((os.path.join(dirName, str("mito_number_per_den_value_counts.csv"))), na_rep='NaN')

# print out the frequencies:
print("\n Frequency of mitochondria counts per dendrite segment in SR and SLM:\n\n", num_mito_per_den_counts)

print("\nFiltering out dendrite segments that don't have any mitochondria...")

# Filter the dendrite data to get only dendrites that have at least one mitochondria in them. Start by removing NaN entries
dendrites_with_mitos = dendrite_data.dropna(subset=['Number of Dendritic mitochondria within Dendrite Object']).reset_index(drop=True)[["Object ID", "Subregion", "Genotype", "Animal", "Layer", "Stub", "Tile", "Number of Dendritic mitochondria within Dendrite Object", "Den_Length_um", "Den_Width_um", "Den_Area_um_sq", "Dendrite Object parent ID"]].rename(columns={'Number of Dendritic mitochondria within Dendrite Object': 'Num_mitos_in_dendrite'})

# Also filter to remove any rows where the number of mitos are 0
dendrites_with_mitos = dendrites_with_mitos.query('Num_mitos_in_dendrite > 0').reset_index()

# Let's take a closer look at the dendrite statistics before and after subsetting to only dendrites with mitochondria

# Get summary statistics for the dendrite objects before subsetting
den_stats_before_sub = dendrite_data[[comparison, "Layer", "Den_Length_um", "Den_Width_um", "Den_Area_um_sq", "Number of Dendritic mitochondria within Dendrite Object"]] .groupby([comparison, "Layer"]).describe().rename(columns={'Number of Dendritic mitochondria within Dendrite Object': 'Num_mitos_in_dendrite'}).T

den_stats_sub = dendrites_with_mitos[[comparison, "Layer", "Den_Length_um", "Den_Width_um", "Den_Area_um_sq", "Num_mitos_in_dendrite"]] .groupby([comparison, "Layer"]).describe().T

# Save the summary stats as a CSV file in case we need to look at that later
den_stats_before_sub.to_csv(os.path.join(stats_path, str("Dendrite_stats_before_sub.csv")))
den_stats_sub.to_csv(os.path.join(stats_path, str("Dendrite_stats_after_sub.csv")))

print("\nPulling the mitochondria data for just the mitochondria in dendrites...")

# We need to pull the mitos from the mito_data that match the mitos in dendrites_with_mitos

# Create a tile ID that includes the subregion, in case there are tiles with the same tile name in both subregions
dendrites_with_mitos["Tile_ID"] = dendrites_with_mitos["Subregion"] + "_" + dendrites_with_mitos["Tile"]
mito_data["Tile_ID"] = mito_data["Subregion"] + "_" + mito_data["Tile"]

# get the IDs of the dendrites in each tile that have at least 1 mito so we can associate them with their children mitochondria in mito_data
den_IDs = dendrites_with_mitos[["Tile_ID","Object ID", "Num_mitos_in_dendrite"]]

# For each tile, get list of dendrite IDs and a sum of how many mitos in dendrites in that tile.
tile_list = dendrites_with_mitos["Tile_ID"].unique()
tile_den_mito_df = pd.DataFrame(index=tile_list) # empty dataframe to hold the tile data

# Get the dendrite IDs for each tile in tile_list
tile_den_mito_df["Den_ID_list"] = [list(den_IDs[den_IDs["Tile_ID"] == x]["Object ID"]) for x in tile_list]
# Also get the number of dendrites that have mitos and total number of mitos in dendrites for each tile in tile_list
tile_den_mito_df["Num_den_with_mitos"] = [len(x) for x in tile_den_mito_df["Den_ID_list"]]
tile_den_mito_df["total_mitos_in_dendrites"] = [int(sum(den_IDs[den_IDs["Tile_ID"] == x]["Num_mitos_in_dendrite"])) for x in tile_list]

# Use the list of dendrite IDs to find the same parent IDs in the mito dataframe for each dendrite

# will house dendrite level data for dendrites that have mitos
den_with_mitos_df = pd.DataFrame()

# will house the individual mito data for only mitos that have parent dendrites
mitos_with_den_df = pd.DataFrame()

for tile in tile_list:
    # list of dendrite IDs for the current tile
    den_IDs = tile_den_mito_df.loc[tile]["Den_ID_list"]

    den_curr_tile = dendrites_with_mitos[dendrites_with_mitos["Tile_ID"] == tile]

    # get all the mito_data for the current tile
    mitos_curr_tile = mito_data[mito_data["Tile_ID"] == tile]

    # For each dendrite, grab the mito IDs, the mito count (for confirmation) and the avg mito area and range
    for den in den_IDs:
        # Get the mitos in the current dendrite
        mitos_curr_den = mitos_curr_tile[pd.to_numeric(mitos_curr_tile["Dendrite Object parent ID"], errors='coerce') == float(den)].reset_index(drop=True).rename(columns={"Object ID": "Mito_ID"})  # will rename the mito "Object ID" to "Mito_ID"

        # Make sure to filter for mitos that have a parent dendrite
        if len(mitos_curr_den) > 0:
            # Grab dendrite ID, area and length of the current dendrite

            # Get the number of mitos in the current dendrite
            num_curr_mitos = len(mitos_curr_den)

            # Get the lengh and area of the current dendrite which will apply to each mito
            mitos_curr_den["Den_Area_um_sq"] = list(den_curr_tile[den_curr_tile["Object ID"] == den]["Den_Area_um_sq"]) * num_curr_mitos
            mitos_curr_den["Den_Length_um"] = list(den_curr_tile[den_curr_tile["Object ID"] == den]["Den_Length_um"]) * num_curr_mitos
            mitos_curr_den["Den_Width_um"] = list(den_curr_tile[den_curr_tile["Object ID"] == den]["Den_Width_um"]) * num_curr_mitos
            mitos_curr_den["Dendrite_ID"] = [den] * num_curr_mitos

            # Add to the mitos_with_den_df
            mitos_with_den_df = pd.concat([mitos_with_den_df, mitos_curr_den], axis=0)

            # Now get some data for the dendrite level dataframe
            # list the mito object IDs
            curr_mito_IDs = list(mitos_curr_den["Mito_ID"])

            # Get the average mito area for this dendrite
            curr_mito_area = mitos_curr_den["Area_um_sq"].mean()

            # Get the average mito diameter for this dendrite
            curr_mito_area = mitos_curr_den["Area_um_sq"].mean()

            # Get the sum of dendrite area
            curr_total_mito_area = mitos_curr_den["Area_um_sq"].sum()
            
            # Get the sum of dendrite length
            curr_total_mito_len = mitos_curr_den["Feret_diam_um"].sum()

        # put this together into a temp series with all the data for the current dendrite
            den_meta = den_curr_tile[den_curr_tile["Object ID"] == den][["Object ID", "Animal", "Subregion", "Layer", "Stub", "Tile", "Den_Length_um", "Den_Width_um", "Den_Area_um_sq"]].reset_index(drop=True)
            # temp dataframe for the new dendrite level data
            temp_den_data = pd.DataFrame({"Dendrite_ID": den, "Mito_IDs": [curr_mito_IDs], "Mean_mito_area": curr_mito_area, "Mito_total_area": curr_total_mito_area, "Mito_total_len": curr_total_mito_len, "N_mitos": num_curr_mitos})
            # combine the two dataframes
            temp_den_data = pd.concat([den_meta, temp_den_data], axis=1)

            # Add to the den_with_mitos_df
            den_with_mitos_df = pd.concat([den_with_mitos_df, temp_den_data], axis=0)

print(" done.")

print("\n Normalizing the mitochondria area and length by dendrite area or length...")

# Normalize the denderite data by total mito area or total mito length

# Get total mito area divided by dendrite area for the dendrite df
den_with_mitos_df["Total_mito_area_per_den_area"] = den_with_mitos_df["Mito_total_area"] / den_with_mitos_df["Den_Area_um_sq"]

# Get total mito area divided by dendrite area for the dendrite df
den_with_mitos_df["Total_mito_len_per_den_len"] = den_with_mitos_df["Mito_total_len"] / den_with_mitos_df["Den_Length_um"]

# Normalize the individual mito data by dendrite

# Using the mito data, normalize mito length by dendrite length
mitos_with_den_df["Mito_len_per_den_len"] = mitos_with_den_df["Feret_diam_um"] / mitos_with_den_df["Den_Length_um"]

# Using the mito data, normalize mito area by dendrite area
mitos_with_den_df["Mito_area_per_den_area"] = mitos_with_den_df["Area_um_sq"] / mitos_with_den_df["Den_Area_um_sq"]

# Save the mitos_with_den_df as a csv file
norm_file1 = "Mitos_with_parent_dendrites_df.csv"
mitos_with_den_df.to_csv(os.path.join(dirName, norm_file1), index=False)

# Save the den_with_mitos_df as a csv file
den_with_mitos_df.to_csv(os.path.join(dirName, str("Dendrites_with_mitos_df.csv")), index=False)

# For the dendrite data, get a grouped dataframe as well
den_with_mitos_df_animal_stats = den_with_mitos_df[["Animal", comparison, "Layer", "Mito_total_area", "Den_Width_um", "Den_Length_um", "Den_Area_um_sq", "Mito_total_area", "Total_mito_area_per_den_area", "Mito_total_len", "Total_mito_len_per_den_len"]].groupby(["Animal", comparison, "Layer"], observed=True).agg(["median", "std", "count"]).T
# Save the den_with_mitos_df as a csv file
den_with_mitos_df_animal_stats.to_csv(os.path.join(dirName, str("Dendrites_with_mitos_df_animal_stats.csv")))

den_with_mitos_df_animal = den_with_mitos_df[["Animal", comparison, "Layer", "Den_Width_um", "Den_Length_um", "Den_Area_um_sq", "Mito_total_area", "Total_mito_area_per_den_area", "Mito_total_len", "Total_mito_len_per_den_len"]].groupby(["Animal", comparison, "Layer"], observed=True).median()

den_with_mitos_df_animal.to_csv(os.path.join(dirName, str("Dendrites_with_mitos_df_animal_medians.csv")))

# Get summary statistics for the dendrite data (subsetted)
den_with_mitos_df_summary = den_with_mitos_df[[comparison, "Layer", "Den_Width_um", "Den_Length_um", "Den_Area_um_sq", "Mito_total_area", "Total_mito_area_per_den_area", "Mito_total_len", "Total_mito_len_per_den_len"]].groupby([comparison, "Layer"]).agg(["mean", "median", "std", "count"]).T

# Save the den_with_mitos_df as a csv file
den_with_mitos_df_summary.to_csv(os.path.join(dirName, str("Dendrites_with_mitos_summary_stats.csv")))

# Now let's make a group dataframe with the medians for the dendrite data in each group and get the number of dendrites in each group, which will be normalized by the number of tiles

# Get a grouped dataframe with the metrics of interest, including the total number of dendrites and number of unique tiles in each group
den_with_mitos_group = den_with_mitos_df[[comparison, "Layer", "Tile", "Dendrite_ID", "Den_Area_um_sq", "Den_Length_um", "Den_Width_um", "Mito_total_area", "Mito_total_len", "Total_mito_area_per_den_area", "Total_mito_len_per_den_len"]].groupby([comparison, "Layer"]).agg({"Den_Area_um_sq": "median", "Den_Length_um": "median", "Den_Width_um": "median", "Mito_total_area": "median", "Mito_total_len" : "median", "Total_mito_area_per_den_area": "median", "Total_mito_len_per_den_len": "median", "Dendrite_ID": "count", "Tile": "nunique"})

# Rename the aggregated columns to more appropriate names
den_with_mitos_group = den_with_mitos_group.rename(columns={"Dendrite_ID": "Den_Count", "Tile": "Tile_Count"})

# Normalize the dendrite count by the number of tiles in each group
den_with_mitos_group["Num_den_per_tile"] = round(den_with_mitos_group["Den_Count"] / den_with_mitos_group["Tile_Count"], 1)

# Add dendrite max length to the dataframe
den_with_mitos_group["Max_Den_Length"] = den_with_mitos_df[[comparison, "Layer", "Den_Length_um"]].groupby([comparison, "Layer"]).max()

# Also get the number of primary dendrites that were removed in each group
den_with_mitos_group["Primaries removed"] = dendrites_primary[[comparison, "Layer", "Object ID"]].groupby([comparison, "Layer"]).count()

# Save the grouped dataframe
den_with_mitos_group.to_csv(os.path.join(dirName, str("Dendrites_with_mitos_group_medians.csv")))

# Subset the mito data to have mito area with and without normalization to the dendrites
mito_norm_den_animal_df = mitos_with_den_df[["Animal", comparison, "Layer", "Area_um_sq", "Mito_area_per_den_area", "Feret_diam_um", "Mito_len_per_den_len"]].groupby(["Animal", comparison, "Layer"]).median()

# get total mito count and total number of tiles for each animal in each group
mito_counts_group = mitos_with_den_df[["Animal", comparison, "Layer", "Tile", "Mito_ID"]].groupby(["Animal", comparison, "Layer"]).agg({"Mito_ID": "count", "Tile": "nunique"})

# Rename the aggregated columns to more appropriate names
mito_counts_group = mito_counts_group.rename(columns={"Mito_ID": "Mito_Count", "Tile": "Tile_Count"})

# Normalize the dendrite count by the number of tiles in each group
mito_counts_group["Num_mitos_per_tile"] = round(mito_counts_group["Mito_Count"] / mito_counts_group["Tile_Count"], 1)

# Get dendrite area and length from the dendrite dataframe
den_norm_stats_animal = den_with_mitos_df[["Animal", comparison, "Layer", "Den_Area_um_sq", "Den_Length_um"]].groupby(["Animal", comparison, "Layer"]).median()

# Combine the dataframes so we can save the stats as a CSV file
mito_norm_den_animal_df = pd.concat([mito_norm_den_animal_df, mito_counts_group, den_norm_stats_animal], axis=1)

# Save as a csv file
norm_file2 = "Norm_mito_data_by_dendrites_animal_med.csv"
mito_norm_den_animal_df.to_csv(os.path.join(dirName, norm_file2))

print(f"\n The normalized mitochondria data by animal has been saved as: '{norm_file2}'. \nYou can find the full data for mitochondria in dendrites in the file: '{norm_file1}'\n")


# %%% Save a few CSV summary files for export into Prism

print("\nFormatting and saving the data for easy input into Prism...")

# path to save the prism files (within the output directory)
prism_path = os.path.join(dirName, "Prism Files")
# Create a subfolder for just the prism files
try:
    os.mkdir(prism_path) # make directory 
except FileExistsError:
    pass # skip if the folder already exists

# Get the median, std and count for mito area, norm mito area, mito length norm mito length as well as dendrite area and dendrite length
mito_norm_stats = mitos_with_den_df[["Layer", comparison, "Area_um_sq", "Mito_area_per_den_area", "Feret_diam_um", "Mito_len_per_den_len"]].groupby([comparison, "Layer"]).agg(["median", "std", "count"])

# Get dendrite area and length from the dendrite dataframe
den_norm_stats = den_with_mitos_df[["Layer", comparison, "Den_Area_um_sq", "Den_Length_um"]].groupby([comparison, "Layer"]).agg(["median", "std", "count"])

# Combine the two so we can save the stats as a CSV file
mito_norm_stats = pd.concat([mito_norm_stats, den_norm_stats], axis=1)

# Save as a csv file
mito_norm_stats.to_csv(os.path.join(dirName, str("Summary_stats_norm_mitos_with_den.csv")))

# First the non-normalized area and length

mito_area_trans = prism_format(data=mito_norm_den_animal_df.reset_index(), data_col="Area_um_sq", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_Mito_Area.csv")))

# Same prism file for the normalized mito Ferets
mito_Feret_trans = prism_format(data=mito_norm_den_animal_df.reset_index(), data_col="Feret_diam_um", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_Mito_Feret.csv")))

# Now the normalized area and length

# We also want to export the animal level data for Prism to run some stats
mito_area_norm_trans = prism_format(data=mito_norm_den_animal_df.reset_index(), data_col="Mito_area_per_den_area", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_Norm_Mito_Area.csv")))

# Same prism file for the normalized mito Ferets
mito_Feret_norm_trans = prism_format(data=mito_norm_den_animal_df.reset_index(), data_col="Mito_len_per_den_len", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_Norm_Mito_Feret.csv")))

# Similar dataframe for the total mito dendrite level data

total_mito_area_trans = prism_format(data=den_with_mitos_df_animal.reset_index(), data_col="Mito_total_area", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_total_Mito_Area_per_den.csv")))

total_mito_len_trans = prism_format(data=den_with_mitos_df_animal.reset_index(), data_col="Mito_total_len", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_total_Mito_len_per_den.csv")))

# Get the same data as above but normalized by dendrite area or dendrite length

total_mito_area_norm_trans = prism_format(data=den_with_mitos_df_animal.reset_index(), data_col="Total_mito_area_per_den_area", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_total_Mito_Area_per_den_area.csv")))

total_mito_len_norm_trans = prism_format(data=den_with_mitos_df_animal.reset_index(), data_col="Total_mito_len_per_den_len", col="Layer", row=[comparison, "Animal"], file_path=os.path.join(prism_path, str("Prism_Trans_Med_total_Mito_len_per_den_len.csv")))

# Also save the individual mito level data for area and ferets (for plotting)

mito_indiv_area = prism_format(data = mitos_with_den_df.reset_index(), data_col = "Area_um_sq", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_Mito_area_um2.csv")))

mito_indiv_ferets = prism_format(data = mitos_with_den_df.reset_index(), data_col = "Feret_diam_um", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_Mito_Ferets.csv")))

# Let's also make a combined CA1 / CA2 table with the animal medians
# This should have SR of each animal then SLM of each animal as columns, subregion as rows

# Combined median mito area per animal
med_mito_area_comb = prism_format(data=mito_norm_den_animal_df.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Area_um_sq", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_Mito_Area.csv")))

# Combined median mito Ferets per animal
med_mito_Ferets_comb = prism_format(data=mito_norm_den_animal_df.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Feret_diam_um", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_Mito_Ferets.csv")))

# Combined median mito area per animal normalized to the dendrite area
med_norm_mito_area_comb = prism_format(data=mito_norm_den_animal_df.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Mito_area_per_den_area", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_Norm_Mito_Area.csv")))

# Combined median mito Ferets per animal normalized to the dendrite length
med_norm_mito_Ferets_comb = prism_format(data=mito_norm_den_animal_df.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Mito_len_per_den_len", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_Norm_Mito_Ferets.csv")))

# Combined total mito area per dendrite
total_mito_area_comb = prism_format(data=den_with_mitos_df_animal.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Mito_total_area", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_total_Mito_Area_per_den.csv")))

# Combined total mito length per dendrite
total_mito_area_comb = prism_format(data=den_with_mitos_df_animal.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Mito_total_len", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_total_Mito_Length_per_den.csv")))

# Combined total mito area per dendrite
total_mito_area_norm_den_comb = prism_format(data=den_with_mitos_df_animal.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Total_mito_area_per_den_area", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_total_Mito_Area_per_den_area.csv")))

# Combined total mito length per dendrite
total_mito_len_norm_den_comb = prism_format(data=den_with_mitos_df_animal.reset_index().sort_values(by=["Layer", "Animal"]), data_col="Total_mito_len_per_den_len", col=["Layer", "Animal"], row=comparison, file_path=os.path.join(prism_path, str("Prism_Comb_Med_total_Mito_Length_per_den_len.csv")))


# individual mito level data normalized by the dendrites

mito_indiv_norm_area = prism_format(data = mitos_with_den_df.reset_index(), data_col = "Mito_area_per_den_area", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_Mito_norm_area_um2.csv")))

mito_indiv_norm_ferets = prism_format(data = mitos_with_den_df.reset_index(), data_col = "Mito_len_per_den_len", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_Mito_norm_Ferets.csv")))

# Total mito area and total mito length of individual dendrites

total_mito_area_indiv_den = prism_format(data = den_with_mitos_df.reset_index(), data_col = "Mito_total_area", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_den_total_mito_area.csv")))

total_mito_len_indiv_den = prism_format(data = den_with_mitos_df.reset_index(), data_col = "Mito_total_len", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_den_total_mito_len.csv")))

# Total mito area by dendrite area and total mito length by length of individual dendrites

total_mito_area_by_den_area = prism_format(data = den_with_mitos_df.reset_index(), data_col = "Total_mito_area_per_den_area", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_den_total_mito_area_by_den_area.csv")))

total_mito_len_by_den_length = prism_format(data = den_with_mitos_df.reset_index(), data_col = "Total_mito_len_per_den_len", col=[comparison, "Layer"], file_path = os.path.join(prism_path, str("Prism_indiv_den_total_mito_len_by_den_len.csv")))


# %%% Let's take a look at the distributions for mitochondrial area in the CTL and the KO. Plot histograms and look at the summary stats.

if comparison == "Subregion":

    # Separate the data by subregion
    CA1_data = mito_data[mito_data["Subregion"] == "CA1"]
    CA2_data = mito_data[mito_data["Subregion"] == "CA2"]

    CA1_data_sub = mitos_with_den_df[mitos_with_den_df["Subregion"] == "CA1"]
    CA2_data_sub = mitos_with_den_df[mitos_with_den_df["Subregion"] == "CA2"]

# plot histograms with the get_stats_summary() custom function
# This will return summary stats and whether the data is normal

# plot the CA1 histogram by layer of the original mito data
stats_CA1_full, norm_CA1_full = get_stats_summary(data=CA1_data, data_col="Area_um_sq", x_label="Area (um2)", save_dir=stats_path, group_name="CA1_Mitos_Full", binsize=0.02, hist_comp="proportion", figsize=(8, 5), xlim=[0, 0.8], ylim=[0, 0.15])

# similar histogram for CA2 for the original mito data
stats_CA2_full, norm_CA2_full = get_stats_summary(data=CA2_data, data_col="Area_um_sq", x_label="Area (um2)", save_dir=stats_path, group_name="CA2_Mitos_Full", binsize=0.02, hist_comp="proportion", figsize=(8, 5), xlim=[0, 0.8], ylim=[0, 0.15])

# Get a similar histogram for the subsetted mito population with dendrites
stats_CA1_sub, norm_CA1_sub = get_stats_summary(data=CA1_data_sub, data_col="Area_um_sq", x_label="Area (um2)", save_dir=stats_path, group_name="CA1_Mitos_Sub", binsize=0.02, hist_comp="proportion", figsize=(8, 5), xlim=[0, 0.8], ylim=[0, 0.15])

# Get a similar histogram for the subsetted mito population with dendrites
stats_CA2_sub, norm_CA2_sub = get_stats_summary(data=CA2_data_sub, data_col="Area_um_sq", x_label="Area (um2)", save_dir=stats_path, group_name="CA2_Mitos_Sub", binsize=0.02, hist_comp="proportion", figsize=(8, 5), xlim=[0, 0.8], ylim=[0, 0.15])

# Histograms for Mito Ferets

# Start with the full dataset, histogram of CA1 across layers
stats_CA1_Feret, norm_CA1_feret = get_stats_summary(data=CA1_data, data_col="Feret_diam_um", x_label="Area (um2)", save_dir=stats_path, group_name="CA1_Mitos_Full", binsize=0.05, hist_comp="proportion", figsize=(8, 5), xlim=[0, 2], ylim=[0, 0.1])

# similar histogram for CA2 for the original mito data
stats_CA2_Feret, norm_CA2_feret = get_stats_summary(data=CA2_data, data_col="Feret_diam_um", x_label="Area (um2)", save_dir=stats_path, group_name="CA2_Mitos_Full", binsize=0.05, hist_comp="proportion", figsize=(8, 5), xlim=[0, 2], ylim=[0, 0.1])

# plot the CA1 histogram by layer for the subset of mitos in dendrites
stats_CA1_Feret_sub, norm_CA1_feret_sub = get_stats_summary(data=CA1_data_sub, data_col="Feret_diam_um", x_label="Feret's Diameter (um)", save_dir=stats_path, group_name="CA1_Mitos_Sub", binsize=0.05, hist_comp="proportion", figsize=(8, 5), xlim=[0, 2], ylim=[0, 0.1])

# plot the CA1 histogram by layer for the subset of mitos in dendrites
stats_CA2_Feret_sub, norm_CA2_feret_sub = get_stats_summary(data=CA2_data_sub, data_col="Feret_diam_um", x_label="Feret's Diameter (um)", save_dir=stats_path, group_name="CA2_Mitos_Sub", binsize=0.05, hist_comp="proportion", figsize=(8, 5), xlim=[0, 2], ylim=[0, 0.1])

# Histograms for dendrite area and dendrite length, comparing subregions for SR and SLM

# filter the subsetted dendrite data to get just SR and SLM

den_data_SR = den_with_mitos_df[den_with_mitos_df["Layer"] == "Proximal"]
den_data_SLM = den_with_mitos_df[den_with_mitos_df["Layer"] == "Distal"]

# Dendrite length in SR (CA1 vs CA2)
stats_SR_den_length = get_stats_summary(data=den_data_SR, data_col="Den_Length_um", x_label="Length (um)", save_dir=stats_path, hue_col = "Subregion", colors = ["dark rose", "ocean"], group_name="SR_Dendrites", binsize=0.1, hist_comp="proportion", figsize=(8, 5), xlim=[0, 7], ylim=[0, 0.1])

# Dendrite length in SLM (CA1 vs CA2)
stats_SLM_den_length = get_stats_summary(data=den_data_SLM, data_col="Den_Length_um", x_label="Length (um)", save_dir=stats_path, hue_col = "Subregion", colors = ["dark rose", "ocean"], group_name="SLM_Dendrites", binsize=0.1, hist_comp="proportion", figsize=(8, 5), xlim=[0, 7], ylim=[0, 0.1])

# Dendrite area in SR (CA1 vs CA2)
stats_SR_den_area = get_stats_summary(data=den_data_SR, data_col="Den_Area_um_sq", x_label="Area (um2)", save_dir=stats_path, hue_col = "Subregion", colors = ["dark rose", "ocean"], group_name="SR_Dendrites", binsize=0.1, hist_comp="proportion", figsize=(8, 5), xlim=[0, 5], ylim=[0, 0.15])

# Dendrite area in SLM (CA1 vs CA2)
stats_SLM_den_area = get_stats_summary(data=den_data_SLM, data_col="Den_Area_um_sq", x_label="Area (um2)", save_dir=stats_path, hue_col = "Subregion", colors = ["dark rose", "ocean"], group_name="SLM_Dendrites", binsize=0.1, hist_comp="proportion", figsize=(8, 5), xlim=[0, 5], ylim=[0, 0.15])


# For mito area, get the summary stats for each layer from each animal
area_summary_table = mitos_with_den_df[["Area_um_sq", "Animal", "Layer"]].groupby(["Animal", "Layer"]).describe()

# save file to CSV
area_summary_table.to_csv(os.path.join(stats_path, str("mito_area_stats_by_animal_sub_den.csv")))

# For mito Ferets, get the summary stats for each layer from each animal
Feret_summary_table = mitos_with_den_df[["Feret_diam_um", "Animal", "Layer"]].groupby(["Animal", "Layer"]).describe()

# save file to CSV
Feret_summary_table.to_csv(os.path.join(stats_path, str("mito_diam_stats_by_animal_sub_den.csv")))


# %% Get some basic metadata information about this code and the analysis to save as a txt file with the data.

# Total number of tiles analyzed in each region
basal_CA1_count = list(mito_avgs["Layer"].where(mito_avgs["Subregion"] == "CA1").dropna()).count("Basal")
prox_CA1_count = list(mito_avgs["Layer"].where(mito_avgs["Subregion"] == "CA1").dropna()).count("Proximal")
dist_CA1_count = list(mito_avgs["Layer"].where(mito_avgs["Subregion"] == "CA1").dropna()).count("Distal")

basal_CA2_count = list(mito_avgs["Layer"].where(mito_avgs["Subregion"] == "CA2").dropna()).count("Basal")
prox_CA2_count = list(mito_avgs["Layer"].where(mito_avgs["Subregion"] == "CA2").dropna()).count("Proximal")
dist_CA2_count = list(mito_avgs["Layer"].where(mito_avgs["Subregion"] == "CA2").dropna()).count("Distal")


total_analyzed = "\nCA1:  " + str(basal_CA1_count) + " Basal; " + str(prox_CA1_count) + " Proximal; " + str(dist_CA1_count) + " Distal" + "\n\nCA2:  " + str(basal_CA2_count) + " Basal; " + str(prox_CA2_count) + " Proximal; " + str(dist_CA2_count) + " Distal"


# For reference, we'll add the name of the code used to analyse the data (this code)
# Note: Keep this updated if you update the name or date on the code
code_name = inspect.getfile(inspect.currentframe()).split("\\")[-1]

# Let's get how many animals and how many stubs are in the dataset
animal_list = list(set(mito_avgs["Animal"]))  # list of animal IDs in the data
num_animals = len(animal_list)  # number of animals analyzed
animals = ", ".join(animal_list)  # animal IDs formatted for text output

groups = mito_data[comparison].unique()
groups_str = ", ".join(groups)

# Get how many stubs per animal
num_stub_per_animal = [len(list(set(mito_avgs["Stub"].where(mito_avgs["Animal"] == x).dropna()))) for x in animal_list]

if len(list(set(num_stub_per_animal))) > 1:  # if animals have different numbers of stubs
    num_stub = str(min(num_stub_per_animal)) + "-" + str(max(num_stub_per_animal)) + " stubs / mouse"
else:  # if all animals have same number of stubs
    num_stub = str(num_stub_per_animal[1]) + " stubs / mouse"

# List of animals and stubs in the data
stub_list = mito_data[["Animal", "Stub"]].drop_duplicates().reset_index(drop=True)

# Once we get WT and KO data we may want to separate out how many WT and KO mice

analyzed = str(num_animals) + " mice analyzed (" + str(num_stub) + ")"

# explain any processing of the data by this code
processing = "\n" + str(mitos_removed) + " mitochondria larger than 2.1 um2 or smaller than 0.01 um2, " + str(den_removed_spines) + " dendritic spine heads (< 0.4 um2 area), and " + str(den_removed_primaries) + " primary dendrites (> 1.9 um width) were removed.\nAll areas and diameters converted from pixels to micron squared. \n\n For mitochondria within segmented dendrites, mitochondria area was normalized to dendrite area and mitochondria Feret's diameter was normalized todendrite Feret's diameter of each mitochondria's parent dendrite."

# Put together the metadata information for the text file
meta_summary = {"Date": date.today().isoformat(), "Dataset": dataset, "Animals": animals, "Analyzed": analyzed, "Groups": groups_str, "Sampling": sampling,
                "Tile Size": tile_size, "N Tiles": total_analyzed, "Processing": processing, "Data CSV": data_loc.split("\\")[-1], "Biodock AI": AI_name, "Analysis Code": code_name}

# String of how many tiles were excluded
excluded = str(num_flag) + " tiles\n" + ", ".join(flagged.index)

if num_flag >= 1:  # If any tiles were excluded due to low counts
    meta_summary["Excluded"] = excluded  # add excluded column to the metadata

if len(notes) > 1:  # if there are notes
    meta_summary["Notes"] = notes  # add notes column to the metadata

# write the metadata dictionary to a text file
meta_file = open(os.path.join(dirName, 'analysis_summary.txt'), 'wt')
for line in meta_summary:
    meta_file.write(str(line) + ":  " + str(meta_summary[line]) + "\n\n")
meta_file.close()

curr_dir = os.path.join(os.getcwd(), dirName)

print(f"\nData has been saved in the directory: \n '{curr_dir}'")

print("\nCode is finished! Go take a nap.")


