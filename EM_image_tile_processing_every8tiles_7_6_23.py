# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 2023 2:11pm
Last updated on Fri Feb 24 2023

Quick code to sample and process SEM tiles (6144 x 4969 pixels) for analysis. The tiles should already be cropped to the desired size and inverted (if needed). For images that were scanned into quadrants, you should run the code "EM_image_tiling_14x20_from_quad.py" first to invert and break the quadrants up into the smaller tiles. This code will sample every 8th tile (by file name), set the proper scale, gaussian blur, rename and save them in a separate directory for analysis.

INPUT: A directory with a folder for each tile set you want to process for analysis. The tile set should contain all tiles from the same image and should be labeled with the tile number first (ex: "Tile_001-001..." to "Tile_020-014..."). Please name the folder with any identifying information you wish to be added to the image file names (ie. Animal, stub, dendritic layer, etc.) in the output. Alternatively, you can specify the identifying information to be used instead of the folder names, if you're just running a single group. NOTE: Make sure to include "Basal", "Proximal" or "Distal" in the folder name so the code knows how to sort it in the output.

If you specify the animal name as the very beginning of the folder name (eg. MCC210), the code will also sort the tiles by animal in the output. For now, it only recognizes MCC mice since we are focused on the MCU KO, but more genotypes can be added as necessary.

OUTPUT: The output will be an analysis directory "For Analysis" containing the sampled, processed tiles ready for analysis. There will be a folder for each dendritic layer (basal, proximal and distal) and then a subfolder for each animal (if provided). If any folders in the input directory do not specify one of those three dendritic layers in their name, the sampled tiles from that tile set will be placed in a separate folder labeled as "Unsorted" to be manually sorted into the appropriate group(s). 

If the analysis folder already exists, sample tiles will be added to the exising directory. A warning will be thrown if there is already a file with the same filename in the analysis folder, and the file will be skipped.

@author: pannoni
"""

# import needed packages
import easygui
import os
from datetime import date
from PIL import Image, ImageFilter
import tifffile

#GUI to have user select the main directory with the EM images. Make sure your directory is organized as described above, with a folder for each tile set to be sampled.
main_folder = easygui.diropenbox(msg='Please select the main directory with your SBFEM images.', title='Select Main Folder', default='*')

#!!! if desired, you can specify any identifying information for the selected data (ie. Animal, stub name, cre, etc). This will be added to the image names in the output directory. Note that this identifying information will be applied to the entire dataset in your main folder. This can be useful if you're only running one group.
# By default, this code will use the folder name as the idenfifying information, unless specified below.
ID = ""

# List the folders within the main directory -> each tile set should have a folder
image_sets = [f for f in os.listdir(main_folder) if os.path.isdir(os.path.join(main_folder, f)) and f != "For Analysis"] # only get folders, ignore the "For Analysis" folder if it exists
# Number of image sets we have to process
num_img_sets = len(image_sets)

animal_list = [] # empty list to put the animal IDs

# Start counters to track how many sampled tiles we have for each dendritic layer
basal_count = 0
prox_count = 0
dist_count = 0
total_N_sampled = 0
tiles_skipped = 0

# make an output directory in the same location as the main folder
out_dir = os.path.join(main_folder, "For Analysis")

try:
    os.mkdir(out_dir) # Create output directory
    print("\nDirectory '" + str(out_dir) +  "' Created.\n") 
except FileExistsError:
    print("\nDirectory '" + str(out_dir) +  "' already exists. \nThe sampled tiles will be added to the existing directory.\n")
    
# Also make folders for each of the three dendritic layers
try:
    os.mkdir(os.path.join(out_dir, "Basal"))
    os.mkdir(os.path.join(out_dir, "Proximal"))
    os.mkdir(os.path.join(out_dir, "Distal"))
except FileExistsError:
    pass # skip if these folders already exist

# loop through each image set to process the image tiles
for s in image_sets:
    
    # group into basal, proximal or distal based on the name
    if "Basal" in s:
        layer = "Basal"
    elif "Proximal" in s:
        layer = "Proximal"
    elif "Distal" in s:
        layer = "Distal"
    else:
        layer = "Unsorted"
        try: # If needed, create a folder for the unsorted tiles in the output directory
            os.mkdir(os.path.join(out_dir, "Unsorted"))
        except FileExistsError:
            pass
    
    # get the full file path for the current set of tiles
    image_path = os.path.join(main_folder, s)
    # get the image files in the folder (tifs only)
    tile_images = [t for t in os.listdir(image_path) if t.endswith('.tif')]
    
    # Check the number of tiles and throw some warnings
    if len(tile_images) == 1:
        print("\nWARNING: The folder '" + str(s) + "' only has one tile image in it. Something has likely gone wrong.")
        
    elif len(tile_images) == 0:
        print("\nWARNING: The folder '" + str(s) + "' does not contain any .tif images. This folder will be skipped.")
        pass # go to next image set
    
    print("\nProcessing '" + str(s) + "'...")
    
    # First let's get the appropriate location to save the tiles
    # Get the animal ID from the folder name, if provided (should be the first 6 characters)
    if "MCC" in s[0:3]:
        animal = s[0:6]
        # path where the tiles will be saved
        out_path = os.path.join(out_dir, layer, animal)
        if animal not in animal_list: # add animal to list if it's not already there
            animal_list.append(animal)
            
        try: # Create a folder for the current animal
            os.mkdir(os.path.join(out_path))
        except FileExistsError:
            pass
    
    else: # if the animal ID is not in the folder name
        animal = []
        # path where the tiles will be saved
        out_path = os.path.join(out_dir, layer)
    
    # Loop through each tile to pull out every 10th tile
    for tile in tile_images:
        
        # only pull the tiles we want to sample (can be changed)
        sample = ["Tile_001-008", "Tile_002-002", "Tile_002-010", "Tile_003-004", "Tile_003-012", "Tile_004-006", "Tile_004-014", "Tile_005-008", "Tile_006-002", "Tile_006-010", "Tile_007-004", "Tile_007-012", "Tile_008-006", "Tile_008-014", "Tile_009-008", "Tile_010-002", "Tile_010-010", "Tile_011-004", "Tile_011-012", "Tile_012-006", "Tile_012-014", "Tile_013-008", "Tile_014-002", "Tile_014-010", "Tile_015-004", "Tile_015-012", "Tile_016-006", "Tile_016-014", "Tile_017-008", "Tile_018-002", "Tile_018-010", "Tile_019-004", "Tile_019-012" , "Tile_020-006", "Tile_020-014"]
        
        # Get the current tile number
        tile_num = tile[:12]

        if tile_num in sample:
            
            # Create a new image file name
            if len(ID) > 1:  # if there is user provided ID
                # create a new image file name using the provided ID
                tile_name = str(ID) + "_" + str(tile_num) + ".tif"
            else:   
                # Create a new image file name using the folder name
                tile_name = str(s) + "_" + str(tile_num) + ".tif"
                
            # Before doing anything else, check whether this file name already exists in the output folder
            if tile_name not in os.listdir(out_path):
                
                # open the image file to get the size
                Image.MAX_IMAGE_PIXELS = None # gets past a limit on the pixel size of the images
                with Image.open(os.path.join(image_path, tile)) as image: # open tile image           
                    width, height = image.size # get image dimensions
                    
                    # throw a warning if the image dimensions are not 6144 x 4969 pixels
                    if width != 6144 or height != 4096:
                        print("Uh oh. " + str(tile) + " appears to be the wrong size. Tiles should be 6144 x 4969. Please double check the dimensions of this tile.")
                        
                # Our scale is 2nm per pixel, we want to gaussian blur with a radius of 1 pixel (2nm)
                twonm = 1
                    # open the image file to gaussian blur it
                with Image.open(os.path.join(image_path, tile)) as image:
                    # get the image mode. For most of the images we're analyzing they will be 16 bit unsigned ("I;16")
                    image_mode = image.mode
                    
                    if image_mode != "L":
                        # The filter doesn't work on 16 bit images, so line below converts the image to "L" which is 8 bit B&W pixels
                        # Note that converting the entire image doesn't work, but converting it point by point seems to work
                        image_gauss = image.point(lambda i:i*(1./256)).convert('L').filter(ImageFilter.GaussianBlur(radius = twonm))

                    elif image_mode == "L": # if image is already 8 bit, no need to convert
                        image_gauss = image.filter(ImageFilter.GaussianBlur(radius = twonm))
                
                    # save the image in the correct output folder
                    save_path = os.path.join(out_path, tile_name)
                    
                    # save the image in the folder
                    tifffile.imsave(
                        save_path,
                        image_gauss)
                        # resolution=(500, 500, 'MICROMETER')) # set resolution of image to 500 px / micron, same as 2nm / pixel
                    
                # Counter to track how many sampled tiles we have for each layer
                if layer == "Basal":
                    basal_count = basal_count + 1
                elif layer == "Proximal":
                    prox_count = prox_count + 1
                elif layer == "Distal":
                    dist_count = dist_count + 1
                    
                total_N_sampled = total_N_sampled + 1
                    
            else: # if tile is already in the output folder, skip it
                print("\n'" + str(tile_name) + "' is already in the 'For Analysis' folder. Tile has been skipped.")
                    
                tiles_skipped = tiles_skipped + 1
                
#%% get some metadata information to save in a text file with the images
# include date run, location of input directory, processing (gaussian blur), scale set (2nm), N of tile sets processed, list of animals and N of stubs per animal

# get the tile size for the metadata
tile_size = "6144 x 4969 pixels"

# tile sampling
sampling = "Every 8 tiles"

total_sampled = str(basal_count) + " Basal; " + str(prox_count) + " Proximal; " + str(dist_count) + " Distal (Total = " + str(total_N_sampled) + ")"

# list any processing steps
filt = "Gaussian blur, radius 2nm"

# scale set on processed images
scale = "500 pixels / micron"

#!!! Which dataset are the images from? (For now this will need to be set manually below)
dataset = "WT and KO CA1"

# For reference, we'll add the name of the code used to analyse the data (this code)
code_name = "EM_image_tile_processing_every8tiles_7_6_23.py"

animals = ", ".join(animal_list)

# Put together the metadata information for the text file
meta_summary = {"Date": date.today().isoformat(), "Dataset":dataset, "Animals": animals, "Images Processed":num_img_sets, "Sampling": sampling, "Tile Size": tile_size, "Sampled Tiles": total_sampled, "Processing": filt, "Scale": scale, "Image Loc": main_folder, "Analyzed with": code_name}

#!!! Space to put notes about the dataset or any changes made since a previous run.
notes = "Adding more CA1 sections for the CA1 mitochondria analysis."

if len(notes) > 1: # if there are notes
    meta_summary["Notes"] = notes # add notes column to the metadata

# write the metadata dictionary to a text file
meta_file = open(os.path.join(out_dir, 'processing_summary.txt'), 'wt')
for line in meta_summary:
    meta_file.write(str(line) + ":  " + str(meta_summary[line]) + "\n\n")
meta_file.close()

# print out a summary to the console

print("\nCode has finished running! Your processed tile images are located in: \n" + str(out_dir))

print("\n---------------------------------------------\n               CODE SUMMARY\n---------------------------------------------\n")

print("A total of " + str(total_N_sampled) + " tiles were sampled from " + str(len(animal_list)) + " mice in " + str(dataset))

print("(" + str(basal_count) + " Basal; " + str(prox_count) + " Proximal; " + str(dist_count) + " Distal)")

if tiles_skipped != 0: # if some tiles were skipped, print to the console
    print("\n" + str(tiles_skipped) + " tiles were skipped because they are already in 'For Analysis'")
                
            

        
