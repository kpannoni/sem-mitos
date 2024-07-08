# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:58:19 2019

Last updated: 3/12/2024

UPDATE 3/6: Added the option to chose the number to sample at the tile and mito level and a test for normality to the bootstrap function, so that a test of different sampling can be done.

UUPDATE 3/12: Split up the code to run the bootstraps separately and to be able to calculate the probabilities or plot the probability plots without having to rerun the entire bootstrap. Created a custom summary function to plot the distributions of the bootstrap population and print out a nicely formatted summary table.

UPDATE 7/1: Cleaned up to be uploaded on GitHub.


@author: pannoni
"""

# import necessary packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
from datetime import date
import inspect
import os
import sys

# First, import the data and set up the bootstrap

#!!! Get the location of the data. This is where we will create the folder for the bootstrap.
# You will need to update this if the data you want to run the bootstrap on is somewhere else.
# loc = "C:\\Users\\Pannoni\\Documents\\Lab Stuff\\EM_mito_bootstrap_V6_9_19"
loc = "C:\\Users\\Pannoni\\Documents\\Lab Stuff\\Python_Scripts\\local_github_repo\\sem-mitos\\csv_output"
    
#!!! set the number of bootstrap repetitions
n_boot = 5

# Import the mito dataframes for the EM data, which should be located in the folder that "loc" points to
#!!! Note that if you ran the code "EM_image_mito_analysis_biodock_6_20_2024_cleaned", this is what the files should be named
# If they are named something else, change the file name here below.
mito_data = pd.read_csv(os.path.join(loc, "SEM_indiv_mito_data.csv"), index_col=0)
tile_data = pd.read_csv(os.path.join(loc, "SEM_mito_tile_avgs.csv"), index_col=0)

# for mito data, we will need to add a column for section ID that combines animal and stub (if there isn't one already), because some stub names are repeats across different animals
if "Section_ID" not in mito_data.columns:
    mito_data["Section_ID"] = mito_data["Animal"] + "_" + mito_data["Stub"]

#%% Define the get_boot_sample function, which will run the bootstrap and return the means of each repetition in a variable.

def get_boot_sample( orig_df, data_col, lev1_col, lev2_col, lev3_col = [], group_col = [], sort_col = [], nboot=10000, show_resample = False, direct = os.getcwd(), tile_samp = 20, mito_samp = 5):
    '''
    This function performs a hierarchical bootstrap on multiple levels of data in "orig_df". 
    It randomly samples the data at each level indicated (with replacement) to get the mean of the resampled population.
    For the purposes of the EM mitochondrial data, the levels will be Animal, Stub and then Tile.
    The inputs:
        -- orig_df: Dataframe with the data for the bootstrap. Should be organized in columns: One column with the data to be bootstrapped, and each level of the hierarchical data should have unique identifiers in a separate column associated with the data in the data column.
        -- data_col: Name of the data column for the bootstrap. Required.
        -- lev1_col: Name of the column with unique identifiers for the first level (highest level, usually animal). Required.
        -- lev2_col: Name of column with unique identifiers for the second level (often section, stub or hippocampus). Required.
        -- lev3_col: Name of the column with unique identifiers for the third level (optional). Default is 0, no resampling at a third level. Use this variable if you want to sampling the EM data at the individual mitochondria level instead of tile. In that case, you should set this to the column with the tile names.
        -- group_col: Name of the column for the grouping of the first level (optional). Usually this column will be used to group animals by genotype or treatment.
        -- sort_col: Name of a column used to sort the data into categories (optional) after resampling level 2. Any stats (eg. mean) will be applied within these categories. Typically, this would be grouping the resampled sections or hippocampi (level 2) into layers or subregions before pooling data 
        -- nboot: Number of times to perform the bootstrap. The default for nboot is 1000.
        -- show_resample: Optional boolean to print out the resampled levels along the way. Default = False. I recommend only using this option if you're running a small nboot, otherwise it's too much to print to the console.
        -- direct: the path location to save the bootstrap results. A folder will be created in this location for the bootstrap. Default is to create the folder in the current working directory.
        -- tile_samp = Number of random samples to pick at the final data level (if lev3_col = []) or at the level of lev3_col if you've set a third level to resample. In which case, we will go down one more level to resample the final data. Normally, tile_samp will be used for sampling the data at the "Tile" level. Default is 20.
        -- mito_samp = Number of random samples to pick at the final data level if resampling at a third level. Typically, this will be sampling the data at the level of individual mitochondria. Default is 5. This number will only be used if lev3_col is set.
        
    The output:
        -- a variable (bootstats) with the mean for each randomly sampled population of your data. 
        If there is no sort_col given, the data will not be sorted and 'bootstats' will have nboot number of means (one for each repetition). This is the simplest condition.
        If you provide a "sort_col", the data will be sorted by the categories in your column and means will be calculated for each category. In this case, bootstats will have a column for each category in sort_col, and a row for each repetition of the bootstrap.
        If you provide a "group_col", bootstats will be a dictionary containing a dataframe for each group in your grouping column with the group names as the keys. Each dataframe will be formatted as described above, including columns for the categories if "sort_col" is provided.
        -- bootstats will be saved as a separate CSV file for each group in group_col, or just a single file if not grouping.
        -- A csv file with the descriptive stats for the population of bootstats medians will also be saved for each group.
        -- The code will also save a text file with some metadata information about the bootstrap (ie. the sampling scheme, the grouping, n_boot, etc)
        
    --------------------------------------
    
    Bootstrap update:
        
        The code will be used to group the EM mito data by 2 genotypes (cre+ and cre-), resample the animals and then resample the stubs for each sampled animal. The code will then group the tiles in each sampled stub by layer (Basal, Prox, Distal) and sample 20 tiles from each layer. If no third level is provided, the code will pool the sampled tiles from each stub and animal together (n = 180 tiles per group and layer) to get a single mean for each layer in each group. If a third level is specified with "lev3_col", the code will sample 5 mitochondria from each of the 20 sampled tiles, pool the sampled mitochondria together by layer, which should also result in a single mean for each layer in each group (n=900 mitos per group and layer). This process will be repeated n=nboot times. 
        
        Keep in mind that the final level the data is sampled at depends on your input dataframe. For example, if you want your final means to be calculated at the tile level, input a dataframe with the tile level data (each row = a tile). If you want your final means to be at the individual mitochondria level, you should input a dataframe with the individual mito data (each row = a mitochondria) and you would likely want to add the third level of resampling at the tile level. 
        
        The grouping by genotype will also be optional, if you don't provide a grouping column for 'group_col' the the code will skip the first grouping and treat all the animals as a single group. The data will still be grouped by the second grouping (layer) before sampling the tiles, if you provide a "sort_col".
        
        The resulting output should be a dictionary with two entries (one for each genotype) with nboot means organized into columns for layer. There should be nboot rows. This will also be saved as .csv files in the folder specified.
        
        UPDATE 3/6: Added the option to chose the number to sample at the tile and mito level, updated the metadata output. Also, save a file with the summary stats for the resampled median population.
        
        UPDATE 3/12: The function will now also plot a histogram of the resampled median populations for each group
    
    '''
    print("\n\nBEGIN BOOTSTRAP OF " + data_col.upper() + "\n----------------------------------------------------------------\n\nOriginal data:\n", orig_df)
    
    # Record start time
    start_time = time.time()
    
    # create the directory to save into
    dir_name = direct
    save_dir = os.path.join(dir_name, "Bootstrap_median_EM_mitos_" + str(n_boot) + "_tile_samp_" + str(tile_samp) + "_mito_samp_" + str(mito_samp))
     
    try:
        # Create target Directory in the location
        os.mkdir(save_dir)
        print("\nDirectory " , save_dir ,  " Created.\n") 
    except FileExistsError:
        pass
    
    # Grab all the data from the data column
    all_data= orig_df[data_col]
    
    # Check whether there is a grouping, get the number of groups
    if len(group_col) > 0: #if there is grouping
        group = orig_df[group_col] # column with the grouping categories
        group_list = list(set(group))
        try:
            group_order = ["MCC Cre -", "MCC Cre +"] #!!! sort so the CTL is first. 
            # If you have different labels for cKO and CTL you should put them above and make sure the CTL is first in the list
            group_list = sorted(group_list, key=lambda x: group_order.index(x)) # to put the layers in order
        except ValueError:
            pass
        group_num = len(group_list) # number of groups in the data
        
        # print to the console what the bootstrap will be doing
        if len(lev3_col) > 0 and len(sort_col) > 0: # if there's a third level for resampling
            print('\nThe data will be grouped by ' + group_col + ', re-sampled at the level of ' + lev1_col + ', ' + lev2_col + ', and ' + lev3_col + ', and sorted by ' + sort_col + ' to get the resampled medians.\n\nRunning boostrap on "' + data_col + '" with an nboot of ' + str(nboot) + ' ...\n')
        elif len(sort_col) > 0: # only 2 levels
            print('\nThe data will be grouped by ' + group_col + ', re-sampled at the level of ' + lev1_col + ' and ' + lev2_col + ', and sorted by ' + sort_col + ' to get the resampled medians.\n\nRunning boostrap on "' + data_col + '" with an nboot of ' + str(nboot) + ' ...\n')
        else: # no sorting
            print('\nThe data will be grouped by ' + group_col + ' and re-sampled at the level of ' + lev1_col + ' and ' + lev2_col + ' to get the resampled medians.\n\nRunning boostrap on "' + data_col + '" with an nboot of ' + str(nboot) + ' ...\n')
            
    else: # if no grouping
        group_num = 1 
        group_list = ["None"]
        
        # print to the console what the bootstrap will be doing
        if len(lev3_col) > 0 and len(sort_col) > 0: # if there's a third level for resampling
            print('\nThe data will be re-sampled at the level of ' + lev1_col + ', ' + lev2_col + ', and ' + lev3_col + ', and sorted by ' + sort_col + ' to get the resampled medians.\n\nRunning boostrap on "' + data_col + '" with an nboot of ' + str(nboot) + ' ...\n')
        elif len(sort_col) > 0: # only 2 levels
            print('\nThe data will be re-sampled at the level of ' + lev1_col + ' and ' + lev2_col + ', and sorted by ' + sort_col + ' to get the resampled medians.\n\nRunning boostrap on "' + data_col + '" with an nboot of ' + str(nboot) + ' ...\n')
        else: # no sorting
            print('\nThe data will be re-sampled at the level of ' + lev1_col + ' and ' + lev2_col + ' to get the resampled medians.\n\nRunning boostrap on "' + data_col + '" with an nboot of ' + str(nboot) + ' ...\n')
    
    # check whether there is a sort category (ie. Layer) and define the output arrays for "bootstats"
    if len(sort_col) > 0: 
        sort = orig_df[sort_col] # column with the sort categories
        sort_list = list(set(sort))
        try:
            layer_order = ["Basal", "Proximal", "Distal"]
            sort_list = sorted(sort_list, key=lambda x: layer_order.index(x)) # to put the layers in order
        except ValueError:
            pass
        sort_num = len(sort_list) # number of categories to sort into
        
        # Define the output variable "bootstrap". Bootstats will be an array of 0s for each group in the group_col (as a list)
        # If no grouping, the list will contain just one array.
        # If sorting, each array will have a column for each sort category (ie. Layer) and nboot rows
        bootstats = [np.zeros((nboot, sort_num))] * group_num
    else: # if no sort column
        bootstats = [np.zeros(nboot)] * group_num # Compiles nboot means together in a single array for each group, no columns
    
    # Now, let's sample the data at the different levels and get the medians. 
        
    for i in np.arange(nboot): # each nboot iteration
        
        group_loc = 0
        # loop through the groups
        for group in group_list:
            
            temp = pd.DataFrame() # temporary variable to hold the current group's data for the current trial
            
            if show_resample: print("\n", group)
            
            # If there is more than one group
            if len(group_list) > 1:
                # get just the current group's data, we want to sample each group separately
                data = orig_df[orig_df[group_col] == group]
                
            elif len(group_list) == 1: # if there's only 1 group
                data = all_data # sample all the data
            
            # Get the number of samples for the first level. In this case, the sample size is equal to the number of elements in each level. For example, if there are 3 in the data, then the program will pick a random sample of 3 animals (with replacement)
            lev1= data[lev1_col]
            lev1_list = list(set(lev1)) #unique list of level 1 identifiers
            sample1_size = len(lev1_list)
            
            # randomly sample the level 1 identifiers, allow repleacement
            rand_lev1 = np.random.choice(lev1_list,sample1_size)
            if show_resample: print("\nResampled Level 1: ", rand_lev1)
            
            for j in rand_lev1:
                # use the random sample of level 1 to pull out the associated level 2 identifiers from the dataframe
                lev2= data[lev2_col]
                lev2_list = list(set(lev2.where(lev1==j).dropna()))
                sample2_size = len(lev2_list) # get the sample size for level 2
                
                # randomly sample the level 2 identifiers with replacement
                rand_lev2 = np.random.choice(lev2_list,sample2_size)
                if show_resample: print("Resampled Level 2: ", rand_lev2)
                
                # pull out the data from the data_col belonging to each element in the randomly sampled level2 identifiers
                for k in rand_lev2: 
                    lev2_data =  data.where(lev2==k).dropna() # subset the data by the current lev2 identifier (likely Tile or Section)
                    
                    if len(sort_col) > 0:
                        temp_sort = pd.DataFrame()  # empty temproary dataframe for the sorted data
                        # sort the data by sort category and randomly sample within each category
                        for cat in sort_list:
                            curr_sort = lev2_data[lev2_data[sort_col] == cat] # get the data for each element in sort_list
                            # E.g. randomly sample 20 tiles from the current region (default)
                            if len(lev3_col) > 0: # if resampling on a third level
                                lev3= curr_sort[lev3_col] # get the level 3 identifiers for the current category
                                lev3_list = list(set(lev3)) # list of unique identifies for Level 3
                                
                                rand_lev3 = np.random.choice(lev3_list,tile_samp) # randomly sample elements (likely Tiles) from Level 3
                                if show_resample: print("Resampled Level 3: ", rand_lev3)
                                
                                for t in rand_lev3: # for each element selected for level 3, resample 5 datapoints (default)
                                    # This is for the individual mitochondria
                                    # get the data for the current element (Tile)
                                    curr_lev3 = curr_sort[curr_sort[lev3_col] == t] # e.g. all the mitos for the current tile
                                    rand_data = np.random.choice(curr_lev3[data_col],mito_samp) # randomly pick 5 datapoints in the data column
                                    
                                    # add the random sample data to "temp_sort" in a separate column (with the category name as column name)
                                    temp_sort[cat] = rand_data
                                    
                            else: # if not resampling at level 3, just get the data from level 2
                                rand_data = np.random.choice(curr_sort[data_col],tile_samp) # randomly select 20 datapoints from data column in the current category
                                # add the random sample data to "temp_sort" in a separate column (with the category name as column name)
                                temp_sort[cat] = rand_data
                                
                        temp = pd.concat([temp, temp_sort], axis = 0).reset_index(drop = True) # now add temp_sort to the main temp variable (in the next row)
                        # in the end, temp should have 180 datapoints if not resampling at level 3, or 900 if you are resampling at level 3
                    else: # if not sorting the data, just add the current data to the temp varaible
                        # get the data from the unsorted data at level 2 or level 3
                        if len(lev3_col) == 0: # no sampling at level 3
                            rand_data = np.random.choice(lev2_data[data_col],tile_samp) # pull the random data from lev 2
                            # add the random data to the temp dataframe
                            temp = pd.concat([temp, rand_data], axis = 0).reset_index(drop = True)
                        else: # if also resampling at level 3
                            lev3 = lev2_data[lev3_col] # get the level 3 identifiers for the current category
                            lev3_list = list(set(lev3)) # list of level 3 identifiers
                            # randomly select 20 elements from level 3
                            rand_lev3 = np.random.choice(lev3_list,tile_samp)
                            
                            for t in rand_lev3: # for each element selected for level 3, resample 5 datapoints
                                lev3_data = lev2_data[lev2_data[lev3_col] == t] # get just the data for the current level 3 identifier
                                # randomly select 5 elements at level 3
                                rand_data = np.random.choice(lev3_data[data_col],mito_samp)
                                temp = pd.concat([temp, rand_data], axis = 0).reset_index(drop = True)
            
            # get the medians(s) for the current repetition and save it to 'bootstats'. If data is sorted into categories, get the median for each category. 
            # If not sorted, just get a single median for the population.
            # Make sure this is saving the data in the right location in bootstats
            if len(sort_col) > 0: 
                bootstats[group_loc] = np.copy(bootstats[group_loc]) # make a copy of the array for the current group
                bootstats[group_loc][i] = np.median(temp, axis=0).copy() # add the means to the appropriate row for each iteration of nboot
            else: bootstats[group_loc][i] = np.median(temp)
            
            group_loc += 1 # keep track of the groups, increase by 1 each iteration. 
            # This will be used to put the group data in the right place in bootstats
            
            
        # track the progress by printing out the current repetition every 100 repetitions. Adjust the denominator in the fraction below based on how often you want to be updated (for nboot = 1000, every 100 is probably good)
        if ((i+1)/100).is_integer():
            print("...", i+1, "...\n")
        
    # return bootstats. Bootstrap has the population means for each category (columns) calculated for each repetition of the bootstrap (rows)
    print('\n...finished!\n')
    if len(sort_col) > 0 and len(group_col) > 0:
        meta = f"\nThe bootstrapped medians for {data_col} were calculated (nboot = {nboot})\n for the groups {group_list}\nThe means are sorted by column in order of {sort_list}"
    elif len(group_col > 0):
        meta = f"\nThe bootstrapped medians for {data_col} were calculated (nboot = {nboot})\n for the groups {group_list}"
    else: # no grouping or sorting
        meta = f"\nThe bootstrapped medians for {data_col} were calculated (nboot = {nboot})\n"
    print(meta)
    if len(lev3_col) > 0:
         resampling = f"Data was re-sampled at the level of {lev1_col} ({sample1_size}), {lev2_col} ({sample2_size}), and {lev3_col} ({tile_samp}).\n Mitos sampled: {mito_samp} \n"
    else:
         resampling = f"Data was re-sampled at the level of {lev1_col} ({sample1_size}) and {lev2_col} ({sample2_size}).\n Tiles sampled: {tile_samp}\n"
    
    # get the current date and time
    today = date.today().isoformat()
    
    # add the date to the summary and save as a .txt file in the save directory
    # get the name of the current code file
    code_name = inspect.getfile(inspect.currentframe())
    meta = resampling + meta + "\n\n" + today + "\nCode: " + code_name  
    meta_file = open(os.path.join(save_dir, f"bootstrap_summary_{data_col}.txt"), 'wt')
    meta_file.write(meta)
    meta_file.close()
        
        
    # print the run time
    # Record end time
    end_time = time.time()
    runtime = round(end_time - start_time, 2)
    print(f"\nRuntime: {runtime} seconds")
    
    # create a new directory to put the data for plotting later
    plot_dir = os.path.join(save_dir, "Plots")
    
    # create a folder to put data for plotting and eventually the histograms (a separate function)
    try:
        # Create target Directory in the location
        os.mkdir(plot_dir)
    except FileExistsError:
        pass
    
    # Save bootstats as a files with labeled columns
    # Each group in bootstats will be a separate file
    i = 0
    
    for g in bootstats:
        bootstats_df = pd.DataFrame(g, columns = sort_list)
        # get the current group name, remove any spaces so this can be used for the filename
        group_name = group_list[i].replace(" ", "_")
        # save the full bootstats array for the current group
        filename = f"{data_col}_bootstrap_medians_{group_name}_nboot{nboot}.csv"
        # save the CSV file
        bootstats_df.to_csv(os.path.join(save_dir, filename), index = False)
        
        # save the summary stats of the bootstats population for each column
        summary_stats = bootstats_df.describe()
        # also save the summary statistics as a csv file
        filename2 = f"{data_col}_bootstrap_stats_{group_name}_nboot{nboot}.csv"
        # save the CSV file
        summary_stats.to_csv(os.path.join(save_dir, filename2))
        
        # Also save the data in the long format for plotting
        # Melt the DataFrame to long format and specify the name of the new column
        bootstats_long = pd.melt(bootstats_df, var_name=sort_col, value_name=data_col)
        
        # add a column for the metric and the group name. This will make plotting easier.
        
        bootstats_long["Group"] = [group_name] * len(bootstats_long)
        bootstats_long["Metric"] = [data_col] * len(bootstats_long)
        
        # save the long format dataframe
        filename3 = f"{data_col}_bootstrap_medians_{group_name}_for_plotting.csv"
        bootstats_long.to_csv(os.path.join(plot_dir, filename3), index=False)
        
        i += 1 # move index over for the next group
    
    
    return bootstats, save_dir

#%% Define get_direct_prob function which will run the statistics to determine the direct probability that sample 2 is larger than sample 1
    # Note that this code was adapted from Saravanan et al.
    
def get_direct_prob(sample1, sample2, bins = 100):
    '''
    get_direct_prob Returns the direct probability of items from sample2 being
    greater than or equal to those from sample1.
       Sample1 and Sample2 are two bootstrapped samples and this function
       directly computes the probability of items from sample 2 being greater
       than or equal to those from sample1. Since the bootstrapped samples are
       themselves posterior distributions, this is a way of computing a
       Bayesian probability. The joint matrix can also be returned to compute
       directly upon.
       
       Edited so you can control the binning if you want. Default is 100.
    '''
    joint_low_val = min([min(sample1),min(sample2)])
    joint_high_val = max([max(sample1),max(sample2)])
    
    p_joint_matrix = np.zeros((bins,bins)) # was (100,100)
    p_axis = np.linspace(joint_low_val,joint_high_val,num=bins) # was 100 before
    edge_shift = (p_axis[2] - p_axis[1])/2
    p_axis_edges = p_axis - edge_shift
    p_axis_edges = np.append(p_axis_edges, (joint_high_val + edge_shift))

    #Calculate probabilities using histcounts for edges.

    p_sample1 = np.histogram(sample1,bins=p_axis_edges)[0]/np.size(sample1)
    p_sample2 = np.histogram(sample2,bins=p_axis_edges)[0]/np.size(sample2)

    #Now, calculate the joint probability matrix:

    for i in np.arange(np.shape(p_joint_matrix)[0]):
        for j in np.arange(np.shape(p_joint_matrix)[1]):
            p_joint_matrix[i,j] = p_sample1[i]*p_sample2[j]
            
    #Normalize the joint probability matrix:
    p_joint_matrix = p_joint_matrix/np.sum(p_joint_matrix)
    
    #Get the volume of the joint probability matrix in the upper triangle:
    p_test = np.sum(np.triu(p_joint_matrix))
    
    return p_test, p_joint_matrix

#%%% Make a custom function to run all the comparisons we want through the get_direct_prob() function above and save the data

def direct_prob_comp(metric_name, group1, group2, groups=["CTL", "cKO"], save_dir = os.getcwd(), bins = 100, alpha = 0.05):
    ''' A custom function to compare resampled populations across genotypes and across layers. Output will be a table comparing groups for each layer, and a separate table with the layer comparisons for each group. These will also be saved as CSV files. When running this, put the control group first as group 1. The input for group1 and group2 should be the array of bootstrapped medians for each group from the bootstrap output. Optionally, you can specify the save directory location and the alpha value for the comparisons.
    '''
    
    # parse out the different groups to compare
    group1_SO = group1[:,0] 
    group2_SO = group2[:,0]
    
    group1_SR = group1[:,1]
    group2_SR = group2[:,1]
    
    group1_SLM = group1[:,2]
    group2_SLM = group2[:,2]
    
    # create a subdirectory to save the matrices all in one folder   
    matrix_dir = os.path.join(save_dir, "probability_matrices")

    # create a folder within the save_folder to put all the bootstrap matrices
    try:
        # Create target Directory in the location
        os.mkdir(matrix_dir)
    except FileExistsError:
        pass
            
    # run the get_direct_prob function to get the stats comparing CTL and cKO in each layer
    # P here is the probability that the KO is greater than or equal to the CTL
    [p_SO, matrix_SO] = get_direct_prob(group1_SO, group2_SO, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_SO.csv"), matrix_SO, delimiter=",")
    [p_SR, matrix_SR] = get_direct_prob(group1_SR, group2_SR, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_SR.csv"), matrix_SR, delimiter=",")
    [p_SLM, matrix_SLM] = get_direct_prob(group1_SLM, group2_SLM, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_SLM.csv"), matrix_SLM, delimiter=",")
    
    # see if there is already a probability table in the save directory, if so we can just add a row to it
    # file name to look for
    prob_file = os.path.join(save_dir, "bootstrap_prob_table.csv")
    # get the files in the current save location
    file_names = [f for f in os.listdir(save_dir) if f.endswith(('.csv', ".CSV"))] #only grab the CSV files in the save folder
    
    # check if filename exists in this location
    if "bootstrap_prob_table.csv" in file_names:
        # open the probability file
        prob_tb = pd.read_csv(prob_file, index_col = 0)
        print("\nGetting probability file...\n")
    else: # if there is no probability file yet
        prob_tb = pd.DataFrame() # create empty dataframe
        
    # add the new p values as a new row to the probability table
    temp = pd.DataFrame({"SO": [p_SO], "SR":[p_SR], "SLM": [p_SLM]})
    temp.index = [metric_name]
    prob_tb = pd.concat([prob_tb, temp], axis=0)
    
    # save the probability table as a CSV file
    filename1 = "bootstrap_prob_table.csv"
    prob_tb.to_csv(os.path.join(save_dir, filename1))
    
    # Let's make the comparisons across SR and SLM for mito area in the WT and the KO
    [p_g1_SR_SLM, matrix_g1_SR_SLM] = get_direct_prob(group1_SR, group1_SLM, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_" + groups[0] + "_SR_SLM.csv"), matrix_g1_SR_SLM, delimiter=",")
    
    [p_g2_SR_SLM, matrix_g2_SR_SLM] = get_direct_prob(group2_SR, group2_SLM, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_" + groups[1] + "_SR_SLM.csv"), matrix_g2_SR_SLM, delimiter=",")
    
    [p_g1_SO_SLM, matrix_g1_SO_SLM] = get_direct_prob(group1_SO, group1_SLM, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_" + groups[0] + "_SO_SLM.csv"), matrix_g1_SO_SLM, delimiter=",")
    
    [p_g2_SO_SLM, matrix_g2_SO_SLM] = get_direct_prob(group2_SO, group2_SLM, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_" + groups[1] + "_SO_SLM.csv"), matrix_g2_SO_SLM, delimiter=",")
    
    [p_g1_SO_SR, matrix_g1_SO_SR] = get_direct_prob(group1_SO, group1_SR, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_" + groups[0] + "_SO_SR.csv"), matrix_g1_SO_SR, delimiter=",")
    
    [p_g2_SO_SR, matrix_g2_SO_SR] = get_direct_prob(group2_SO, group2_SR, bins = bins)
    np.savetxt(os.path.join(matrix_dir, "boot_matrix_" + metric_name + "_" + groups[1] + "_SO_SR.csv"), matrix_g2_SO_SR, delimiter=",")
    
    # save a prob table comparing the mito area across layers in cKO and CTL
    
    prob_tb_layer = pd.DataFrame({groups[0]: [p_g1_SR_SLM, p_g1_SO_SLM, p_g1_SO_SR], groups[1]:[p_g2_SR_SLM, p_g2_SO_SLM, p_g2_SO_SR] })
    prob_tb_layer.index = ["SR vs SLM", "SO vs SLM", "SO vs SR"]
    # save the probability table
    prob_tb_layer_clean = prob_tb_layer.reset_index()
    filename2 = f"bootstrap_prob_table_{metric_name}_layer.csv"
    prob_tb_layer_clean.to_csv(os.path.join(save_dir, filename2), encoding='utf-8', index=False)
     
    # Using alpha, calculate the low and high cutoffs for significance. 
    # Anything in between these values is not statistically significant.
    sig_cutoff_H = 1-(alpha/2) # p value higher than this number indicates the second group is statistically larger than the first
    sig_cutoff_L = alpha/2 # p value less than this number indicates the first group is statistically larger than the second
    
    
    # print the probabilities to the console
    print("\nBootstrap probabilities comparing " + groups[0] + " to " + groups[1] + ":")
    print(prob_tb)
    
    print("\nAny value above " + str(sig_cutoff_H) + " means " + groups[1] + " is greater than " + groups[0] + ".")
    print("Any value below " + str(sig_cutoff_L) + " means " + groups[0] + " is greater than " + groups[1] + ".")
    
    print("\n Comparing SR and SLM in both groups:")
    print(groups[0] + ":", round(p_g1_SR_SLM, 2))
    print(groups[1] + ":", round(p_g2_SR_SLM,2))
    
    # print whether any of the comparisons are significant
    
    print("\nAny value above " + str(sig_cutoff_H) + " means SLM is greater than the SR.")
    print("Any value below " + str(sig_cutoff_L) + " means SR is greater than the SLM.")
    
    return prob_tb, prob_tb_layer

#%%% Define function to plot a probability matrix

def plot_prob_dist( matrix, p, group_labels, title, subplots = 1, save_dir = os.getcwd(), figsize = (6,4), bins = 100):
    
    ''' Function to plot a probability distribution comparing two groups. Each probability distribution will show the p value in the text. By default, there is one matrix plotted on a single plot. However, you can plot multiple matrices on different subplots if you specify the number of subplots. If you have more than one subplot, you need to include the same number of matrices, p values and group labels as you have subplots. The group labels should be provided in the order of Y first, then X. The title you provide will be used in the plot filename, so make it informative.
    
    Example group_labels = [["SR", "SLM"], ["SO, "SLM"], ["SO", "SR"]]
    This will plot SR vs SLM, SO vs SLM and then SO vs SR on 3 different subplots.
                                            
    '''
    
    fig, axes = plt.subplots(1, subplots, figsize=figsize)
    fig.tight_layout(pad=1)
    
    if subplots == 1:
        plt.title(f"Probability Distribution:\n{title}", y=20, pad=40, fontsize=14, fontweight="bold")
        
        axes.matshow(matrix, cmap=plt.cm.Blues)
        axes.grid(False)
        axes.plot([0,bins],[0,bins], color='Grey', linestyle='--', linewidth=2)
        axes.set_ylabel(group_labels[0], labelpad=3, fontweight="bold")
        axes.set_xlabel(group_labels[1], labelpad=10, fontweight="bold")
        axes.xaxis.set_label_position('top')
        axes.xaxis.set_ticks_position('top') 
        # plt.tick_params(axis='x', which='both', bottom=False) # removes ticks on bottom X axis
        text_label = "P_boot = " + str(round(p,3))
        axes.text(bins/20,bins-(bins/50), text_label, fontsize=10) # was 5,bins-2
        
    else:       
        for ax in range(subplots):
        
            axes[ax].matshow(matrix[ax], cmap=plt.cm.Blues)
            axes[ax].grid(False)
            axes[ax].plot([0,bins],[0,bins], color='Grey', linestyle='--', linewidth=2)
            axes[ax].set_ylabel(group_labels[ax][0], labelpad=-5, fontweight="bold", fontsize=12)
            axes[ax].set_xlabel(group_labels[ax][1], labelpad=6, fontweight="bold", fontsize=12)
            axes[ax].xaxis.set_label_position('top')
            axes[ax].xaxis.set_ticks_position('top') 
            # plt.tick_params(axis='x', which='both', bottom=False) # removes ticks on bottom X axis
            text_label = "P_boot = " + str(round(p[ax],3))
            axes[ax].text(bins/20,bins-(bins/50), text_label, fontsize=10)
            
        plt.suptitle(f"Probability Distribution:\n{title}", y=0.98, fontsize=14, fontweight="bold")
    
    # save the plot as a file in the given save directory
    filename = os.path.join(save_dir, "prob_dist_" + title + ".png")
    plt.savefig(filename, bbox_inches='tight')
    
#%% Function to plot the histograms of the bootstrapped data and get a summary table of the means + standard error for all metrics

def boot_hist(data, data_col, x_label, save_dir = os.getcwd(), hue_col = "Layer", group_name = [], binsize = 0.005, figsize = (8,5), xlim = [0,0.35], ylim = [0,120]):
    ''' Function to create a nice histogram of the resampled median populations comparing across the "hue" column, usually Layer. Each value of hue will be plotted as a separate line on the histogram. The histogram image will be saved in the save_directory.
    
 Required inputs:
    data = long format dataframe with a column for the data to plot and a column for the hue values.
    data_col = column name for the column to be plotted in the histogram
    x_label = Specify the x axis label for the histogram (include units)
    
 Optional inputs:
     save_dir = The location to save the plot as a PNG file. Must be a full path.
     hue_column = The column name to separate the hues. Each hue will be plotted as a separate line. By default, the hue column is "Layer"
     group_name = If desired, specify a group name for the histogram (ie. CTL or cKO). This will be added to the plot title and the filename to distinguish it from other groups. Recommended to do this if you'll have multiple histograms of the same metric in the same locaion, otherwise if the filename is the same this code may save over the previous file.
     binsize = specify the size of the bins. Default is 0.005.
     figsize = Specify the figure size in two dimensions. Width then height. Default is 8x5 inches.
     xlim and ylim: Specify the x or y range of the plot. Provide two numbers. Default: xlim = [0,0.35], ylim = [0,120]
     
    '''
    
    # Plotting the histogram
    sns.set(style="ticks")
    
    plt.figure(figsize=figsize)
    plt.tight_layout()
    sns.histplot(data=data, x=data_col, element="step", hue=hue_col, stat='count', common_norm=False, binwidth=binsize)
    
    # Manually set legend labels
    legend_labels = data[hue_col].unique()
    plt.legend(title=hue_col, labels=legend_labels)
    
    if len(group_name) > 0:
        plt.title(f"Bootstrap medians of {data_col} in {group_name}")
    else:
        plt.title(f"Bootstrap medians of {data_col}")
    plt.xlabel(x_label)
    plt.ylabel('Frequency')
    plt.xlim(xlim)
    plt.ylim(ylim)
    
    # save the plot
    if len(group_name) > 0:
        plt.savefig(os.path.join(save_dir, f"hist_bootstrap_medians_{data_col}_{group_name}.png"))
    else:
        plt.savefig(os.path.join(save_dir, f"hist_bootstrap_medians_{data_col}.png"))

#%%% Create a custom function to summarize the bootstrap results and plot a histogram for all the data in the "Plot" folder. Run this at the end after running the bootstrap on all your metrics of interest. That way all the data will be in the Plot folder and the overall means for all metrics can be added to the means summary table.

def get_boot_summary( boot_loc, data_cols, groups = ["CTL", "cKO"], overall_means = True, hist = True, perc_plot = True, n_boot = 1000):
    '''
    Function to plot a histogram for each long-form dataframe in the "Plot" folder after running the bootstrap. This function will also get the overall means for each metric in the bootstrap main folder and summarize the results in a nice tidy dataframe of means with the standard error in parenthases. For data_cols, put in a list of the data columns that the bootstrap was run on (ie. ["Size", "Count"]). If you don't want to get a table of overall means, change overall_means to False. If you don't want to plot the histograms, change hist to False.
    
    This function is created to accompany the custom "get_boot_sample() function, so it's important that the files are structured and named appropriately. The data you want to plot histograms of should be in separate long form dataframes in a folder called "Plots". If you want to create the overall means table, make sure you have csv files for the summary stats in the main bootstrap folder.
    '''
    if hist:
        
        # get the folder with the plotting data
        plot_dir = os.path.join(boot_loc, "Plots")
        
        if "Plots" in os.listdir(boot_loc): # check if the plot directory exists
            files = [f for f in os.listdir(plot_dir) if f.endswith(('.csv', ".CSV"))]
            
            for f in files:
                # import the long form dataframe
                long_df = pd.read_csv(os.path.join(plot_dir, f))
                
                # Get the metric name from the "Metric" column
                try:
                    metric_name = long_df["Metric"].unique()[0]
                except NameError:
                    print("WARNING: The name of the metric could not be obtained from the file. There should be a columm named ''Metric'' that contains the name of the metric analyzed.")
                    metric_name = "unknown"
                
                # look for the group name in the "Group" column
                try:
                    group_name = long_df["Group"].unique()[0]
                except NameError:
                    print("WARNING: The group name could not be obtained from the file. There should be a columm named ''Group'' that contains the name of the metric analyzed.")
                    group_name = "unknown"
                    
                # define any necessary variables to run boot_hist()
                
                # specify parameters if the metric is count
                if metric_name == "Count":
                    x_label = "Count"
                    bins = 1
                    xlim = [0,30]
                    if n_boot == 1000: # adjust the ylim with the number of repetitions
                        ylim = [0, 500]
                    elif n_boot == 10000:
                        ylim = [0,3000]
                    else:
                        ylim = [1, 1000]
                elif metric_name == "Area_um_sq":
                    x_label = "Area (um2)"
                    bins = 0.005
                    xlim = [0,0.3]
                    if n_boot == 1000: # adjust the ylim with the number of repetitions
                        ylim = [0, 200]
                    elif n_boot == 10000:
                        ylim = [0,1400]
                    else:
                        ylim = [0,1000]
                    
                # plot the histogram for the current data 
                boot_hist(data = long_df, data_col = metric_name, x_label = x_label, save_dir = plot_dir, group_name = group_name, binsize = bins, xlim = xlim, ylim = ylim)
                    
                
    if overall_means: # if you run out of time, just plot the histograms and do this part manually 
        
        # find the right CSV file
        csv_list = os.listdir(boot_loc)
        
        # get only the stats files
        stats_csv = [s for s in csv_list if "stats" in s]
            
        # import the data for each metric and each group (these are the files spit out by the get_boot_sample() function)
        # find the file for mito size control
        CTL_size_file = [s for s in stats_csv if "Area" in s and "Cre_-" in s][0]
        CTL_size = pd.read_csv(os.path.join(boot_loc, CTL_size_file), index_col=0)
        
        # find the file for mito size in the cKO
        cKO_size_file = [s for s in stats_csv if "Area" in s and "Cre_+" in s][0]
        cKO_size = pd.read_csv(os.path.join(boot_loc, cKO_size_file), index_col=0)
        
        # find the file for mito size in the cKO
        CTL_count_file = [s for s in stats_csv if "Count" in s and "Cre_-" in s][0]
        CTL_count = pd.read_csv(os.path.join(boot_loc, CTL_count_file), index_col=0)
        
        # find the file for mito size in the cKO
        cKO_count_file = [s for s in stats_csv if "Count" in s and "Cre_+" in s][0]
        cKO_count = pd.read_csv(os.path.join(boot_loc, cKO_count_file), index_col=0)
        
        # properly format the means and the st dev
        size_WT_SO = f"{round(CTL_size.loc['mean', 'Basal'],2)} (\u00B1{round(CTL_size.loc['std', 'Basal'],2)})"
        size_WT_SR = f"{round(CTL_size.loc['mean', 'Proximal'],2)} (\u00B1{round(CTL_size.loc['std', 'Proximal'],2)})"
        size_WT_SLM = f"{round(CTL_size.loc['mean', 'Distal'],2)} (\u00B1{round(CTL_size.loc['std', 'Distal'],2)})"
        
        size_KO_SO = f"{round(cKO_size.loc['mean', 'Basal'],2)} (\u00B1{round(cKO_size.loc['std', 'Basal'],2)})"
        size_KO_SR = f"{round(cKO_size.loc['mean', 'Proximal'],2)} (\u00B1{round(cKO_size.loc['std', 'Proximal'],2)})"
        size_KO_SLM = f"{round(cKO_size.loc['mean', 'Distal'],2)} (\u00B1{round(cKO_size.loc['std', 'Distal'],2)})"
        
        # same for count
        count_WT_SO = f"{round(CTL_count.loc['mean', 'Basal'],2)} (\u00B1{round(CTL_count.loc['std', 'Basal'],2)})"
        count_WT_SR = f"{round(CTL_count.loc['mean', 'Proximal'],2)} (\u00B1{round(CTL_count.loc['std', 'Proximal'],2)})"
        count_WT_SLM = f"{round(CTL_count.loc['mean', 'Distal'],2)} (\u00B1{round(CTL_count.loc['std', 'Distal'],2)})"
        
        count_KO_SO = f"{round(cKO_count.loc['mean', 'Basal'],2)} (\u00B1{round(cKO_count.loc['std', 'Basal'],2)})"
        count_KO_SR = f"{round(cKO_count.loc['mean', 'Proximal'],2)} (\u00B1{round(cKO_count.loc['std', 'Proximal'],2)})"
        count_KO_SLM = f"{round(cKO_count.loc['mean', 'Distal'],2)} (\u00B1{round(cKO_count.loc['std', 'Distal'],2)})"

            
        overall_mean_tb = pd.DataFrame({"CTL SO": [size_WT_SO, count_WT_SO], "CTL SR": [size_WT_SR, count_WT_SR], "CTL SLM": [size_WT_SLM, count_WT_SLM], "cKO SO": [size_KO_SO, count_KO_SO], "cKO SR": [size_KO_SR, count_KO_SR], "cKO SLM": [size_KO_SLM, count_KO_SLM]})
        overall_mean_tb.index = ["Size", "Count"]
        
        # Save the formatted table as a CSV file
        overall_mean_tb.to_csv(os.path.join(boot_loc, "overall_means_tb.csv"), encoding='utf-16')
        

#%%% Now that we have defined all of our functions, run bootstrap on mito count per tile, resampling at the level of Animal and Stub. We can run this multiple times with different tile sampling. Plan to run this sampling 20, 32, and 50 tiles.

BootCount, save_loc = get_boot_sample(tile_data, data_col="Count", lev1_col="Animal", lev2_col="Section_ID", group_col="Genotype", sort_col="Layer", nboot=n_boot, direct = loc, tile_samp = 50, mito_samp = 100)

# split the data into the groups
# In our case, the CTL is the first array (0), then the cKO (1)
# It's a good idea to check that the order of these groups is correct
count_CTL = BootCount[0] 
count_cKO = BootCount[1]

# Run the direct_prob_comp() function to get the stats comparing CTL and cKO in each layer
# Make sure group1 and group2 are in the same order as the group names listed in "groups"
# P here is the probability that the KO is greater than or equal to the CTL
prob_tb_count, prob_tb_count_layer = direct_prob_comp(metric_name = "Count", group1 = count_CTL, group2 = count_cKO, groups=["CTL", "cKO"], save_dir = save_loc, bins = 10)

    
#%% Also run the bootstrap on individual mito area, resampliing at the Animal, Stub, and Tile level.

BootSize, save_loc = get_boot_sample(mito_data, data_col="Area_um_sq", lev1_col="Animal", lev2_col="Section_ID", lev3_col = "Tile", group_col="Genotype", sort_col="Layer", nboot=n_boot, direct = loc, tile_samp = 50, mito_samp = 100)

# Split the data into the groups
# In our case, the CTL is the first array (0), then the cKO (1)
# It's a good idea to check that the order of these groups is correct
size_CTL = BootSize[0] 
size_cKO = BootSize[1]

# Run the direct_prob_comp() function to get the stats comparing CTL and cKO in each layer
# Make sure group1 and group2 are in the same order as the group names listed in "groups"
# P here is the probability that the KO is greater than or equal to the CTL
(prob_tb_size, prob_tb_size_layer) = direct_prob_comp(metric_name = "Area", group1 = size_CTL, group2 = size_cKO, groups=["CTL", "cKO"], save_dir = save_loc)

# After running all the metrics of interest through the bootstrap, run the function below to plot bootstrap distributions as histograms and get the overall means and standard error. The metrics will be collated together into one table.

get_boot_summary(boot_loc = save_loc, data_cols = ["Area", "Count"], groups = ["CTL", "cKO"], n_boot = n_boot)
    
    
#%%% Now plot the probability distrubutions using the plot_prob_dist() function.

# Here you can import the data if it's not already there, so the following sections can be run without running the entire bootstrap again
#!!! You can manually set the path to your bootstrap directory below (change path name and uncomment the line below)
# save_loc = 'C:\\Users\\example_path\\Bootstrap_median_EM_mitos_10000_tile_samp_50_mito_samp_100'

# Note that the full code should run without this, since save_loc is output by the bootstrap function in the previous section

# First, plot the layer comparisons for mito area in the CTL

matrix_save = os.path.join(save_loc, "probability_matrices")

# import the probability matrices for each layer comparison in the CTL

area_SR_SLM_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_CTL_SR_SLM.csv"), delimiter=",")

area_SO_SLM_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_CTL_SO_SLM.csv"), delimiter=",")

area_SO_SR_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_CTL_SO_SR.csv"), delimiter=",")

WT_matrix_list = [area_SO_SR_matrix, area_SO_SLM_matrix, area_SR_SLM_matrix]


# import the probability table with the p value for CTL SR vs SLM
prob_tb_layer = pd.read_csv(os.path.join(save_loc, "bootstrap_prob_table_Area_layer.csv"), index_col=0)

area_SR_SLM_P = prob_tb_layer.loc["SR vs SLM", "CTL"]
area_SO_SLM_P = prob_tb_layer.loc["SO vs SLM", "CTL"]
area_SO_SR_P = prob_tb_layer.loc["SO vs SR", "CTL"]

WT_prob_list = [area_SO_SR_P, area_SO_SLM_P, area_SR_SLM_P]

group_labels = [["CTL SO", "CTL SR"], ["CTL SO", "CTL SLM"], ["CTL SR", "CTL SLM"]]

plot_title = "Mito_Area_CTL"

plot_prob_dist(matrix = WT_matrix_list, p = WT_prob_list, group_labels = group_labels, title = plot_title, subplots = 3, figsize = (7,4), save_dir = matrix_save)

#%%  Now let's plot the comparison of mitochondrial count in the CTL across layers

count_SR_SLM_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_CTL_SR_SLM.csv"), delimiter=",")

count_SO_SLM_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_CTL_SO_SLM.csv"), delimiter=",")

count_SO_SR_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_CTL_SO_SR.csv"), delimiter=",")

WT_matrix_list_count = [count_SO_SR_matrix, count_SO_SLM_matrix, count_SR_SLM_matrix]

# get the p value for SR vs SLM CTL from the count probability table
prob_tb_layer_count = pd.read_csv(os.path.join(save_loc, "bootstrap_prob_table_Count_layer.csv"), index_col=0)

count_SR_SLM_P = prob_tb_layer_count.loc["SR vs SLM", "CTL"]
count_SO_SLM_P = prob_tb_layer_count.loc["SO vs SLM", "CTL"]
count_SO_SR_P = prob_tb_layer_count.loc["SO vs SR", "CTL"]

WT_prob_list_count = [count_SO_SR_P, count_SO_SLM_P, count_SR_SLM_P]

group_labels = [["CTL SO", "CTL SR"], ["CTL SO", "CTL SLM"], ["CTL SR", "CTL SLM"]]

plot_title2 = "Mito_Count_CTL"

# plot the prob distributions with the plot_prob_dist() function
plot_prob_dist(matrix = WT_matrix_list_count, p = WT_prob_list_count, group_labels = group_labels, title = plot_title2, subplots = 3, figsize = (7,4), save_dir = matrix_save, bins = 10)

#%% plot the probability distributions for mito area in the cKO across layers

# import the probability matrices for each layer comparison
area_SR_SLM_KO_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_cKO_SR_SLM.csv"), delimiter=",")

area_SO_SLM_KO_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_cKO_SO_SLM.csv"), delimiter=",")

area_SO_SR_KO_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_cKO_SO_SR.csv"), delimiter=",")

KO_matrix_list = [area_SO_SR_KO_matrix, area_SO_SLM_KO_matrix, area_SR_SLM_KO_matrix]


# import the probability table with the p value for CTL SR vs SLM
prob_tb_layer = pd.read_csv(os.path.join(save_loc, "bootstrap_prob_table_Area_layer.csv"), index_col=0)

area_KO_SR_SLM_P = prob_tb_layer.loc["SR vs SLM", "cKO"]
area_KO_SO_SLM_P = prob_tb_layer.loc["SO vs SLM", "cKO"]
area_KO_SO_SR_P = prob_tb_layer.loc["SO vs SR", "cKO"]

KO_prob_list = [area_KO_SO_SR_P, area_KO_SO_SLM_P, area_KO_SR_SLM_P]

group_labels = [["cKO SO", "cKO SR"], ["cKO SO", "cKO SLM"], ["cKO SR", "cKO SLM"]]

plot_title3 = "Mito_Area_cKO"

plot_prob_dist(matrix = KO_matrix_list, p = KO_prob_list, group_labels = group_labels, title = plot_title3, subplots = 3, figsize = (7,4), save_dir = matrix_save)

#%%  Now let's plot the comparison of mitochondrial count in the cKO across layers

count_SR_SLM_KO_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_cKO_SR_SLM.csv"), delimiter=",")

count_SO_SLM_KO_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_cKO_SO_SLM.csv"), delimiter=",")

count_SO_SR_KO_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_cKO_SO_SR.csv"), delimiter=",")

KO_matrix_list_count = [count_SO_SR_KO_matrix, count_SO_SLM_KO_matrix, count_SR_SLM_KO_matrix]

# get the p value for SR vs SLM CTL from the count probability table
prob_tb_layer_count = pd.read_csv(os.path.join(save_loc, "bootstrap_prob_table_Count_layer.csv"), index_col=0)

count_KO_SR_SLM_P = prob_tb_layer_count.loc["SR vs SLM", "cKO"]
count_KO_SO_SLM_P = prob_tb_layer_count.loc["SO vs SLM", "cKO"]
count_KO_SO_SR_P = prob_tb_layer_count.loc["SO vs SR", "cKO"]

KO_prob_list_count = [count_KO_SO_SR_P, count_KO_SO_SLM_P, count_KO_SR_SLM_P]

group_labels = [["cKO SO", "cKO SR"], ["cKO SO", "cKO SLM"], ["cKO SR", "cKO SLM"]]

plot_title4 = "Mito_Count_cKO"

# plot the prob distributions with the plot_prob_dist() function
plot_prob_dist(matrix = KO_matrix_list_count, p = KO_prob_list_count, group_labels = group_labels, title = plot_title4, subplots = 3, figsize = (7,4), save_dir = matrix_save, bins = 10)


#%% Plotting mito area CTL to cKO comparisons
# Import the area matrices
gen_SO_area_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_SO.csv"), delimiter=",")
gen_SR_area_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_SR.csv"), delimiter=",")
gen_SLM_area_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Area_SLM.csv"), delimiter=",")

#  Set up the parameters to run the area probability distributions
gen_matrix_list_area = [gen_SO_area_matrix, gen_SR_area_matrix, gen_SLM_area_matrix] # list of matrices

# import the table with the prob comparisons by genotype
prob_tb = pd.read_csv(os.path.join(save_loc, "bootstrap_prob_table.csv"), index_col=0)

# pull out the p values for each comparison for area
gen_SO_area_P = prob_tb.loc["Area", "SO"]
gen_SR_area_P = prob_tb.loc["Area", "SR"]
gen_SLM_area_P = prob_tb.loc["Area", "SLM"]

gen_group_labels = [["CTL SO", "cKO SO"], ["CTL SR", "cKO SR"], ["CTL SLM", "cKO SLM"]]

plot_title5 = "Mito_area_CTL_vs_cKO"

gen_prob_list_area = [gen_SO_area_P, gen_SR_area_P, gen_SLM_area_P] # list of probabilities

# plot the prob distributions with the plot_prob_dist() function for area
plot_prob_dist(matrix = gen_matrix_list_area, p = gen_prob_list_area, group_labels = gen_group_labels, title = plot_title5, subplots = 3, figsize = (7,4), save_dir = matrix_save)

#%% Plotting mito area CTL to cKO comparisons

# now import the count matrices
gen_SO_count_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_SO.csv"), delimiter=",")
gen_SR_count_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_SR.csv"), delimiter=",")
gen_SLM_count_matrix = np.genfromtxt(os.path.join(save_loc, "probability_matrices", "boot_matrix_Count_SLM.csv"), delimiter=",")

gen_matrix_list_count = [gen_SO_count_matrix, gen_SR_count_matrix, gen_SLM_count_matrix] # list of matrices

# pull out the p values for each comparison for count
gen_SO_count_P = prob_tb.loc["Count", "SO"]
gen_SR_count_P = prob_tb.loc["Count", "SR"]
gen_SLM_count_P = prob_tb.loc["Count", "SLM"]

gen_group_labels = [["CTL SO", "cKO SO"], ["CTL SR", "cKO SR"], ["CTL SLM", "cKO SLM"]]

plot_title6 = "Mito_count_CTL_vs_cKO"

gen_prob_list_count = [gen_SO_count_P, gen_SR_count_P, gen_SLM_count_P] # list of probabilities

# plot the prob distributions with the plot_prob_dist() function for area
plot_prob_dist(matrix = gen_matrix_list_count, p = gen_prob_list_count, group_labels = gen_group_labels, title = plot_title6, subplots = 3, figsize = (7,4), save_dir = matrix_save, bins = 10)

print(f"Bootstrap has finished! Probability matrices have been plotted for each comparison across layers and genotypes.\n\nData was saved in the location:\n{save_loc}")
    
    




