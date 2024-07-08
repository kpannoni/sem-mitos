# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 13:37:38 2022

A file containing common functions for the mitochondrial analysis that are repeated in multiple analysis codes.

The first function listed is to subtract background fluorescence from the mean ROI fluorescence. The function input will be the dataframe with the negative control data and the main dataframe with the name of the column with ROI fluorescence. The output will be the same main dataframe with a new column containing the background subtracted fluorescence.

The second function is to rank data from largest to smallest and separate the data into quartiles, which will be plotted as violin plots for each group in the data. Input is the dataframe and data column name as well as the desired save path. This function will save CSV files with the data for each quartile and a PNG image of each group violin plot in the save path, and returns a list with the quartile data for each group.

@author: pannoni
"""

def sub_BG_fl(data, fl_col, neg_ctrl, bg_col = 0, group = 0):
    '''    
    Function to subtract background fluorescence from mean ROI fluorescence separately for each layer and subregion for each IHC cohort in your dataframe. There should a backround fluorescence for each layer, subregion and IHC cohort in the "neg_ctrl" dataframe. If there is only one cohort, you should still have a column for IHC Cohort. If you want group your controls by an additional factor (ie. treatment or animal), add the column name for your grouping factor under "group =". By default, the code assumes the name of the fluorescence column in your neg_ctrl dataframe is the same as your main dataframe (since they are pulled from the same excel sheet). If they are named differently, you can include the background fluorescence column name under "bg_col ="
    
    Example syntax: 
    output_df = sub_BG_fl(data = main_df, fl_col = "ROI fluor", neg_ctrl = neg_ctrl_df, group = "Treatment")
    
    In the case above, you would need to have a column called "ROI fluor", "Layer", "Subregion", "IHC Cohort" and a column "Treatment" in both your main dataframe and your negative control dataframe.
    
    This code will output the same dataframe with an added column for the background fluorescence ("bg_fl") and a column for the fluorescence without the background ("Fl_no_bg").
    '''
    import pandas as pd
    print("\nSubtracting background fluorescence from the mean ROI fluorescence...")
    
    # Check that the dataframe provided has the correct columms.
    
    
    # empty dataframe to put the new fluorescence data
    fluor_df = pd.DataFrame()
    
    # iterate through the neg controls in the given neg_ctrl dataframe
    for neg in neg_ctrl.iterrows():
        # get the current layer, subregion and cohort
        curr_layer = neg[1]["Layer"]
        curr_sub = neg[1]["Subregion"]
        curr_cohort = neg[1]["IHC Cohort"]
        
        # if a third grouping is specified
        if group != 0:
            curr_group = neg[1][group] # get the current group
            # pull the data for the current layer, subregion and group
            curr_data = data.where((data["Layer"] == curr_layer) & (data["Subregion"] == curr_sub) & (data[group] == curr_group) & (data["IHC Cohort"] == curr_cohort)).dropna(how="all")
        else: # if no other grouping specified
            # pull the data for the current layer and subregion
            curr_data = data.where((data["Layer"] == curr_layer) & (data["Subregion"] == curr_sub) & (data["IHC Cohort"] == curr_cohort)).dropna(how="all")
        
        # Get the background fluor from the current neg control
        bg_fl = neg[1][fl_col]
        
        # subtract the bg fluor from the current whole ROI fluor data
        fl_no_bg = curr_data[fl_col] - bg_fl
        
        # create temp dataframe with the data
        temp = pd.concat([curr_data, pd.DataFrame({"bg_fl": bg_fl, "Fl_no_bg": fl_no_bg})], axis=1)
        
        #append the temp dataframe to norm_fl_df
        fluor_df = fluor_df.append(temp)
    
    # Re-sort the index in the same order as the original dataframe
    fluor_df = fluor_df.sort_index()
    
    return fluor_df

#%%% Short function to perform two-way ANOVA and post hoc

def ANOVA_posthoc(data, data_col, savepath, group="Group", para=False):
    '''This function will run a standard two-way ANOVA (genotype and layer) with pairwise posthoc tests using the pingouin package. By default, the ANOVA will be done with the factors "Genotype" and "Layer" and the posthoc tests will be done across the group names in a column called "Group", but a different grouping column can be specified. Make sure these columns exist and are named correctly. Every group you want to compare in the posthoc should be in the grouping column (You will probably need to combine Genotype and Layer to make the group column). The posthoc tests will be non-parametric with Sidak's adjustment unless otherwise specified.
    '''
    import pingouin as pg
    import os
    
    # run a two way ANOVA
    ANOVA = pg.anova(dv=data_col, between=["Layer","Genotype"], data=data)
    
    # run the posthoc tests
    posthoc = pg.pairwise_tests(data=data, dv=data_col, between=group, parametric=para, padjust='sidak')
    # Save posthoc table as .CSV file
    filename_ANOVA = data_col + "_ANOVA.csv"
    ANOVA.to_csv(os.path.join(savepath, filename_ANOVA))
    filename_posthoc = data_col + "_posthoc.csv"
    posthoc.to_csv(os.path.join(savepath, filename_posthoc))
    return ANOVA, posthoc


#%%% Make a function to rank the data by a column, separate the data into quartiles and then plot each quartile as a violin plot on the same plot within each dendritic layer.

def data_quartiles(data, data_col, save_path, group_col="Layer", output_cols = 0, b_color = "#1e488f", p_color = "#de7e5d", d_color = "#0a888a", colors = 0, y_range = [[0, 0.07], [0.02, 0.11], [0.05, 0.18], [0.1, 2.5]], units = 'Area (\u00B5m\u00B2)'):
    '''This function will take the data in the specified "data_col" column, rank it from largest to smallest, and then separate the data into quartiles within each group in the "group_col" column. The default group used is dendritic layer, in a column named "Layer". If you would like to group the data differently, or you've used a different name for your layer column, please specify the column name of the grouping column with "group_col= ". If another group column is not specified and the "Layer" column doesn't exist, the code will run into an error.
    
    Within each group, the quartiles will be plotted as four violin plots on the same plot for each group. The plots will be saved as .PNG images in the specified "save_path" location, along with CSV files with the ranked data for each group. If grouping by Layer, you can specify which colors you want the plot for each layer to be using "b_color" (basal), p_color (proximal), and/or "d_color" (distal). If no other color is specified, the layer plots will be colored dark blue, rusty orange, and dark teal, respectively. If grouping by some other factor, you can specify a list of colors that will be used for each group plot (in order).
    
    Required input:
        data = dataframe with the data
        data_col = column name of the data
        save_path = full path location to save the images and CSV files.
        
    Optional inputs:
        group_col = Column name of the grouping column. Default (if left blank) is "Layer".
        
        b_color, p_color, d_color = Optional variables to specify a color hex code to use for plotting the three different layers (basal, proximal or distal) instead of the default layer colors (see above). You can specify the color for each layer separately. These colors will only be used if you are grouping by Layer (default). If using another grouping, use "colors" below instead.
        
        colors = If not grouping by layer, use this to specify a list of colors (hex codes or RGB) for plotting your groups. The number of colors in the list should equal the number of groups in your data. Colors will be plotted in the order of your groups in your grouped column. If no color list is specified here, the code will use colors from the default colorblind seaborn palette.
        
        y_range = Set the range of the Y axis for each of the 4 violin plots. Should be a list of 4 sets of ranges (in brackets). Example: [[1,2], [2,4], [4-6], [5-10]] will range from 1-2 for Q1, 2-4 for Q2, etc. A default range is set that works well for mitochondria area, but new ranges will probably need to be set if you're plotting other metrics. Note that all the layer plots will be plotted using the same Y-axis ranges for each quartile. Pick a range that includes the largest and smallest values across all the groups in each quartile.
        
        output_cols = optional list of columns to include in the output dataframe and saved csv file. Default is just to save the rank column and the speficied data column. Columns listed must be in the "data" dataframe.
        
        units = A string with the data units that will be displayed on the Y-axis of the violin plots. Default is Area in um^2.
    
    '''
    
    import pandas as pd
    import os
    import seaborn as sns
    import matplotlib.pyplot as plt

    # Get a list of the groups in the data
    group_list = list(set(data[group_col]))
    
    # Start counter for how many groups there are
    g = 0
    
    # Create an empty list for the output
    output_df = []
    
    # If the grouping is by layer, set the plot color for each layer
    if group_col == "Layer":
        layer_colors = {"Basal":b_color, "Proximal":p_color, "Distal":d_color}
        # sort the group list in order of Basal, Prox then Distal
        layer_order = ["Basal", "Proximal", "Distal"]
        try: # avoids an error if the layer names in the data don't match layer_order
            group_list.sort(key=lambda l:layer_order.index(l))
        except ValueError:
            pass # will not sort the layers
        
    # Create an empty dataframe for the summary statistics with the groups as columns
    summary_stats = pd.DataFrame(columns = group_list, index = ["Mean", "Median", "N_mitos", "Min", "Max"])
    
    # Set the Y axis range for each quartile. We want the same range across layers for better comparison
    Q_cut = y_range # use the specified Y-axes ranges
        
    # Catch if the user entered y_range isn't formatted correctly or doesn't contain enough elements
    if len(y_range) != 4:
        print ("\nERROR: User defined y_range should include a list of four elements separated by brackets. If you wish to specify the Y-axis ranges, please provide a range for each quartile plot formatted as in the following example: \ny_range = [[min_Q1, max_Q1], [min_Q2, max_Q2], [min_Q3, max_Q3], [min_Q4, max_Q4]]\n\nThe default Y-axis range has been used instead.")
        Q_cut = [[0, 0.07], [0.02, 0.11], [0.05, 0.18], [0.1, 2.5]] # use the default Y-axis range instead of y_range
    
    for group in group_list:
        
        # increase group counter
        g = g + 1
        
        # get the data for the current group
        group_data = data.where(data[group_col] == group).dropna(how="all").reset_index(drop=True)
        # order the data in the specified data column, from largest to smallest
        rank_df = group_data.sort_values(by=data_col, ascending=False).reset_index(drop=True)
        
        # add a columm with the rank
        rank_df = rank_df.rename_axis("Rank").reset_index()
        rank_df["Rank"] = rank_df["Rank"] + 1 # Add 1 so the rank starts at 1 and not zero
        
        try: # Change the data type of the Object ID column to "int" so there's no decimals
            rank_df["Object ID"] = rank_df["Object ID"].astype(int)
        except KeyError:
            pass
        
        # separate the ranked data into quartiles, which we will save as CSV files
    
        # get the quartiles of the data
        quartiles = pd.DataFrame(rank_df[data_col]).quantile(q=[.25, .5, .75], interpolation='midpoint')
        
        # Get the borders of the quartiles 
        Q1_cut = quartiles.loc[0.25, data_col]
        med = quartiles.loc[0.5, data_col]
        Q3_cut = quartiles.loc[0.75, data_col]
        
        # columns to include in the output dataframe and csv files
        if output_cols == 0: # no output columns specified
            cols = ["Rank", data_col, group_col] # default, just include the "rank" and the data column
        else:
            cols = output_cols # use the user specified columns instead
            if data_col not in cols: # in case the user left out the data column in their list, add it
                cols.append(data_col)
        
        # get the mito data in each quartile
        Q1 = rank_df[cols][rank_df[data_col]<Q1_cut]
        Q2 = rank_df[cols][(rank_df[data_col]>=Q1_cut) & (rank_df[data_col]<med)]
        Q3 = rank_df[cols][(rank_df[data_col]>=med) & (rank_df[data_col]<Q3_cut)]
        Q4 = rank_df[cols][rank_df[data_col]>=Q3_cut] 
        
        # For easy plotting, let's add a column for quadrant to each quartile dataframe
        Q1["Quartile"] = ["Q1"] * len(Q1["Rank"])
        Q2["Quartile"] = ["Q2"] * len(Q2["Rank"])
        Q3["Quartile"] = ["Q3"] * len(Q3["Rank"])
        Q4["Quartile"] = ["Q4"] * len(Q4["Rank"])
        
        # # Save each quartile as a separate CSV file
        Q1.to_csv(os.path.join(save_path, str("Rank_" + group + "_" + data_col + "_Q1.csv")), index=False)
        Q2.to_csv(os.path.join(save_path, str("Rank_" + group + "_" + data_col + "_Q2.csv")), index=False)
        Q3.to_csv(os.path.join(save_path, str("Rank_" + group + "_" + data_col + "_Q3.csv")), index=False)
        Q4.to_csv(os.path.join(save_path, str("Rank_" + group + "_" + data_col + "_Q4.csv")), index=False)
        
        # Now let's prepare to plot the quartiles
        quart_list = [Q1, Q2, Q3, Q4]
        
        # Get the color for plotting the current group
        if group_col == "Layer": # if grouping by layer
            group_color = layer_colors[group] # use the layer-specific color
        elif colors != 0: # if other colors specified
            # Use the user specified colors (in order)
            group_color = colors[g]
        elif colors == 0: # no color specified by user
        # use the default colorblind seaborn pallet
            group_color = next(iter(sns.color_palette("colorblind")))
        
        # Create a violin plot with all 4 quartiles
        
        fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(10,4)) # plot with 4 subplots
        fig.subplots_adjust(hspace=1, wspace=0.7) # spacing between subplots
        sns.despine(fig=fig, top=True, right=True) # remove spines for aethetics
        
        title = str(group) + " " + str(data_col) + " By Quartile"
        plt.suptitle(title, y = 1.15, fontsize=24, fontweight="bold") # overall plot title
        # plt.ylabel(units, labelpad=10, fontsize=28, fontweight="bold")
        
        # counter to keep track of subplots
        i = 0
        # plot a violin plot for each quarant of data
        for quart in quart_list:
            # Plot the violin plot on current subplot
            sns.violinplot(y=data_col, data=quart, ax=axes[i], color = group_color)
            
            # Set plot and axes titles
            axes[i].set_title("Q" + str(i+1), y=1.05, fontsize=20, fontweight='bold') # subplot title
            axes[i].set_ylabel(None) # remove Y-axis label from the subplots. We'll add it back to the first plot later
            plt.setp(axes[i].collections, linewidth=3, edgecolor="k") # color of the violin border
            plt.setp(axes[i].lines, color='k', alpha=0.5) # color of inside lines
            
            # set the Y axis range based on the default or provided y_range
            axes[i].set_ylim(Q_cut[i])
            axes[i].tick_params(axis='y', labelsize=14) # set tick label size for Y-axis
            
            # increase subplot counter
            i = i+1
        
        axes[0].set_ylabel(units, labelpad=10, fontsize=18, fontweight="bold") # only need the first Y axis label
        
        # Save the figure as a PNG file
        plt.savefig(os.path.join(save_path, str(group + "_" + data_col + "_quart.png")), format='png', dpi=300, bbox_inches='tight')
        
        plt.close() # don't show the figure, but it will be saved
        
        # add the quartile data to a large dataframe with all three layers
        output_df.append(quart_list)
        
        # Add summary stats to the table
        summary_stats.loc["Mean", group] = round(rank_df[data_col].mean(), 3)
        summary_stats.loc["Median", group] = round(rank_df[data_col].median(),3)
        summary_stats.loc["N_mitos", group] = len(rank_df)
        summary_stats.loc["Min", group] = round(rank_df[data_col].min(),3)
        summary_stats.loc["Max", group] = round(rank_df[data_col].max(),3)
        
    # print the summary stats to the console
        
    print("\n", summary_stats)
        
    # Also save the summary stats as a CSV file
    summary_stats.to_csv(os.path.join(save_path, str("Summary_stats_" + data_col + ".csv")))
        
    return output_df # returns the quartile dataframes for each layer in a list
    # Layers should be returned in the order Basal, Proximal and Distal
    # If you want to capture each layer in a separate variables, call the function with 3 variables on the left
    # i.e. (basal, prox, distal) = data_quartiles()
    # You can call the function with just one variable, which will combine the three layers into a list, but it gets a little messy to unpack.
    
#%%% Function to create violin plots comparing WT and MCU KO in each layer. 
    
# Updated to perform the correct posthoc comparison
    
def WT_KO_violin(data, data_col, save_path, plot_type = "violin", group = "Genotype", y_range = [0, 0.5], units = 'Area (\u00B5m\u00B2)', title = 0, colors = ["#7F7F7F", "#FFFFFF"], cut = 1, stats = 0, SO=True, ctl_line = False):
    
    ''' This function will produce and save a violin plot for the given data comparing the groups you provide under "group" by layer. By default, the groups are Genotype, which would be the KO compared to control. If there is no grouping and you just want to plot the 3 layers, set "group" to be 0 when you call the function. The violins will be scaled by area so that all violins are the same total area.
    
    Required inputs are the dataframe with the data, the name of the data column, and the path to save the violin plot. Note that your dataframe should contain a column called "Layer", as that is what is going to be plotted on the X with your data column on the Y. 
    
    Optionally, you can change the group column you are comparing across, the range of the Y axis, the units of the Y-axis, the title of the plot and the colors of the groups. The default units are area in microns and the default range is from 0 to 0.5 microns. If you provide colors, make sure the number of colors matches the number of groups. If no grouping, the colors will apply to the three layers. I may add more customization later on.
    
    Added 8.31.23: This function will perform a standard 2-way ANOVA with Layer and Genotype as factors, then a post hoc comparison of cre + vs cre - within each layer. The significant results will be printed and corresponding significance bars will be added to the plots. The full statistics table will be saved as a CSV with the plot. To run these stats, set stats = 1 (default is to not run any stats and just print the plot).
    
    Added 9.20.23: Added a "ctl_line" option, which can be used to add a dotted line at 1 on the X-axis. Useful for normalized data. By default, no line is plotted. To add the line, set ctl_line = True.
    
    Updated 12.4.23: Updated to use the stats function ANOVA_posthoc() above with the corrected posthoc tests. 
    
    Updated 4.16.24: This now has the option to do a boxplot instead of violin plot. Set the variable plot_type to "box" to plot a box plot instead of violin. Violin is default.
    
    '''
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import pingouin as pg
    import pandas as pd
    import os
    from numpy import mean
    import numpy as np
    import warnings
    
    # Now lets make a violin plot comparing WT and KO mitochondria in each layer
    # This will be individual mito data
    
    # Make a palette out of the colors
    my_pal = sns.color_palette(colors)
    
    # First, if SO is False, remove it from the plot data
    # Note, SO will still be included in the stats
    if SO == False:
        plot_data = data.where(data["Layer"] != "Basal").dropna(how="all")
        plot_data["Layer"] = plot_data["Layer"].cat.remove_unused_categories()
    else:
        plot_data = data
        
    # First, if stat is set to 1 and the group is "Genotype", do the 2way ANOVA
    if stats == 1 and group == "Genotype":
        
        # Create a new column for grouping the given group and layer
        plot_data["Group"] = plot_data['Genotype'] + "_" + plot_data['Layer'].astype('str')
        
        # test if the data is normal
        if len(data[data_col]) < 5000:
            test = "shapiro" # will run a Shapiro-Wilk test for normality
        else:
            test = "normaltest" # will run the omnibus test for normality (better for large N)
            
            
        norm_test = pg.normality(data[data_col], method = test)
        norm = norm_test["normal"][0] # boolean for normality, true if normal
        
        # custom function to perform two-way standard ANOVA across genotype and layer with posthoc tests
        # Will run parametric posthocs if normal data, non-parametric if not
        ANOVA, posthoc = ANOVA_posthoc(data, data_col, save_path, para = norm)
    
        # Pick out which statistics are significant so we can add the stats to the plot
        # First check the overall ANOVA. If the overall ANOVA is significant, we can run the posthoc
        overall_ANOVA_sig = pd.DataFrame(columns = ["Factor", "P-val", "Sig"]) # empty dataframe to add the overall significance
        overall_ANOVA_sig["Sig"] = overall_ANOVA_sig["Sig"].astype("bool") # change datatype to boolean
        posthoc_sig = pd.DataFrame(columns = ["Group A", "Group B", "P-val", "Sig", "Test", "N Group A", "N Group B"]) # empty dataframe to add the multiple comparison sig
        posthoc_sig["Sig"] = posthoc_sig["Sig"].astype("bool")
    
        # determine whether overall anova is significant for each factor
        for row in ANOVA.iterrows():
            factor = row[1]["Source"] # the factor being compared
            p_val = row[1]["p-unc"] # the overall ANOVA p-value
            
            if p_val <= 0.05: # only pull out the stats that are significant
                overall_sig = True
                curr_stats = pd.DataFrame({"Factor": [factor], "P-val": [p_val], "Sig": [overall_sig]}) # pull relevant stats from overall ANOVA
                overall_ANOVA_sig = pd.concat(objs=[overall_ANOVA_sig, curr_stats], ignore_index=True) # add to the overall_ANOVA_sig dataframe
        
        # print out the overall ANOVA results
        print("\n----------------------------------\n  " + data_col + " Two-way ANOVA\n----------------------------------\n")
        
        if overall_ANOVA_sig.empty:   # if there's no significant results     
            print("Overall ANOVA was not significant for either factor.")
        else:
            print(overall_ANOVA_sig)
        
        # if genotype is significant, look at the post hoc comparing genotype across each layer
        if "Genotype" in list(overall_ANOVA_sig["Factor"]):
            for i in range(len(posthoc)): # loop through all the comparisons from the post hoc 
                groupA = posthoc.loc[i]["A"] # first group
                groupB = posthoc.loc[i]["B"] # second group
                # get the N of group A and group B from the data
                N_A = data["Group"][data["Group"]==groupA].count()
                N_B = data["Group"][data["Group"]==groupB].count()
                p_val_corr = posthoc.loc[i]["p-corr"] # corrected p-value for the comparison
                if norm:
                    test = "Unpaired t-test"
                else:
                    test = "Mann-Whit"
                
                # determine if any comparisons on the plot are significant
                # narrow down just our comparisons of same layer across genotype
                if groupA.split("_")[1] == groupB.split("_")[1]: # only comparing same layer
                    if p_val_corr <= 0.05: # if significant
                        sig_comp = True
                        comp_stats = pd.DataFrame({"Group A": [groupA], "Group B": [groupB], "P-val": [p_val_corr], "Sig": [sig_comp], "Test": [test], "N Group A":N_A, "N Group B":N_B}) # pull relevant stats from any significant posthoc comparisons
                        posthoc_sig = pd.concat(objs=[posthoc_sig, comp_stats], ignore_index=True) # add to posthoc dataframe
    
            # print out the posthoc results
            if posthoc_sig.empty:
                print("\nNone of the post hoc comparisons between genotype were significant.\n")
            else:
                print("\nMultiple comparisons across genotype:\n")
                print(posthoc_sig)
                
        
    # Now we can make the plot!
    if SO and group != 0: # need a bigger figure size
        fig_size = (5,3)
    else:
        fig_size = (3,4) # plot will be shorter because only 2 layers or there's no grouping
    fig, ax = plt.subplots(figsize=fig_size) # size of plot
    # fig.subplots_adjust(hspace=1, wspace=0.7) # spacing between subplots
    sns.despine(fig=fig, top=True, right=True) # remove spines for aethetics

    if title == 0: # if no title provided
        title = ""
    else: 
        title = title + ": " + str(data_col)
        
    if stats == 1:
        plt.suptitle(title, y = 1.2, fontsize=24, fontweight="bold") # title needs to be moved up a little to fit the stats
    else: 
        plt.suptitle(title, y = 1.05, fontsize=24, fontweight="bold")

    if group == 0 and plot_type == "violin": # no grouping, plot violin plot
        try: # catch incase columns are missing or named improperly
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sns.violinplot(y=data_col, x = "Layer", hue="Layer", data=plot_data, palette = my_pal, density_norm = "area", inner = "quartile", cut = cut, bw_method=0.1, alpha=0.95) # create a violin plot with just layer
        except ValueError:
            print("Warning: The specified columns '" + str(data_col) + "' and/or 'Layer' are not in the dataframe provided. Please make sure the columns in your dataframe are present and named as shown.")
            pass
    elif group != 0 and plot_type == "violin": # if the data is grouped
        # Create a new column for grouping the given group and layer
        try: # catch
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sns.violinplot(y=data_col, x = "Layer", data=plot_data, hue = group, palette = my_pal, density_norm = "area", inner = "quartile", cut = cut, bw_method=0.1, alpha=0.95) # create the violin plot of layer grouped by the specified group
        except ValueError:
            print("Warning: The columns '" + str(data_col) + "', 'Layer' and '" + str(group) + "' must be in your dataframe. Please make sure your columns and inputs are named correctly. By default, the violin plot is grouped by Layer and by Genotype. To group by something other than genotype, please specify the column name under the 'group' variable. If you only want to group the plot by layer and not any additional grouping, please set 'group = 0' when you call the function.")
            pass
    elif group == 0 and plot_type == "box":
        try: # catch incase columns are missing or named improperly
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sns.boxplot(y=data_col, x = "Layer", hue = "Layer", data=plot_data, palette = my_pal, linewidth=2, linecolor="k", whis=(0, 100)) # create a box plot by layer
        except ValueError:
            print("Warning: The specified columns '" + str(data_col) + "' and/or 'Layer' are not in the dataframe provided. Please make sure the columns in your dataframe are present and named as shown.")
            pass
    elif group != 0 and plot_type == "box": # if the data is grouped
        try: # catch
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sns.boxplot(y=data_col, x = "Layer", data=plot_data, hue = group, palette = my_pal, linewidth=2, linecolor="k", whis=(0, 100)) # create the box plot of layer grouped by the specified group
        except ValueError:
            print("Warning: The columns '" + str(data_col) + "', 'Layer' and '" + str(group) + "' must be in your dataframe. Please make sure your columns and inputs are named correctly. By default, the box plot is grouped by Layer and by Genotype. To group by something other than genotype, please specify the grouping column name under the 'group' variable. If you only want to group the plot by layer and not any additional grouping, please set 'group = 0' when you call the function.")
            pass
    
    # Format the plot and axes
    if plot_type == "violin":
        plt.setp(ax.collections, linewidth=2, edgecolor="k", alpha=.75) # color of the violin border
        plt.setp(ax.lines, color='k', linestyle=(0, (4, 1)), linewidth=1.5) # color and style of inside lines
        plt.setp(ax.lines[1:17:3], color='k', linestyle="solid", linewidth=3)
        
    # set the Y axis range based on the default or provided y_range
    # if stats == 1:
    #     y_range[1] += 0.04 # extend y range to account for adding stats bars
    ax.set_ylim(y_range)
    ax.tick_params(axis='y', labelsize=16) # set tick label size for Y-axis
    ax.tick_params(axis='x', labelsize=20)
    ax.set_ylabel(units, labelpad=10, fontsize=18, fontweight="bold") # set the y axis label with the given units
    ax.set_xlabel(None) # set the x axis label with the given units
    ax.set_xticks(ticks = np.arange(3), labels=["SO", "SR", "SLM"], fontweight="bold") # remove x-axis label
    
    # if an ANOVA was run and any post hocs were significant, add sig bars to the plot  
    if stats == 1 and not posthoc_sig.empty:
        if SO: # if including SO
        # find which violins to draw significance between on X-axis (0-1 Basal, 2-3 Prox, 4-5 Distal)
            for comp in posthoc_sig["Group A"]: # loop over each sig posthoc comparison to get the layer
                layer = comp.split("_")[1] # get the layer being compared
                if layer == "Basal":
                    x1 = -0.25
                    x2 = 0.25
                if layer == "Proximal":
                    x1 = 0.75
                    x2 = 1.25
                if layer == "Distal":
                    x1 = 1.75
                    x2 = 2.25
                    
                # definte the placement of the bar at +2 above the data max with tails of length 2
                y, h = max(y_range) + .02, .02
                
                # plot the bar for the current statistical comparison
                ymin,ymax = plt.ylim()
                xmin,xmax = plt.xlim()
                # line=Line2D([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                # ax.add_line(line)
                # Get number of stars based on significance
                layer_p = posthoc_sig[posthoc_sig["Group A"].str.contains(layer)]["P-val"].reset_index(drop=True)[0]
                
                if layer_p < 0.0001:
                    sig_text = "****"
                elif layer_p < 0.001:
                    sig_text = "***"
                elif layer_p < 0.01:
                    sig_text = "**"
                elif layer_p < 0.05:
                    sig_text = "*"
                    
                if plot_type == "violin":
                    ax.text((x1+x2)/2, y+h, sig_text, ha='center', va='bottom', color="k", fontsize=30)
                else:
                    ax.text((x1+x2)/2, y-h, sig_text, ha='center', va='bottom', color="k", fontsize=30)
                    
        elif SO == False:
        # find which violins to draw significance between on X-axis (0-1 Basal, 2-3 Prox, 4-5 Distal)
            for comp in posthoc_sig["Group A"]: # loop over each significant posthoc comparison
                layer = comp.split("_")[1] # get the layer being compared
                if layer == "Basal":
                    pass
                if layer == "Proximal":
                    x1 = -0.25
                    x2 = 0.25
                if layer == "Distal":
                    x1 = 0.75
                    x2 = 1.25
                    
                try: # definte the placement of the bar at +2 above the data max with tails of length 2
                    y, h = max(y_range) + .02, .02
                    
                    # plot the bar for the current statistical comparison
                    ymin,ymax = plt.ylim()
                    xmin,xmax = plt.xlim()
                    # Create a new axis to plot the bars on
                    # ax2 = plt.axes([xmin,ymin,xmax,ymax+0.04], facecolor=(1,1,1,0))
                    line=Line2D([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                    ax.add_line(line)
                    
                    # Determine the level of significance to print on the plot
                    layer_p = posthoc_sig[posthoc_sig["Group A"].str.contains(layer)]["P-val"].reset_index(drop=True)[0]
                    
                    if layer_p < 0.0001:
                        sig_text = "****"
                    elif layer_p < 0.001:
                        sig_text = "***"
                    elif layer_p < 0.01:
                        sig_text = "**"
                    elif layer_p < 0.05:
                        sig_text = "*"
                        
                    ax.text((x1+x2)/2, y+h, sig_text, ha='center', va='bottom', color="k", fontsize=18)
                except UnboundLocalError:
                    pass
                
    # plot a line at Y= 1 if ctl_line is true
    if ctl_line:
        
        # remove the median lines from the violins
        # ax.lines[1].remove()
        # ax.lines[3].remove()
        # ax.lines[5].remove()
        
        # Add line at Y=1 for the control
        ax.axhline(y=1, color='k', linestyle=':', lw=2, zorder=-1)
        
        # add lines for the means
        sns.pointplot(x = "Layer", y=data_col, data=plot_data, estimator=mean, color ="k", join=False, scale=0.75)
    
    
    # Remove legend if there's no grouping. The X-axis will say the layers
    if group == 0:
        ax.legend("", frameon=False) # remove legend
        ax.set_xlabel(None) # remove x-axis label
        ax.set_xticks(ticks = np.arange(3), labels=["SO", "SR", "SLM"]) # remove x-axis label
        ax.set_ylabel(units, labelpad=10, fontsize=18, fontweight="bold") # set the y-axis with the units provided. Default units is "Area (um2)"
        
    else: # move and resize the legend
    
        from matplotlib.collections import PolyCollection
    
        plt.legend(title=group, fontsize='12', title_fontsize='14', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        
        # individually change violin colors in the order of the colors specified in the function call
        if len(colors) > 2:
            for ind, violin in enumerate(ax.findobj(PolyCollection)):
                rgb = colors[ind]
                violin.set_facecolor(rgb)

    # Save the figure as a PNG file
    if group == 0: # if no further grouping
        plt.savefig(os.path.join(save_path, str(data_col + "_by_layer_violin.png")), format='png', dpi=300, bbox_inches='tight')
    else: # if a group was used, include it in the filename
        plt.savefig(os.path.join(save_path, str(data_col + "_by_" + str(group) + "_violin.png")), format='png', dpi=300, bbox_inches='tight')
    
        
        
#%% Prism formatting function prism_format()

# Function to format a table for prism, remove any NaN values, and save the CSV
def prism_format(data, data_col, file_path, col = "Layer", row = 0, sort = 0):
    '''Quick little function to format a table for prism with the image averages for your metric of interest. Unless otherwise specified, the columns of the table will be the dendritic layers and the rows will be each tile or ROI. The function will remove NaNs from the table and save the table as a csv file using the filename you input into the function.
    
Please input into the function the dataframe with the data, the column with the metric of interest, and the full file path where you wish to save the file, including the file name (with .csv at the end). All three of these inputs are required.

Optionally, if you want the columns of the table to be something other than dendritic layer, please specify this as "col =" and then the column name of the factor you would like for the columns. You can also specify what you want the rows

You can specify sort = "Sort Col" if you would like to sort the table by a specific column. Otherwise, the table will be sorted in the order that the dataframe is in. '''

    # import relevant packages
    import pandas as pd
    
    # Format the table by pivoting the data using the provided column. 
    # If a column name is provided for "col", this will be used instead of "layer" as the columns of the formatted table.
    # If a sort list is provided for "sort", then use that to sort the resulting table
    if row == 0:
        data_tb = pd.pivot(data, values=data_col, columns=col).reset_index(drop=True)
    else: 
        data_tb = pd.pivot(data, values=data_col, columns=col, index=row).reset_index()
        
    # drop the NANs in the middle of each column
    for c in data_tb:
         column = data_tb[c].dropna().reset_index(drop=True) # remove the NaNs
         data_tb[c] = column
                
    # # trim any trailing NaNs
    data_tb = data_tb.dropna(how="all")
    
    if sort != 0: # may not actually work to sort
        data_tb.sort_index(axis= 1, level = sort) # sort the columns by the specified column data
    
    # Save the table as a CSV file using the provided filepath
    data_tb.to_csv(file_path, index=False)
    
    return data_tb

#%%% Statistical summary function 

import os

def get_stats_summary(data, data_col, x_label, save_dir = os.getcwd(), hue_col = "Layer", colors = ['twilight', 'teal', 'dark peach'], group_name = [], hist_comp = "count", binsize = 0.005, figsize = (8,5), xlim = [0,0.35], ylim = [0,120]):
    
    ''' Function to create a nice histogram of the resampled median populations comparing across the "hue" column, usually Layer. Each value of hue will be plotted as a separate line on the histogram. The histogram image will be saved in the save_directory
    
 Required inputs:
    data = long format dataframe with a column for the data to plot and a column for the hue values.
    data_col = column name for the column to be plotted in the histogram
    x_label = Specify the x axis label for the histogram(s) (include units)
    
 Optional inputs:
     save_dir = The location to save the plot as a PNG file. Must be a full path.
     hue_column = The column name to separate the hues. Each hue will be plotted as a separate line. By default, the hue column is "Layer"
     colors = Specify a list of colors the length of hues in hue_col (use XKCD color names: https://klaash.github.io/xkcdcolorpicker)
     group_name = If desired, specify a group name for the histogram (ie. CTL or cKO). This will be added to the plot title and the filename to distinguish it from other groups. Recommended to do this if you'll have multiple histograms of the same metric in the same locaion, otherwise if the filename is the same this code may save over the previous file.
     hist_comp = define how the histogram is calculated. Default is "count" but you can change this to "proportion" if you want to normalize the bars to sum to 1. (ie. if count is different between groups)
     binsize = specify the size of the bins. Default is 0.005.
     figsize = Specify the figure size in two dimensions. Width then height. Default is 8x5 inches.
     xlim and ylim: Specify the x or y range of the plot. Provide two numbers. Default: xlim = [0,0.35], ylim = [0,120]
     
    '''
    import seaborn as sns
    import matplotlib.pyplot as plt
    import os
    from scipy import stats
    import pandas as pd
        
    # run normality test on the data
    if len(data[data_col]) < 5000:
        (norm_stat, p_val) = stats.shapiro(data[data_col])
        # define normality
        if p_val > 0.05:
            normal = True
        elif p_val <= 0.05:
            normal = False
        test = "Shapiro-Wilk"
    else:
        result = stats.anderson(data[data_col])
        norm_stat = result.statistic
        crit_value = result.critical_values[2] # critical value at 5% significance level
        # define normality
        if norm_stat > crit_value:
            normal = False
        else:
            normal = True
        test = "Anderson-Darling"
    
    # Plotting the histogram
    sns.set(style="ticks")
    
    # set the palette of colors
    pal= sns.xkcd_palette(colors)
    
    # list of hues in proper order:
    try:
        hue_list = list(data[hue_col].cat.categories) # keep the order of the categories
    except AttributeError: # to catch the error if the hue column is not categorical
        hue_list = list(data[hue_col].unique()) # keep the original order in hue col from the dataframe
    
    plt.figure(figsize=figsize)
    plt.tight_layout()
    sns.histplot(data=data, x=data_col, element="bars", fill=True, kde=True, hue=hue_col, hue_order= hue_list, stat=hist_comp, palette = pal, alpha=0.2, common_norm=False, binwidth=binsize)
    
    # Manually set legend labels
    # legend_labels = data[hue_col].unique()
    # plt.legend(title=hue_col, labels=legend_labels)
    
    if len(group_name) > 0:
        plt.title(f"Distribution of {data_col} in {group_name}")
    else:
        plt.title(f"Distribution of {data_col}")
    plt.xlabel(x_label)
    plt.ylabel(hist_comp)
    plt.xlim(xlim)
    plt.ylim(ylim)
    
    # add a horizontal line for the mean or median of each group
    # get the list of hues in the hue column
    
    for hue in hue_list:
        ind = hue_list.index(hue)
        # get the mean or median of the curren group
        if normal:
            hue_avg = data[data[hue_col] == hue][data_col].mean()
            avg = "mean"
            # maybe add line to legend stating mean or median
        else:
            hue_avg = data[data[hue_col] == hue][data_col].median()
            avg = "median"
            
        # plot a line at the mean or median
        plt.axvline(hue_avg, linestyle="--", color = pal[ind])
        
    
    # save the plot
    if len(group_name) > 0:
        plt.savefig(os.path.join(save_dir, f"hist_{data_col}_{avg}_{group_name}.png"))
    else:
        plt.savefig(os.path.join(save_dir, f"hist_{data_col}_{avg}.png"))
        
    # print out the results of the Shapiro Wilk test to the console
    if normal: 
        print(f"\nDistribution of {data_col} is normal.")
    else: 
        print(f"\nDistribution of {data_col} is not normal.")
        
    if test == "Shapiro-Wilk":
        print(f"{test}, N = {len(data[data_col])}, P = {round(p_val,3)}")
    elif test == "Anderson-Darling":
        print(f"{test}, N = {len(data[data_col])}, Statitic = {round(norm_stat,1)}, Critical Value = {crit_value}")
                
    # now let's create a summary stats table with the mean or median, std and N for each group in "hue" 
    summary_tb = round(data[[hue_col,data_col]].groupby(hue_col).describe(),3)
    summary_tb_full = summary_tb[data_col].T
    
    # subset just the most important stats to return
    if normal:
        summary_tb_sub = round(summary_tb[data_col][["mean", "std", "count"]], 3)
    else:
        summary_tb_sub = round(summary_tb[data_col][["50%", "std", "count"]], 3).rename(columns={"50%":"median"})
    # change count to integer
    summary_tb_sub["count"] = summary_tb_sub["count"].astype("int")
        
    print(f"\n Summary of {data_col} by {hue_col} in {group_name}:\n\n", summary_tb_sub)

    # save the full stats table as a CSV file
    if len(group_name) > 0:
        summary_tb_full.reset_index().to_csv(os.path.join(save_dir, f"{group_name}_{data_col}_summary_stats.csv"), index=False)
    else:
        summary_tb_full.reset_index().to_csv(os.path.join(save_dir, f"{data_col}_summary_stats.csv"), index=False)
        
    return summary_tb_sub.reset_index(), normal
    