import pandas as pd
"""Script to combine the region data to variable data and display intersection."""
import timeit
import math

start_time = timeit.default_timer()


def conv_to_string(gene_list, index):
    """Pull out all the values in gene_list at the index
        'index' and joins them with a comma. The resulting string is
        returned.

        gene_list: list
        index: int
    """
    return ",".join([str(item[index]) for item in gene_list])


def conv_to_string_NVR(gene_list, index, start_val, stop_val):
    """Pull out all the values in gene_list at the index
        'index', finds the percentage voerlap and joins
        them with a comma. The resulting string is returned.

        gene_list: list
        index: int
        start_val: int
        stop_val: int
    """

    rng = stop_val - start_val
    temp_list = [item[index] / float(rng) * 100 for item in gene_list]
    return_list = [0 if math.isnan(x) else x for x in temp_list]
    return ','.join(str(e) for e in return_list)

def conv_to_string_gene(gene_list, index, index2):
    """Pull out all the values in gene_list at the index
        'index', finds the percentage voerlap and joins 
        them with a comma. The resulting string is returned.
        
        gene_list: list
        index: int
    """
    temp_list = [item[index] / float(item[index2]) * 100 for item in gene_list]
    return_list =  [0 if math.isnan(x) else x for x in temp_list]
    return ','.join(str(e) for e in return_list)



def add_to_region(list_of_genes, data_frame, start, stop, type_overlap):
    """
    Append dataframe results ot the list in the format below
    """
    for i in xrange(len(data_frame)):
        diff_gene = data_frame.iloc[i,2] -  data_frame.iloc[i,1] + 1
        if type_overlap in [1,4]:
            diff_region = data_frame.iloc[i,2] - start + 1
        elif type_overlap in [2,3,5,6]:
            diff_region = stop - start + 1
        elif type_overlap in [8,9]:
            diff_region = stop - data_frame.iloc[i,1] + 1
        else:
            diff_region = data_frame.iloc[i,2]- data_frame.iloc[i,1] + 1
        list_of_genes.append([data_frame.iloc[i,1], data_frame.iloc[i,2],
                              type_overlap, diff_region, diff_gene])
#         print type(diff_region)
#         print type(diff_gene)
    return list_of_genes

# read the two files (the first is the variable file and the second is the variable file)
df = pd.read_csv('variable_data.txt', sep='\t') # Variable data
df1 = pd.read_csv('region_data.txt', sep='\t') # Region data

# list of choromozomes and thier sizes (sizes irrlevant for this task)
all_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
           '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

# iterate through each region data range and incrementally add the variables that match
with open('results_file.txt', 'ab') as myfile:

    # add the headers to the file
    string = 'chr' + '\t' + 'Region_start' + '\t' + 'Region_stop' + '\t' + 'Overlapping_Variable_Start' + '\t' + \
        'Overlapping_Variable_Stop' + '\t' + 'Overlap_Type' + '\t' + \
        'Number_of_Overlapping_Variable' + '\n'
    myfile.write(string)

    for chrm in all_chr:

        # choose chromosomes only in all_chr list since there maybe others
        to_keep = [chrm]

        # only keep rows that contain the string in to_keep
        df_variable = df[df['chr'].isin(to_keep)]
#         df_variable = df_variable.sort_values(by='start')

        df_nvr = df1[df1['chr'].isin(to_keep)]
        df_nvr = df_nvr.sort_values(by='start')

        # initalize variable_overlapping_regions to store start, stop and the overlap type
        variable_overlapping_regions = []
        
        for i in xrange(len(df_nvr)):

            # for each row get the start and stop point for each line in region_data 

            region_start = df_nvr.iloc[i, 1] # region data start point
            region_stop = df_nvr.iloc[i, 2] #  region data stop point
            
            df_result = df_variable[(df_variable['start'] < region_start) & (df_variable['stop'] < region_stop) & 
                                (df_variable['stop'] > region_start)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 1)
            df_result = df_variable[(df_variable['start'] < region_start) & (df_variable['stop'] == region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 2)
            df_result = df_variable[(df_variable['start'] < region_start) & (df_variable['stop'] > region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 3)
            df_result = df_variable[(df_variable['start'] == region_start) & (df_variable['stop'] < region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 4)
            df_result = df_variable[(df_variable['start'] == region_start) & (df_variable['stop'] == region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 5)
            df_result = df_variable[(df_variable['start'] == region_start) & (df_variable['stop'] > region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 6)
            df_result = df_variable[(df_variable['start'] > region_start) & (df_variable['stop'] < region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 7)
            df_result = df_variable[(df_variable['start'] > region_start) & (df_variable['stop'] == region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 8)
            df_result = df_variable[(df_variable['start'] > region_start) & (df_variable['stop'] > region_stop) &
                                (df_variable['start'] < region_stop)]
            variable_overlapping_regions = add_to_region(variable_overlapping_regions, df_result, region_start, region_stop, 9)

            # if there is an overlap use conv_to_string to convert each item in the
            # list to a comma separated string, else add '-' in the respective
            # columns
            if len(variable_overlapping_regions) > 0:
                string = chrm + '\t' + str(region_start) + '\t' + str(region_stop) + '\t' + \
                conv_to_string(variable_overlapping_regions, 0) + '\t' + \
                conv_to_string(variable_overlapping_regions, 1) + '\t' + \
                conv_to_string(variable_overlapping_regions, 2) + '\t' + \
                str(len(variable_overlapping_regions)) + '\n'
            else:
                string = chrm + '\t' + str(region_start) + '\t' + str(region_stop) + '\t' + '-' + '\t' + '-' + '\t' + \
                '-' + '\t' + str(0) +  '\n'
            myfile.write(string)

            variable_overlapping_regions = []

        print "finished " + str(chrm)

stop_time = timeit.default_timer()
print "TIME  finished"
print stop_time - start_time

