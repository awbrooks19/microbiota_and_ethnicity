#!/usr/bin/env python

""" 
Author: Andrew W. Brooks
Affiliation: Vanderbilt Genetics Institute
Date: October 27, 2016

Released under MIT License
Copyright (c) 2017 Andrew W. Brooks
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, 
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to 
permit persons to whom the Software is furnished to do so, subject to the following 
conditions: The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software. The software is provided "as is",
without warranty of any kind, express or implied, including but not limited to the 
warranties of merchantability, fitness for a particular purpose and noninfringement. 
In no event shall the authors or copyright holders be liable for any claim, damages 
or other liability, whether in an action of contract, tort or otherwise, arising from,
out of or in connection with the software or the use or other dealings in the software.
"""

### IMPORT REQUIRED ###
import argparse # Parse Command Line Arguments to Script

### FUNCTIONALITY ###
import os # Tool for terminal and operating system type calls
import glob # Tool to Regex Search for Files
import copy
import random       # Generate Random Values
random.seed(54321)  # Set Random Seed for Reproducibility
import itertools # Iterate through data

### DATA STRUCTURES ###
import pandas as pd # Pandas Dataframes
import biom # Biom Format Used for Microbiome Analyses
from skbio.tree import TreeNode # SciKit - Bio TreeNode object to Store Phylogeny
from skbio.stats.composition import ancom # Import ANCOM Test of Observation Differential Abundance 
import ete3 as ete # Ete3 tree viewing tools
from skbio.stats.distance import DistanceMatrix # SciKit-Bio DistanceMatrix Object

### TOOLKITS ###
import numpy as np # Numpy Numerical Toolkit
import skbio as sk # SciKit-Bio Toolkit
import scipy as sp  # Scipy Scientific Toolkit
from scipy.cluster.hierarchy import linkage # Scipy Tree Linkage Function
import statsmodels.api as sm # StatsModels Statsistics Toolkit
from statsmodels.sandbox.stats.multicomp import multipletests
import patsy # Patsy format Dataframe for Regression
import vcf # Import PyVCF Tools for working with Variant Call Format (VCF) SNP Data
from IPython.display import display

### PLOTTING ###
import matplotlib.pyplot as plt # PyPlot for General Plotting
import seaborn as sns # Seaborn plotting with Pandas
#%matplotlib inline # Run Matplotlib PyPlot in Jupyter Notebook

##################################################################################
######################### GENERAL & UI FUNCTIONS #################################
##################################################################################

##################################################################################
### FUNCTION - Write String to File ###
def file_write(string, path): f = open(path,'w'); f.write(string); f.close()
        
##################################################################################
### FUNCTION - Append String to File ###
def file_append(string, path): f = open(path,'a'); f.write(string); f.close()
    
##################################################################################
### FUNCTION - Print a List of Strings Line By Line to the Terminal ###
def list_print(listIn): 
    for listValIn in listIn: print(listValIn)
    
##################################################################################
### FUNCTION - Write 'w' or Append 'a' A List of Strings Line By Line to a File ###
def list_write(listIn, filePathIn, addBlankTail=False, writeAppend='a'):
    fOut = open(filePathIn, writeAppend)
    for listValIn in listIn: fOut.write((listValIn+"\n"))
    if addBlankTail == True: fOut.write("\n")
    fOut.close()

##################################################################################
### FUNCTION - Make Directory if it Doesn't Exist ###
def dir_make(dirPath):
    dirPath = dirPath.replace('//','/') # Replace // with / in path
    if not os.path.isdir(dirPath): os.makedirs(dirPath)

##################################################################################
### FUNCTION - Convert Number to String with Specified Decimal Places ###
def convert_num_to_str(val, decimals='0'): return ('{:.'+decimals+'f}').format(val)

##################################################################################
### FUNCTION - Convert Iterable of Numbers to List of Strings with Specified Decimal Places ###
def convert_nums_to_strs(vals, decimals='0'): return [('{:.'+decimals+'f}').format(x) for x in vals]

##################################################################################
### Function: Query User for List of Files or Directories Specified by Regex (*.txt or */)
def ui_input():
    ### Get User Input ###
    idx = input('   - Enter Value: ')
    # Try Coercing to Int and Float #
    if '.' in idx:
        try: idx = float(idx); print('   - Coerced to Float: '+str(idx))
        except: pass
    else:
        try: idx = int(idx); print('   - Coerced to Int: '+str(idx))
        except: pass
    print('   - Returning: '+str(idx))
    return idx
#ui_input()

##################################################################################
### FUNCTION - Ask User for True False Value ###
def ui_tf(yesno=False):
    # Print Prompt and Query User #
    if yesno == False:
        print('\n   - Enter True or False -')
        print('    (t/f, T/F, true/false, True/False)')
    elif yesno == True: 
        print('   - Enter Yes or No -')
        print('    (y/n, Y/N, yes/no, Yes/No)')
    else: 
        print('   - Error: yesno value invalid, must be True or False -')
        return None
    x = input('   - Entry: ')
    # Evaluate Input #
    trueVals = ['t','T','true','True','y','yes','Y','Yes']
    falseVals = ['f','F','false','False','n','no','N','No']
    if x in trueVals: 
        if yesno == False: print('   - Returning True - ') 
        elif yesno == True: print('   - Returning Yes - ')
        return True
    elif x in falseVals: 
        if yesno == False: print('   - Returning False - ') 
        elif yesno == True: print('   - Returning No - ')
        return False    
    else:
        print('   - Error: Invalid Value, Please Try Again...\n')
        return ui_tf(yesno=yesno)
#ui_tf(yesno=True)

##################################################################################
### Function: Query User for List of Files or Directories Specified by Regex (*.txt or */)
def ui_choose_list(listIn):
    print('   - Choose Index of Desired Location -')
    for curIdx, cur in enumerate(listIn):
        print('     '+str(curIdx)+' : '+str(cur))
    idx = int(input('   - Enter Integer Index: '))
    return listIn[idx]
#ui_choose_list([1,2,3,4])

##################################################################################
### Function: Query User for List of Files or Directories Specified by Regex (*.txt or */)
def ui_choose_regex(regexStr):
    print('   - Choose Index of Desired Location -')
    for curIdx, cur in enumerate(glob.glob(regexStr)):
        print('     '+str(curIdx)+' '+cur)
    idx = int(input('   - Enter Integer Index: '))
    return glob.glob(regexStr)[idx]
#ui_choose_regex('*.ipynb')


##################################################################################
##################### PANDAS & STATISTICS FUNCTIONS ##############################
##################################################################################


##################################################################################
### FUNCTION - Scale a list of numbers proportionally to a new min and max range ###
def stats_rescale(listIn, newMin, newMax, oldMin=None, oldMax=None):
    if (oldMin == None): oldMin = min(listIn)
    if (oldMax == None): oldMax = max(listIn)
    OldRange = (oldMax - oldMin)
    NewRange = (newMax - newMin)
    return [(((OldValue - min(listIn) * NewRange) / OldRange) + newMin) for OldValue in listIn]

##################################################################################
### FUNCTION - Input a List of P values and return list of FDR Benjamin corrected P values ###
def stats_fdr(pVals, alpha=0.10):
    return multipletests(pVals, alpha=0.05, method='fdr_bh')[1]

##################################################################################
### FUNCTION - Calculate Kruskal Wallis ###
# Pass in a Pandas Dataframe, and the names of a continuous and categorical (nominal) column
def stats_kw(df, conCol, catCol, fullStats=True, printOut=True):
    kwStat, kwP = sp.stats.mstats.kruskalwallis(*[np.array(df[df[catCol] == x][conCol]) for x in df[catCol].unique()])
    if printOut == True:
        print('   - Kruskal Wallis Test of '+conCol+' by '+catCol+' -')
        print('     - Significance P   : '+str(kwP))
        print('     - Test Statistic   : '+str(kwStat))
    # Compute Full Stats #
    if fullStats == True:
        kwDegFreedom = len(df[catCol].unique())-1
        kwRanks = sp.stats.rankdata(df[conCol])
        kwCorrection = sp.stats.tiecorrect(kwRanks)
        kwUncorrectedStat = (1 - kwCorrection)*kwStat
        kwUncorrectedP = 1 - sp.stats.chi2.cdf(kwUncorrectedStat, 3)
        if printOut == True:
            print('     - Degree Freedom   : '+str(kwDegFreedom))
            print('     - Correction Div   : '+str(1 - kwCorrection))
            print('     - Uncorrected P    : '+str(kwUncorrectedP))
            print('     - Uncorrected Stat : '+str(kwUncorrectedStat))
        # Return complete stats
        return kwP, kwStat, kwDegFreedom, kwUncorrectedP, kwUncorrectedStat 
    # Return Minimum Stats #
    return kwStat, kwP
#kwRes = stats_kw(varDict['mp'], conCol='bmi', catCol='race')

##################################################################################
### FUNCTION - Calculate Welch's T-test ###
# Pass in a Pandas Dataframe, and the names of a continuous and categorical (nominal) column
def stats_welch(list1, list2):
    welchStat, welchP = sp.stats.ttest_ind(list1, list2, axis=0, equal_var=False, nan_policy='propagate')
    return welchP, welchStat
#stats_welch(list1=xOut['Mean_decrease_in_accuracy'], list2=xOut['Standard_deviation'])

##################################################################################
### FUNCTION - Calculate Mann-Whitney-U ###
# Pass in a Pandas Dataframe, and the names of a continuous and categorical (nominal) column
def stats_mwu(df, conCol, catCol, altTest='two-sided', printOut=False):
    allOut = pd.DataFrame(columns=['g0', 'g1', 'p', 'stat', 'enough_obs', 'g0_count', 'g1_count'])
    ### For Each Combination Perform Mann-Whitney-U ###
    for groupIdx, groupCombination in enumerate(itertools.combinations(df[catCol].unique(),2)):
        if printOut == True:
            print('   - Performing Mann-Whitney-U between '+str(groupCombination[0])+' and '+str(groupCombination[1]))
        try:
            # Get Lists of Values for Each Group #
            list1 = df.loc[df[catCol]==groupCombination[0],[conCol]]
            list2 = df.loc[df[catCol]==groupCombination[1],[conCol]]    
            # Calculate MWU #
            mannStat, mannP = sp.stats.mannwhitneyu(list1, list2, use_continuity=True, alternative=altTest)
            # Check for Enough Values - minimum 20#
            if (len(list1)<20) or (len(list2)<20):
                if printOut == True: 
                    if len(list1)<20: print('     - Warning: Less than 20 observations in '+groupCombination[0])
                    if len(list2)<20: print('     - Warning: Less than 20 observations in '+groupCombination[1])
                allOut.loc[groupIdx] = groupCombination[0],groupCombination[1], mannP, mannStat, False, len(list1), len(list2)
            else:
                allOut.loc[groupIdx] = groupCombination[0],groupCombination[1], mannP, mannStat, True, len(list1), len(list2)         
                if printOut == True:
                    print('     - Significance P   :'+str(mannP))
                    print('     - Test Statistic   : '+str(mannStat))
        # If MWU Errors out 
        except NameError:
            if printOut == True: print('     - Error Mann-Whitney-U: Could not compare '+str(groupCombination[0])+' and '+str(groupCombination[1]))
            allOut.loc[groupIdx] = groupCombination[0],groupCombination[1], np.NaN, np.NaN, False, np.NaN, np.NaN 
            
    return allOut
#mwuRes = stats_mwu(varDict['mp'], conCol='bmi', catCol='race', printOut=False)

##################################################################################
### FUNCTION - LINEAR REGRESSION - ORDINARY LEAST SQUARES
# regResults, regDependent, regPredictors = stats_regression(mapDF, regEquation=" age_years ~ bmi + race + race:sex ")
def stats_regression(dfIn, regEquation=" age_years ~ bmi + sex:race ", printOut=True):
    ### GET VARIABLES INTO X[n,p] predictors and y[n,1] outcome
    y, X = patsy.dmatrices(regEquation, dfIn, return_type='dataframe')
    ### GENERATE OLS REGRESSION MODEL ###
    statsModel = sm.OLS(y, X, disp=0)
    ### FIT DATA TO A LINE USING OLS MINIMIZATION ###
    statsModelFit = statsModel.fit(disp=0)
    ### PRINT RESULTS ###
    if printOut == True: print(statsModelFit.summary())
    ### RETURN: resultsObject, y, X
    return statsModelFit, y, X

"""
### DISPLAY IN NOTEBOOK ###
def display(self): display(self.df)
### DISPLAY PANDAS - SCIENTIFIC NOTATION ###
def display_scientific(self, dispPos=3): pd.set_option('display.float_format', '{:.'+str(dispPos)+'g}'.format)
### DISPLAY PANDAS - FLOAT NOTATION ###
def display_float(self, dispPos=3): pd.set_option('display.float_format', '{:.5f}'.format)
### DISPLAY PANDAS - MAX COLUMN WIDTH (-1 = No Wrap) ###
def display_colwidth(self, widthIn=-1): pd.set_option('display.max_colwidth', widthIn)
### Set Float Display ###
#pd.set_option('display.float_format', lambda x: '%.3f' % x)
"""
##################################################################################
### FUNCTION - Input Delimited File into Pandas Dataframe ###
def pd_in(path, sep='\t', index_col=None, skiprows=0, nrows=None, comment=None): 
    return pd.read_csv(path, sep=sep, index_col=index_col, skiprows=skiprows, nrows=nrows, low_memory=False, comment=comment, na_values=None)

##################################################################################
### FUNCTION - Read Excel File into Pandas Dataframe ###
def pd_in_excel(filePath, sheetIn=None):
    print(' - Loading Excel File: '+filePath+' - ')
    ### Load Available Sheets ###
    fileLoad = pd.ExcelFile(filePath)
    ### Check for Multiple Sheets ###
    if len(fileLoad.sheet_names) > 1: 
        if isinstance(sheetIn, str):
            sheetVal = sheetIn
        else:
            print('   - More than One Sheet, Please Choose One: ')
            sheetVal = ui_choose_list(fileLoad.sheet_names)
    elif len(fileLoad.sheet_names) == 1: fileLoad.sheet_names[0]
    else: print(' - Error: No Sheets Available in Excel File - '); return
    ### Load Specified Sheet ###
    print('   - Loading Sheet '+sheetVal+' - ')
    curTable = fileLoad.parse(sheetVal)
    print('   - Table Loaded - ')
    return curTable
#pd_in_excel(filePath, sheetIn=None)

##################################################################################
### FUNCTION - Check if Passed Item is Functional Pandas Dataframe ###
def pd_check(pdIn):
    # Check if Pandas Dataframe #
    if not isinstance(pdIn, pd.DataFrame): 
        print('   - Error: Not a Pandas Dataframe Object - '); return False
    # Check for Columns and Index #
    elif len(pdIn.columns) < 1:  
        print('   - Error: No Columns in the Pandas Dataframe - '); return False
    elif len(pdIn.index) < 1:  
        print('   - Error: No Indices in the Pandas Dataframe - '); return False
    # Otherwise Return True #
    else: return True

##################################################################################
### FUNCTION - Set Index of Pandas Dataframe ###
def pd_setindex(df, col): return df.set_index(col, drop=True, append=False, inplace=False)

##################################################################################
### FUNCTION - Write Pandas Dataframe to Delimited File ###
def pd_out(df, path, sep='\t', append=False, header=True): 
    if append == False: modeOut = 'w'
    else: modeOut = 'a'
    df.to_csv(path, sep=sep, na_rep='', mode=modeOut, header=header)

##################################################################################
### FUNCTION - Remove Dataframe Indices not in the Provided List ###
def pd_filter_indices(df, vals): return df[df.index.map(lambda x: x in vals)]

##################################################################################
### FUNCTION - MANN-U ON PANDAS DATAFRAME #####
### TAKES IN A PANDAS DATAFRAME AND CALCULATES KRUSKAL WALLIS AND PAIRWISE MANN WHITNEY U's ###
def pd_stats_mannu(dfIn, continuousColumn, categoricalColumn, performMWU=True):
    ### Setup Output Table ###
    dfOut = pd.DataFrame(columns=['Test','OneGroup','OneCount','OneMean','TwoGroup','TwoCount','TwoMean','pValue','TestStatistic','pFDR','pBonferroni','NumBonferroni']);tNum=0
    ### Perform Kruskal Wallis ###
    kruskalStat, kruskalP = sp.stats.mstats.kruskalwallis(*[np.array(dfIn[dfIn[categoricalColumn] == x][continuousColumn]) for x in dfIn[categoricalColumn].unique()])
    dfOut.loc[str(tNum)] = ['KruskalWallis', 'All',len(dfIn),np.mean(dfIn[continuousColumn]),np.NaN,np.NaN,np.NaN, kruskalP, kruskalStat,np.NaN,np.NaN,np.NaN]; tNum+=1
    if performMWU == True:
        ### For Each Combination Perform Mann-Whitney-U ###
        numberOfComparisons = ((len(dfIn[categoricalColumn].unique()))*((len(dfIn[categoricalColumn].unique()))-1))/2
        for groupCombination in itertools.combinations(dfIn[categoricalColumn].unique(),2): 
            try:
                mannStat, mannP = sp.stats.mannwhitneyu(dfIn.loc[dfIn[categoricalColumn]==groupCombination[0],[continuousColumn]], dfIn.loc[dfIn[categoricalColumn]==groupCombination[1],[continuousColumn]],use_continuity=True, alternative='two-sided')
                dfOut.loc[str(tNum)] = ['MannWhitneyU',groupCombination[0],len(dfIn[dfIn[categoricalColumn]==groupCombination[0]]),np.mean(dfIn[dfIn[categoricalColumn]==groupCombination[0]][continuousColumn]),groupCombination[1],len(dfIn[dfIn[categoricalColumn]==groupCombination[1]]),np.mean(dfIn[dfIn[categoricalColumn]==groupCombination[1]][continuousColumn]),mannP, mannStat,np.NaN,mannP*numberOfComparisons,numberOfComparisons]
                tNum+=1
            except:
                dfOut.loc[str(tNum)] = ['MannWhitneyU',groupCombination[0],len(dfIn[dfIn[categoricalColumn]==groupCombination[0]]),np.mean(dfIn[dfIn[categoricalColumn]==groupCombination[0]][continuousColumn]),groupCombination[1],len(dfIn[dfIn[categoricalColumn]==groupCombination[1]]),np.mean(dfIn[dfIn[categoricalColumn]==groupCombination[1]][continuousColumn]),np.NaN, np.NaN,np.NaN,np.NaN,np.NaN]
        ### Perform FDR Correction ###
        dfOut.loc[convert_nums_to_strs(np.arange(1,tNum)),'pFDR'] = multipletests(dfOut.loc[convert_nums_to_strs(np.arange(1,tNum)),'pValue'],alpha=0.05,method='fdr_bh')[1]
        ### Set Corrected P-values > 1.0 to 1.0 ###
        dfOut.loc[dfOut['pBonferroni']>1.0, ['pBonferroni']]=1.0
    return dfOut



##################################################################################
############################ BIOM FUNCTIONS ######################################
##################################################################################
##################################################################################
############################## BIOM I/O ##########################################

##################################################################################
### FUNCTION - Input a Single BIOM Table Given a Path
def biom_in(path): return biom.load_table(path)

##################################################################################
### Function: Import Directory of BIOM Tables
def biom_indir(dir_path, print_out=True):
    if print_out==True: print(' - Loading Directory of BIOM Tables - '); print('   '+dir_path)
    biomIn = []
    for curIdx, curBiom in enumerate(glob.glob(dir_path+'*.biom')):
        if (curIdx%10 == 0) and (print_out==True): print('   - Progress: '+str(curIdx))
        biomIn.append(biom_in(curBiom))
    return biomIn

##################################################################################
### FUNCTION - Write BIOM to TSV File
def biom_tsv(table, path): f = open(path,'w'); f.write(table.to_tsv(header_key='str', header_value='str2')); f.close()

##################################################################################
### FUNCTION - Write BIOM to JSON Format
def biom_json(table, path): f = open(path,'w'); f.write(table.to_json("biom", direct_io=None)); f.close()

##################################################################################
### FUNCTION - Return Biom List of Observations (i.e. OTUs)
def biom_obs(table): return table.ids(axis='observation')

##################################################################################
### FUNCTION - Return List of Observation Counts
def biom_obs_counts(bt): return bt.sum('observation')

##################################################################################
### FUNCTION - Return Biom List of Samples
def biom_sam(table): return table.ids(axis='sample')

##################################################################################
### FUNCTION - Return List of Sample Counts
def biom_sam_counts(bt): return bt.sum('sample')

##################################################################################
### Yield: Biom Samples one at a time by returning id, values, metadata
def biom_yieldsam(table):
    iterSamples=table.iter(axis='sample')
    for values, id, metadata in iterSamples: yield id, values, metadata

##################################################################################
### Yield: Biom Observations one at a time by returning id, values, metadata
def biom_yieldobs(table): 
    iterSamples=table.iter(axis='observation')
    for values, id, metadata in iterSamples: yield id, values, metadata

##################################################################################
### Yield: Rarefied BIOM Tables Filtered to just Observations and Samples in Passed BIOM Table ###
def biom_yield_filtered_rarefactions(table, rarefied_tables):
    ### Get Remaining Observations ###
    taxaOTUs = biom_obs(table)
    ### Get Remaining Observations ###
    taxaSamples = biom_sam(table)
    ### For each Rarefied Table...
    print(' - Trimming Rarefied tables -')
    for curIdx, curRare in enumerate(rarefied_tables):
        if curIdx%10==0: print('   - Progress: '+str(curIdx)+' - ')
        ### Trim Rarefied Table to Just OTUs and Samples from Table ###
        yield biom_filtersam_keep(biom_filterobs_keep(table=copy.deepcopy(curRare), listKeep=taxaOTUs), listKeep=taxaSamples)        

##################################################################################
### Yield: Ubiquity Fraction, Table Filtered of Observations by Minimumum Ubiquity
def biom_yield_minubiquity(varDict, ubiqFractions=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]):
    for ubiqCur in ubiqFractions:
        print( ' - Yielding: BIOM Filtered of OTUs in Less than '+str(ubiqCur*100)+' Percent of Samples -')
        yield ubiqCur, biom_filterobs_minubiq(copy.deepcopy(varDict['bt']), minubiq=ubiqCur)

##################################################################################
### Yield: Ubiquity Fraction, Table Filtered of Observations by Maximum Ubiquity
def biom_yield_maxubiquity(varDict, ubiqFractions=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]):
    for ubiqCur in ubiqFractions:
        print( ' - Yielding: BIOM Filtered of OTUs in More than '+str(ubiqCur*100)+' Percent of Samples -')
        yield ubiqCur, biom_filterobs_maxubiq(copy.deepcopy(varDict['bt']), maxubiq=ubiqCur)
        
##################################################################################
### FUNCTION - Add Sample Data from Map to Biom Table
def biom_metadata(table, df): return table.add_metadata(df.to_dict(orient='dict'))

##################################################################################
### FUNCTION - Summarize Key Stats of Biom Table
def biom_summarize(bt, description=None):
    outText = []; 
    if description == None: outText.append(' - BIOM Summary -')
    else: outText.append(' - '+str(description)+' -')
    outText.append('   - Samples         : '+str(len(biom_sam(bt)))+'   (i.e.)  '+str(biom_sam(bt)[0:3])+'...')
    outText.append('   - Observations    : '+str(len(biom_obs(bt)))+'   (i.e.)  '+str(biom_obs(bt)[0:3])+'...')
    outText.append('   - Total Counts    : '+str(float(bt.sum())))
    outText.append('   - Number Non-Zero : '+str(bt.nnz))
    #if bt.metadata() == None: outText.append('   - No Metadata Available')
    outText.append(''); list_print(outText)

##################################################################################
############################ BIOM FILTERING ######################################

##################################################################################
### FUNCTION - Filter Biom Observations with less than mincount
def biom_filterobs_mincount(table, mincount):
    filter_func = lambda values, id, md: sum(values) >= mincount
    return table.filter(filter_func, axis='observation', inplace=False)
    
##################################################################################
### FUNCTION - Filter Biom Observations with more than maxcount 
def biom_filterobs_maxcount(table, maxcount):
    filter_func = lambda values, id, md: sum(values) <= maxcount
    return table.filter(filter_func, axis='observation', inplace=False)
    
##################################################################################
### FUNCTION - Filter Biom Observations in less than minubiq fraction of samples (0.0 - 1.0)
def biom_filterobs_minubiq(table, minubiq):
    filter_func = lambda val, id_, md: sum(val>0)/len(val) >= minubiq
    return table.filter(filter_func, axis='observation', inplace=False)

##################################################################################
### FUNCTION - Filter Biom Observations in more than maxubiq fraction of samples (0.0 - 1.0)
def biom_filterobs_maxubiq(table, maxubiq):
    filter_func = lambda val, id_, md: sum(val>0)/len(val) <= maxubiq
    return table.filter(filter_func, axis='observation', inplace=False)

##################################################################################
### FUNCTION - Filter Biom to only Observations in listKeep
def biom_filterobs_keep(table, listKeep):
    filter_func = lambda values, id, md: id in listKeep
    return table.filter(filter_func, axis='observation', inplace=False)

##################################################################################
### FUNCTION - Filter Biom Observations in listRemove
def biom_filterobs_remove(table, listRemove):
    filter_func = lambda values, id, md: id not in listRemove
    return table.filter(filter_func, axis='observation', inplace=False)

##################################################################################
### FUNCTION - Filter Biom Samples with less than mincount
def biom_filtersam_mincount(table, mincount):
    filter_func = lambda values, id, md: sum(values) >= mincount
    return table.filter(filter_func, axis='sample', inplace=False)

##################################################################################
### FUNCTION - Filter Biom Samples with more than maxcount
def biom_filtersam_maxcount(table, maxcount):
    filter_func = lambda values, id, md: sum(values) <= maxcount
    return table.filter(filter_func, axis='sample', inplace=False)

##################################################################################
### FUNCTION - Filter Biom to only Samples in listKeep
def biom_filtersam_keep(table, listKeep): 
    filter_func = lambda values, id, md: id in listKeep 
    return table.filter(filter_func, axis='sample', inplace=False)

##################################################################################
### FUNCTION - Filter Biom Samples in listRemove
def biom_filtersam_remove(table, listRemove):
    filter_func = lambda values, id, md: id not in listRemove
    return table.filter(filter_func, axis='sample', inplace=False)

##################################################################################
### FUNCTION - Filter Biom samples not in metadata group
def biom_filtersam_mdkeep(table, metadata_category, metadata_group):
    filter_f = lambda values, id_, md: md[metadata_category] == metadata_group
    return table.filter(filter_f, axis='sample', inplace=False)

##################################################################################
### FUNCTION - Filter Biom samples in metadata group
def biom_filtersam_mdremove(table, metadata_category, metadata_group):
    filter_f = lambda values, id_, md: md[metadata_category] != metadata_group
    return table.filter(filter_f, axis='sample', inplace=False)

##################################################################################
### FUNCTION - Filter Biom Observations that are not in list of taxonomies
# Must specify the taxLevel if not to species level (i.e. 4:order, 5:family, 6:genus)
def biom_filtertaxa_keep(table, taxaKeep, taxLevel):
    allOTU = []
    for curOTU in biom_yieldobs(table):
        if (';'.join(curOTU[2]['taxonomy'][0:taxLevel])) in taxaKeep:
            allOTU.append(curOTU[0])
    return biom_filterobs_keep(table=table, listKeep=allOTU)

##################################################################################
### FUNCTION - Filter Biom Observations that are in list of taxonomies
# Must specify the taxLevel if not to species level (i.e. 4:order, 5:family, 6:genus)
def biom_filtertaxa_remove(table, taxaKeep, taxLevel):
    allOTU = []
    for curOTU in biom_yieldobs(table):
        if (';'.join(curOTU[2]['taxonomy'][0:taxLevel])) in taxaKeep:
            allOTU.append(curOTU[0])
    return biom_filterobs_remove(table=table, listKeep=allOTU)

##################################################################################
############################ BIOM CONVERSION #####################################

##################################################################################
### FUNCTION - Return Biom Counts as Pandas Dataframe
def biom_to_pd(table): return pd.DataFrame(table.matrix_data.todense().T.astype('float'), index=table.ids(axis='sample'), columns=table.ids(axis='observation')).T

##################################################################################
### FUNCTION - Return Biom Counts as Pandas Dataframe Collapsed by a Mapping Category (mean=False - returns sum)
def biom_to_pd_categorize(table, mapIn, mapCol, mean=True):
    # Convert Table to Dataframe #
    tablePD = biom_to_pd(table)
    # Add Mapping File and Collapse by Mean #
    if mean == True:
        btPD = pd.concat([tablePD.T, mapIn[mapCol]], axis=1).groupby([mapCol], axis=0).mean().T
    else:
        btPD = pd.concat([tablePD.T, mapIn[mapCol]], axis=1).groupby([mapCol], axis=0).sum().T
    # Return Collapsed Frame #
    return btPD

##################################################################################
### FUNCTION - Return Biom Counts as Numpy Array of Arrays (loop: [samples outer [observations inner]])
def biom_to_array(table): return table.matrix_data.todense().T

##################################################################################
### FUNCTION - Return Biom Table as Relative Abundance
def biom_relative(table): return table.norm(axis='sample', inplace=False)

##################################################################################
### FUNCTION - Return Biom Table as Presence/Absence (1 if present, 0 if not)
def biom_presence_absence(table): return table.pa()

##################################################################################
### FUNCTION - Repeatedly Rarefy Samples rareReps Times from Biom Table to rareDepth -> List of Rarefied Tables w/o taxonomy or metadata
def biom_rare(table, rareDepth, rareReps): 
    dfCounts = biom_to_pd(table); rareTables=[] # Get Biom Counts as Pandas Dataframe; empty list to store rarefied tables #
    rare_ind = lambda x: sk.stats.subsample_counts(x.astype(int), rareDepth, replace=False) # Subsample Function
    for i in np.arange(rareReps): rareTables.append( biom.Table( (dfCounts.apply(rare_ind, axis=0).as_matrix()) , dfCounts.index, dfCounts.columns) ) # For the number of rarefactions; rarefy df inside Biom Table Call and add table to rarefied List #
    return rareTables

##################################################################################
############################# BIOM TAXONOMY ######################################

##################################################################################
### FUNCTION - Collapse Observation Counts at Taxonomic Level
# tax_level: 0 = Kingdom | 1 = Phylum | 2 = Class | 3 = Order | 4 = Family | 5 = Genus | 6 = Species
def biom_taxonomy_collapse(table, taxLevel):
    collapse_f = lambda id_, md: '; '.join(md['taxonomy'][:taxLevel + 1])
    return table.collapse(collapse_f, axis='observation',norm=False)

### FUNCTION - Create Dataframe of Bacterial Taxonomy
# tax_level: 0 = Kingdom | 1 = Phylum | 2 = Class | 3 = Order | 4 = Family | 5 = Genus | 6 = Species
def biom_taxonomy_df(inputTable):
    print(' - Generating Dataframe of Taxonomy - ')
    xOut = pd.DataFrame(index=(biom_obs(inputTable)), columns=['kingdom','phylum','class','order','family','genus','species'])
    for curTax in biom_yieldobs(table=inputTable):
        for iX in np.arange(1,8): 
            xOut.loc[curTax[0]].iloc[(iX-1)]  = ';'.join(curTax[2]['taxonomy'][0:(iX)])
    return xOut

##################################################################################
### FUNCTION - Examine Observation Structuring by Categorical Mapping Column ###
# Takes in a BIOM Table, mapping dataframe and mapping column
# Returns a dataframe of each observations taxonomy and stats by group
def biom_group_stats(table, mapIn, mapCat,
                     nnz_group_count=[5,10]):
    print(' - Entering Group Statistics for Observations Pipeline - ')
    ### Get Dataframe of Relative Abundance Counts ###
    print('   - Getting Relative Abundance as Dataframe - ')
    btrDF = biom_to_pd( biom_relative(table)).T
    ### Add the Mapping Column ###
    print('   - Adding Mapping Column - ')
    btrDF[mapCat] = mapIn[mapCat]
    ### Get Unique Groupings ###
    mapUnique = mapIn[mapCat].unique()
    print('   - Getting Unique Categories: '+str(mapUnique))
    ### Get Dataframe of Bacterial Taxonomies ###
    print('   - Getting Dataframe of Bacterial Taxonomies - ')
    taxDf = biom_taxonomy_df(table)
    ### Add Columns for Group Stats ###
    for curGroup in mapUnique:
        taxDf[curGroup+'_nnz'] = 0
        taxDf[curGroup+'_ubiquity'] = 0.0
        taxDf[curGroup+'_sum'] = 0.0
        taxDf[curGroup+'_mean'] = 0.0
        taxDf[curGroup+'_median'] = 0.0
    ### For each Observation...
    print('   - Looping Through '+str(len(btrDF.columns))+' Observations - ')
    for curIdx, curOtu in enumerate(btrDF.columns):
        ### Check not map column ###
        if curOtu == mapCat: continue
        elif curIdx % 100==0: print('     - Progress: '+str(curIdx))
        ### For each unique group...
        for curGroup in mapUnique:
            ### Get Counts ###
            curCounts = btrDF[btrDF[mapCat] == curGroup][curOtu]
            ### Get Statistics ###
            taxDf.loc[curOtu, curGroup+'_nnz'] = len(curCounts[curCounts > 0])
            taxDf.loc[curOtu, curGroup+'_ubiquity'] = (len(curCounts[curCounts > 0]))/(len(curCounts))
            taxDf.loc[curOtu, curGroup+'_sum'] = curCounts.sum()
            taxDf.loc[curOtu, curGroup+'_mean'] = curCounts.mean()
            taxDf.loc[curOtu, curGroup+'_median'] = curCounts.median()
    ### Find the Number of Groups with At Least X Observation Counts ###        
    if (nnz_group_count != None) and (len(nnz_group_count)>0):
        for curNNZ in nnz_group_count:
            for curObs in taxDf.filter(regex=("._nnz")).index:
                taxDf.loc[curObs,'nnz_group_'+str(curNNZ)] = sum(taxDf.ix[curObs].filter(regex=("._nnz")) >= curNNZ)   
    return taxDf

#gsDf = biom_group_stats(table=varDict['bt'], mapIn=varDict['mp'], mapCat='race')

##################################################################################
### FUNCTION - Perform ANCOM Test of Microbial Differential Abundance ###
# Takes in a BIOM Table, mapping dataframe and mapping column
# Returns a dataframe of each observations taxonomy and stats by group
def biom_ancom(table, mapIn, mapCat, alpha=0.05, tau=0.02, theta=0.10):
    ### Get Dataframe of Relative Abundance Counts ###
    btrDf = biom_to_pd(biom_relative(table)).T
    ### Get Dataframe of Counts with Multiplicative Replacement ###
    btrDFmr = pd.DataFrame(data=sk.stats.composition.multiplicative_replacement(btrDf), 
                           columns=btrDf.columns, index=btrDf.index)
    ### Perform ANCOM Test ###
    ancomRes = ancom(table=btrDFmr, grouping=mapIn[mapCat],
                     alpha=alpha, tau=tau, theta=theta,
                     significance_test=sp.stats.mstats.kruskalwallis)
    return ancomRes

##################################################################################
### FUNCTION - Taxonomy Barchart ###
# Overview: Calculate Taxon Abundance Profiles
# Input: Biom Table (relative or absolute abundance table) & Taxonomic Level (taxLevel) 1=Phylum:6=Genus
### Plot Samples or Groups Abundance Profiles at Taxonomic Level ###
def biom_taxa_barchart(table, mp=None, collapseCat=None, outPath=None):
    print(' --- Entering Taxonomy Barchart Pipeline --- ')
    # Convert Table to Dataframe #
    bttdf = biom_to_pd(table)
    # Collapse by Mapping Category if Passed #
    if collapseCat != None: 
        if not (isinstance(mp, pd.DataFrame)): 
            print(' - Error: Must Pass Mapping File if Passing Collapse Category -')
            return
        else:
            print(' - Collapsing Mean Abundance by '+str(collapseCat)+' - ')
            bttdf = pd.concat([bttdf.T, mp[collapseCat]], axis=1).groupby([collapseCat], axis=0).mean().T
    # Create Figure #
    f = plt.figure(figsize=(len(bttdf.columns)*0.15,8)); ax=f.gca()
    bttdf.sort_index(inplace=True) # Sort so taxa in alphabetical order
    bttdf.T.plot.bar(stacked=True, cmap=plt.cm.get_cmap('rainbow'), ax=ax, width=1)
    ax.legend(loc='center left', bbox_to_anchor=(1.1, -0.1))
    # Check if Relative Abundance and Set ylim if So #
    if biom_sam_counts(table).max() <= 1.1: 
        print(' - Looks Like Relative Abundance: Setting Y-lim 0 to 1 -')
        plt.ylim(0,1)
    # Output Figure to File or Display #
    print(' - Outputting Figure - ')
    if outPath != None: 
        plt.savefig(outPath, bbox_inches='tight');  plt.close()
    else: 
        plt.show(); plt.close()
    print(' --- Exiting Taxonomy Barchart Pipeline --- \n')
    return

##################################################################################
### FUNCTION - Taxonomy Correlation ###
# Overview: Calculate Spearman Correlation between Bacterial Taxon Relative Abundance
# Input: Biom Table (table) & Taxonomic Level (taxLevel) 1=Phylum:6=Genus
# Output: Returns Dataframes of
    # Collapsed Taxon Relative Abundance
    # Spearman R Correlation by Taxa
    # Spearman p-values 
    # Mask matrix (True/False) of correlations passing bonferroni correction
# Author: Andrew Brooks - License: MIT - Date: 
def biom_taxa_correlation(table, outPath=None):
    print(' --- Entering Taxonomy Correlation Pipeline --- ')
    # Convert Table to Dataframe #
    bttdf = biom_to_pd(table)
    
    ### Calculate Spearman Correlation and Convert to Dataframes ###
    print(' - Calculating Correlation - ')
    spearmanR, spearmanP = sp.stats.spearmanr(bttdf, axis=1)
    spearmanR = pd.DataFrame(data=spearmanR, columns=bttdf.index, index=bttdf.index)
    spearmanP = pd.DataFrame(data=spearmanP, columns=bttdf.index, index=bttdf.index)
    
    ### Create Mask of Bonferroni Corrected p-values < 0.05 ###
    print(' - Creating Mask - ')
    mask = np.zeros_like(spearmanP)
    mask[np.where((spearmanP*((len(spearmanP)*len(spearmanP))-len(spearmanP))/2) > 0.05)] = True
    if outPath != None:
        
        ### Plot Spearman Correlation as Heatmaps and Clustermaps + masked versions ##
        plt.figure(figsize=(len(spearmanR)*1.2, len(spearmanR)))
        kwargs={'annot':True}
        # Plot ClusterMap #
        fig = sns.clustermap(data=spearmanR, method='average', figsize=(len(spearmanR)*1.1,len(spearmanR)*1.1), **kwargs)
        plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        plt.title('Clustermap of Spearman Correlation')
        plt.savefig(outPath+'.pdf'); plt.close()
        # Plot Masked Clustermap #
        fig = sns.clustermap(data=spearmanR, method='average', figsize=(len(spearmanR)*1.1,len(spearmanR)*1.1), mask=mask, **kwargs)
        plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        plt.title('Clustermap of Spearman Correlation Masked Bonferroni')
        plt.savefig(outPath+'_masked.pdf'); plt.close()
    print(' --- Exiting Taxonomy Correlation Pipeline --- \n')
    return spearmanR,spearmanP,pd.DataFrame(data=mask, columns=spearmanP.index, index=spearmanP.index)

##################################################################################
### FUNCTION - Taxonomy Network ###
# Overview: Input Taxonomy Dataframe and Correlation Stats and Generate Correlation Network File for Cytoscape
# Author: Andrew Brooks - License: MIT - Date: 
def biom_taxa_network(table, r, p, m, outPath=None):
    
    print(' --- Entering Taxonomy Network Pipeline --- ')
    
    # Convert Table to Dataframe #
    bttdf = biom_to_pd(table)
    
    # Check if Relative Abundance and Set ylim if So #
    if biom_sam_counts(table).max() <= 1.1: 
        print(' - Looks Like Relative Abundance: Setting Y-lim 0 to 1 -')
        plt.ylim(0,1)
        
    ### Create Network for Cytoscpe based on Correlation ###
    print('   - Masking with Bonferroni Significance - ')
    mR = r[(m == 0.0)]
    mP = p[(m == 0.0)]
    ### For Each Taxa... Write Node File ###
    file_write(path=outPath+'_nodes.txt',string='')
    print(' - Writing Node File - ')
    for tIdx, tId in enumerate(mR.index):
        ### Setup Empty Node String ###
        nodeStr = ''+tId+'\t'
        ### Get Terminal Taxa Name ###
        nodeStr += tId.split(';')[-1].strip()+'\t'
        ### Get Additional Taxonomic Levels ###
        for tRange in np.arange(2,10):
            try: nodeStr += tId.split(';')[-(tRange)]+'\t'
            except: continue
        nodeStr += str(bttdf.loc[tId,].mean())+'\n'
        file_append(path=outPath+'_nodes.txt', string=nodeStr)

    ### For Each Significant Correlation ###
    file_write(path=outPath+'_edges.txt',string='')
    print(' - Writing Edge File - ')
    for corrIdx, corrPair in enumerate(itertools.combinations(mR.index, 2)):
        edgeStr = ''
        ### Check if Correlation Bonferroni Significant ###
        if mP.loc[corrPair[0],corrPair[1]] <= 0.05:

            ### Get Terminal Taxa Name First Taxon ###
            edgeStr += corrPair[0]+'\t'+corrPair[0].split(';')[-1].strip()+'\t'
            ### Get Additional Taxonomic Levels ###
            for tRange in np.arange(2,10):
                try: edgeStr += corrPair[0].split(';')[-(tRange)]+'\t'
                except: continue

            ### Get Terminal Taxa Name Second Taxon ###
            edgeStr += corrPair[1]+'\t'+corrPair[1].split(';')[-1].strip()+'\t'
            ### Get Additional Taxonomic Levels ###
            for tRange in np.arange(2,10):
                try: edgeStr += corrPair[1].split(';')[-(tRange)]+'\t'
                except: continue
            ### Add Correlation and Significance ###
            edgeStr += str(mR.loc[corrPair[0],corrPair[1]])+'\t'+str(mP.loc[corrPair[0],corrPair[1]])+'\n'
            file_append(path=outPath+'_edges.txt', string=edgeStr)
    return

##################################################################################
######################### BIOM ALPHA DIVERSITY ###################################

##################################################################################
### FUNCTION - Compute Alpha Diversity for Samples Across Rarefied Biom Tables ### 'chao1','observed_otus','shannon','simpson'
def biom_alpha(rareTables, alphaMetric):
    for curIdx, curTable in enumerate(rareTables):
        if curIdx == 0: alphaAvg = sk.diversity.alpha_diversity(alphaMetric, biom_to_array(curTable).astype(int), ids=curTable.ids(axis='sample'), validate=True)
        else: alphaAvg += sk.diversity.alpha_diversity(alphaMetric, biom_to_array(curTable).astype(int), ids=curTable.ids(axis='sample'), validate=True)
    return alphaAvg / (len(rareTables))

##################################################################################
########################## BIOM BETA DIVERSITY ###################################

##################################################################################
### FUNCTION - Calculate Beta Diversity for Provided Metric
def biom_beta(table, betaMetric, treeIn=None): 
    if betaMetric == 'braycurtis': return biom_beta_bc(table)
    elif betaMetric == 'jaccard': return biom_beta_bj(table)
    elif betaMetric == 'unweighted_unifrac': return biom_beta_uu(table, treeIn)
    elif betaMetric == 'weighted_unifrac': return biom_beta_wu(table, treeIn)
    else: print('Error: Invalid Beta Diversity Metric '+betaMetric)

##################################################################################
### FUNCTION - Calculate Beta Diversity Distance Matrix - Bray Curtis
def biom_beta_bc(table): return sk.diversity.beta_diversity('braycurtis', (table.matrix_data.todense().astype(int).T), ids=table.ids(axis='sample'))

##################################################################################
### FUNCTION - Calculate Beta Diversity Distance Matrix - Binary Jaccard
def biom_beta_bj(table): return sk.diversity.beta_diversity('jaccard', (table.matrix_data.todense().astype(int).T), ids=table.ids(axis='sample'))

##################################################################################
### FUNCTION - Calculate Beta Diversity Distance Matrix - Unweighted Unifrac 
def biom_beta_uu(table, tree): return sk.diversity.beta_diversity('unweighted_unifrac', (table.matrix_data.todense().astype(int).T), ids=table.ids(axis='sample'), otu_ids=table.ids(axis='observation'), tree=tree)

##################################################################################
### FUNCTION - Calculate Beta Diversity Distance Matrix - Weighted Unifrac 
def biom_beta_wu(table, tree): return sk.diversity.beta_diversity('weighted_unifrac', (table.matrix_data.todense().astype(int).T), ids=table.ids(axis='sample'), otu_ids=table.ids(axis='observation'), tree=tree)

##################################################################################
### FUNCTION - Calculate Consensus Distance Matrix from List of Distance Matrices (converted through numpy matrices)
def biom_beta_consensus(betaTables):
    dmAvg = betaTables[0].data
    for curTab in np.arange(1,len(betaTables)): dmAvg = dmAvg + betaTables[curTab].data
    return DistanceMatrix((dmAvg / len(betaTables)), betaTables[0].ids)

##################################################################################
### FUNCTION - UPGMA Cluster Beta Diversity Distance Matrix 
def biom_beta_upgma(betaTable):
    betaLinkage = linkage(betaTable.condensed_form(), method='average')
    return TreeNode.from_linkage_matrix(betaLinkage, betaTable.ids)

##################################################################################
### FUNCTION - CALCULATE CONSENSUS TREE FROM LIST OF UPGMA TREES
def biom_beta_upgma_consensus(upgmaTrees): return sk.tree.majority_rule(upgmaTrees, weights=None, cutoff=0.5)[0]

##################################################################################
### FUNCTION - Calculate Average ANOSIM Distinguishability across Beta Distance Matrices for Specified Mapping Category ###
def biom_beta_anosim(betaTables, mapDF, category, permutations=99):
    print(' - Performing ANOSIM Analysis - ')
    ### If Passed a List of Distance Matrices ###
    if isinstance(betaTables, list):
        print(' - Assessing for Average Across List of Distance Matrices - ')
        anosimAvg = sk.stats.distance.anosim(betaTables[0], mapDF, column=category, permutations=permutations)
        for betaIdx in np.arange(1,len(betaTables)): 
            if betaIdx%10==0: print('   - Progress: '+str(betaIdx)+' - ')
            anosimAvg += sk.stats.distance.anosim(betaTables[betaIdx], mapDF, column=category, permutations=permutations)
        return (anosimAvg[2:])/len(betaTables)
    ### Elif Passed Single Distance Matrix ###
    elif isinstance(betaTables, DistanceMatrix):
        print(' - Assessing for Single Distance Matrix - ')
        anosimAvg = sk.stats.distance.anosim(betaTables, mapDF, column=category, permutations=permutations)
        return anosimAvg
    ### Else Return Input Error ###
    else: 
        print(' - Error: Must Pass List of Distance Matrices or Single Distance Matrix - ')
        return None

##################################################################################
### FUNCTION - CORRELATION BETWEEN DISTANCE MATRIX AND COVARIATE DATAFRAME FROM REGRESSION STYLE EQUATION
# TEST MODIFIED TO MEASURE CORRELATION OF ONLY SPECIFIED EQUATION, NOT ALL
# POSSIBLE SUBSETS OF VARIABLES
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with scikit-bio.
# ----------------------------------------------------------------------------
def biom_beta_bioenv(dmIn, dfIn, regEquation='age + sex:race + bmi + week', factorOfInterest = 'Week', permutations=10, recursive=False):
    if recursive == False: 
        print(' --- Performing BioENV Test --- ')
        print('   - DM ~ '+regEquation)
    # Store Factor of Interest #
    factorOfInterestList =  list(copy.deepcopy(dfIn[factorOfInterest]))
    ### EXTRACT REGRESSION EQUATION FROM COVARIATES DATAFRAME ###
    mapDfTemp = copy.deepcopy(dfIn) # Make a Copy of the Mapping File
    mapDfTemp['temp'] = np.arange(len(mapDfTemp)) # Add Filler Numeric Column
    updateRegEq = ' temp ~ '+regEquation
    y, mapDfIn = patsy.dmatrices(updateRegEq, data=copy.deepcopy(mapDfTemp), eval_env=False,  return_type='dataframe')
    del mapDfIn['Intercept']
    ### MAKE SURE DISTACE MATRIX IS OBJECT ###
    if not isinstance(dmIn, DistanceMatrix): raise TypeError("Must provide a DistanceMatrix as input.")
    ### MAKE SURE MAP IS PANDAS DF ###
    if not isinstance(mapDfIn, pd.DataFrame): raise TypeError("Must provide a pandas.DataFrame as input.")  
    ### IF NO LIST OF COLUMNS TO USE THEN USE ALL ###
    columns = mapDfIn.columns.values.tolist()
    ### CHECK FOR DUPLICATE COLUMNS ###
    if len(set(columns)) != len(columns): raise ValueError("Duplicate column names are not supported.")
    ### CHECK TO MAKE SURE AT LEAST ON COLUMN FACTOR ###
    if len(columns) < 1: raise ValueError("Must provide at least one column.")
    ### MAKE SURE COLUMNS ARE IN DATAFRAME ###
    for column in columns:
        if column not in mapDfIn: raise ValueError("Column '%s' not in data frame." % column)
    ### ORDER COVARIATES TO SAMPLES IN DISTANCE MATRIX ###
    variablesDF = mapDfIn.loc[dmIn.ids, columns]
    ### CHECK FOR MISSING IDS IN DF FROM DM ###
    if variablesDF.isnull().any().any(): raise ValueError("One or more IDs in the distance matrix are not in the data frame, or there is missing data in the data frame.")
    ### CHECK ALL COVARIATES ARE CASTABLE AS FLOATS ##
    try: variablesDF = variablesDF.astype(float)
    except ValueError: raise TypeError("All specified columns in the data frame must be numeric.")
    ##### CENTER DATAFRAME AROUND MEAN OF EACH COVARIATE AND COMPUTE AS STD #####
    if recursive == False: print('   - Transforming the Data - ')
    # Modified from http://stackoverflow.com/a/18005745
    dfCenterScale = variablesDF.copy()
    dfCenterScale -= dfCenterScale.mean()
    dfCenterScale /= dfCenterScale.std()
    ### MAKE SURE COULD BE CENTERED AROUND STD ###
    if dfCenterScale.isnull().any().any(): raise ValueError("Column(s) in the data frame could not be scaled, likely because the column(s) had no variance.")
    ### CAST AS ARRAY ###
    dfCovariatesArray = dfCenterScale.values
    ##### FLATTEN THE DISTANCE MATRIX #####
    dmInFlattened = dmIn.condensed_form()
    ### COMPUTE PREDICTOR DISTANCE MATRIX ###
    dfEuclidianDistance = sp.spatial.distance.pdist(dfCovariatesArray,metric='euclidean')
    ### CALCULATE CORRELATION BETWEEN DISTANCE MATRIX AND COVARIATE DISTANCES ###
    if recursive == False: print('   - Calculating Correlation - ')
    rhoOut = sp.stats.spearmanr(dmInFlattened, dfEuclidianDistance)[0]
    if recursive == False:  print('   - Correlation: '+str(rhoOut))
    ### RANDOMLY PERMUTE PREDICTOR A NUMBER OF TIMES ###
    if (factorOfInterest != None) and (recursive == False):
        print('   - Randomly Permuting Factor of interest: '+factorOfInterest+' - ')
        dfOutF = pd.DataFrame(columns=['equation','permute','rho','p'])
        repPs = []
        ### STORE CORRELATION WHEN PREDICTOR IS RANDOMLY PERMUTED ###
        for repSig in np.arange(permutations):
            if repSig%10==0:print('     - Progress: '+str(repSig)+' - ')
            ### PERMUTE FACTOR OF INTEREST ###
            np.random.shuffle(factorOfInterestList)
            mapDfTemp[factorOfInterest] = factorOfInterestList
            ### PERFORM REGRESSION AND STORE P-VALUES ###
            repPs.append(biom_beta_bioenv(dmIn, dfIn=mapDfTemp, regEquation=regEquation, factorOfInterest=factorOfInterest, permutations=permutations, recursive=True)[0])
        sigOut = (int(len(repPs))-int((sum(rhoOut > repPs)))) /  int(len(repPs))
        print('   - Significance: '+ str(sigOut)+' ('+str(len(repPs))+'-'+str((sum(rhoOut > repPs)))+')/'+str(len(repPs)))
    else: sigOut = None
    if recursive == False:  print(' --- Finished BioENV Test --- ')
    return rhoOut, sigOut

##################################################################################
########################## PHYLOGENY FUNCTIONS ###################################
##################################################################################

##################################################################################
### FUNCTION - Input Phylogeny into SciKit-Bio TreeNode Object and Validate ###
def tree_in(trePath, ete=False):
    # Try to open as TreeNode(skbio) #
    if ete == True:
        return tree_ete_in(trePath)
    try:
        treeIn = TreeNode.read(trePath)
        print('   - Read TreeNode(skbio) from: '+trePath)
        # Make Name File Source #
        if isinstance(trePath, str): treeIn.name = trePath
        # Validate the Tree and Return #
        return tree_validate(treeIn)
    except:
        print('   - Error: Could not Read as TreeNode(skbio) - ')
        return None
        
##################################################################################
### FUNCTION - Validate Tree Object - Includes Steps Below ###
def tree_validate(treeIn):
    # 1. If any branch lengths are None then set to 0.0 - Necessary for Error Correction #
    for idx, e in enumerate(treeIn.traverse()):
        if e.length == None: e.length = 0.0
    # 2. Root tree at midpoint - Necessary for Error Correction #
    return treeIn.root_at_midpoint()

##################################################################################
### FUNCTION - Return List of Tree Tips ###
def tree_tips(treeIn): return [x.name for x in treeIn.tips()]

##################################################################################
### FUNCTION - Trim TreeNode to Only Observations in Biom Table (Table can be missing observations in Tree - filter by tree_tips(treeIn) if necessary!) ###
def tree_match_tips(table, treeIn):
    # Adapted from the GNEISS package. Morton JT... Knight R. 2017. Balance trees reveal microbial niche differentiation. mSystems 2:e00162-16. https://doi.org/10.1128/mSystems.00162-16. DOWNLOADED FROM: https://github.com/biocore/gneiss/
    sharedTips = list(set(tree_tips(treeIn)) & set(table.ids(axis='observation'))) # Get list of shared observations
    treeOut = treeIn.shear(names=sharedTips) # Trim to only list of shared
    treeOut.bifurcate(); treeOut.prune(); treeOut = tree_validate(treeOut); return treeOut # Clean the trimmed tree and return

##################################################################################
### FUNCTION - Write Tree Object to Newick File ###
def tree_out(treeIn, path): treeIn.write(path, format='newick')

##################################################################################
############################# VCF FUNCTIONS ######################################
##################################################################################

##################################################################################
### FUNCTION - Open VCF File and Yield Variants One-by-One ###
# verifySNP=True: Check to make sure the variant is a SNP and not an Indel
def vcf_yield(vcfPath, verifySNP=True):
    vcfReader = vcf.Reader(open(vcfPath)) # Open VCF File
    for SNP in vcfReader: # For Each SNP...
        if (verifySNP==True) and not (SNP.is_snp): continue # Verify it is a SNP... or else skip
        yield SNP

##################################################################################
### FUNCTION - Print Basic SNP Info ###
def vcf_info(SNP): print(SNP); MAF = (len([x.gt_type for x in SNP.samples if x.gt_type > 0]) / (len(SNP.samples)*2)); print('   MAF: '+str(MAF))

##################################################################################
### FUNCTION - Loop through Samples in SNP and Add SNP Genotype to Sample Indexed Pandas Dataframe ###
def vcf_populations(SNP ,sampleInfo):
    sampleInfo[SNP.ID] = 0 # Add SNP Column as RSID
    for curSample in SNP.samples: # For each sample call in the SNP...
        sampleInfo.loc[curSample.sample, SNP.ID] = curSample.gt_type # Set Sample to genotype (0,1,2)
    return sampleInfo # Return Updated Sample Info        

##################################################################################
##################################################################################
##################################################################################
############################## PIPELINES #########################################
##################################################################################
##################################################################################
##################################################################################

##################################################################################
############################ PIPELINE INPUT ######################################
##################################################################################
# rare_dir and dist_matrices must be REGEX strings, i.e. beta/dm_*.txt
def pipe_input(dir_path=None, tsv_path=None, csv_path=None, biom_path=None, map_path=None, rare_dir=None, dist_matrix=None, dist_matrices=None, tre_path=None, vocalPipe=True):
    
    ### Start Pipeline and Create Output Object ###
    varDict = {}
    
    ### Set Output Directory ###
    if dir_path != None: 
        print(" - Setting Output Directory -o: \n   "+str(dir_path)+'\n')
        varDict['dir'] = dir_path
        dir_make(dirPath=dir_path)
    
    ### Input TSV File ###
    if tsv_path != None: 
        print(" - Loading TSV File -x: \n   "+str(tsv_path))
        varDict['tsv'] = pd_in(path=tsv_path,sep='\t')
        print("   - Loaded TSV File -x\n")
    
    ### Input CSV File ###
    if csv_path != None: 
        print(" - Loading CSV File -c: \n   "+str(csv_path))
        varDict['csv'] = pd_in(path=csv_path,sep=',')
        print("   - Loaded CSV File -c\n")
        
    ### Input Biom Table ###
    if biom_path != None: 
        print(" - Loading BIOM Table -b: \n   "+str(biom_path))
        varDict['bt'] = biom_in(biom_path)
        print("   - Loaded BIOM Table -b")
        biom_summarize(varDict['bt'])
        try: 
            varDict['rt'] = biom_relative(table=varDict['bt'])
            print('   - Relative Abundance Table Computed - \n')
        except:
            print('   - Warning: Could not Calculate Relative Abundance - \n')
    ### Input Mapping File ###
    if map_path != None: 
        print(" - Loading Mapping File -m: \n   "+str(map_path))
        mp = pd_in(map_path)
        varDict['mp'] = pd_setindex(mp, '#SampleID')
        print("   - Loaded Mapping File -m\n")
        ### If the Map Index is Numeric Recast as Int (for HMP sample indices) ###
        if isinstance(varDict['mp'].index, pd.indexes.numeric.Int64Index):
            print('   - Recasting Map Index as Numeric -')
            varDict['mp'].index = convert_nums_to_strs(varDict['mp'].index)
    
    ### Input Directory of Rarefied Tables ###
    if rare_dir != None:
        varDict['rts'] = []
        print(" - Loading Directory of Rarefied Tables -r: \n   "+str(rare_dir))
        for curIdx, curTable in enumerate(glob.glob(rare_dir)): 
            varDict['rts'].append(biom_in(curTable))
            if curIdx % 10 == 0: 
                print("     - Completed: "+str(curIdx)+" -")
        print("   - Loaded Directory of Rarefied Tables -r\n")
    
    ### Input Distance Matrix ###
    if dist_matrix != None: 
        print(" - Loading Distance Matrix -d: \n   "+str(dist_matrix))
        varDict['dm'] = DistanceMatrix.read(dist_matrix)
        print("   - Loaded Distance Matrix -d\n")
    
    ### Input Directory of Distance Matrices ###
    if dist_matrices != None:
        varDict['dms'] = []
        print(" - Loading Distance Matrices from Folder -f: \n   "+str(dist_matrices))
        for curIdx, curDM in enumerate(glob.glob(dist_matrices)):
            varDict['dms'].append(DistanceMatrix.read(curDM))
            if curIdx % 10 == 0: 
                print("     - Completed: "+str(curIdx)+" -")
        print("   - Loaded Distance Matrices from Folder -f\n")
    
    ### Input Phylogeny or Dendrogram ###
    if tre_path != None: 
        print(" - Loading Tree Dendrogram -t: \n   "+str(tre_path))
        varDict['tre'] = tree_in(tre_path)
        print("   - Loaded Tree Dendrogram -t\n")
    
    ### Return Input Objects ###
    return varDict

##################################################################################
######################## PIPELINE BIOM PROCESS ###################################
##################################################################################
### Function: Flexibly Import BIOM Objects (bt,mp,tre), Perform QC, Generate Realtive Abundance, Write Results, Return varDict
def pipe_biom_process(bt,  # Must provide bt (biom.Table Object or Path to BIOM table as TSV or JSON) 
                      mp=None,  # Must be a pd.Dataframe Object or Path to TSV File in QIIME Mapping Format
                      tre=None, # Must be a TreeNode Object or Path to Newick File
                      dir_path=None, # Must be a String Path to Output Directory (will make if doesn't exist)
                      skip_qc=False,
                      filter_by_map=True, # Filter to Only Samples in Mapping File
                      filter_by_tre=False, # Filter to Only Observations in Tree
                      filter_sam_keep=None,    filter_sam_remove=None, # Filter List of Samples
                      filter_sam_min=None,     filter_sam_max=None, # Filter Samples by Abundance
                      filter_obs_keep=None,    filter_obs_remove=None, # Filter List of Observations
                      filter_obs_min=None,     filter_obs_max=None, # Filter Observations by Abundance
                      filter_obs_minubiq=None, filter_obs_maxubix=None, # Filter Observations by Ubiquity
                      gen_relative=True, # If true will generate relative abundance table
                      write_biom=True, write_tsv=True): # Write BIOM Table into output folder
    ##################################################################################
    print(' --- Entering BIOM Process Pipeline --- \n'); 
    varDict = {} # Initialize the Empty Dictionary Object #

    ### Import BIOM Table ### bt ###
    # If passed a BIOM Table #
    if isinstance(bt, biom.Table): 
        varDict['bt'] = copy.deepcopy(bt)
    # If passed a string path #
    elif isinstance(bt, str): 
        varDict['bt'] = pipe_input(biom_path=bt)['bt']
    # If passed a varDict object #
    elif isinstance(bt, dict):
        if 'bt' in bt.keys():
            varDict['bt'] = copy.deepcopy(bt['bt'])
        else: 
            print(" - Error: Passed Dictionary Doesn't Contain Biom Table varDict['bt'] - ")
            return None
        if ('mp' in bt.keys()) and (mp == None): 
            mp = copy.deepcopy(bt['mp'])
        if ('tre' in bt.keys()) and (tre == None): 
            tre = copy.deepcopy(bt['tre'])
        if ('dir' in bt.keys()) and (dir_path == None): 
            dir_path = copy.deepcopy(bt['dir'])
    else: 
        print(' - Error: bt Object Must be biom.Table or Path to Table or formatted varDict - ')
        return None
    
    ### Import Map ### mp ###
    if isinstance(mp, pd.DataFrame): 
        varDict['mp'] = mp
    elif isinstance(mp, str): 
        varDict['mp'] = pipe_input(map_path=mp)['mp']
    elif 'mp' in varDict.keys(): 
        pass
    else: 
        print(' - Warning: No Mapping File Loaded -\n')
    
    ### Import Tree ### tre ###
    if isinstance(tre, TreeNode): 
        varDict['tre'] = tre
    elif isinstance(tre, str): 
        varDict['tre'] = pipe_input(tre_path=tre)['tre']
    elif 'tre' in varDict.keys(): 
        pass
    else: 
        print(' - Warning: No Tree File Loaded -\n')
    
    ### Set Output Directory ### dir ###
    if isinstance(dir_path, str): 
        varDict['dir'] = dir_path
    elif 'dir' in varDict.keys(): 
        dir_path = varDict['dir']
    else: 
        print(' - Warning: No Output Directory Set\n')
    
    if skip_qc == True: print(' --- Exiting BIOM Process Without QC ---\n'); return varDict
    ##################################################################################
    ### Summarize Input Table ###
    biom_summarize(varDict['bt'], description='BIOM Input Table')
    
    ### Filter Biom Table to only Samples in Map Index ###
    if 'mp' in varDict.keys():
        if (filter_by_map == True):
            print(' - Filtering Table to Only Samples in Map -')
            varDict['bt'] = biom_filtersam_keep(varDict['bt'], varDict['mp'].index)

    ### Trim the Biom Table to Only Observations in the Tree ###
    if 'tre' in varDict.keys():
        if (isinstance(varDict['tre'], TreeNode)) & (filter_by_tre == True):
            print(' - Filtering Table to Only Observations in Tree')
            varDict['bt'] = biom_filterobs_keep(varDict['bt'], tree_tips(varDict['tre']))
        
    ### Filter Samples by List ###
    if filter_sam_keep != None: 
        print(' - Filtering to Only Samples in Passed List -')
        varDict['bt'] = biom_filtersam_keep(varDict['bt'], filter_sam_keep)
    if filter_sam_remove != None: 
        print(' - Filtering Samples Not in Passed List -')
        varDict['bt'] = biom_filtersam_remove(varDict['bt'], filter_sam_remove)
    
    ### Filter Observations by List ###
    if filter_obs_keep != None: 
        print(' - Filtering to Only Observations in Passed List -')
        varDict['bt'] = biom_filterobs_keep(varDict['bt'], filter_obs_keep)
    if filter_obs_remove != None: 
        print(' - Filtering Observations Not in Passed List -')
        varDict['bt'] = biom_filterobs_remove(varDict['bt'], filter_obs_remove)

    ### Filter Observations by min and max count ###
    if filter_obs_min != None: 
        print(' - Filtering Observations with Less than '+str(filter_obs_min)+' Counts -')
        varDict['bt'] = biom_filterobs_mincount(varDict['bt'], filter_obs_min)
    if filter_obs_max != None: 
        print(' - Filtering Observations with More than '+str(filter_obs_max)+' Counts -')
        varDict['bt'] = biom_filterobs_maxcount(varDict['bt'], filter_obs_max)

    ### Filter Observations by Ubiquity ###
    if filter_obs_minubiq != None:
        print(' - Filtering Observations with Less than '+str(filter_obs_minubiq)+' Ubiquity -')
        varDict['bt'] = biom_filterobs_minubiq(varDict['bt'], filter_obs_minubiq)
    if filter_obs_maxubix != None: 
        print(' - Filtering Observations with More than '+str(filter_obs_maxubix)+' Ubiquity -')
        varDict['bt'] = biom_filterobs_maxubiq(varDict['bt'], filter_obs_maxubix)     

    ### Filter Samples by min and max count ###
    if filter_sam_min != None:
        print(' - Filtering Samples with Less than '+str(filter_sam_min)+' Counts -')
        varDict['bt'] = biom_filtersam_mincount(varDict['bt'], filter_sam_min)
    if filter_sam_max != None: 
        print(' - Filtering Samples with More than '+str(filter_sam_max)+' Counts -')
        varDict['bt'] = biom_filtersam_maxcount(varDict['bt'], filter_sam_max)
        
    ### Trim the Tree to Only Observations Shared with Biom Table ###
    if 'tre' in varDict.keys():
        print(' - Trimming Tree to Only Observations Shared with Biom Table -')
        varDict['tre'] = tree_match_tips(varDict['bt'], varDict['tre'])
    
    ### Filter Map to Only Samples in Biom Table ###
    if 'mp' in varDict.keys():
        print(' - Filtering Map to Only Samples in Table -')
        varDict['mp'] = pd_filter_indices(varDict['mp'], biom_sam(varDict['bt']))
        
        ### Add Sample Data from Map to Biom Table ###
        print(' - Adding Mapping Data to Table Object -')
        varDict['bt'].add_metadata(varDict['mp'].to_dict(orient='dict'))
        
    ### Summarize Table ###
    print(); biom_summarize(varDict['bt'], description='BIOM Post-QC Filtering')
    
    ##################################################################################
    ##### 1_0_Output QC Results #####
    if 'dir' in varDict.keys(): # Check if output directory specified
        if ((write_tsv == True) or (write_biom  == True)): # If set to write to TSV or BIOM
            ### Make Output Directory ###
            dir_make(dir_path); print('\n - Making Output Directory and Saving QC Files - ')
            ### Write Quality Controlled Biom Table to TSV (Doesn't contain taxonomy) ###
            if write_tsv  == True: 
                print('   - Writing BIOM as TSV: '+varDict['dir']+'/1_0_biom.txt - ')
                biom_tsv(varDict['bt'], varDict['dir']+'/1_0_biom.txt')
            ### Write Quality Controlled Biom Table to JSON (Contains Taxonomy) ###
            if write_biom == True: 
                print('   - Writing BIOM as JSON: '+varDict['dir']+'/1_0_biom.biom - ')
                biom_json(varDict['bt'], varDict['dir']+'/1_0_biom.biom')
            ### Write Quality Controlled Mapping File ###
            if isinstance(varDict['mp'],pd.DataFrame): 
                print('   - Writing Mapping File as TSV: '+varDict['dir']+'/1_0_map.txt - ')
                pd_out(varDict['mp'], varDict['dir']+'/1_0_map.txt')
            ### Write Quality Controlled Tree to Newick File ###
            if isinstance(varDict['tre'],TreeNode): 
                print('   - Writing Tree Phylogeny as Newick: '+varDict['dir']+'/1_0_tree.tre -\n')
                tree_out(varDict['tre'], varDict['dir']+'/1_0_tree.tre')
    
    ##################################################################################
    ### 1_1_Generate: BIOM Table with Relative Abundance ###
    if gen_relative==True: # If generate relative abundance tables is specified
        print(" - Generating Relative Abundance BIOM Table - varDict['btr']")
        varDict['btr'] = biom_relative(varDict['bt']) # Make relaitve abundance table -> varDict['btr']
        if 'dir' in varDict.keys(): # Check if output directory specified
            if ((write_tsv == True) or (write_biom  == True)):
                ### Write Relative Biom Table to TSV (Doesn't contain taxonomy) ###
                outPath=varDict['dir']+'/1_1_biom_relative'
                if write_tsv  == True: 
                    print('   - Writing Relative Abundance Table as TSV -'); print('     - '+outPath+'.txt')
                    biom_tsv( varDict['btr'], outPath+'.txt')
                if write_biom == True: 
                    print('   - Writing Relative Abundance Table as JSON  -'); print('     - '+outPath+'.biom')
                    biom_json(varDict['btr'], outPath+'.biom')  
                    
    # Print Dictionary Keys #
    print('\n - Dictionary Keys -')
    print(varDict.keys())
    print('\n --- Exiting BIOM Process Pipeline --- \n')
    ### Return Dictionary ###
    return varDict


##################################################################################
###################### PIPELINE BIOM RAREFACTION #################################
##################################################################################
def pipe_biom_rarefaction(btIn, dirPath=None, minSam=1000, rareReps=100, tsvTrue=True, biomTrue=True):
    ### Get List of Rarefied Biom Tables (Note: Don't have taxonomy or metadata to save space) ###
    rareTables = biom_rare(table=btIn, rareDepth=minSam, rareReps=rareReps)
    ### Output rarefied tables in tsv format ###
    if (tsvTrue == True) and (dirPath != None):
        dir_make(dirPath+'/');dir_make(dirPath+'/rarefaction_tsv_'+str(minSam)+'/')
        for curIdx, curTable in enumerate(rareTables): 
            biom_tsv(curTable, dirPath+'/rarefaction_tsv_'+str(minSam)+'/rarefaction_'+str(curIdx)+'.txt')
    ### Output rarefied tables in biom format ###
    if biomTrue == True and (dirPath != None):
        dir_make(dirPath+'/')
        dir_make(dirPath+'/rarefaction_biom_'+str(minSam)+'/')
        for curIdx, curTable in enumerate(rareTables): 
            biom_json(curTable, dirPath+'/rarefaction_biom_'+str(minSam)+'/rarefaction_'+str(curIdx)+'.biom')
    return rareTables # Return Rarefied Tables
    
#pipe_rarefaction(outDict['bt'], dirPath='rarefaction/',minSam=200,rareReps=10)

##################################################################################
######################## PIPELINE BIOM DIVERSITY #################################
##################################################################################
### Function: Compute Alpha and Beta Diversity 
def pipe_biom_diversity(varDict, rarefaction='', dir_path=None, 
                        alpha=['shannon','pielou_e','observed_otus','simpson','chao1'], 
                        beta=['braycurtis', 'jaccard', 'unweighted_unifrac', 'weighted_unifrac'], consensus=True, upgma=False):
    print(' --- Entering BIOM Diversity Pipeline --- ')
    
    ### Check for Rarefied and Relative Abundance Tables ###
    if 'rts' not in varDict.keys(): 
        print(" - Error: Rarefied Tables not Computed in varDict['rts'] -"); return
    if 'mp' not in varDict.keys():
        print(" - Error: Mapping File not in varDict['mp] -")
        
    ##### 1_3_Alpha_Diversity #####
    if len(alpha) > 0:
        ### Calculate Alpha Diversity for Each Sample ### 
        print('\n - Computing Alpha Diversity - ')
        alphaDiversity=pd.DataFrame(index=varDict['mp'].index) # Make a copy of the mapping file
        for alphaMetric in alpha: # For Each Alpha Diversity Metric
            print('   - '+alphaMetric+' - ')
            alphaDiversity[alphaMetric] = biom_alpha(varDict['rts'], alphaMetric) # Calculate Alpha Diversity for Metric
        ### Add Alpha Diversity to Mapping File ###
        varDict['mp'] = pd.concat([varDict['mp'], alphaDiversity], axis=1)
        ### Write Alpha Diversity to File ###
        if dir_path != None: 
            alphaDiversity.to_csv(dir_path+'/1_3_alpha_'+str(rarefaction)+'.txt', sep='\t',header=True, index=True)
            ### Write Quality Controlled and Alpha Diversity Mapping File ###
            pd_out(varDict['mp'], dir_path+'/1_3_alpha_map_'+str(rarefaction)+'.txt')
    
    ##### 1_4_Beta_Diversity #####
    if len(beta) > 0:
        print('\n - Computing Beta Diversity - ')
        ### Initialize Storage ###
        varDict['dms']={}; varDict['dm']={}; 
        if upgma == True: varDict['upgmas']={}; varDict['upgma']={}
        ### For Each Beta Diversity Metric ###
        for betaMetric in beta:
            ### Calculate and Store Beta Diversity Distance Matrices in List ###
            print('   - Computing '+betaMetric+' Beta Diversity - ')
            varDict['dms'][betaMetric] = []
            for curIdx, curTable in enumerate(varDict['rts']):
                if curIdx %10==0: print('     - Progress: '+str(curIdx))
                varDict['dms'][betaMetric].append(biom_beta(curTable, betaMetric, treeIn=varDict['tre']))
            ### Calculate Consensus Distance Matrix ###
            if consensus==True: varDict['dm'][betaMetric] = biom_beta_consensus(varDict['dms'][betaMetric])
            
            ### If Output Directory Specified ###
            if dir_path != None:
                print(' - Writing '+betaMetric+' Beta Diversity -> /1_4_beta/'+betaMetric+'_'+str(rarefaction)+'/')
                ### Make Folder to Store Distance Matrices ###
                dir_make(dir_path+'/1_4_beta/')
                dir_make(dir_path+'/1_4_beta/'+betaMetric+'_'+str(rarefaction)+'/')
                ### Save Each Beta Distance Matrix ###
                for curIdx, curTable in enumerate(varDict['dms'][betaMetric]): 
                    if curIdx %10==0: 
                        print('     - Progress '+str(curIdx)+' - ')
                    curTable.write(dir_path+'/1_4_beta/'+betaMetric+'_'+str(rarefaction)+'/dm_'+str(curIdx)+'.txt',format='lsmat')
                ### Save Consensus Distance Matrix ###
                if consensus==True: catchOut = varDict['dm'][betaMetric].write((dir_path+'/1_4_beta/'+betaMetric+'_'+str(rarefaction)+'/consensus.txt'),format='lsmat')
            
            ### If upgma is True then Calculate Relationship Dendrograms ###
            if upgma == True:
                print(' - Computing UPGMA Dendrograms for '+betaMetric+' Beta Diversity -')
                varDict['upgmas'][betaMetric] = []
                ### Calculate UPGMA Trees from Distance Matrices ###
                for curTable in varDict['dms'][betaMetric]: 
                    varDict['upgmas'][betaMetric].append(biom_beta_upgma(curTable)) 
                ### Calculate Consensus Tree from UPGMA List of Individuals Trees ###
                varDict['upgma'][betaMetric] = biom_beta_upgma_consensus(varDict['upgmas'][betaMetric])
                ### Save Each UPGMA Tree ###
                if dir_path != None:
                    for curIdx, curTree in enumerate(varDict['upgmas'][betaMetric]): tree_out(curTree, dir_path+'/1_4_beta/'+betaMetric+'_'+str(rarefaction)+'/upgma_'+str(curIdx)+'.tre')
                    ### Save Consensus UPGMA Tree ###
                    tree_out(varDict['upgma'][betaMetric], dir_path+'/1_4_beta/'+betaMetric+'_'+str(rarefaction)+'/consensus_upgma.tre')

    print(' --- Exiting BIOM Diversity Pipeline --- ')
    return varDict
#varDict = pipe_biom_diversity(varDict=varDict, rarefaction='', dir_path=varDict['dir'], alpha=['shannon','pielou_e','observed_otus','simpson','chao1'],beta=['braycurtis', 'jaccard', 'unweighted_unifrac', 'weighted_unifrac'])

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
####### PANDAS & STATISTICAL TOOLS (DataFrame - pandas ) #########################
# Class to contain a pandas dataframe and additional info
# Tools include: Statistics, Plotting
class DF:

    ### CLASS INITIALIZATION ##################################################
    def __init__(self, objIn, sepIn='\t', verboseIn=True, columnsIn=None, indexIn=None):
        
        ##### DATA STRUCTURES #####
        ### VERBOSE BOOL .verbose ###
        self.verbose = verboseIn
        ### SEPARATOR BOOL .sep ### Default '\t'
        self.sep=sepIn
        ##### INPUT PATH .path ###
        self.source=objIn
        ### SET PALETTE ###
        self.palette = None
        ### SET COLUMNS AND INDEX IF PROVIDED ###
        if columnsIn != None: self.columns = columnsIn
        if indexIn   != None: self.index = indexIn
        
        ##### LOAD FROM PROVIDED SOURCE #####
        self = self.load()
        ### UPDATE COLUMNS ###
        self.columns = list(self.df.columns)
        ### UPDATE INDEX ###
        self.index = list(self.df.index)         
    
    ### SET ACTION IF CALLED BY PRINT() ######################################
    def __str__(self):
        self.display()
        return ""

    ### SET ACTION IF CALLED WITH NO COMMAND #################################
    def __repr__(self):
        self.display()
        return ""
    
    ##########################################################################
    ########################## COPY ##########################################
    
    ### COPY .copy() Return unconnected copy of class object.
    def copy(self): return copy.deepcopy(self)

    ##########################################################################
    ########################## DISPLAY #######################################
    
    ### DISPLAY IN NOTEBOOK ###
    def display(self): display(self.df)
    ### DISPLAY PANDAS - SCIENTIFIC NOTATION ###
    def display_scientific(self, dispPos=3): pd.set_option('display.float_format', '{:.'+str(dispPos)+'g}'.format)
    ### DISPLAY PANDAS - FLOAT NOTATION ###
    def display_float(self, dispPos=3): pd.set_option('display.float_format', '{:.5f}'.format)
    ### DISPLAY PANDAS - MAX COLUMN WIDTH (-1 = No Wrap) ###
    def display_colwidth(self, widthIn=-1): pd.set_option('display.max_colwidth', widthIn)
 
    ##########################################################################
    ########################### I/O ##########################################    
    
    ### MAIN LOAD FUNCTION ###
    def load(self):
        ##### DATAFRAME .df #########
        iB=False # BOOL INPUT SUCCESS #
        ### IF DF() CLASS THEN SET AS SELF ###
        if isinstance(self.source, DF): self = copy.deepcopy(self.source); iB=True
        ### IF STRING IMPORT AS TSV PATH ###
        elif isinstance(self.source, str): self.load_tsv(); iB=True 
        ### IF NONE THEN MAKE EMPTY DATAFRAME
        elif self.source is None: 
            self.df = pd.DataFrame(index=self.index, columns=self.columns, dtype=None, copy=False)
            self.df.fillna(0); iB = True
        ### IF INT OR FLOAT THEN FILL DATARAME WITH VALUE ###
        elif (isinstance(self.source, int)) or (isinstance(self.source,float)):
            self.df = pd.DataFrame(index=self.index, columns=self.columns, dtype=None, copy=False)
            self.df.fillna(self.source); iB = True
        ### IF PANDAS DF THEN COPY ###
        elif isinstance(self.source, pd.DataFrame): self.df = self.source.deepcopy(); iB=True 
        ### IF DICT THEN FORMAT AS DF ###
        elif isinstance(self.source, dict): self.load_dict(self.source); iB=True
        ### ELSE TRY IMPORTING AS OBJECT ### 
        else:  self.load_obj(self.source);iB=True
        ### ELSE ISSUE ERROR AND INITIALIZE AS EMPTY
        if not iB: 
            print("  WARNING: Could not import DF from provided source. Initialized empty.")
            self.df = pd.DataFrame(index=indexIn, columns=columnsIn, dtype=None, copy=False);self.df.fillna(0)
        return self
        
    ### DF FROM TSV ### -> self
    def load_tsv(self):
        if self.verbose is True: print("  LOADING: File " + self.source);
        self.df = pd.read_csv(self.source, sep=self.sep, index_col=None, skiprows=0, verbose=False, low_memory=False) 
        if self.verbose is True: print("  SUCCESS: Loaded DF.")
        return self
    
    ### DF FROM DICTIONARY ### -> self
    def load_dict(self, dictIn):
        if self.verbose is True: print("  CREATING: DF from Dictionary.")
        self.df = pd.DataFrame.from_dict(dictIn, orient='columns', dtype=None)
        if self.verbose is True: print("  SUCCESS: Created DF from Dictionary.")
        return self
    
    ### DF FROM RECORD ###  -> self (objIn i.e. ndarray(i.e.structured dtype), list of tuples, dict, or DataFrame
    def load_obj(self, objIn):
        if self.verbose is True: print("  CREATING: DF from Object.")
        self.df = pd.DataFrame.from_records(data=objIn, columns=self.columns, index=self.index, coerce_float=True)
        if self.verbose is True: print("  SUCCESS: Created DF from Object.")
        return self   
    
    ### DF TO TSV FILE ### -> self
    def write(self, outPath):
        if self.verbose is True: print("  WRITING: DF to File "+outPath)
        self.df.to_csv(outPath, sep=self.sep)
        if self.verbose is True: print("  SUCCESS: DF Written to File "+outPath)
        return self
    
    ##########################################################################
    ####################### TABLE MANIPULATION ###############################
    
    ### SET INDEX COLUMN OR COLUMNS ###
    def set_index(self, colNameOrList, drop=True, append=False): self.df.set_index(colNameOrList, drop=drop, append=append, inplace=True, verify_integrity=False)
    
    ### YIELD BY INDICES (rows) ###
    def yield_index(self): 
        for iCur in self.i: yield self.df.loc[[iCur]]
    ### YIELD BY COLUMNS (cols) ###
    def yield_columns(self): 
        for cCur in self.c: yield self.df[[cCur]]
    
    ### REPLACE IN DATAFRAME ###
    def replace(self, toReplace, replaceWith): self.df.replace(to_replace=toReplace, value=replaceWith, inplace=True, limit=None, regex=False, method='pad')
    ### REPLACE STRING WITH REGEX ###
    def replace_regex(self, regexToReplace, replaceWithStr): self.df.replace(to_replace=regexToReplace, value=replaceWithStr, inplace=True, limit=None, regex=True)    
    
    ### RESET .reset()
    def reset(self): return self.__init__(self.source)
    
    ##########################################################################
    ############################ PALETTE #####################################
    ### DISPLAY A PREVIEW OF THE PALETTE ...ooooh pretty ###
    def palette_preview(self): sns.palplot(self.palette)
    ### SET .palette FOR A RAINBOW OF COLORS ###
    def palette_heatmap(self, numColors): self.palette = sns.palplot(sns.color_palette("coolwarm", numColors))
    ### SET .palette FOR A HEATMAP OF COLORS ###
    def palette_rainbow(self, numColors, lightIn=0.5, saturationIn=0.8, previewIn=False): self.palette = sns.hls_palette(numColors, l=lightIn, s=saturationIn)

    ##########################################################################
    ########################### PLOTTING #####################################
    
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
    
    ### LINEAR MODEL PLOT ### continuousXcontinuous, 95%CI, can color and split plots by categories.
    def plot_lm(self, continuousColX, continuousColY, categoricalColor=None, categoricalRow=None, categoricalColumn=None):
        plt.clf(); f = plt.figure(figsize=(5,5))
        ### SET STYLE ###
        sns.set(style="ticks", context="talk"); pal=None
        ### GET COLOR PALETTE ###
        if categoricalColor != None: self.palette_rainbow(len(self.df[categoricalColor].unique()), lightIn=0.7, saturationIn=1.0 )
        ### MAKE LINEAR MODEL PLOT ###
        sns.lmplot(x=continuousColX,y=continuousColY,row=categoricalRow,col=categoricalColumn,hue=categoricalColor,data=self.df, size=7, palette=self.palette, ci=95)
        plt.show(); plt.clf()
    
    ### VIOLIN PLOT ### categoricalXcontinuous, color ideal with two groups.
    def plot_violin(self, categoricalColX, continuousColY, categoricalColor=None):
        f = plt.figure(figsize=(16,10))
        ### IF CATEGORICAL COLOR ON EACH HALF OF VIOLIN ###
        splitIn = False;
        if categoricalColor != None:
            if len(self.df[categoricalColor].unique()) == 2: splitIn = True
            ### GET COLOR PALETTE ###
            self.palette_rainbow(len(self.df[categoricalColor].unique()), lightIn=0.7, saturationIn=1.0 )
        ### MAKE VIOLIN PLOT ###
        return sns.violinplot(x=categoricalColX, y=continuousColY, hue=categoricalColor, data=self.df, split=splitIn, palette=self.palette)
        #sns.despine(left=True)

    ##########################################################################
    ########################## STATS ####################################
    ##### .stats_X STATISTICS FUNCTIONS #####################
    
    ### LINEAR REGRESSION - ORDINARY LEAST SQUARES
    # regResults, regDependent, regPredictors = stats_regression(mapDF, regEquation=" age_years ~ bmi + race + race:sex ")
    def stats_regression(self, regEquation=" age_years ~ bmi + sex:race ", printOut=True):
        ### GET VARIABLES INTO X[n,p] predictors and y[n,1] outcome
        y, X = patsy.dmatrices(regEquation, self.df, return_type='dataframe')
        ### GENERATE OLS REGRESSION MODEL ###
        statsModel = sm.OLS(y, X, disp=0)
        ### FIT DATA TO A LINE USING OLS MINIMIZATION ###
        statsModelFit = statsModel.fit(disp=0)
        ### PRINT RESULTS ###
        if printOut == True: print(statsModelFit.summary())
        ### RETURN: resultsObject, y, X
        return statsModelFit, y, X
    
    """
    ### LOGISTIC REGRESSION
    # regResults, regDependent, regPredictors = stats_regression(mapDF, regEquation=" age_years ~ bmi + race + race:sex ")
    def stats_logistic(self, regEquation=" age_years ~ bmi + sex:race ", printOut=True):
        ### GET VARIABLES INTO X[n,p] predictors and y[n,1] outcome
        y, X = patsy.dmatrices(regEquation, self.df, return_type='matrix')
        ### GENERATE OLS REGRESSION MODEL ###
        statsModel = sm.Logit(y[:, 0], X, disp=0)
        ### FIT DATA TO A LINE USING OLS MINIMIZATION ###
        statsModelFit = statsModel.fit(method='bfgs', disp=0)
        ### PRINT RESULTS ###
        if printOut == True: print(statsModelFit.summary())
        ### RETURN: resultsObject, y, X
        return statsModelFit, y, X
    """
    
    ### MANN-WHITNEY-U ON TWO COLUMNS ###
    def stats_mwu(self,continuousColOne,continuousColTwo):
        testStat, pValue = sp.stats.mannwhitneyu(self.df[continuousColOne], self.df[continuousColTwo], use_continuity=True, alternative='two-sided')
        if self.verbose is True: print(" - Mann-Whitney-U: P="+str('%.2E' % Decimal(pValue))+" TS="+str(testStat)+" [ "+continuousColOne+" : "+continuousColTwo+" ] - ")
        return pValue, testStat 
    
    ### MANN-WHITNEY-U PAIRWISE CONTINUOUS VARIABLE BY CATEGORICAL COLUMN ###
    def stats_mwu_pairwise(self,continuousColumn,categoricalColumn):
        if self.verbose is True: print(" - Comparing groups of "+categoricalColumn+" by "+continuousColumn+" - ")
        resDF = DF(None, columnsIn=['group1','group2','p','t'], indexIn=[])
        ### GET ALL PAIRWISE COMBINATIONS OF DISTANCE MATRICES ### 
        for catIDX, (cat1,cat2) in enumerate(itertools.combinations(self.df[categoricalColumn].unique(),2)):
            ### CALCULATE MANN-WHITNEY-U ###
            testStat, pValue = sp.stats.mannwhitneyu(self.df[self.df[categoricalColumn]==cat1][continuousColumn] , self.df[self.df[categoricalColumn]==cat2][continuousColumn], use_continuity=True, alternative='two-sided')
            if self.verbose is True: 
                if pValue > 0.00001: pPrint = str(pValue)
                else: pPrint = str('%.2E' % Decimal(pValue))
                print(" - Mann-Whitney-U: p="+ pPrint+" | t="+str(testStat)+" [ "+cat1+"(n="+str(len(self.df[self.df[categoricalColumn]==cat1][continuousColumn]))+") : "+cat2+"(n="+str(len(self.df[self.df[categoricalColumn]==cat2][continuousColumn]))+")] - ")
            resDF.df.loc[catIDX] = [cat1,cat2,pValue,testStat]
        return resDF
    
    
    ##########################################################################
    ############################## PLAYPEN ###################################
    
    #### ADD COLUMN TO DATAFRAME ####
    #def add_c(self, newColumn, name=None, overwrite=True):
    #    newColumn = pd.Series(data=newColumn, index=self.df.index, dtype=None, name=name, copy=False)
    #    self.df = pd.concat([self.df,newColumn], axis=0, join='outer', join_axes=useIndex, ignore_index=overwrite,
    #                            keys=None, levels=None, names=None, verify_integrity=False,copy=True)
    #### ADD ROW TO DATAFRAME ###
    #def add_r(self, newRow, name=None, overwrite=True):    
    #    newRow = pd.Series(data=newRow, index=self.df.columns, dtype=None, name=name, copy=False)
    #    self.df = pd.concat([self.df,newRow], axis=0, join='outer', join_axes=None, ignore_index=overwrite,
    #                            keys=None, levels=None, names=None, verify_integrity=False,copy=True)
    
    #### MERGE DATAFRAME CLASS OBJECT INTO CURRENT OBJECT ###
    #def merge(self, inDF, useIndex=None, unionIn=False, replaceColnames=False):
    #    if unionIn is true: joinIn = 'inner'
    #    else: joinIn = 'outer'
    #    return DF(pd.concat([self.df,inDF], axis=1, join=joinIn, join_axes=useIndex, ignore_index=overwrite,
     
    ##########################################################################
    ############################ TODO ########################################
    
    ### MISSINGNESS TOOLS - BEAUTIFUL FIGURES 
    #import missingno as msno
    #msno.matrix(collisions.sample(250))
        
### EXAMPLE USAGE ###
#X = DF("test_data/ag_analysis/1_1_qc_1000_map.txt")    
#X = DF(None, columnsIn=["Ace","Duce","Cinco"], indexIn=np.arange(6))
### PERFORM MANN-WHITNEY-U ON PAIRSWISE COMBINATIONS OF AGE BY RACE ###
# X.stats_mwu_pairwise('age_years','race')
#X.df
### PLOT A LINEAR MODEL ###
#f = X.plot_violin(categoricalColX='race',continuousColY='age_years',categoricalColor='economic_region')
#g = X.plot_lm(continuousColX='age_years',continuousColY='bmi',categoricalColor='economic_region')#,categoricalColumn='race',categoricalRow='sex')
#X.stats_regression("bmi ~ age_years + race + sex")[0]
#X.set_c(colName='winning',newCol=np.arange(len(X.i))).d()
#X.set_r(index='winning',newCol=).d()
#X.set_r(index=6, newRow = ['Nan','Nan','Nan'])
#X.update()

### LOAD MAPPING FILE & SET SAMPLID TO INDEX ###
#mapX = DF(mapPath,sepIn='\t',verboseIn=True)
#mapX.set_index("#SampleID")

### EXAMPLE USAGE ###
# SET DISPLAY TO FLOAT #
#mapX.display_float()

# PERFORM PAIRWISE MANN-WHITNEY-U #
#mwuRes = mapX.stats_mwu_pairwise(categoricalColumn='race',continuousColumn='bmi')
#mwuRes.out_tsv('t98.txt'); mwuRes.display()

# PERFORM REGRESSION #
#mapX.stats_regression("age_years ~ race")

# PLOT VIOLIN #
#mapX.plot_violin(categoricalColX='sex',continuousColY='bmi',categoricalColor='race')

# PLOT LINEAR MODEL #
#xLM = mapX.plot_lm(continuousColX='age_years',continuousColY='bmi',categoricalColor='economic_region')
#xLM.savefig('t99.pdf')

# x.s(0, index=[2,3,6],column=['Duce','Ave']).df

##################################################################################
################################### MAIN #########################################
##################################################################################
def main():
    return
if __name__ == '__main__': main()
