#!/bin/bash

# Andrew W. Brooks
# Vanderbilt Genetics Institute

### MIT License ###########################################################################

# Copyright (c) 2016, 2017 Andrew W. Brooks

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

###########################################################################################

###########################################################################################
##### Data/1_Stata_Filtering - STATA QUALITY CONTROL ###############################################################

### INITIAL QUALITY CONTROL OF THE MAPPING FILE WAS CARRIED OUT USING STATA
### WHICH CAN BE REPRODUCED BY RUNNING THE 10_metadata_analysis.do STATA .DO SCRIPT
### THE ONLY CHANGE WHICH MAY NEED TO BE EXECUTED IS APPENDING THE PATH TO THE 
### RAW AMERICAN GUT MAPPING FILE ON LINE 9

### THE 10_metadata_analysis.do FILE IN 1_Stata_Filtering FOLDER CONTAINS
### ALL COMMANDS USED TO EXAMINE THE RAW AMERICAN GUT MAPPING FILE 00_raw_ag-cleaned.txt
### AND SHOWS ALL OUTPUT WHICH I HAVE INCLUDED BUT LEFT COMMENTED TO ALLOW THE SCRIPT TO REMAIN EXECUTABLE

### THIS STATA SCRIPT WILL OUTPUT THE 12_map_trim_basic.csv FILE WHICH HAS HAD SAMPLES
### EXCLUDED REPRESENTING THE FOLLOWING QUALITY CONTROL FROM THE MANUSCRIPT:
### remove samples (Raw N=9,475): with BMI more than 60 (-988 [8,487]) or less than 10 (-68 [8,419]), 
### missing age (-661 [7,758]), 
### with age greater than 55 years old (-2,777 [4,981]) or less than 18 years old (-582 [4,399]), 
### and blank samples or those not appearing in the mapping file (-482 [3,917])

### MICROSOFT EXCEL WAS USED TO EXPORT 12_map_trim_basic.csv AS A TAB DELIMITED TEXT FILE
### ALTERNATIVELY ONE COULD CONVERT TO A TSV USING THE FOLLOWING FUNCTION
### CSV TO TSV ###################################
# This simple function will replace commas with tabs in csv (or any) file
# $1 = input csv file
# $2 = output tsv file
csvtotsv(){
    cat $1 | tr "," "\\t" > $2
}
### SIMPLY BY EXECUTING THE FOLLOWING COMMAND IN THE TERMINAL:
# csvtotsv 12_map_trim_basic.csv 12_map_trim_basic.txt

### JUST THE SAMPLE COLUMN WAS EXTRACTED FROM 12_map_trim_basic.csv AND USED TO GENERATE
### THE FILE 11_samples_to_keep.txt WHICH WAS USED IN THE FOLLOWING STEPS...

###########################################################################################
##### Data/2_Qiime_Filtering - QIIME QUALITY CONTROL ###############################################################

### ALL SUBSEQUENT STEPS CAN BE RUN IN THE 2_Qiime_Filtering FOLDER IN THE QIIME 1.9 ENVIRONMENT

### FILTER SAMPLES FROM OTU TABLE USING LIST OF SAMPLES REMAINING AFTER QUALITY CONTROL IN STATA ###
filter_samples_from_otu_table.py -i 01_otu_table.biom -o 11_filter_samples_from_otu_table.biom --sample_id_fp=11_samples_to_keep.txt

### SAVE RESOURCE PATHS ###
bt=11_filter_samples_from_otu_table.biom
map=12_map_trim_basic.txt
tre=00_raw_97_otus.tree

### FILTER UNKNOWN / OTHER != RACE ###
filter_samples_from_otu_table.py -i $bt -o 20_0_filter_race.biom -m $map -s 'race:*,!Unknown,!Other'
biom summarize-table -i 20_0_filter_race.biom -o 20_0_filter_race_summary.txt &
biom summarize-table -i 20_0_filter_race.biom -o 20_0_filter_race_summary_qualitative.txt --qualitative &
biom convert -i 20_0_filter_race.biom -o 20_0_filter_race.txt --to-tsv --header-key taxonomy &

### FILTER ONLY FECAL BODY SITE ###
filter_samples_from_otu_table.py -i 20_0_filter_race.biom -o 20_0_filter_simple_body_site.biom -m $map -s 'simple_body_site:FECAL'
biom summarize-table -i 20_0_filter_simple_body_site.biom -o 20_0_filter_simple_body_site_summary.txt &
biom summarize-table -i 20_0_filter_simple_body_site.biom -o 20_0_filter_simple_body_site__summary_qualitative.txt --qualitative &
biom convert -i 20_0_filter_simple_body_site.biom -o 20_0_filter_simple_body_site.txt --to-tsv --header-key taxonomy &

### FILTER UNKNOWN / OTHER != SEX ###
filter_samples_from_otu_table.py -i 20_0_filter_simple_body_site.biom -o 20_1_filter_sex.biom -m $map -s 'sex:*,!Unknown,!Other,!other'
biom summarize-table -i 20_1_filter_sex.biom -o 20_1_filter_sex_summary.txt &
biom summarize-table -i 20_1_filter_sex.biom -o 20_1_filter_sex_summary_qualitative.txt --qualitative &
biom convert -i 20_1_filter_sex.biom -o 20_1_filter_sex.txt --to-tsv --header-key taxonomy &

### FILTER USA == COUNTRY ###
filter_samples_from_otu_table.py -i 20_1_filter_sex.biom -o 20_2_filter_country.biom -m $map -s 'country:USA'
biom summarize-table -i 20_2_filter_country.biom -o 20_2_filter_country_summary.txt &
biom summarize-table -i 20_2_filter_country.biom -o 20_2_filter_country_summary_qualitative.txt --qualitative &
biom convert -i 20_2_filter_country.biom -o 20_2_filter_country.txt --to-tsv --header-key taxonomy &

### CONVERT TO TXT WITH TAXONOMY DATA ###
biom convert -i 20_2_filter_country.biom -o 20_3_filter_final_taxon.txt --to-tsv --header-key taxonomy

### ADD METADATA TO BIOM TABLE AND CONVERT TO TXT ###
biom add-metadata -i 20_2_filter_country.biom -o 20_4_filter_with_metadata.biom --sample-metadata-fp 12_map_trim_basic.txt 
biom convert -i 20_4_filter_with_metadata.biom -o 20_4_filter_with_metadata.txt --to-tsv --header-key taxonomy

###########################################################################################
##### 3_Custom_Filtering - nQIIME QUALITY CONTROL #########################################

######################################################################
##### CUSTOM FILTERING IN SCRIPTS REPRODUCED IN QIIME #####
# THIS WAS BOTH FOR VALIDATION, AND TAXA SUMMARIES WERE USED IN THE 
# CHRISTENSENELLACEAE ANALYSIS.

##### DEPTH 1,000 FOR SAMPLES #####
### FILTER OTUS FROM OTU TABLE WITH < 10 COUNTS ###
filter_otus_from_otu_table.py -n 10 -i 20_2_filter_country.biom -o 20_3_filter_otus_10_1000.biom

### FILTER SAMPLES FROM OTU TABLE WITH < 1,000 COUNTS ###
filter_samples_from_otu_table.py -n 1000 -i 20_3_filter_otus_10_1000.biom -o 20_4_filter_samples_1000.biom

### REFILTER OTUS FROM OTU TABLE WITH < 10 COUNTS ###
filter_otus_from_otu_table.py -n 10 -i 20_4_filter_samples_1000.biom -o 20_5_filter_otus_10_1000.biom

### FILTER SAMPLES FROM OTU TABLE WITH < 1,000 COUNTS ###
filter_samples_from_otu_table.py -n 1000 -i 20_5_filter_otus_10_1000.biom -o 20_6_final_table_1000.biom


##### DEPTH 10,000 FOR SAMPLES #####
### FILTER OTUS FROM OTU TABLE WITH < 10 COUNTS ###
filter_otus_from_otu_table.py -n 10 -i 20_2_filter_country.biom -o 20_3_filter_otus_10_10000.biom

### FILTER SAMPLES FROM OTU TABLE WITH < 1,000 COUNTS ###
filter_samples_from_otu_table.py -n 10000 -i 20_3_filter_otus_10_10000.biom -o 20_4_filter_samples_10000.biom

### REFILTER OTUS FROM OTU TABLE WITH < 10 COUNTS ###
filter_otus_from_otu_table.py -n 10 -i 20_4_filter_samples_10000.biom -o 20_5_filter_otus_10_10000.biom

### FILTER SAMPLES FROM OTU TABLE WITH < 1,000 COUNTS ###
filter_samples_from_otu_table.py -n 10000 -i 20_5_filter_otus_10_10000.biom -o 20_6_final_table_10000.biom


##### SUMMARIZE TAXONOMY #####
summarize_taxa.py -i 20_6_final_table_1000.biom -o 20_7_summarize_taxa_1000_relative/ &
summarize_taxa.py -i 20_6_final_table_1000.biom -o 20_7_summarize_taxa_1000_absolute/ -a &
summarize_taxa.py -i 20_6_final_table_10000.biom -o 20_7_summarize_taxa_10000_relative/ &
summarize_taxa.py -i 20_6_final_table_10000.biom -o 20_7_summarize_taxa_10000_absolute/ -a &


##### GROUP SIGNIFICANCE #####
# TEST IF TAXA VARY BY METADATA CATEGORIES #

for taxL in {2..6}
do
    group_significance.py -i "20_7_summarize_taxa_1000_relative/20_6_final_table_1000_L"$taxL".biom" -o "21_0_group_significance_1000_relative_L"$taxL"_race.txt" -m 12_map.txt -c race &
done

###########################################################################################
##### SUMMARIZE TAXA ######################################################################

### SUMMARIZE TAXA THROUGH PLOTS - GET RELATIVE ABUNDANCE PLOTS BY RACE ###
summarize_taxa_through_plots.py -i 20_2_filter_country.biom -o 20_5_summarize_taxa_through_plots -m 12_map_trim_basic.txt -c race -s

###########################################################################################
##### BETA DIVERSITY ######################################################################

### PERFORM 100 RAREFACTIONS OF THE BIOM TABLE AT A DEPTH OF 1,000 & 10,000 COUNTS ###
parallel_multiple_rarefactions.py -i 20_2_filter_country.biom -o 23_0_multiple_rarefactions_1000 -m 1000 -x 1000 -n 100 -O 8 -Z 500 -s 5000
parallel_multiple_rarefactions.py -i 20_2_filter_country.biom -o 23_0_multiple_rarefactions_10000 -m 10000 -x 10000 -n 100 -O 8 -Z 500 -s 5000

##### COMPUTE BRAY CURTIS BETA DIVERISTY #####
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_1000/ -o 24_0_beta_diversity_1000_bray_curtis/ -m bray_curtis -t $tre -O 8 &
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_10000/ -o 24_0_beta_diversity_10000_bray_curtis/ -m bray_curtis -t $tre -O 8 &

##### COMPUTE UNWEIGHTED UNIFRAC BETA DIVERISTY #####
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_1000/ -o 24_0_beta_diversity_1000_unweighted_unifrac/ -m unweighted_unifrac -t $tre -O 8 &
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_10000/ -o 24_0_beta_diversity_10000_unweighted_unifrac/ -m unweighted_unifrac -t $tre -O 8 &

##### COMPUTE WEIGHTED UNIFRAC BETA DIVERISTY #####
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_1000/ -o 24_0_beta_diversity_1000_weighted_unifrac/ -m weighted_unifrac -t $tre -O 8 &
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_10000/ -o 24_0_beta_diversity_10000_weighted_unifrac/ -m weighted_unifrac -t $tre -O 8 &

##### COMPUTE BINARY JACCARD BETA DIVERISTY #####
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_1000/ -o 24_0_beta_diversity_1000_binary_jaccard/ -m binary_jaccard -t $tre -O 8 &
parallel_beta_diversity.py -i 23_0_multiple_rarefactions_10000/ -o 24_0_beta_diversity_10000_binary_jaccard/ -m binary_jaccard -t $tre -O 8 &

###########################################################################################
##### ANOSIM & CLUSTER QUALITY ############################################################

### CALCULATE ANOSIM TESTS FOR EACH RAREFIED BIOM TABLE - BRAY CURTIS 1,000 ###
for i in 24_0_beta_diversity_1000_bray_curtis/*
do
 echo $i
 pathName=`dirname $i`; fileName=`basename $i`; noExtension="${fileName%.*}";
 compare_categories.py -i $i -m $map -c race -o "24_2_compare_categories_1000_bray_curtis_race/"$noExtension"/" --method=anosim
done

### CALCULATE CLUSTER QUALITY ON EACH BIOM TABLE - BRAY CURTIS 1,000 ###
mkdir 24_1_cluster_quality_1000_bray_curtis
for i in 24_0_beta_diversity_1000_bray_curtis/*
do
 echo $i
 pathName=`dirname $i`
 fileName=`basename $i`
 noExtension="${fileName%.*}"
 cluster_quality.py -i $i -m $map -c race -o "24_1_cluster_quality_1000_bray_curtis/"$noExtension".txt"
done


### REMOVE EACH RACE AND PERFORM ANOSIM ###

map=12_map_trim_basic.txt
mkdir "24_7_remove_races_"$rare"_"$metric"/"
mkdir "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander"
mkdir "24_7_remove_races_"$rare"_"$metric"/Caucasian"
mkdir "24_7_remove_races_"$rare"_"$metric"/Hispanic"
mkdir "24_7_remove_races_"$rare"_"$metric"/African_American"

mkdir "24_8_compare_categories_remove_race_"$rare"_"$metric"/"
mkdir "24_8_compare_categories_remove_race_"$rare"_"$metric"/Asian_or_Pacific_Islander"
mkdir "24_8_compare_categories_remove_race_"$rare"_"$metric"/Caucasian"
mkdir "24_8_compare_categories_remove_race_"$rare"_"$metric"/Hispanic"
mkdir "24_8_compare_categories_remove_race_"$rare"_"$metric"/African_American"

rare=10000
metric=bray_curtis
for i in {0..99}
do
 inPath="24_0_beta_diversity_"$rare"_"$metric"/"$metric"_rarefaction_"$rare"_"$i".txt"
 echo $inPath
 # FILTER ASIAN PACIFIC ISLANDER #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Asian or Pacific Islander'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Asian_or_Pacific_Islander/"$i --method=anosim
 # FILTER CAUCASIAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Caucasian'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Caucasian/"$i --method=anosim
 # FILTER HISPANIC #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Hispanic'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Hispanic/"$i --method=anosim
 # FILTER AFRICAN AMERICAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!African American'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/African_American/"$i --method=anosim
done

rare=10000
metric=unweighted_unifrac
for i in {0..99}
do
 inPath="24_0_beta_diversity_"$rare"_"$metric"/"$metric"_rarefaction_"$rare"_"$i".txt"
 echo $inPath
 # FILTER ASIAN PACIFIC ISLANDER #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Asian or Pacific Islander'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Asian_or_Pacific_Islander/"$i --method=anosim
 # FILTER CAUCASIAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Caucasian'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Caucasian/"$i --method=anosim
 # FILTER HISPANIC #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Hispanic'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Hispanic/"$i --method=anosim
 # FILTER AFRICAN AMERICAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!African American'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/African_American/"$i --method=anosim
done

rare=10000
metric=weighted_unifrac
for i in {0..99}
do
 inPath="24_0_beta_diversity_"$rare"_"$metric"/"$metric"_rarefaction_"$rare"_"$i".txt"
 echo $inPath
 # FILTER ASIAN PACIFIC ISLANDER #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Asian or Pacific Islander'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Asian_or_Pacific_Islander/"$i --method=anosim
 # FILTER CAUCASIAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Caucasian'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Caucasian/"$i --method=anosim
 # FILTER HISPANIC #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Hispanic'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Hispanic/"$i --method=anosim
 # FILTER AFRICAN AMERICAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!African American'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/African_American/"$i --method=anosim
done

rare=1000
metric=bray_curtis
for i in {0..99}
do
 inPath="24_0_beta_diversity_"$rare"_"$metric"/"$metric"_rarefaction_"$rare"_"$i".txt"
 echo $inPath
 # FILTER ASIAN PACIFIC ISLANDER #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Asian or Pacific Islander'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Asian_or_Pacific_Islander/"$i --method=anosim
 # FILTER CAUCASIAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Caucasian'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Caucasian/"$i --method=anosim
 # FILTER HISPANIC #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Hispanic'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Hispanic/"$i --method=anosim
 # FILTER AFRICAN AMERICAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!African American'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/African_American/"$i --method=anosim
done

rare=1000
metric=unweighted_unifrac
for i in {0..99}
do
 inPath="24_0_beta_diversity_"$rare"_"$metric"/"$metric"_rarefaction_"$rare"_"$i".txt"
 echo $inPath
 # FILTER ASIAN PACIFIC ISLANDER #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Asian or Pacific Islander'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Asian_or_Pacific_Islander/"$i --method=anosim
 # FILTER CAUCASIAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Caucasian'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Caucasian/"$i --method=anosim
 # FILTER HISPANIC #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Hispanic'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Hispanic/"$i --method=anosim
 # FILTER AFRICAN AMERICAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!African American'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/African_American/"$i --method=anosim
done

rare=1000
metric=weighted_unifrac
for i in {0..99}
do
 inPath="24_0_beta_diversity_"$rare"_"$metric"/"$metric"_rarefaction_"$rare"_"$i".txt"
 echo $inPath
 # FILTER ASIAN PACIFIC ISLANDER #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Asian or Pacific Islander'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Asian_or_Pacific_Islander/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Asian_or_Pacific_Islander/"$i --method=anosim
 # FILTER CAUCASIAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Caucasian'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Caucasian/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Caucasian/"$i --method=anosim
 # FILTER HISPANIC #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!Hispanic'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/Hispanic/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/Hispanic/"$i --method=anosim
 # FILTER AFRICAN AMERICAN #
 filter_distance_matrix.py -i $inPath -o "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -s 'race:*,!African American'
 compare_categories.py -i "24_7_remove_races_"$rare"_"$metric"/African_American/filtered_"$rare"_"$metric"_"$i".txt" -m $map -c race -o "24_8_compare_categories_remove_race_"$rare"_"$metric"/African_American/"$i --method=anosim
done

###########################################################################################
##### PCoA ################################################################################

### PERFORM PCoA ORDINATION ON RAREFIED TABLES ###
principal_coordinates.py -i 24_0_beta_diversity_1000_bray_curtis/ -o 24_3_principle_coordinates_1000_bray_curtis/ &

### MAKE EMPEROR PCoA PLOTS ###
make_emperor.py -i 24_3_principle_coordinates_1000_bray_curtis/ -m 12_map_trim_basic.txt -b race,sex --ignore_missing_samples --number_of_segments=14 -o 24_3_make_emperor_1000_bray_curtis &

### PCOA WITHOUT EACH RACE ###
principal_coordinates.py -i 24_7_remove_races_1000_bray_curtis/Asian_or_Pacific_Islander/ -o 24_9_remove_races_PCoA_1000_bray_curtis_Asian_or_Pacific_Islander/ &
principal_coordinates.py -i 24_7_remove_races_1000_bray_curtis/African_American/ -o 24_9_remove_races_PCoA_1000_bray_curtis_African_American/ &
principal_coordinates.py -i 24_7_remove_races_1000_bray_curtis/Caucasian/ -o 24_9_remove_races_PCoA_1000_bray_curtis_Caucasian/ &
principal_coordinates.py -i 24_7_remove_races_1000_bray_curtis/Hispanic/ -o 24_9_remove_races_PCoA_1000_bray_curtis_Hispanic/ &

principal_coordinates.py -i 24_7_remove_races_10000_bray_curtis/Asian_or_Pacific_Islander/ -o 24_9_remove_races_PCoA_10000_bray_curtis_Asian_or_Pacific_Islander/ &
principal_coordinates.py -i 24_7_remove_races_10000_bray_curtis/African_American/ -o 24_9_remove_races_PCoA_10000_bray_curtis_African_American/ &
principal_coordinates.py -i 24_7_remove_races_10000_bray_curtis/Caucasian/ -o 24_9_remove_races_PCoA_10000_bray_curtis_Caucasian/ &
principal_coordinates.py -i 24_7_remove_races_10000_bray_curtis/Hispanic/ -o 24_9_remove_races_PCoA_10000_bray_curtis_Hispanic/ &

principal_coordinates.py -i 24_7_remove_races_1000_unweighted_unifrac/Asian_or_Pacific_Islander/ -o 24_9_remove_races_PCoA_1000_unweighted_unifrac_Asian_or_Pacific_Islander/ &
principal_coordinates.py -i 24_7_remove_races_1000_unweighted_unifrac/African_American/ -o 24_9_remove_races_PCoA_1000_unweighted_unifrac_African_American/ &
principal_coordinates.py -i 24_7_remove_races_1000_unweighted_unifrac/Caucasian/ -o 24_9_remove_races_PCoA_1000_unweighted_unifrac_Caucasian/ &
principal_coordinates.py -i 24_7_remove_races_1000_unweighted_unifrac/Hispanic/ -o 24_9_remove_races_PCoA_1000_unweighted_unifrac_Hispanic/ &

principal_coordinates.py -i 24_7_remove_races_10000_unweighted_unifrac/Asian_or_Pacific_Islander/ -o 24_9_remove_races_PCoA_10000_unweighted_unifrac_Asian_or_Pacific_Islander/ &
principal_coordinates.py -i 24_7_remove_races_10000_unweighted_unifrac/African_American/ -o 24_9_remove_races_PCoA_10000_unweighted_unifrac_African_American/ &
principal_coordinates.py -i 24_7_remove_races_10000_unweighted_unifrac/Caucasian/ -o 24_9_remove_races_PCoA_10000_unweighted_unifrac_Caucasian/ &
principal_coordinates.py -i 24_7_remove_races_10000_unweighted_unifrac/Hispanic/ -o 24_9_remove_races_PCoA_10000_unweighted_unifrac_Hispanic/ &

principal_coordinates.py -i 24_7_remove_races_1000_weighted_unifrac/Asian_or_Pacific_Islander/ -o 24_9_remove_races_PCoA_1000_weighted_unifrac_Asian_or_Pacific_Islander/ &
principal_coordinates.py -i 24_7_remove_races_1000_weighted_unifrac/African_American/ -o 24_9_remove_races_PCoA_1000_weighted_unifrac_African_American/ &
principal_coordinates.py -i 24_7_remove_races_1000_weighted_unifrac/Caucasian/ -o 24_9_remove_races_PCoA_1000_weighted_unifrac_Caucasian/ &
principal_coordinates.py -i 24_7_remove_races_1000_weighted_unifrac/Hispanic/ -o 24_9_remove_races_PCoA_1000_weighted_unifrac_Hispanic/ &

principal_coordinates.py -i 24_7_remove_races_10000_weighted_unifrac/Asian_or_Pacific_Islander/ -o 24_9_remove_races_PCoA_10000_weighted_unifrac_Asian_or_Pacific_Islander/ &
principal_coordinates.py -i 24_7_remove_races_10000_weighted_unifrac/African_American/ -o 24_9_remove_races_PCoA_10000_weighted_unifrac_African_American/ &
principal_coordinates.py -i 24_7_remove_races_10000_weighted_unifrac/Caucasian/ -o 24_9_remove_races_PCoA_10000_weighted_unifrac_Caucasian/ &
principal_coordinates.py -i 24_7_remove_races_10000_weighted_unifrac/Hispanic/ -o 24_9_remove_races_PCoA_10000_weighted_unifrac_Hispanic/ &

###########################################################################################
##### ALPHA DIVERSITY #####################################################################

### PERFORM RAREFACTIONS AT MULTIPLE DEPTHS FOR ALPHA DIVERSITY ANALAYSES ###
bt=20_2_filter_country_json.biom
parallel_multiple_rarefactions.py -i 20_2_filter_country.biom -o 22_0_multiple_rarefactions_1000_81000_5000 -m 1000 -x 81000 -n 100 -O 8 -Z 500 -s 5000

### CALCULATE ALPHA DIVERSITY WITH A VARIETY OF METRICS ###
parallel_alpha_diversity.py -i 22_0_multiple_rarefactions_1000_81000_5000/ -o 21_1_alpha_diversity_1000_81000_5000 -O 8 -Z 100 -T -m chao1,shannon,observed_otus,observed_species,simpson -t 00_raw_97_otus.tree

### COLLATE ALPHA DIVERSITY RESULTS ACROSS RAREFACTIONS ###
collate_alpha.py -i 22_1_alpha_diversity_1000_81000_5000 -o 22_2_collate_alpha_1000_81000_5000

# COMPARE RACE AND SEX ALPHA DIVERSITY #
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/shannon.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_shannon_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/chao1.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_chao1_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_otus.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_observed_otus_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_species.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_observed_species_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/simpson.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_simpson_1000 -d 1000 &

compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/shannon.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_shannon_11000 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/chao1.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_chao1_11000 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_otus.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_observed_otus_11000 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_species.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_observed_species_11000 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/simpson.txt -m 12_map_trim_basic.txt -c race,sex -o 22_3_compare_alpha_diversity_simpson_11000 -d 11000 &

# COMPARE AGE AND BMI ALPHA DIVERSITY #
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/shannon.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_shannon_1000_age_bmi -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/chao1.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_chao1_1000_age_bmi -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_otus.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_observed_otus_1000_age_bmi -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_species.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_observed_species_1000_age_bmi -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/simpson.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_simpson_1000_age_bmi -d 1000 &

compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/shannon.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_shannon_11000_age_bmi -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/chao1.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_chao1_11000_age_bmi -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_otus.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_observed_otus_11000_age_bmi -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_species.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_observed_species_11000_age_bmi -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/simpson.txt -m 12_map_age_categories_custom.txt -c bmi_cat,age_cat,age_group_custom -o 22_3_compare_alpha_diversity_simpson_11000_age_bmi -d 11000 &

# COMPARE AGE ALTERNATIVE 2 (5 YEAR INCREMENTS) #
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/shannon.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_shannon_1000_age2 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/chao1.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_chao1_1000_age2 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_otus.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_observed_otus_1000_age2 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_species.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_observed_species_1000_age2 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/simpson.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_simpson_1000_age2 -d 1000 &

compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/shannon.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_shannon_11000_age2 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/chao1.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_chao1_11000_age2 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_otus.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_observed_otus_11000_age2 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/observed_species.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_observed_species_11000_age2 -d 11000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000_81000_5000/simpson.txt -m 12_map_age_categories_custom.txt -c age_group_custom2 -o 22_3_compare_alpha_diversity_simpson_11000_age2 -d 11000 &
