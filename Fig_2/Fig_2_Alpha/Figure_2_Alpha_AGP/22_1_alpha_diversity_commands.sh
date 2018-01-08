#!/bin/bash

### Calculate Alpha Diversity ###
parallel_alpha_diversity.py -i 23_0_multiple_rarefactions_1000/ -o 22_1_alpha_diversity_1000/ -O 8 -m chao1,shannon,observed_otus,simpson,equitability -t 00_raw_97_otus.tree &
parallel_alpha_diversity.py -i 23_0_multiple_rarefactions_10000/ -o 22_1_alpha_diversity_10000/ -O 8 -m chao1,shannon,observed_otus,simpson,equitability -t 00_raw_97_otus.tree &

### Collate Alpha Diversity ###
collate_alpha.py -i 22_1_alpha_diversity_1000 -o 22_2_collate_alpha_1000 &
collate_alpha.py -i 22_1_alpha_diversity_10000 -o 22_2_collate_alpha_10000 &

### Compare Alpha Diversity Across Groups ###
mkdir 22_3_compare_alpha_diversity_1000/
compare_alpha_diversity.py -i 22_2_collate_alpha_1000/shannon.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_1000/22_3_compare_alpha_diversity_shannon_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000/chao1.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_1000/22_3_compare_alpha_diversity_chao1_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000/observed_otus.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_1000/22_3_compare_alpha_diversity_observed_otus_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000/equitability.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_1000/22_3_compare_alpha_diversity_equitability_1000 -d 1000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_1000/simpson.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_1000/22_3_compare_alpha_diversity_simpson_1000 -d 1000 &

mkdir 22_3_compare_alpha_diversity_10000/
compare_alpha_diversity.py -i 22_2_collate_alpha_10000/shannon.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_10000/22_3_compare_alpha_diversity_shannon_10000 -d 10000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_10000/chao1.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_10000/22_3_compare_alpha_diversity_chao1_10000 -d 10000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_10000/observed_otus.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_10000/22_3_compare_alpha_diversity_observed_otus_10000 -d 10000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_10000/equitability.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_10000/22_3_compare_alpha_diversity_equitability_10000 -d 10000 &
compare_alpha_diversity.py -i 22_2_collate_alpha_10000/simpson.txt -m 12_map_age_categories_custom.txt -c race,sex,bmi_cat,age_cat,age_group_custom,age_group_custom2 -o 22_3_compare_alpha_diversity_10000/22_3_compare_alpha_diversity_simpson_10000 -d 10000 &

