clear
set more off
log using "10_stata_log.smcl",replace

*** NAVIGATE TO FOLDER WITH AMERICAN GUT MAPPING FILE ***
* cd /Users/brooks/Work/16_6_22_American_Gut/data/analyses/

*** IMPORT TAB DELIMITED MAPPING FILE ***
import delimited 00_raw_ag-cleaned.txt
* (459 vars, 8,735 obs)

*******************************************************************************
********************** EXAMINING KEY VARIABLES ********************************
*******************************************************************************
drop vioscreen_*
order race, before(non_food_allergies_beestings)
order sex, before(non_food_allergies_beestings)
order age_years, before(non_food_allergies_beestings)
order bmi, before(non_food_allergies_beestings)

*******************************************************************************
*** CLEAN BMI *** deleted observations are done for factor on raw dataset (not all in order)
destring bmi, replace ignore("Unknown")
drop if bmi > 60
* (988 observations deleted)
drop if bmi < 10
* (68 observations deleted)
sum bmi, d
/*
                             BMI
-------------------------------------------------------------
      Percentiles      Smallest
 1%        13.98          11.09
 5%        17.12          11.15
10%        18.73          11.46       Obs               7,679
25%        20.73          11.49       Sum of Wgt.       7,679

50%        23.04                      Mean           23.83447
                        Largest       Std. Dev.      5.146806
75%        25.95          56.38
90%        29.85          58.48       Variance       26.48961
95%        33.28          58.63       Skewness       1.493022
99%        42.94          59.88       Kurtosis       7.942063
*/

hist bmi, frequency
graph save Graph "11_BMI_Histogram_10_60.gph"

*******************************************************************************
*** CLEAN age_years ***
destring age_years, gen(ageF) ignore("Unknown")
* age_years: characters U n k o w removed; ageF generated as int
* (661 missing values generated)
drop if ageF==.
* (661 observations deleted)

hist ageF, frequency

graph save Graph "11_Age_Histogram.gph"
*(file /Users/brooks/Work/16_6_22_American_Gut/data/analyses/11_Age_Histogram.gph saved)

sum ageF , d
/*
                          AGE_YEARS
-------------------------------------------------------------
      Percentiles      Smallest
 1%            2              0
 5%            9              0
10%           22              1       Obs               8,074
25%           34              1       Sum of Wgt.       8,074

50%           48                      Mean           46.07035
                        Largest       Std. Dev.      17.91236
75%           60             94
90%           67             94       Variance       320.8527
95%           71            101       Skewness      -.4825545
99%           79            101       Kurtosis       2.815875
*/

drop if ageF>65
*(1,063 observations deleted)

drop if ageF>55
*(1,714 observations deleted)

drop if ageF<10
*(434 observations deleted)

drop if ageF<18
*(148 observations deleted)






*** TABULATE RACE ***
tabulate race
/*
                     RACE |      Freq.     Percent        Cum.
--------------------------+-----------------------------------
         African American |         70        0.84        0.84
Asian or Pacific Islander |        343        4.10        4.93
                Caucasian |      7,439       88.87       93.80
                 Hispanic |        156        1.86       95.66
                    Other |        199        2.38       98.04
                  Unknown |        164        1.96      100.00
--------------------------+-----------------------------------
                    Total |      8,371      100.00
*/
tabulate sex
/*
        SEX |      Freq.     Percent        Cum.
------------+-----------------------------------
    Unknown |        511        6.10        6.10
     female |      4,036       48.21       54.32
       male |      3,818       45.61       99.93
      other |          6        0.07      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tabulate age_cat
/*
    AGE_CAT |      Freq.     Percent        Cum.
------------+-----------------------------------
        20s |        689        8.23        8.23
        30s |      1,444       17.25       25.48
        40s |      1,479       17.67       43.15
        50s |      1,673       19.99       63.13
        60s |      1,555       18.58       81.71
        70+ |        560        6.69       88.40
    Unknown |        443        5.29       93.69
       baby |          7        0.08       93.78
      child |        357        4.26       98.04
       teen |        164        1.96      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab age_years 
/*
  AGE_YEARS |      Freq.     Percent        Cum.
------------+-----------------------------------
        0.0 |          2        0.02        0.02
        1.0 |         47        0.56        0.59
       10.0 |         25        0.30        0.88
      101.0 |          2        0.02        0.91
       11.0 |         22        0.26        1.17
       12.0 |         11        0.13        1.30
       13.0 |         23        0.27        1.58
       14.0 |         19        0.23        1.80
       15.0 |         17        0.20        2.01
       16.0 |         19        0.23        2.23
       17.0 |         12        0.14        2.38
       18.0 |         27        0.32        2.70
       19.0 |         50        0.60        3.30
        2.0 |         84        1.00        4.30
       20.0 |         39        0.47        4.77
       21.0 |         67        0.80        5.57
       22.0 |         47        0.56        6.13
       23.0 |         55        0.66        6.79
       24.0 |         56        0.67        7.45
       25.0 |         82        0.98        8.43
       26.0 |         77        0.92        9.35
       27.0 |         89        1.06       10.42
       28.0 |         90        1.08       11.49
       29.0 |         91        1.09       12.58
        3.0 |         27        0.32       12.90
       30.0 |         92        1.10       14.00
       31.0 |        126        1.51       15.51
       32.0 |        127        1.52       17.02
       33.0 |        160        1.91       18.93
       34.0 |        172        2.05       20.99
       35.0 |        157        1.88       22.86
       36.0 |        147        1.76       24.62
       37.0 |        170        2.03       26.65
       38.0 |        141        1.68       28.34
       39.0 |        156        1.86       30.20
        4.0 |         49        0.59       30.78
       40.0 |        143        1.71       32.49
       41.0 |        133        1.59       34.08
       42.0 |        147        1.76       35.84
       43.0 |        149        1.78       37.62
       44.0 |        135        1.61       39.23
       45.0 |        148        1.77       41.00
       46.0 |        164        1.96       42.96
       47.0 |        153        1.83       44.79
       48.0 |        182        2.17       46.96
       49.0 |        126        1.51       48.46
        5.0 |         61        0.73       49.19
       50.0 |        148        1.77       50.96
       51.0 |        160        1.91       52.87
       52.0 |        168        2.01       54.88
       53.0 |        192        2.29       57.17
       54.0 |        172        2.05       59.23
       55.0 |        177        2.11       61.34
       56.0 |        152        1.82       63.16
       57.0 |        183        2.19       65.34
       58.0 |        146        1.74       67.09
       59.0 |        178        2.13       69.22
        6.0 |         47        0.56       69.78
       60.0 |        212        2.53       72.31
       61.0 |        163        1.95       74.26
       62.0 |        188        2.25       76.50
       63.0 |        187        2.23       78.74
       64.0 |        126        1.51       80.24
       65.0 |        179        2.14       82.38
       66.0 |        151        1.80       84.18
       67.0 |        144        1.72       85.90
       68.0 |        105        1.25       87.16
       69.0 |        102        1.22       88.38
        7.0 |         37        0.44       88.82
       70.0 |         87        1.04       89.86
       71.0 |         80        0.96       90.81
       72.0 |         75        0.90       91.71
       73.0 |         76        0.91       92.62
       74.0 |         50        0.60       93.21
       75.0 |         37        0.44       93.66
       76.0 |         27        0.32       93.98
       77.0 |         27        0.32       94.30
       78.0 |         14        0.17       94.47
       79.0 |         12        0.14       94.61
        8.0 |         28        0.33       94.95
       80.0 |         28        0.33       95.28
       81.0 |          5        0.06       95.34
       82.0 |          7        0.08       95.42
       83.0 |          4        0.05       95.47
       84.0 |          8        0.10       95.57
       85.0 |          6        0.07       95.64
       86.0 |          4        0.05       95.69
       87.0 |          1        0.01       95.70
       88.0 |          2        0.02       95.72
       89.0 |          2        0.02       95.75
        9.0 |         52        0.62       96.37
       91.0 |          1        0.01       96.38
       92.0 |          2        0.02       96.40
       94.0 |          4        0.05       96.45
    Unknown |        297        3.55      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/

*******************************************************************************
*** TABULATE RACE BY SEX ***
 tabulate race sex
/*
                      |                     SEX
                 RACE |   Unknown     female       male      other |     Total
----------------------+--------------------------------------------+----------
     African American |         4         42         24          0 |        70 
Asian or Pacific Is.. |         5        155        182          1 |       343 
            Caucasian |       395      3,643      3,397          4 |     7,439 
             Hispanic |         5         67         83          1 |       156 
                Other |         9         91         99          0 |       199 
              Unknown |        93         38         33          0 |       164 
----------------------+--------------------------------------------+----------
                Total |       511      4,036      3,818          6 |     8,371 
*/
tab bmi_cat 
/*
    BMI_CAT |      Freq.     Percent        Cum.
------------+-----------------------------------
     Normal |      4,533       54.15       54.15
      Obese |        786        9.39       63.54
 Overweight |      1,706       20.38       83.92
Underweight |        758        9.06       92.98
    Unknown |        588        7.02      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab body_habitat 
/*
      BODY_HABITAT |      Freq.     Percent        Cum.
-------------------+-----------------------------------
      UBERON:feces |      7,481       89.37       89.37
       UBERON:nose |         10        0.12       89.49
UBERON:oral cavity |        543        6.49       95.97
       UBERON:skin |        337        4.03      100.00
-------------------+-----------------------------------
             Total |      8,371      100.00
*/
tab body_product 
/*
 BODY_PRODUCT |      Freq.     Percent        Cum.
--------------+-----------------------------------
 UBERON:feces |      7,794       89.23       89.23
 UBERON:mucus |         10        0.11       89.34
UBERON:saliva |        572        6.55       95.89
 UBERON:sebum |        359        4.11      100.00
--------------+-----------------------------------
        Total |      8,735      100.00
*/
tab body_site
/*
          BODY_SITE |      Freq.     Percent        Cum.
--------------------+-----------------------------------
       UBERON:feces |      7,481       89.37       89.37
     UBERON:nostril |         10        0.12       89.49
UBERON:skin of hand |        174        2.08       91.57
UBERON:skin of head |        163        1.95       93.51
      UBERON:tongue |        543        6.49      100.00
--------------------+-----------------------------------
              Total |      8,371      100.00
*/





*******************************************************************************
******************** EXAMINING SECONDARY VARIABLES ****************************
*******************************************************************************
 tabulate acid_reflux 
/*
                            ACID_REFLUX |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |         80        0.96        0.96
           I do not have this condition |        420        5.02        5.97
                         Self-diagnosed |         29        0.35        6.32
                                Unknown |        842       10.06       16.38
                                no_data |      7,000       83.62      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tabulate acne_medication 
/*
ACNE_MEDICA |
       TION |      Freq.     Percent        Cum.
------------+-----------------------------------
    Unknown |        297        3.55        3.55
      false |      7,877       94.10       97.65
       true |        197        2.35      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tabulate acne_medication_otc 
/*
ACNE_MEDICA |
   TION_OTC |      Freq.     Percent        Cum.
------------+-----------------------------------
    Unknown |        264        3.15        3.15
      false |      7,482       89.38       92.53
       true |        625        7.47      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
 tabulate add_adhd 
/*
                               ADD_ADHD |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |        186        2.22        2.22
Diagnosed by an alternative medicine .. |         10        0.12        2.34
           I do not have this condition |      3,738       44.65       47.00
                         Self-diagnosed |        148        1.77       48.76
                                Unknown |      4,289       51.24      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab alcohol_consumption 
/*
ALCOHOL_CON |
   SUMPTION |      Freq.     Percent        Cum.
------------+-----------------------------------
    Unknown |        196        2.34        2.34
      false |      1,978       23.63       25.97
       true |      6,197       74.03      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab alcohol_consumption race
/*
ALCOHOL_CO |                               RACE
 NSUMPTION | African..  Asian o..  Caucasian   Hispanic      Other    Unknown |     Total
-----------+------------------------------------------------------------------+----------
   Unknown |         2          3         66          3          4        118 |       196 
     false |        32        150      1,671         45         66         14 |     1,978 
      true |        36        190      5,702        108        129         32 |     6,197 
-----------+------------------------------------------------------------------+----------
     Total |        70        343      7,439        156        199        164 |     8,371 
*/
tab alcohol_frequency 
/*
            ALCOHOL_FREQUENCY |      Freq.     Percent        Cum.
------------------------------+-----------------------------------
                        Daily |        913       10.91       10.91
                        Never |      1,978       23.63       34.54
Occasionally (1-2 times/week) |      1,635       19.53       54.07
   Rarely (a few times/month) |      2,012       24.04       78.10
   Regularly (3-5 times/week) |      1,637       19.56       97.66
                      Unknown |        196        2.34      100.00
------------------------------+-----------------------------------
                        Total |      8,371      100.00
*/
tab alcohol_frequency race 
/*
                      |                               RACE
    ALCOHOL_FREQUENCY | African..  Asian o..  Caucasian   Hispanic      Other    Unknown |     Total
----------------------+------------------------------------------------------------------+----------
                Daily |         2         10        881         10          7          3 |       913 
                Never |        32        150      1,671         45         66         14 |     1,978 
Occasionally (1-2 t.. |         6         53      1,501         27         40          8 |     1,635 
Rarely (a few times.. |        25        104      1,755         49         63         16 |     2,012 
Regularly (3-5 time.. |         3         23      1,565         22         19          5 |     1,637 
              Unknown |         2          3         66          3          4        118 |       196 
----------------------+------------------------------------------------------------------+----------
                Total |        70        343      7,439        156        199        164 |     8,371 
*/
tab allergic_to_i_have_no_food_aller 
/*
ALLERGIC_TO |
_I_HAVE_NO_ |
FOOD_ALLERG |
IES_THAT_I_ |
    KNOW_OF |      Freq.     Percent        Cum.
------------+-----------------------------------
      false |      5,468       65.32       65.32
       true |      2,903       34.68      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab allergic_to_other 
/*
ALLERGIC_TO |
     _OTHER |      Freq.     Percent        Cum.
------------+-----------------------------------
      false |      7,802       93.20       93.20
       true |        569        6.80      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
	  tab allergic_to_peanuts 
/*
ALLERGIC_TO |
   _PEANUTS |      Freq.     Percent        Cum.
------------+-----------------------------------
      false |      8,267       98.76       98.76
       true |        104        1.24      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab allergic_to_shellfish 
/*
ALLERGIC_TO |
 _SHELLFISH |      Freq.     Percent        Cum.
------------+-----------------------------------
      false |      8,259       98.66       98.66
       true |        112        1.34      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab allergic_to_tree_nuts 
/*
ALLERGIC_TO |
 _TREE_NUTS |      Freq.     Percent        Cum.
------------+-----------------------------------
      false |      8,274       98.84       98.84
       true |         97        1.16      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab allergic_to_unspecified 
/*
ALLERGIC_TO |
_UNSPECIFIE |
          D |      Freq.     Percent        Cum.
------------+-----------------------------------
      false |      3,590       42.89       42.89
       true |      4,781       57.11      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab alzheimers 
/*
                             ALZHEIMERS |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |          4        0.05        0.05
Diagnosed by an alternative medicine .. |          1        0.01        0.06
           I do not have this condition |      4,135       49.40       49.46
                         Self-diagnosed |          1        0.01       49.47
                                Unknown |      4,230       50.53      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab antibiotic_history 
/*
                     ANTIBIOTIC_HISTORY |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
                               6 months |      1,051       12.56       12.56
I have not taken antibiotics in the p.. |      5,485       65.52       78.08
                                  Month |        228        2.72       80.80
                                Unknown |        238        2.84       83.65
                                   Week |        161        1.92       85.57
                                   Year |      1,208       14.43      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab antibiotic_history race
/*
                      |                               RACE
   ANTIBIOTIC_HISTORY | African..  Asian o..  Caucasian   Hispanic      Other    Unknown |     Total
----------------------+------------------------------------------------------------------+----------
             6 months |        11         40        954         17         26          3 |     1,051 
I have not taken an.. |        50        212      4,962        103        132         26 |     5,485 
                Month |         0         16        202          5          4          1 |       228 
              Unknown |         1         11         86          2          9        129 |       238 
                 Week |         0          8        146          3          3          1 |       161 
                 Year |         8         56      1,089         26         25          4 |     1,208 
----------------------+------------------------------------------------------------------+----------
                Total |        70        343      7,439        156        199        164 |     8,371 
*/
tab appendix_removed 
/*
APPENDIX_RE |
      MOVED |      Freq.     Percent        Cum.
------------+-----------------------------------
         No |      6,910       82.55       82.55
   Not sure |         22        0.26       82.81
    Unknown |        289        3.45       86.26
        Yes |        774        9.25       95.51
      false |        315        3.76       99.27
       true |         61        0.73      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
 tab artificial_sweeteners 
/*
        ARTIFICIAL_SWEETENERS |      Freq.     Percent        Cum.
------------------------------+-----------------------------------
                        Daily |         27        0.31        0.31
                        Never |        414        4.74        5.05
Occasionally (1-2 times/week) |         32        0.37        5.41
   Rarely (a few times/month) |         85        0.97        6.39
   Regularly (3-5 times/week) |         28        0.32        6.71
                      Unknown |        867        9.93       16.63
                      no_data |      7,282       83.37      100.00
------------------------------+-----------------------------------
                        Total |      8,735      100.00
*/
tab asd
/*
                                    ASD |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |        118        1.41        1.41
Diagnosed by an alternative medicine .. |          6        0.07        1.48
           I do not have this condition |      3,949       47.17       48.66
                         Self-diagnosed |         34        0.41       49.06
                                Unknown |      4,264       50.94      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab asd race
/*
                      |                               RACE
                  ASD | African..  Asian o..  Caucasian   Hispanic      Other    Unknown |     Total
----------------------+------------------------------------------------------------------+----------
Diagnosed by a medi.. |         0          8         82         13         14          1 |       118 
Diagnosed by an alt.. |         0          0          5          0          1          0 |         6 
I do not have this .. |        30        180      3,546         70         98         25 |     3,949 
       Self-diagnosed |         0          1         31          0          2          0 |        34 
              Unknown |        40        154      3,775         73         84        138 |     4,264 
----------------------+------------------------------------------------------------------+----------
                Total |        70        343      7,439        156        199        164 |     8,371 
*/
tab autoimmune 
/*
                             AUTOIMMUNE |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |        465        5.55        5.55
Diagnosed by an alternative medicine .. |         26        0.31        5.87
           I do not have this condition |      3,562       42.55       48.42
                         Self-diagnosed |         30        0.36       48.78
                                Unknown |      4,288       51.22      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
 tab autoimmune race
/*
                      |                               RACE
           AUTOIMMUNE | African..  Asian o..  Caucasian   Hispanic      Other    Unknown |     Total
----------------------+------------------------------------------------------------------+----------
Diagnosed by a medi.. |         1          7        438          5         10          4 |       465 
Diagnosed by an alt.. |         0          2         20          1          2          1 |        26 
I do not have this .. |        29        178      3,161         73        102         19 |     3,562 
       Self-diagnosed |         0          1         27          1          1          0 |        30 
              Unknown |        40        155      3,793         76         84        140 |     4,288 
----------------------+------------------------------------------------------------------+----------
                Total |        70        343      7,439        156        199        164 |     8,371 
*/
tab bowel_movement_frequency 
/*
BOWEL_MOVEMEN |
  T_FREQUENCY |      Freq.     Percent        Cum.
--------------+-----------------------------------
 Five or more |         66        0.79        0.79
         Four |         99        1.18        1.97
Less than one |        450        5.38        7.35
          One |      2,006       23.96       31.31
        Three |        328        3.92       35.23
          Two |      1,092       13.05       48.27
      Unknown |      4,330       51.73      100.00
--------------+-----------------------------------
        Total |      8,371      100.00
*/
tab bowel_movement_quality 
/*
                 BOWEL_MOVEMENT_QUALITY |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
I don't know, I do not have a point o.. |        136        1.62        1.62
I tend to be constipated (have diffic.. |        424        5.07        6.69
I tend to be constipated (have diffic.. |        151        1.80        8.49
 I tend to have diarrhea (watery stool) |        410        4.90       13.39
I tend to have diarrhea (watery stool.. |        172        2.05       15.45
     I tend to have normal formed stool |      1,878       22.43       37.88
I tend to have normal formed stool - .. |        907       10.84       48.72
                                Unknown |      4,293       51.28      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab breastmilk_formula_ensure 
/*
              BREASTMILK_FORMULA_ENSURE |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
I eat both solid food and formula/bre.. |         40        0.48        0.48
                                     No |      3,687       44.04       44.52
                                Unknown |      4,434       52.97       97.49
                                    Yes |         51        0.61       98.10
                                  false |        158        1.89       99.99
                                   true |          1        0.01      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab cancer 
/*
                                 CANCER |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |         41        0.49        0.49
Diagnosed by an alternative medicine .. |          1        0.01        0.50
           I do not have this condition |        492        5.88        6.38
                                Unknown |        837       10.00       16.38
                                no_data |      7,000       83.62      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab cancer race
/*
                      |                               RACE
               CANCER | African..  Asian o..  Caucasian   Hispanic      Other    Unknown |     Total
----------------------+------------------------------------------------------------------+----------
Diagnosed by a medi.. |         0          1         40          0          0          0 |        41 
Diagnosed by an alt.. |         0          0          1          0          0          0 |         1 
I do not have this .. |         4         25        423          8         30          2 |       492 
              Unknown |         9         71        675         20         18         44 |       837 
              no_data |        57        246      6,300        128        151        118 |     7,000 
----------------------+------------------------------------------------------------------+----------
                Total |        70        343      7,439        156        199        164 |     8,371 
*/
tab cancer_treatment 
/*
 CANCER_TREATMENT |      Freq.     Percent        Cum.
------------------+-----------------------------------
     Chemotherapy |         10        0.11        0.11
     No treatment |          2        0.02        0.14
Radiation therapy |          8        0.09        0.23
     Surgery only |         22        0.25        0.48
          Unknown |      1,411       16.15       16.63
          no_data |      7,282       83.37      100.00
------------------+-----------------------------------
            Total |      8,735      100.00
*/
tab cardiovascular_disease 
/*
                 CARDIOVASCULAR_DISEASE |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |        149        1.78        1.78
Diagnosed by an alternative medicine .. |          1        0.01        1.79
           I do not have this condition |      3,973       47.46       49.25
                         Self-diagnosed |          8        0.10       49.35
                                Unknown |      4,240       50.65      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab cat
/*
        CAT |      Freq.     Percent        Cum.
------------+-----------------------------------
    Unknown |        274        3.27        3.27
      false |      5,836       69.72       72.99
       true |      2,261       27.01      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab cdiff 
/*
                                  CDIFF |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |         80        0.96        0.96
Diagnosed by an alternative medicine .. |          9        0.11        1.06
           I do not have this condition |      3,953       47.22       48.29
                         Self-diagnosed |         10        0.12       48.41
                                Unknown |      4,319       51.59      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/
tab census_region 
/*
CENSUS_REGI |
         ON |      Freq.     Percent        Cum.
------------+-----------------------------------
    Midwest |        846       10.11       10.11
  Northeast |      1,437       17.17       27.27
      South |      1,483       17.72       44.99
    Unknown |      1,955       23.35       68.34
       West |      2,650       31.66      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab chickenpox 
/*
 CHICKENPOX |      Freq.     Percent        Cum.
------------+-----------------------------------
         No |      1,272       15.20       15.20
   Not sure |        262        3.13       18.33
    Unknown |        278        3.32       21.65
        Yes |      6,559       78.35      100.00
------------+-----------------------------------
      Total |      8,371      100.00
*/
tab chickenpox race
/*
           |                               RACE
CHICKENPOX | African..  Asian o..  Caucasian   Hispanic      Other    Unknown |     Total
-----------+------------------------------------------------------------------+----------
        No |        27        124      1,017         52         46          6 |     1,272 
  Not sure |         2         14        229          6          9          2 |       262 
   Unknown |         1          7        129          5          8        128 |       278 
       Yes |        40        198      6,064         93        136         28 |     6,559 
-----------+------------------------------------------------------------------+----------
     Total |        70        343      7,439        156        199        164 |     8,371 
*/
tab clinical_condition 
/*
                     CLINICAL_CONDITION |      Freq.     Percent        Cum.
----------------------------------------+-----------------------------------
Diagnosed by a medical professional (.. |        973       11.62       11.62
Diagnosed by an alternative medicine .. |         37        0.44       12.07
           I do not have this condition |      2,899       34.63       46.70
                         Self-diagnosed |         44        0.53       47.22
                                Unknown |      4,418       52.78      100.00
----------------------------------------+-----------------------------------
                                  Total |      8,371      100.00
*/









*******************************************************************************
********************* REMOVE UNNECESSARY VARIABLES ****************************
*******************************************************************************
drop alcohol_types_beercider alcohol_types_red_wine alcohol_types_sour_beers alcohol_types_spiritshard_alcoho alcohol_types_unspecified alcohol_types_white_wine
drop altitude
drop artificial_sweeteners
drop assigned_from_geo
drop breastmilk_formula_ensure
drop barcodesequence 
drop cancer_treatment

drop vioscreen_*
*******************************************************************************
export delimited using "12_map_trim_basic.csv", delimiter(tab) replace
save "10_metadata_console.dta", replace
log close
*******************************************************************************

*******************************************************************************

