
##### AG TESTING IF CLUSTER IS ENRICHED FOR ETHNIC ASSOCIATIONS #####

### CLUSTER FDR P VALUES FOR ETHNIC STRUCTURING AG ###
inCluster = [1.55E-09,
9.29E-09,
3.38E-06,
3.93E-06,
4.10E-06,
7.36E-06,
0.000845582,
0.000869401,
0.001378597,
0.001416892,
0.001707589,
0.002807231,
0.003363839,
0.00342074,
0.005172632,
0.006058413,
0.015516904,
0.245608353
0.603480075
0.749635875]

### NONCLUSTER FDR P VALUES FOR ETHNIC STRUCTURING AG ###
nonCluster = [2.03E-10,
4.10E-06,
2.11E-05,
3.54E-05,
5.90E-05,
0.00066741,
0.000962969,
0.005172632,
0.006058413,
0.006058413,
0.012581671,
0.012581671,
0.012581671,
0.012581671,
0.012581671,
0.015074454,
0.034737156,
0.044956657,
0.050019333,
0.063723536,
0.097242944,
0.097242944,
0.097242944,
0.097242944,
0.097242944,
0.111017855,
0.122383354,
0.129432358,
0.167013908,
0.170191384,
0.18492238,
0.201363621,
0.20231235,
0.2044017,
0.21068217,
0.219176901,
0.222553233,
0.24789991,
0.254527568,
0.259303602,
0.261284489,
0.262890396,
0.278476793,
0.286988055,
0.352165394,
0.35627965,
0.369445273,
0.386892761,
0.386892761,
0.436487185,
0.486235712,
0.537953432,
0.542354644,
0.61176106,
0.61176106,
0.61176106,
0.61176106,
0.61176106,
0.61176106,
0.626407804,
0.627370539,
0.640880172,
0.643620982,
0.646858921,
0.646858921,
0.646858921,
0.65798878,
0.679723351,
0.679723351,
0.689988229,
0.700822915,
0.700822915,
0.72183505,
0.749635875,
0.751398206,
0.788544353,
0.789275074,
0.840966001,
0.842907268,
0.876960728,
0.876960728,
0.898117819,
0.903890933,
0.903890933,
0.903890933,
0.903890933,
0.903890933,
0.903890933,
0.905944047,
0.932955479,
0.947308217,
0.947308217,
0.963303905,
0.985075304,
0.985075304,
0.985075304,
0.985574026,
0.985574026,
0.989195518,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952,
0.990414952]

mannStat, mannP = sp.stats.mannwhitneyu(inCluster, nonCluster, use_continuity=True, alternative='less')
mannStat
mannP
# 257.5
# 2.3131606387082258e-10


##### HMP TESTING IF CLUSTER IS ENRICHED FOR ETHNIC ASSOCIATIONS #####

### CLUSTER FDR P VALUES FOR ETHNIC STRUCTURING HMP ###
inCluster = [0.001929418,
0.004633401,
0.007179918,
0.023362802,
0.088356111,
0.187890492]

### NONCLUSTER FDR P VALUES FOR ETHNIC STRUCTURING HMP ###
nonCluster = [0.000128281,
0.001929418,
0.003760497,
0.006294281,
0.006294281,
0.00659445,
0.015294782,
0.022998127,
0.025581646,
0.036075007,
0.047229589,
0.064139032,
0.064139032,
0.088356111,
0.088356111,
0.11277003,
0.118344007,
0.15278534,
0.182943438,
0.195580563,
0.235921196,
0.235921196,
0.235921196,
0.235921196,
0.310717138,
0.310717138,
0.310717138,
0.310717138,
0.368904048,
0.368904048,
0.513146761,
0.513146761,
0.520192473,
0.520192473,
0.559648812,
0.579534218,
0.599644117,
0.698419824,
0.729606279,
0.729606279,
0.781824035,
0.903010109,
0.903010109,
0.911048862,
0.95607139,
0.974804846]

mannStat, mannP = sp.stats.mannwhitneyu(inCluster, nonCluster, use_continuity=True, alternative='less')
mannStat
mannP
# 51.5
# 0.0068532129682198949