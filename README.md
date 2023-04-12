# capstone_rmi
Master's of Science in Business Analytics (MSBA) Capstone work in conjunction with RMI

### Overview
A pilot project to develop an interactive tool mapping California’s most harmful methane flaring and the socially vulnerable communities least able to cope with that hazard.

###
Research Question 1: which flares are most impactful? 
RQ2: which block groups are being the most impacted? 
01_data_clean_load.ipynb
02.1_bg_flare_impact.ipynb
02.2_inter_adj_bg_df.ipynb
04_flare_vol_by_year.ipynb

I calculated the extent that block groups (Census spatial units) intersect with methane flaring in California at a variety of buffer radius. I then use those spatial calculations to visualize the social impact of methane flaring on vulnerable communities, using Tableau. 

###
RQ3: Is there evidence of disproportionate impact of flaring on minorities? 
03.1_bg_buffer_analysis.ipynb
03.2_bg_buffer_analysis_refs.ipynb
03.3_permutation_testing.R
03.4_permutation_testing_refs.R

I explored whether evidence exists of disporportionate impact on BIPOC communities from methane flaring relative to white communities. Put simply, are minorities more likely to live within a certain buffer zone of a flare than whites?

I did the spatial calculations in Python, then transferred to R to complete the permutation testing itself. 


