
# covid-19-uk-seroprevalence-model

This document provides replication code and data for our UK model of covid-19 infections. The R script estimates the level of infections in the country in any given week, based on data from the ONS and the study "Age-specific mortality and immunity patterns of SARS-CoV-2" by Megan O’Driscoll et al (2020). 

To replicate our analysis, clone this repository to your local computer and run "DC_uk_seropositivity_model_rep.R".

This analysis is contingent on the assumptions we lay out, particularly on the accuracy of the relation between registered deaths and underlying associated infections, which may change over time. This approach is not perfect. If the IFR in reality is higher or more deaths are counted than previously, infections are overestimated -- if the IFR in reality is lower, infections are underestimated. And perhaps more importantly: if treatments have improved since the investigations of O’Driscoll et al, it would mean the current outbreak is even larger than the one we suggest. 

We hope future work can improve and expand upon the approach presented here, and thank the research community and data providers for providing the data and estimates we make use of. 
