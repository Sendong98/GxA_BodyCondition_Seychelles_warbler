# GxA in Body Condition in the Seychelles warbler
R code and dataset for all analyses related to: No evidence of age-dependent changes in the additive genetic variance of body condition in a wild avian population.
The 'dataset' and 'figures' folders contain all the raw data and figures used in the paper.

## File information:
'1.data_extraction_GxA_bodycondition.R' is the R script that extracts data from the Seychelles warbler database and formats the dataset for analysis.
'2.dataset_distribution_description.R' contains the code for a basic distribution description of age, body mass, tarsus length, sex differences, etc.
'3.senescence_pattern_sw.R' describes how we estimate the senescence pattern of body condition in the Seychelles warbler (Describes the structure and dataset of Model 1 to Model 6).
'4.heritability_estimates_of_body_condition.R' describes the structure of pedigree- and GRM-based animal models and also estimates for different age classes (Describes the structure and dataset of Model 7 to Model 10).
'5.IxA_models.R' describes the random regression models for testing individual-by-age interactions with polynomial terms of age in the random part (Describes the structure and dataset of Model 11S1 to Model 11S3).
'6.GxA_models.R' describes the random regression animal models for testing genotype-by-age interactions with polynomial terms of age and permanent environment in the random part (Describes the structure and dataset of Model 12S1 to Model 12S11).
'7.character-state_model.R' contains the character state model that tests the genetic correlation between different age classes (Describes the structure and dataset of Model 13).
'8.reacnrom_model.R' In this script, the reacnorm package was used to separate the additive genetic variance into Vg and Vgxa (Model 14).
