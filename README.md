# GxA of Body Condition in the Seychelles warbler
R code and dataset for all analyses related to: No evidence of age-dependent changes in the additive genetic variance of body condition in a wild avian population.

The 'dataset' and 'figures' folders contain all the raw data and figures used in the paper.

## Abstract
Age-related declines in performance, commonly referred to as senescence, are theorised to arise due to reduced selection in later life. This relaxation of selection can allow the accumulation of deleterious mutations (mutation accumulation) and the persistence of genes with age-specific antagonistic effects (antagonistic pleiotropy). Both processes are predicted to generate genotype-by-age interactions for fitness-related traits, highlighting the need to investigate age-dependent additive genetic variation. Here, we use longitudinal data collected over 27 years from a pedigreed wild population of Seychelles warblers (Acrocephalus sechellensis) to explore sources of variation in, and the quantitative genetic basis of, senescence in a key fitness-related morphological trait: body condition. Our analyses, based on 1,443 individuals, reveal a clear pattern of senescence in body condition, and that this trait has a heritability estimated at 0.11 (95% CrI: 0.08–0.15). Furthermore, we find evidence for individual-by-age interactions, but not genotype-by-age interactions, in body condition trajectories. These results offer new insights into the genetic architecture of senescence in natural populations and suggest that age-related changes in morphological traits may be influenced and potentially limited by evolutionary processes acting on underlying genetic variance and the genetic architecture of the traits.


## File information
`1.data_extraction_GxA_bodycondition.R` is the R script that extracts data from the Seychelles warbler database and formats the dataset for analysis.

`2.dataset_distribution_description.R` contains the code for a basic distribution description of age, body mass, tarsus length, sex differences, etc.

`3.senescence_pattern_sw.R` describes how we estimate the senescence pattern of body condition in the Seychelles warbler (describes the structure and dataset for Models 1-6).

`4.heritability_estimates_of_body_condition_stepwise.R` describes the structure of heritability estimates - stepwise and using pedigree data. (Model 7S1-7S5)

`5.heritability_GRM.R` describes the structure of GRM-based animal models (Model 8).

`6.heritability_age_classes.R` describes the structure of pedigree- and GRM-based animal models for different age classes (Models 9-10) and '7.plot_heritability_age_class.R' for plotting the age classes h2

`8.IxA_models.R` describes random regression models for testing individual-by-age interactions, with polynomial age terms in the random part (describes the structure and dataset for Models 11S1 to 11S3).

`9.GxA_models.R` describes the random regression animal models for testing genotype-by-age interactions with polynomial terms of age and permanent environment in the random part (Describes the structure and dataset of Model 12S1 to Model 12S11). `10.IxA_GxA_checks_plots.R` contains the code to diagnose the models and plots. 

`11.character-state_model_plots.R` contains the character state model that tests the genetic correlation between different age classes and plots (Describes the structure and dataset of Model 13).

`12.submit_gxa_array.sh` describes how to submit the R scripts / heavy jobs to the HPC, including all the gxa models and character-state_model.

`13.Reactnorm_models.R` In this script, the reacnorm package was used to separate the additive genetic variance into Vg and Vgxa (Model 14).

## Having issues
If you have any trouble, please file an issue in the GitHub repository.

## License
MIT

<img src="https://github.com/Sendong98/GxA_BodyCondition_Seychelles_warbler/blob/main/figures/IMG_1053_seychelles_warbler.jpg" width="300" height="200">  <img src="https://github.com/Sendong98/GxA_BodyCondition_Seychelles_warbler/blob/main/figures/DRG_logo_small.2400x2400.jpeg" width="200" height="200">


