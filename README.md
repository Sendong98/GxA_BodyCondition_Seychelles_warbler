# GxA of Body Condition in the Seychelles warbler
R code and dataset for all analyses related to: No evidence of age-dependent changes in the additive genetic variance of body condition in a wild avian population.

The 'dataset' and 'figures' folders contain all the raw data and figures used in the paper.

## Abstract
Age-related declines in performance, commonly referred to as senescence, are theorised to arise due to reduced selection in later life. This relaxation of selection can allow the accumulation of deleterious mutations (mutation accumulation) and the persistence of genes with age-specific antagonistic effects (antagonistic pleiotropy). Both processes are predicted to generate genotype-by-age interactions for fitness-related traits, highlighting the need to investigate age-dependent additive genetic variation. Here, we use longitudinal data collected over 27 years from a pedigreed wild population of Seychelles warblers (Acrocephalus sechellensis) to explore sources of variation in, and the quantitative genetic basis of, senescence in a key fitness-related morphological trait: body condition. Our analyses, based on 1,443 individuals, reveal a clear pattern of senescence in body condition, and that this trait has a heritability estimated at 0.11 (95% CrI: 0.08–0.15). Furthermore, we find evidence for individual-by-age interactions, but not genotype-by-age interactions, in body condition trajectories. These results offer new insights into the genetic architecture of senescence in natural populations and suggest that age-related changes in morphological traits may be influenced and potentially limited by evolutionary processes acting on underlying genetic variance and the genetic architecture of the traits.


## File information
`1.data_extraction_GxA_bodycondition.R` is the R script that extracts data from the Seychelles warbler database and formats the dataset for analysis.

`2.dataset_distribution_description.R` contains the code for a basic distribution description of age, body mass, tarsus length, sex differences, etc.

`3.senescence_pattern_sw.R` describes how we estimate the senescence pattern of body condition in the Seychelles warbler (Describes the structure and dataset of Model 1 to Model 6).

`4.heritability_estimates_of_body_condition.R` describes the structure of pedigree- and GRM-based animal models and also estimates for different age classes (Describes the structure and dataset of Model 7 to Model 10).

`5.IxA_models.R` describes the random regression models for testing individual-by-age interactions with polynomial terms of age in the random part (Describes the structure and dataset of Model 11S1 to Model 11S3).

`6.GxA_models.R` describes the random regression animal models for testing genotype-by-age interactions with polynomial terms of age and permanent environment in the random part (Describes the structure and dataset of Model 12S1 to Model 12S11).

`7.character-state_model.R` contains the character state model that tests the genetic correlation between different age classes (Describes the structure and dataset of Model 13).

`8.reacnrom_model.R` In this script, the reacnorm package was used to separate the additive genetic variance into Vg and Vgxa (Model 14).

## Having issues
If you have any trouble, please file an issue in the GitHub repository.

## License
MIT

<img src="https://github.com/Sendong98/GxA_BodyCondition_Seychelles_warbler/blob/main/figures/IMG_1053_seychelles_warbler.jpg" width="300" height="200"> | <img src="[https://github.com/Sendong98/GxA_BodyCondition_Seychelles_warbler/blob/main/figures/IMG_1053_seychelles_warbler.jpg](https://lh7-us.googleusercontent.com/sitesv-images-rt/AMxu72v2xvltOOCe-9rR7Nutgc8i-NQTYbnqXw_8MZp8O5II_oznRYXOgc5zxSosn4gquMKvezNLNxaoaS6CPCTEkAUKflpXVinfnogfIP-B5Ny6n9ktv01k0J6SCCADo36IQLZcZX-Q3VBm1Gixe20cVFGkFiPXmYnkqiKjOWxxSx0Dg4M2ZlbZy3qj4EgsNofX_d2apizEGo9ESxgwQz0XjjJx8CM66dbsALdHYiXEk9I=w1280)" width="300" height="200">

<img src="https://github.com/Sendong98/GxA_BodyCondition_Seychelles_warbler/blob/main/figures/IMG_1053_seychelles_warbler.jpg" width="300" height="200">

