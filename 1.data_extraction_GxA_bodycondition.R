# Load packages
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(lubridate)
library(RODBC)
library(MasterBayes)
library(pedantics)

#set up database paths
#here you need to load the warbler database

"%!in%"<-Negate("%in%")
DRIVERINFO <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)};"
MDBPATH <- "E:\\download_12092025\\SW_database_2024_summer\\SeychellesWarbler1.11.accdb"
PATH <- paste0(DRIVERINFO, "DBQ=", MDBPATH)
swdb<-odbcDriverConnect(PATH)

# Load data from the datasets
body_mass_all <- read_xlsx("catch_data_new.xlsx")
pedigree_all <- read_xlsx("sys_PedigreeCombined.xlsx")
field_period <- read_xlsx("sys_FieldPeriodBreedingSeason.xlsx")
catch_detail <- read_xlsx("catch_detail.xlsx")
sex <- read_xlsx("sys_SexEstimates.xlsx")
sex_birthdate <- read_xlsx("sex_birthdate.xlsx")

#extract from database
body_mass_data_status <- sqlFetch(swdb,"sys_CatchPlusStatus",stringsAsFactors = F)
birdid <- sqlFetch(swdb,"tblBirdID",stringsAsFactors = F)
bugs <- sqlFetch(swdb,"usys_qTerrQualityInsectsPerDM2",stringsAsFactors = F)
GroupSize <- sqlFetch(swdb,"usys_qBreedStatusByFieldPeriodPlusBreedGroup",stringsAsFactors = F)
BirdIDPeriodTerritory <- sqlFetch(swdb,"usys_qBreedStatusBirdIDPeriodTerritory",stringsAsFactors = F)
write.csv(BirdIDPeriodTerritory,file="BirdIDPeriodTerritory.csv")

#add periodid and sex
body_mass_data <- body_mass_all %>% drop_na(BodyMass) %>% filter(Island=="CN")
body_mass_data_period <- left_join(body_mass_data,field_period,by="FieldPeriodID")
body_mass_data_period_sex <- left_join(body_mass_data_period,sex[,1:2],by="BirdID")

#add observer
body_mass_data_period_sex_info <- left_join(body_mass_data_period_sex,catch_detail[,c("CatchID","Observer")],by="CatchID")

#add season
body_mass_data_period_sex_info_season <- body_mass_data_period_sex_info %>% mutate(summer=if_else(month(OccasionDate)>3 & month(OccasionDate)<10,1,0))
#add age, filter >= 3months (post-fledging)
body_mass_data_period_sex_info_season_age <- left_join(body_mass_data_period_sex_info_season,sex_birthdate[,c(1,3)],by="BirdID")
body_mass_data_period_sex_info_season_age <- body_mass_data_period_sex_info_season_age %>% mutate(
  age_year = time_length(interval(BirthDate, OccasionDate), "years"),
  age_month = time_length(interval(BirthDate, OccasionDate), "months"))
body_mass_data_period_sex_info_season_age_3mon <- body_mass_data_period_sex_info_season_age %>% drop_na(SexEstimate) %>% filter(age_month>=3)

#birth year
body_mass_data_period_sex_info_season_age_3mon$BirthYear <- year(body_mass_data_period_sex_info_season_age_3mon$BirthDate)

#minutes (consider quadratic effect)- catch time from 6am
catch_time <- hm(body_mass_data_period_sex_info_season_age_3mon$CatchTime)-hm("6:00")
body_mass_data_period_sex_info_season_age_3mon_minutes <- body_mass_data_period_sex_info_season_age_3mon %>% mutate(minutes=(hour(catch_time)*60+minute(catch_time)))

#Catch year
body_mass_data_period_sex_info_season_age_3mon_minutes$CatchYear <- year(body_mass_data_period_sex_info_season_age_3mon_minutes$OccasionDate)

#insect abundance
bugs<-bugs[bugs$Island=="CN",]
bugs$invertsum<-rowSums(bugs[,c(29:47)], na.rm=T)
total_invert<-bugs[,c("Island","Year","FieldPeriodID","Location","invertsum")]
total_invert<-total_invert%>%
  group_by(Year, FieldPeriodID)%>%
  mutate(avg_invert=mean(invertsum))
total_invert$occasionyear<-total_invert$Year
total_invert<-total_invert[,c("occasionyear", "FieldPeriodID","avg_invert")]
total_invert<-unique(total_invert)
write.csv(total_invert, "mean_insect_new.csv")

#Add insects
body_mass_data_period_sex_info_season_age_3mon_minutes_insects <- left_join(body_mass_data_period_sex_info_season_age_3mon_minutes,total_invert[,2:3],by="FieldPeriodID")

#Group Size
birdid_FPid_groupid <- GroupSize %>% dplyr::select("BirdID","FieldPeriodID","BreedGroupID")
GroupSize_groupid <- GroupSize %>% 
  filter(Status %in% c("BrM","BrF","BrU","H","AB","ABX","FL","OFL")) %>% 
  group_by(FieldPeriodID, BreedGroupID) %>% 
  summarise(group_size = n_distinct(BirdID), .groups = "drop")
birdid_FPid_groupid <- birdid_FPid_groupid %>% 
  left_join(GroupSize_groupid, 
            by = c("FieldPeriodID", "BreedGroupID"))

body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize <- left_join(body_mass_data_period_sex_info_season_age_3mon_minutes_insects,birdid_FPid_groupid,by=c("BirdID","FieldPeriodID"))

#convert type
body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$BodyMass <- as.numeric(body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$BodyMass)
body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$LeftTarsus <- as.numeric(body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$LeftTarsus)
body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$RightTarsus <- as.numeric(body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$RightTarsus)
body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$age_year <- as.numeric(body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$age_year)
body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$minutes <- as.numeric(body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$minutes)
body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$group_size <- as.numeric(body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$group_size)
body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$SexEstimate <- factor(body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize$SexEstimate)


#Left tarsus as Right if not measured
df <- body_mass_data_period_sex_info_season_age_3mon_minutes_insects_groupsize %>%
  mutate(RightTarsus = if_else(is.na(RightTarsus), LeftTarsus, RightTarsus))
nrow(df)

#drop na
df_clean <- df %>% drop_na(BodyMass,age_year,RightTarsus,SexEstimate,summer,minutes,
                           avg_invert,BirthYear,CatchYear,Observer,BirdID,group_size)
nrow(df_clean)
write.csv(df_clean,file = "dataset_mass_20260502.csv")

#Pedigree dataset for RRMs
pedigree_Gen <- pedigree_all %>% mutate(Sire=GeneticFather,Dam=GeneticMother)

# rename and reorder pedigree data
pedigree_reorder <- pedigree_Gen[,c("BirdID","Dam","Sire")]
pedigree_reorder <- as.data.frame(pedigree_reorder)
pedigree_reorder <- orderPed(pedigree_reorder)
pedigree_reorder <- pedigree_reorder %>% distinct(BirdID,.keep_all = TRUE)

# check missing individual
pedigree_check <- pedigree_reorder[,c("BirdID","Dam","Sire")]
pedigree_check[!pedigree_check$Dam%in%pedigree_check$BirdID&!is.na(pedigree_check$Dam),]
pedigree_check[!pedigree_check$Sire%in%pedigree_check$BirdID&!is.na(pedigree_check$Sire),]

# pedigree stats
pedigreeStats <- pedStatSummary(pedigreeStats(pedigree_reorder, retain = "informative", graphicalReport = "n"))
write_rds(pedigreeStats,"pedigreeStats_20260501.rds")
pedigreeStats <- read_rds("pedigreeStats.rds")

PrunedPed_bodymass <- prunePed(pedigree_reorder, keep = df_clean$BirdID, make.base=TRUE)
write.csv(PrunedPed_bodymass, file = "PrunedPed_bodymass_20260502.csv", row.names = FALSE)
PrunedPed_bodymass <- read.csv("PrunedPed_bodymass_20260502.csv")

Measured <- ifelse(pedigree_reorder$BirdID %in% df_clean$BirdID, 1,0) #Phenotypic data recorded?
PrunedPedStats_measured <- pedStatSummary(pedigreeStats(pedigree_reorder, dat = Measured, retain = "informative", graphicalReport = "n"))
PrunedPedStats_measured
write_rds(PrunedPedStats_measured,"PrunedPedStats_measured_20260503.rds")
PrunedPedStats_measured <- read_rds("PrunedPedStats_measured_20260503.rds")
pdf("PrunedPedGraph_20260502.pdf", width = 10, height = 10)

drawPedigree(
  pedigree_reorder,
  dat = Measured,
  plotfull = "y",
  dots = "y",
  sexColours = c('lightpink2','royalblue')
)

dev.off()

