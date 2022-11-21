# Load RCT information
tcList <- read.csv("C:/git/Troy/TroyCohortDiagnostics/inst/settings/TargetComparatorList.csv") #Should be change path
trials <- "ENGAGE_AF-TIMI_48"
idx <- which(tcList[,"trial"] == trials)
RCTbaseline <- read.csv(file.path("C:/git/Troy/TroyCohortDiagnostics/inst/csv", paste0(trials[1], ".csv")))
RCTbaseline <- RCTbaseline %>% filter(.data$isNa != 'Y')
RCTbaseline$target <- as.numeric(RCTbaseline$target)
RCTbaseline$comparator <- as.numeric(RCTbaseline$comparator)
RCTbaseline$targetsd <- as.numeric(RCTbaseline$targetsd)
RCTbaseline$comparatorsd <- as.numeric(RCTbaseline$comparatorsd)
RCTbaseline$targetsize <- as.numeric(RCTbaseline$targetsize)
RCTbaseline$comparatorsize <- as.numeric(RCTbaseline$comparatorsize)

cohort <- Andromeda::andromeda()  
sql <- "select * from @cohortDatabaseSchema.@cohortTable"
sql <- SqlRender::render(sql, cohortDatabaseSchema=cohortDatabaseSchema, cohortTable=cohortTable)
DatabaseConnector::querySqlToAndromeda(connection, sql, andromeda = cohort, andromedaTableName = "cohort")


#########################
### Eligible criteria ###
#########################
tc <- tcList[idx[1],c(1,2)]
targetNumber <- tc[1,1]
comparatorNumber <- tc[1,2]
targetCohort <- cohort$cohort %>% filter(.data$COHORT_DEFINITION_ID == targetNumber) %>% mutate(treatment = 1)
comparatorCohort <- cohort$cohort %>% filter(COHORT_DEFINITION_ID == comparatorNumber) %>% mutate(treatment = 0)
population <- targetCohort %>% union(comparatorCohort)
population <- rename(population, rowId = SUBJECT_ID)

sizeOE <- data.frame(treatment = -1, n = population %>% summarise(n = n_distinct(.data$rowId)))
sizeOE <- rbind(sizeOE, data.frame(population %>% group_by(.data$treatment) %>% tally()))

covariateSettings <- createCovariateSettings(useDemographicsGender = TRUE,
                                             useDemographicsAge = TRUE,
                                             useDemographicsRace = TRUE,
                                             useMeasurementValueAnyTimePrior = TRUE,
                                             useDrugGroupEraLongTerm = TRUE,
                                             useConditionGroupEraLongTerm = TRUE,
                                             useProcedureOccurrenceLongTerm = TRUE,
                                             useChads2 = TRUE)

covariates <- getDbCovariateData(connectionDetails = connectionDetails,
                                 cdmDatabaseSchema = cdmDatabaseSchema,
                                 cohortDatabaseSchema = cohortDatabaseSchema,
                                 cohortTable = cohortTable,
                                 cohortId = c(targetNumber, comparatorNumber),
                                 covariateSettings = covariateSettings)

covariates$covariates <- covariates$covariates %>% left_join(population, by=("rowId"), copy = TRUE)
covariates$covariates <- covariates$covariates %>% left_join(covariates$covariateRef)

statisticsPooled <- covariates$covariates %>%
  group_by(.data$covariateId) %>%
  summarise(sum = sum(as.numeric(.data$covariateValue), na.rm = TRUE),
            mean = mean(as.numeric(.data$covariateValue), na.rm = TRUE),
            sumSqr = sum(as.numeric(.data$covariateValue)^2, na.rm = TRUE),
            median = median(as.numeric(.data$covariateValue), na.rm = TRUE),
            n = n(),
            min = min(as.numeric(.data$covariateValue), na.rm = TRUE),
            max = max(as.numeric(.data$covariateValue), na.rm = TRUE)) %>%
  mutate(sd = sqrt(abs(.data$sumSqr - .data$mean^2))) 

statisticsPooled <- statisticsPooled %>% left_join(covariates$covariateRef) 
statisticsPooled <- statisticsPooled %>% left_join(select(covariates$analysisRef, analysisId, isBinary), by=("analysisId"))

# Result for ENGAGE AF-TIMI 48 trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

covariateSettings <- createCovariateSettings(useConditionOccurrenceLongTerm = TRUE)

covariatesAfib <- getDbCovariateData(connectionDetails = connectionDetails,
                                     cdmDatabaseSchema = cdmDatabaseSchema,
                                     cohortDatabaseSchema = cohortDatabaseSchema,
                                     cohortTable = cohortTable,
                                     cohortId = c(targetNumber, comparatorNumber),
                                     covariateSettings = covariateSettings)

# paroxysmal 
paroxysmal <- c(4154290102)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesAfib$covariates %>% filter(covariateId %in% paroxysmal) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "paroxysmal"

# age>=75yr
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(.data$covariateId == 1002, .data$covariateValue >= 75) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "age >= 75yr"

# Prior stroke, TIA, or systemic embolism
sts <- c(198446,312337,314965,372435,372924,373503,374055,375557,376713,377254,379778,432656,432923,433545,433832,437060,439847,441874,442053,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3171253,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4043734,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4048785,4077086,4108356,4108360,4108375,4108376,4108377,4108380,4108381,4110189,4110190,4110192,4110194,4110195,4110331,4110333,4111714,4111848,4111850,4111852,4111853,4112020,4112161,4112162,4112165,4119140,4121618,4129534,4131383,4138327,4139517,4141405,4142739,4145897,4146185,4148906,4170620,4264734,4298750,4319146,4338523,35610084,35610085,36712812,36712813,36717605,37109920,37109921,37109924,37109925,37110678,37110679,37116421,37116473,37395562,40479572,40485430,40489292,42535465,42535466,42535523,42535524,42535829,42539410,42599802,42599894,43530683,43530727,43531607,44782773,44784293,45767658,45768887,45768888,45771016,45772786,46270031,46270380,46270381,46272244,46273649)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% sts) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Prior stroke, TIA, or systemic embolism"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442, 316139, 443587, 45766164, 443580)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Hypertension requiring treatment
hptdrugs <- c(1319998,1332418,1314002,1335471,1322081,1338005,1351557,1340128,950370,1346823,19049145,1395058,19050216,19089969,1328165,1341927,1346686,19063575,1353776,1363749,974166,19122327,978555,1347384,1326012,1386957,19004539,19015802,1308216,1367500,19071995,19102106,907013,1307046,1310756,1314577,1318137,1318853,19113063,1319133,1319880,19020061,40226742,19024904,1327978,1373225,1345858,1353766,1331235,1334456,1317640,1342439,1308842,1307863,19010493,19102107)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hptdrugs) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Hypertension requiring treatment"

# # CHADS2 score mean
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1903)))

# CHADS2 - <= 3
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue <= 3) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - <= 3"

# CHADS2 - 4-6
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue > 3) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 4-6"


#Dose reduction at randomization
quinidine <- c(1307863, 1360421, 21600249, 21600250)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- length(unique(c(data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue <= 50))$rowId, data.frame(covariates$covariates %>% filter(conceptId %like% 3025315, covariateValue <= 60))$rowId, data.frame(covariates$covariates %>% filter(conceptId %in% quinidine))$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Dose reduction at randomization"

#creatinine clearance <= 50
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue <= 50) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "creatinine clearance <= 50"

# Weight ≤60 kg
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3025315, covariateValue <= 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Weight ≤60 kg"

#Use of verapamil or quinidine

resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% quinidine) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Verapamil or quinidine"

#aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

# Thienopyridine
thienopyridine <- c(1322184,40163718,40241186,40252640,1302398, 1310149)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% thienopyridine)  %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "thienopyridine"

# amiodarone
amiodarone <- c(1309944)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% amiodarone) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "amiodarone"

# digoxin
digoxin <- c(1326303)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% digoxin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "digoxin"

resultTablePooledOE <- resultTablePooled

#########################
### Indication cohrot ###
#########################
tc <- tcList[idx[2],c(1,2)]
targetNumber <- tc[1,1]
comparatorNumber <- tc[1,2]
targetCohort <- cohort$cohort %>% filter(.data$COHORT_DEFINITION_ID == targetNumber) %>% mutate(treatment = 1)
comparatorCohort <- cohort$cohort %>% filter(COHORT_DEFINITION_ID == comparatorNumber) %>% mutate(treatment = 0)
population <- targetCohort %>% union(comparatorCohort)
population <- rename(population, rowId = SUBJECT_ID)

sizeOI <- data.frame(treatment = -1, n = population %>% summarise(n = n_distinct(.data$rowId)))
sizeOI <- rbind(sizeOI, data.frame(population %>% group_by(.data$treatment) %>% tally()))

covariateSettings <- createCovariateSettings(useDemographicsGender = TRUE,
                                             useDemographicsAge = TRUE,
                                             useDemographicsRace = TRUE,
                                             useMeasurementValueAnyTimePrior = TRUE,
                                             useDrugGroupEraLongTerm = TRUE,
                                             useConditionGroupEraLongTerm = TRUE,
                                             useProcedureOccurrenceLongTerm = TRUE,
                                             useChads2 = TRUE)

covariates <- getDbCovariateData(connectionDetails = connectionDetails,
                                 cdmDatabaseSchema = cdmDatabaseSchema,
                                 cohortDatabaseSchema = cohortDatabaseSchema,
                                 cohortTable = cohortTable,
                                 cohortId = c(targetNumber, comparatorNumber),
                                 covariateSettings = covariateSettings)

covariates$covariates <- covariates$covariates %>% left_join(population, by=("rowId"), copy = TRUE)
covariates$covariates <- covariates$covariates %>% left_join(covariates$covariateRef)

statisticsPooled <- covariates$covariates %>%
  group_by(.data$covariateId) %>%
  summarise(sum = sum(as.numeric(.data$covariateValue), na.rm = TRUE),
            mean = mean(as.numeric(.data$covariateValue), na.rm = TRUE),
            sumSqr = sum(as.numeric(.data$covariateValue)^2, na.rm = TRUE),
            median = median(as.numeric(.data$covariateValue), na.rm = TRUE),
            n = n(),
            min = min(as.numeric(.data$covariateValue), na.rm = TRUE),
            max = max(as.numeric(.data$covariateValue), na.rm = TRUE)) %>%
  mutate(sd = sqrt(abs(.data$sumSqr - .data$mean^2))) 

statisticsPooled <- statisticsPooled %>% left_join(covariates$covariateRef) 
statisticsPooled <- statisticsPooled %>% left_join(select(covariates$analysisRef, analysisId, isBinary), by=("analysisId"))

# Result for ENGAGE AF-TIMI 48 trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

covariateSettings <- createCovariateSettings(useConditionOccurrenceLongTerm = TRUE)

covariatesAfib <- getDbCovariateData(connectionDetails = connectionDetails,
                                     cdmDatabaseSchema = cdmDatabaseSchema,
                                     cohortDatabaseSchema = cohortDatabaseSchema,
                                     cohortTable = cohortTable,
                                     cohortId = c(targetNumber, comparatorNumber),
                                     covariateSettings = covariateSettings)

# paroxysmal 
paroxysmal <- c(4154290102)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesAfib$covariates %>% filter(covariateId %in% paroxysmal) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "paroxysmal"

# age>=75yr
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(.data$covariateId == 1002, .data$covariateValue >= 75) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "age >= 75yr"

# Prior stroke, TIA, or systemic embolism
sts <- c(198446,312337,314965,372435,372924,373503,374055,375557,376713,377254,379778,432656,432923,433545,433832,437060,439847,441874,442053,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3171253,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4043734,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4048785,4077086,4108356,4108360,4108375,4108376,4108377,4108380,4108381,4110189,4110190,4110192,4110194,4110195,4110331,4110333,4111714,4111848,4111850,4111852,4111853,4112020,4112161,4112162,4112165,4119140,4121618,4129534,4131383,4138327,4139517,4141405,4142739,4145897,4146185,4148906,4170620,4264734,4298750,4319146,4338523,35610084,35610085,36712812,36712813,36717605,37109920,37109921,37109924,37109925,37110678,37110679,37116421,37116473,37395562,40479572,40485430,40489292,42535465,42535466,42535523,42535524,42535829,42539410,42599802,42599894,43530683,43530727,43531607,44782773,44784293,45767658,45768887,45768888,45771016,45772786,46270031,46270380,46270381,46272244,46273649)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% sts) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Prior stroke, TIA, or systemic embolism"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442, 316139, 443587, 45766164, 443580)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Hypertension requiring treatment
hptdrugs <- c(1319998,1332418,1314002,1335471,1322081,1338005,1351557,1340128,950370,1346823,19049145,1395058,19050216,19089969,1328165,1341927,1346686,19063575,1353776,1363749,974166,19122327,978555,1347384,1326012,1386957,19004539,19015802,1308216,1367500,19071995,19102106,907013,1307046,1310756,1314577,1318137,1318853,19113063,1319133,1319880,19020061,40226742,19024904,1327978,1373225,1345858,1353766,1331235,1334456,1317640,1342439,1308842,1307863,19010493,19102107)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hptdrugs) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Hypertension requiring treatment"

# # CHADS2 score mean
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1903)))

# CHADS2 - <= 3
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue <= 3) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - <= 3"

# CHADS2 - 4-6
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue > 3) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 4-6"


#Dose reduction at randomization
quinidine <- c(1307863, 1360421, 21600249, 21600250)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- length(unique(c(data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue <= 50))$rowId, data.frame(covariates$covariates %>% filter(conceptId %like% 3025315, covariateValue <= 60))$rowId, data.frame(covariates$covariates %>% filter(conceptId %in% quinidine))$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Dose reduction at randomization"

#creatinine clearance <= 50
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue <= 50) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "creatinine clearance <= 50"

# Weight ≤60 kg
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3025315, covariateValue <= 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Weight ≤60 kg"

#Use of verapamil or quinidine

resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% quinidine) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Verapamil or quinidine"

#aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

# Thienopyridine
thienopyridine <- c(1322184,40163718,40241186,40252640,1302398, 1310149)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% thienopyridine)  %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "thienopyridine"

# amiodarone
amiodarone <- c(1309944)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% amiodarone) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "amiodarone"

# digoxin
digoxin <- c(1326303)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% digoxin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "digoxin"

resultTablePooledOI <- resultTablePooled

############## compute value
RCTbaseline$OEvalue_RCT <- 0
RCTbaseline$OIvalue_RCT <- 0

#mean
idx <- which(RCTbaseline$typeofstatistics=="mean")
RCTbaseline$OEvalue_RCT[idx] <- (resultTablePooledOE[idx,]$mean - (RCTbaseline[idx,]$target * (RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + RCTbaseline[idx,]$comparator * (RCTbaseline$comparatorsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])))) / ((RCTbaseline[idx,]$targetsd * (RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + RCTbaseline[idx,]$comparatorsd * (RCTbaseline$comparatorsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1]))) + resultTablePooledOE[idx,]$sd)/2
RCTbaseline$OIvalue_RCT[idx] <- (resultTablePooledOI[idx,]$mean - (RCTbaseline[idx,]$target * (RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + RCTbaseline[idx,]$comparator * (RCTbaseline$comparatorsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])))) / ((RCTbaseline[idx,]$targetsd * (RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + RCTbaseline[idx,]$comparatorsd * (RCTbaseline$comparatorsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1]))) + resultTablePooledOI[idx,]$sd)/2

#median

#category 
idx <- which(RCTbaseline$typeofstatistics=="category")
RCTbaseline$OEvalue_RCT[idx] <- (resultTablePooledOE$n[idx] / sizeOE[1,2]) - (RCTbaseline[idx,]$target + RCTbaseline[idx,]$comparator) / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])
RCTbaseline$OIvalue_RCT[idx] <- (resultTablePooledOI$n[idx] / sizeOI[1,2]) - (RCTbaseline[idx,]$target + RCTbaseline[idx,]$comparator) / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])

RCTbaseline$OEvalue <- NA
RCTbaseline$OIvalue <- NA
RCTbaseline$OEvalue[idx] <- resultTablePooledOE$n[idx] / sizeOE[1,2]
RCTbaseline$OIvalue[idx] <- resultTablePooledOI$n[idx] / sizeOI[1,2]

#percent
idx <- which(RCTbaseline$typeofstatistics=="percent")
RCTbaseline$OEvalue_RCT[idx] <- (resultTablePooledOE$n[idx] / sizeOE[1,2]) - (((RCTbaseline[idx,]$target * RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + ((RCTbaseline[idx,]$comparator * RCTbaseline$comparatorsize[1]) / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1]))) / 100)
RCTbaseline$OIvalue_RCT[idx] <- (resultTablePooledOI$n[idx] / sizeOI[1,2]) - (((RCTbaseline[idx,]$target * RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + ((RCTbaseline[idx,]$comparator * RCTbaseline$comparatorsize[1]) / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1]))) / 100)

filename <- paste0(trials, '.csv')
filename <- file.path("C:/output/TroyOhdsi", filename)
write.csv(RCTbaseline, filename)
