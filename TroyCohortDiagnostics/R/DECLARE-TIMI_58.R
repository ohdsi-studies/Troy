# Load RCT information
tcList <- read.csv("C:/git/Troy/TroyCohortDiagnostics/inst/settings/TargetComparatorList.csv") #Should be change path
trials <- "DECLARE-TIMI_58"
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
                                             useProcedureOccurrenceLongTerm = TRUE)

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

# Result for DECLARE TIMI trial
resultTablePooled <- data.frame()

# Mean age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# BMI 
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# glycated hemoglobin
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004410)))

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

#EGFR
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3029859))[1,])

# Medical history: History of atherosclerotic vascular disease: any
acsangina <- c(312327,315296,315830,315832,319039,321318,434376,436706,438170,438438,438447,441579,444406,761735,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4068938,4078531,4108669,4108670,4116486,4119455,4119456,4119457,4119942,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4155008,4155009,4155963,4161456,4161457,4161973,4161974,4178129,4184827,4198141,4201629,4231426,4243372,4262446,4264145,4267568,4270024,4275436,4296653,4303359,4310270,4324413,4324893,35610091,35610093,35611570,35611571,35615052,35615053,36712982,36712983,36712984,37209632,37309713,43020460,43531588,44782712,44782769,45766075,45766076,45766115,45766116,45766150,45766151,45771322,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
ich <- c(372435,377254,379778,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4077086,4108356,4110189,4110190,4110192,4111714,4119140,4129534,4131383,4138327,4141405,4142739,4145897,4146185,4298750,4319146,35610084,35610085,36717605,37110678,37110679,37116473,37395562,40479572,42535465,42535466,42535523,42535524,43530683,43531607,44782773,45767658,45772786,46270031,46270380,46270381,46273649)
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
anyathero <- c(acsangina, ich, pvd)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% anyathero) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: any"

# Medical history: ACS + angina
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% acsangina) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: coronary "

# Medical history: peripheral arterial disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: peripheral"

# Medical history: Nonhemorrhagic stroke 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ich) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: cerebrovascular"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442, 316139, 443587, 45766164, 443580)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"


# insulin 
insulin  <- c(21600713,35605670,1531601,1567198,35602717,1516976,1502905,1544838,1588986,1550023,1513876,19078608,1590165,1596977,1586346,1513843,1513849,1562586,1586369)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% insulin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "insulin "

# metformin 
metformin <- c(21600747, 1503297)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "metformin"

# SU 
SU <- c(1594973,1597756,1560171,19097821,1559684,1502809,1502855)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% SU) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "SU "

# dpp4
dpp4  <- c(43013884,40239216,40166035,1580747,19122137)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dpp4) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "dpp4"

# glp1 
glp1 <- c(44816332,45774435,1583722,40170911,44506754,793143, 1123618)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% glp1) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "glp1"

#antiplatelets
antiplatelets <- c(1322184,40163718,1302398,40252640,40241186)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% antiplatelets) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "antiplatelets"

# ACEi or ARB
aceiarb <- c(21601782)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aceiarb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ACEi or ARB"

# bb
bb <- c(21601664)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% bb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "betablocker"


# statin and ezetimibe
statins <- c(21601853, 1526475) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins and ezetimibe"

# Diuretics
Diuretics <- c(21601461)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% Diuretics) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Diuretics"



##
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
                                             useProcedureOccurrenceLongTerm = TRUE)

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

# Result for DECLARE TIMI trial
resultTablePooled <- data.frame()

# Mean age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# BMI 
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# glycated hemoglobin
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004410)))

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

#EGFR
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3029859))[1,])

# Medical history: History of atherosclerotic vascular disease: any
acsangina <- c(312327,315296,315830,315832,319039,321318,434376,436706,438170,438438,438447,441579,444406,761735,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4068938,4078531,4108669,4108670,4116486,4119455,4119456,4119457,4119942,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4155008,4155009,4155963,4161456,4161457,4161973,4161974,4178129,4184827,4198141,4201629,4231426,4243372,4262446,4264145,4267568,4270024,4275436,4296653,4303359,4310270,4324413,4324893,35610091,35610093,35611570,35611571,35615052,35615053,36712982,36712983,36712984,37209632,37309713,43020460,43531588,44782712,44782769,45766075,45766076,45766115,45766116,45766150,45766151,45771322,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
ich <- c(372435,377254,379778,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4077086,4108356,4110189,4110190,4110192,4111714,4119140,4129534,4131383,4138327,4141405,4142739,4145897,4146185,4298750,4319146,35610084,35610085,36717605,37110678,37110679,37116473,37395562,40479572,42535465,42535466,42535523,42535524,43530683,43531607,44782773,45767658,45772786,46270031,46270380,46270381,46273649)
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
anyathero <- c(acsangina, ich, pvd)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% anyathero) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: any"

# Medical history: ACS + angina
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% acsangina) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: coronary "

# Medical history: peripheral arterial disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: peripheral"

# Medical history: Nonhemorrhagic stroke 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ich) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: cerebrovascular"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442, 316139, 443587, 45766164, 443580)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"


# insulin 
insulin  <- c(21600713,35605670,1531601,1567198,35602717,1516976,1502905,1544838,1588986,1550023,1513876,19078608,1590165,1596977,1586346,1513843,1513849,1562586,1586369)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% insulin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "insulin "

# metformin 
metformin <- c(21600747, 1503297)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "metformin"

# SU 
SU <- c(1594973,1597756,1560171,19097821,1559684,1502809,1502855)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% SU) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "SU "

# dpp4
dpp4  <- c(43013884,40239216,40166035,1580747,19122137)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dpp4) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "dpp4"

# glp1 
glp1 <- c(44816332,45774435,1583722,40170911,44506754,793143, 1123618)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% glp1) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "glp1"

#antiplatelets
antiplatelets <- c(1322184,40163718,1302398,40252640,40241186)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% antiplatelets) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "antiplatelets"

# ACEi or ARB
aceiarb <- c(21601782)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aceiarb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ACEi or ARB"

# bb
bb <- c(21601664)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% bb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "betablocker"


# statin and ezetimibe
statins <- c(21601853, 1526475) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins and ezetimibe"

# Diuretics
Diuretics <- c(21601461)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% Diuretics) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Diuretics"


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

#percent
idx <- which(RCTbaseline$typeofstatistics=="percent")
RCTbaseline$OEvalue_RCT[idx] <- (resultTablePooledOE$n[idx] / sizeOE[1,2]) - (((RCTbaseline[idx,]$target * RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + ((RCTbaseline[idx,]$comparator * RCTbaseline$comparatorsize[1]) / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1]))) / 100)
RCTbaseline$OIvalue_RCT[idx] <- (resultTablePooledOI$n[idx] / sizeOI[1,2]) - (((RCTbaseline[idx,]$target * RCTbaseline$targetsize[1] / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1])) + ((RCTbaseline[idx,]$comparator * RCTbaseline$comparatorsize[1]) / (RCTbaseline$targetsize[1] + RCTbaseline$comparatorsize[1]))) / 100)

filename <- paste0(trials, '.csv')
filename <- file.path("C:/output/TroyOhdsi", filename)
write.csv(RCTbaseline, filename)
