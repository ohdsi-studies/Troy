# Load RCT information
tcList <- read.csv("C:/git/Troy/TroyCohortDiagnostics/inst/settings/TargetComparatorList.csv") #Should be change path
trials <- "ROCKET-AF"
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

# Result for ROCKET-AF trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# Median  BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

# Diastolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3012888)))

covariateSettings <- createCovariateSettings(useConditionOccurrenceLongTerm = TRUE)

covariatesAfib <- getDbCovariateData(connectionDetails = connectionDetails,
                                     cdmDatabaseSchema = cdmDatabaseSchema,
                                     cohortDatabaseSchema = cohortDatabaseSchema,
                                     cohortTable = cohortTable,
                                     cohortId = c(targetNumber, comparatorNumber),
                                     covariateSettings = covariateSettings)

# atrial fibrillation 
afib <- c(313217102)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesAfib$covariates %>% filter(covariateId %in% afib) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "afib"

# paroxysmal 
paroxysmal <- c(4154290102)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesAfib$covariates %>% filter(covariateId %in% paroxysmal) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "paroxysmal"

#aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

#vitamin k antagonists
vka <- c(21600962, 1310149)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% vka) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "vitamin k antagonists"

# # CHADS2 score mean
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1903)))

# CHADS2 - 2
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 2) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 2"

# CHADS2 - 3
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 3) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 3"

# CHADS2 - 4
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 4) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 4"

# CHADS2 - 5
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 5) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 5"

# CHADS2 - 6
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 6) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 6"

# Prior stroke, TIA, or systemic embolism
sts <- c(198446,312337,314965,372435,372924,373503,374055,375557,376713,377254,379778,432656,432923,433545,433832,437060,439847,441874,442053,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3171253,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4043734,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4048785,4077086,4108356,4108360,4108375,4108376,4108377,4108380,4108381,4110189,4110190,4110192,4110194,4110195,4110331,4110333,4111714,4111848,4111850,4111852,4111853,4112020,4112161,4112162,4112165,4119140,4121618,4129534,4131383,4138327,4139517,4141405,4142739,4145897,4146185,4148906,4170620,4264734,4298750,4319146,4338523,35610084,35610085,36712812,36712813,36717605,37109920,37109921,37109924,37109925,37110678,37110679,37116421,37116473,37395562,40479572,40485430,40489292,42535465,42535466,42535523,42535524,42535829,42539410,42599802,42599894,43530683,43530727,43531607,44782773,44784293,45767658,45768887,45768888,45771016,45772786,46270031,46270380,46270381,46272244,46273649)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% sts) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Prior stroke, TIA, or systemic embolism"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442, 316139, 443587, 45766164, 443580)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# Hypertension
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Medical history: MI
mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "MI"

# Medical history: peripheral arterial disease
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: peripheral"

# Medical history: COPD
copd <- c(255573,257004,259043,261325,261895,440748,3185549,3654571,3654836,3654837,4046986,4050732,4050733,4050734,4050961,4056405,4083395,4110048,4110056,4110635,4112828,4112836,4115044,4136683,4145496,4148124,4166508,4166517,4177944,4193588,4196712,4209097,4246105,4281815,4286497,4315386,4321845,42574216,42598711,43530693,44791725,44807895,45769389,46269701,46274062)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% copd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: COPD"

#creatinine clearance
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027108))[1,])

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

# Result for ROCKET-AF trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# Median  BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

# Diastolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3012888)))

covariateSettings <- createCovariateSettings(useConditionOccurrenceLongTerm = TRUE)

covariatesAfib <- getDbCovariateData(connectionDetails = connectionDetails,
                                     cdmDatabaseSchema = cdmDatabaseSchema,
                                     cohortDatabaseSchema = cohortDatabaseSchema,
                                     cohortTable = cohortTable,
                                     cohortId = c(targetNumber, comparatorNumber),
                                     covariateSettings = covariateSettings)

# atrial fibrillation 
afib <- c(313217102)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesAfib$covariates %>% filter(covariateId %in% afib) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "afib"

# paroxysmal 
paroxysmal <- c(4154290102)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesAfib$covariates %>% filter(covariateId %in% paroxysmal) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "paroxysmal"

#aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

#vitamin k antagonists
vka <- c(21600962, 1310149)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% vka) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "vitamin k antagonists"

# # CHADS2 score mean
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1903)))

# CHADS2 - 2
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 2) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 2"

# CHADS2 - 3
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 3) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 3"

# CHADS2 - 4
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 4) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 4"

# CHADS2 - 5
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 5) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 5"

# CHADS2 - 6
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(covariateId == 1903, covariateValue == 6) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "CHADS2 - 6"

# Prior stroke, TIA, or systemic embolism
sts <- c(198446,312337,314965,372435,372924,373503,374055,375557,376713,377254,379778,432656,432923,433545,433832,437060,439847,441874,442053,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3171253,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4043734,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4048785,4077086,4108356,4108360,4108375,4108376,4108377,4108380,4108381,4110189,4110190,4110192,4110194,4110195,4110331,4110333,4111714,4111848,4111850,4111852,4111853,4112020,4112161,4112162,4112165,4119140,4121618,4129534,4131383,4138327,4139517,4141405,4142739,4145897,4146185,4148906,4170620,4264734,4298750,4319146,4338523,35610084,35610085,36712812,36712813,36717605,37109920,37109921,37109924,37109925,37110678,37110679,37116421,37116473,37395562,40479572,40485430,40489292,42535465,42535466,42535523,42535524,42535829,42539410,42599802,42599894,43530683,43530727,43531607,44782773,44784293,45767658,45768887,45768888,45771016,45772786,46270031,46270380,46270381,46272244,46273649)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% sts) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Prior stroke, TIA, or systemic embolism"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442, 316139, 443587, 45766164, 443580)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# Hypertension
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Medical history: MI
mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "MI"

# Medical history: peripheral arterial disease
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: peripheral"

# Medical history: COPD
copd <- c(255573,257004,259043,261325,261895,440748,3185549,3654571,3654836,3654837,4046986,4050732,4050733,4050734,4050961,4056405,4083395,4110048,4110056,4110635,4112828,4112836,4115044,4136683,4145496,4148124,4166508,4166517,4177944,4193588,4196712,4209097,4246105,4281815,4286497,4315386,4321845,42574216,42598711,43530693,44791725,44807895,45769389,46269701,46274062)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% copd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: COPD"

#creatinine clearance
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027108))[1,])

##
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
