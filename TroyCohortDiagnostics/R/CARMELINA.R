# Load RCT information
tcList <- read.csv("C:/git/Troy/TroyCohortDiagnostics/inst/settings/TargetComparatorList.csv") #Should be change path
trials <- "CARMELINA"
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

# Result for CARMELINA trial
resultTablePooled <- data.frame()

# Mean age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Male
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8507001)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# Medical history: Heart failure
hf <- c(312927,314378,316139,316994,319835,439694,439696,439698,439846,442310,443580,443587,444031,444101,762002,762003,764871,764872,764873,764874,764876,764877,3184320,3656094,4004279,4009047,4014159,4023479,4030258,4071869,4079296,4079695,4103448,4108244,4108245,4111554,4124705,4138307,4139864,4141124,4142561,4172864,4177493,4185565,4193236,4195785,4195892,4199500,4205558,4206009,4215446,4215802,4229440,4233224,4233424,4242669,4259490,4264636,4267800,4273632,4284562,4307356,4311437,4327205,35615055,36712927,36712928,36712929,36713488,36716182,36716748,36717359,37110330,37309625,37311948,40479192,40479576,40480602,40480603,40481042,40481043,40482857,40486933,42598803,43020421,43020657,43021735,43021736,43021825,43021826,43021840,43021841,43021842,43022054,43022068,43530642,43530643,43530961,44782428,44782655,44782713,44782718,44782719,44782728,44782733,44784345,44784442,45766164,45766165,45766166,45766167,45766964,45773075)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Heart failure"

# Medical history: Ischemic heart disease
ihd <- c(4131824,4108677,40483752,4155008,37309626,3661520,4126801,4264145,4317287,4119949,434376,4069185,46273495,441579,4119944,4173632,42536628,40489421,44825431,4161456,315296,43022035,4161973,45766151,4215140,4119946,3661641,4185932,4219755,4100397,4186397,4324893,4124686,4270024,4119455,35615119,45766150,3661645,315830,37118984,4184827,761736,35610091,40487039,764123,37108686,4119943,4068741,37016755,4206867,765132,43021610,36712982,4119457,45766075,4119950,4030582,3661502,4303359,3183233,43021955,42536629,4121467,4124685,35615053,46270164,4201629,4329847,35622329,4145721,37017177,4108217,436706,3661642,37115756,4263712,4121464,35615052,4273462,4173171,3661524,4100132,3661503,764149,35610089,4310270,46270163,35610087,319039,43021066,37109910,4108218,4101319,4182190,45766117,4078531,3654466,45771322,4236169,3654467,37309713,761735,438170,4119945,4108670,45766114,4100871,4161974,45766116,4170094,315286,45766115,4119456,3661547,45766113,43021064,314666,46270160,4185302,37209632,4119942,35611571,4121466,4189939,36712984,43020460,35611570,4068938,40487573,4121465,4155009,46270161,4200113,4124683,45766076,4207921,4124684,4243372,4209541,4011131,37110242,35610093,3189643,439693,40481919,315832,4155963,4198141,438438,3661646,45773170,4119947,4119948,4178129,4108669,4108722,4138833,44782712,444406,43021065,321318,4051874,44783791,438168,316427,4102852,36712983,4267568,36714444,46270159,45766241,319844,46270162,4172865,4096252,4323202,4275436,312327,4116486,4296653,46270158,46274044,3661644,3661643,3654465,4324413,4151046,3661504,761737,4121468,4215259,3655133,4161457,44782769,4231426,4262446,43531588,438447)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ihd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Ischemic heart disease"

# hpt
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# afib
afib <- c(313217,4108832,4117112,4119601,4119602,4141360,4154290,4199501,4232691,4232697,37395821,44782442,45768480)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% afib) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Atrial fibrillation"

# EGFR
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3029859))[1,])

# eGFR >= 90
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 90) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 90"

# eGFR >= 60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 60"

# eGFR >= 45 to < 60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 45, covariateValue <60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 45 to < 60"

# eGFR >= 30 to < 45
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 30, covariateValue <45) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 30 to < 45"

# eGFR <60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue <60 ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR < 60"

#Median albumin-to-creatinine ratio
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3020682))[1,])

# UACR < 30
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue <30 ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "UACR < 30"

# UACR 30-300
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue >= 30, covariateValue <= 300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "UACR 30-300"

# UACR > 300
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue > 300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "UACR > 300"

# Mean BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# glycated hemoglobin
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004410)))

# Fasting plasma glucose
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027997)))

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

# Diastolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3012888)))

# Heart rate / min
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027018)))

# 	Cholesterol [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027114)))

# 	Cholesterol in LDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3028437)))

# 	Cholesterol in HDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3007070)))

#	Triglyceride [Mass/volume] in Serum or Plasma
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3022192)))

# Antidiabetics
diabetesdrugs<- c(21600712)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Antidiabetics"

# metformin 
metformin <- c(21600747, 1503297)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "metformin"

# SU 
SU <- c(1594973,1597756,1560171,19097821,1559684,1502809,1502855)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% SU) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "SU "

# insulin 
insulin  <- c(21600713,35605670,1531601,1567198,35602717,1516976,1502905,1544838,1588986,1550023,1513876,19078608,1590165,1596977,1586346,1513843,1513849,1562586,1586369)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% insulin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "insulin "

# Antihypertensive drugs
hptdrugs <- c(1319998,1332418,1314002,1335471,1322081,1338005,1351557,1340128,950370,1346823,19049145,1395058,19050216,19089969,1328165,1341927,1346686,19063575,1353776,1363749,974166,19122327,978555,1347384,1326012,1386957,19004539,19015802,1308216,1367500,19071995,19102106,907013,1307046,1310756,1314577,1318137,1318853,19113063,1319133,1319880,19020061,40226742,19024904,1327978,1373225,1345858,1353766,1331235,1334456,1317640,1342439,1308842,1307863,19010493,19102107)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hptdrugs) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Antihypertensive drugs"

# ACEi or ARB
aceiarb <- c(21601782)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aceiarb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ACEi or ARB"

# bb
bb <- c(21601664)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% bb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "betablocker"

# Diuretics
Diuretics <- c(21601461)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% Diuretics) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Diuretics"

# CCB
ccb <- c(21601744)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ccb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ccb"

# aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

# statin
statins <- c(21601853) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins"


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

# Result for CARMELINA trial
resultTablePooled <- data.frame()

# Mean age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Male
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8507001)))

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# Medical history: Heart failure
hf <- c(312927,314378,316139,316994,319835,439694,439696,439698,439846,442310,443580,443587,444031,444101,762002,762003,764871,764872,764873,764874,764876,764877,3184320,3656094,4004279,4009047,4014159,4023479,4030258,4071869,4079296,4079695,4103448,4108244,4108245,4111554,4124705,4138307,4139864,4141124,4142561,4172864,4177493,4185565,4193236,4195785,4195892,4199500,4205558,4206009,4215446,4215802,4229440,4233224,4233424,4242669,4259490,4264636,4267800,4273632,4284562,4307356,4311437,4327205,35615055,36712927,36712928,36712929,36713488,36716182,36716748,36717359,37110330,37309625,37311948,40479192,40479576,40480602,40480603,40481042,40481043,40482857,40486933,42598803,43020421,43020657,43021735,43021736,43021825,43021826,43021840,43021841,43021842,43022054,43022068,43530642,43530643,43530961,44782428,44782655,44782713,44782718,44782719,44782728,44782733,44784345,44784442,45766164,45766165,45766166,45766167,45766964,45773075)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Heart failure"

# Medical history: Ischemic heart disease
ihd <- c(4131824,4108677,40483752,4155008,37309626,3661520,4126801,4264145,4317287,4119949,434376,4069185,46273495,441579,4119944,4173632,42536628,40489421,44825431,4161456,315296,43022035,4161973,45766151,4215140,4119946,3661641,4185932,4219755,4100397,4186397,4324893,4124686,4270024,4119455,35615119,45766150,3661645,315830,37118984,4184827,761736,35610091,40487039,764123,37108686,4119943,4068741,37016755,4206867,765132,43021610,36712982,4119457,45766075,4119950,4030582,3661502,4303359,3183233,43021955,42536629,4121467,4124685,35615053,46270164,4201629,4329847,35622329,4145721,37017177,4108217,436706,3661642,37115756,4263712,4121464,35615052,4273462,4173171,3661524,4100132,3661503,764149,35610089,4310270,46270163,35610087,319039,43021066,37109910,4108218,4101319,4182190,45766117,4078531,3654466,45771322,4236169,3654467,37309713,761735,438170,4119945,4108670,45766114,4100871,4161974,45766116,4170094,315286,45766115,4119456,3661547,45766113,43021064,314666,46270160,4185302,37209632,4119942,35611571,4121466,4189939,36712984,43020460,35611570,4068938,40487573,4121465,4155009,46270161,4200113,4124683,45766076,4207921,4124684,4243372,4209541,4011131,37110242,35610093,3189643,439693,40481919,315832,4155963,4198141,438438,3661646,45773170,4119947,4119948,4178129,4108669,4108722,4138833,44782712,444406,43021065,321318,4051874,44783791,438168,316427,4102852,36712983,4267568,36714444,46270159,45766241,319844,46270162,4172865,4096252,4323202,4275436,312327,4116486,4296653,46270158,46274044,3661644,3661643,3654465,4324413,4151046,3661504,761737,4121468,4215259,3655133,4161457,44782769,4231426,4262446,43531588,438447)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ihd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Ischemic heart disease"

# hpt
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# afib
afib <- c(313217,4108832,4117112,4119601,4119602,4141360,4154290,4199501,4232691,4232697,37395821,44782442,45768480)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% afib) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Atrial fibrillation"

# EGFR
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3029859))[1,])

# eGFR >= 90
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 90) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 90"

# eGFR >= 60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 60"

# eGFR >= 45 to < 60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 45, covariateValue <60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 45 to < 60"

# eGFR >= 30 to < 45
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue >= 30, covariateValue <45) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR >= 30 to < 45"

# eGFR <60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue <60 ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "eGFR < 60"

#Median albumin-to-creatinine ratio
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3020682))[1,])

# UACR < 30
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue <30 ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "UACR < 30"

# UACR 30-300
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue >= 30, covariateValue <= 300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "UACR 30-300"

# UACR > 300
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue > 300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "UACR > 300"

# Mean BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# glycated hemoglobin
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004410)))

# Fasting plasma glucose
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027997)))

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

# Diastolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3012888)))

# Heart rate / min
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027018)))

# 	Cholesterol [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027114)))

# 	Cholesterol in LDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3028437)))

# 	Cholesterol in HDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3007070)))

#	Triglyceride [Mass/volume] in Serum or Plasma
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3022192)))

# Antidiabetics
diabetesdrugs<- c(21600712)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Antidiabetics"

# metformin 
metformin <- c(21600747, 1503297)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "metformin"

# SU 
SU <- c(1594973,1597756,1560171,19097821,1559684,1502809,1502855)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% SU) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "SU "

# insulin 
insulin  <- c(21600713,35605670,1531601,1567198,35602717,1516976,1502905,1544838,1588986,1550023,1513876,19078608,1590165,1596977,1586346,1513843,1513849,1562586,1586369)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% insulin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "insulin "

# Antihypertensive drugs
hptdrugs <- c(1319998,1332418,1314002,1335471,1322081,1338005,1351557,1340128,950370,1346823,19049145,1395058,19050216,19089969,1328165,1341927,1346686,19063575,1353776,1363749,974166,19122327,978555,1347384,1326012,1386957,19004539,19015802,1308216,1367500,19071995,19102106,907013,1307046,1310756,1314577,1318137,1318853,19113063,1319133,1319880,19020061,40226742,19024904,1327978,1373225,1345858,1353766,1331235,1334456,1317640,1342439,1308842,1307863,19010493,19102107)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hptdrugs) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Antihypertensive drugs"

# ACEi or ARB
aceiarb <- c(21601782)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aceiarb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ACEi or ARB"

# bb
bb <- c(21601664)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% bb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "betablocker"

# Diuretics
Diuretics <- c(21601461)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% Diuretics) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Diuretics"

# CCB
ccb <- c(21601744)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ccb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ccb"

# aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

# statin
statins <- c(21601853) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins"


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
filename <- paste0(trials, '.csv')
filename <- file.path("C:/output/TroyOhdsi", filename)
write.csv(RCTbaseline, filename)
