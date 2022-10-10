# Load RCT information
tcList <- read.csv("C:/git/Troy/TroyCohortDiagnostics/inst/settings/TargetComparatorList.csv") #Should be change path
trials <- "TRITON-TIMI_38"
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


# history feature
covariateSettingsHistory <- createCovariateSettings(useConditionOccurrenceLongTerm = TRUE,
                                                    useProcedureOccurrenceLongTerm = TRUE,
                                                    endDays = -7)

covariatesHistory <- getDbCovariateData(connectionDetails = connectionDetails,
                                        cdmDatabaseSchema = cdmDatabaseSchema,
                                        cohortDatabaseSchema = cohortDatabaseSchema,
                                        cohortTable = cohortTable,
                                        cohortId = c(targetNumber, comparatorNumber),
                                        covariateSettings = covariateSettingsHistory)
covariatesHistory$covariates <- covariatesHistory$covariates %>% left_join(covariatesHistory$covariateRef) 

# final diagnosis feature
covariateSettingsFinalDiagnosis <- createCovariateSettings(useConditionOccurrenceShortTerm = TRUE,
                                                           useProcedureOccurrenceShortTerm = TRUE,
                                                           useDeviceExposureShortTerm = TRUE,
                                                           useDrugGroupEraShortTerm = TRUE,
                                                           shortTermStartDays = -7)

covariatesFinalDiagnosis <- getDbCovariateData(connectionDetails = connectionDetails,
                                               cdmDatabaseSchema = cdmDatabaseSchema,
                                               cohortDatabaseSchema = cohortDatabaseSchema,
                                               cohortTable = cohortTable,
                                               cohortId = c(targetNumber, comparatorNumber),
                                               covariateSettings = covariateSettingsFinalDiagnosis)

covariatesFinalDiagnosis$covariates <- covariatesFinalDiagnosis$covariates %>% left_join(covariatesFinalDiagnosis$covariateRef)

# Result for TRITON-TIMI 38 trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))
resultTablePooled <- resultTablePooled[-1,]

# Final diagnosis of NSTEMI
nstemi <- c(4270024, 315296)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% nstemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Non-ST-elevation MI and unstable angina"

# Final diagnosis of STEMI
stemi <- c(312327,319039,434376,436706,438170,438438,438447,441579,444406,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4178129,4243372,4267568,4275436,4296653,4303359,4324413,35610091,35610093,35611570,35611571,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% stemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: ST-elevation MI"

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# age>=75yr
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(.data$covariateId == 1002, .data$covariateValue >= 75) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "age >= 75yr"

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# Median BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Hypertension
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# dyslipidemia
dyl <- c(432867,437521,437827,438720,440360,4029258,4029259,4029260,4029261,4029263,4029305,4029890,4029891,4029892,4030586,4030618,4030619,4031945,4031947,4079876,4079885,4079886,4079887,4096215,4104485,4120314,4134862,4142496,4143177,4144326,4144529,4159131,4220010,4223495,4270878,4291436,4292672,4294296,4294297,4295608,4295609,4298010,4298723,4298733,4299409,4300461,4301409,35608140,36674388,36676683,36715325,37016144,37016353,40482885,43530660,43531564,43531651,45757265,45757280,45757432,45757500,45770880)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dyl) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: dyslipidemia"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Medical history: MI
mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "MI"

# Medical history: CABG
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: CABG"

#creatinine clearance <= 50
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue <= 50) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "creatinine clearance <= 60"

# Index procedure: PCI
pci <- c(43531440,2617369,4225903,4264286,43527999,4265293,45770795,4214516,43527998,4020653,44511273,43527997,4171077,44511532,3190422,44784573,44789455,4175997,2313811,43527994,43531438,37111313,4283892,44512256,43527996,43527909,3180695,44511269,3173558,44511272,4216356,3181809,2313802,43533247,44511270,43533352,4328103,4181025,4006788,44511131,4337738,4329263,44511130,2313803,2313804,4139198,44511268,35607959,40756929,43531439,4264285,2313801,44511133,3170433,43527908,4178148,2617370,43527995,2001505,43533248,44511271,4238755,2000064,43533353,2313810,2001506,3186322)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% pci) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: PCI"

# Index procedure: CABG
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: CABG"

# Index procedure: Stent #DES 45772824 #BMS 45771110
stent <- c(45771110, 45772824)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% stent) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: stent"

# Index procedure: Stent #DES 45772824 #BMS 45771110
bms <- c(45771110)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% bms) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: BMS"

# Index procedure: Stent #DES 45772824 #BMS 45771110
des <- c(45772824)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% des) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: DES"

# Index drug: heparin 1367571 21600973
heparin <- c(1367571, 21600973)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% heparin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: heparin"

# Index drug: LMWH 35807344 1301025
lmwh <- c(35807344, 1301025)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% lmwh) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: LMWH"

# Index drug: Bivalirudin 19084670
bivalirudin <- c(19084670)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% bivalirudin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: bivalirudin"

# ACEi or ARB
aceiarb <- c(21601782)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aceiarb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ACEi or ARB"

# bb
bb <- c(21601664)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% bb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "betablocker"

# statin
statins <- c(21601853) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins"

# CCB
ccb <- c(21601744)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ccb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ccb"

#aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

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


# history feature
covariateSettingsHistory <- createCovariateSettings(useConditionOccurrenceLongTerm = TRUE,
                                                    useProcedureOccurrenceLongTerm = TRUE,
                                                    endDays = -7)

covariatesHistory <- getDbCovariateData(connectionDetails = connectionDetails,
                                        cdmDatabaseSchema = cdmDatabaseSchema,
                                        cohortDatabaseSchema = cohortDatabaseSchema,
                                        cohortTable = cohortTable,
                                        cohortId = c(targetNumber, comparatorNumber),
                                        covariateSettings = covariateSettingsHistory)
covariatesHistory$covariates <- covariatesHistory$covariates %>% left_join(covariatesHistory$covariateRef) 

# final diagnosis feature
covariateSettingsFinalDiagnosis <- createCovariateSettings(useConditionOccurrenceShortTerm = TRUE,
                                                           useProcedureOccurrenceShortTerm = TRUE,
                                                           useDeviceExposureShortTerm = TRUE,
                                                           useDrugGroupEraShortTerm = TRUE,
                                                           shortTermStartDays = -7)

covariatesFinalDiagnosis <- getDbCovariateData(connectionDetails = connectionDetails,
                                               cdmDatabaseSchema = cdmDatabaseSchema,
                                               cohortDatabaseSchema = cohortDatabaseSchema,
                                               cohortTable = cohortTable,
                                               cohortId = c(targetNumber, comparatorNumber),
                                               covariateSettings = covariateSettingsFinalDiagnosis)

covariatesFinalDiagnosis$covariates <- covariatesFinalDiagnosis$covariates %>% left_join(covariatesFinalDiagnosis$covariateRef)

# Result for TRITON-TIMI 38 trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))
resultTablePooled <- resultTablePooled[-1,]

# Final diagnosis of NSTEMI
nstemi <- c(4270024, 315296)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% nstemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Non-ST-elevation MI and unstable angina"

# Final diagnosis of STEMI
stemi <- c(312327,319039,434376,436706,438170,438438,438447,441579,444406,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4178129,4243372,4267568,4275436,4296653,4303359,4324413,35610091,35610093,35611570,35611571,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% stemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: ST-elevation MI"

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# age>=75yr
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(.data$covariateId == 1002, .data$covariateValue >= 75) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "age >= 75yr"

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# Median BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Hypertension
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# dyslipidemia
dyl <- c(432867,437521,437827,438720,440360,4029258,4029259,4029260,4029261,4029263,4029305,4029890,4029891,4029892,4030586,4030618,4030619,4031945,4031947,4079876,4079885,4079886,4079887,4096215,4104485,4120314,4134862,4142496,4143177,4144326,4144529,4159131,4220010,4223495,4270878,4291436,4292672,4294296,4294297,4295608,4295609,4298010,4298723,4298733,4299409,4300461,4301409,35608140,36674388,36676683,36715325,37016144,37016353,40482885,43530660,43531564,43531651,45757265,45757280,45757432,45757500,45770880)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% dyl) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: dyslipidemia"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Medical history: MI
mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "MI"

# Medical history: CABG
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: CABG"

#creatinine clearance <= 50
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue <= 50) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "creatinine clearance <= 60"

# Index procedure: PCI
pci <- c(43531440,2617369,4225903,4264286,43527999,4265293,45770795,4214516,43527998,4020653,44511273,43527997,4171077,44511532,3190422,44784573,44789455,4175997,2313811,43527994,43531438,37111313,4283892,44512256,43527996,43527909,3180695,44511269,3173558,44511272,4216356,3181809,2313802,43533247,44511270,43533352,4328103,4181025,4006788,44511131,4337738,4329263,44511130,2313803,2313804,4139198,44511268,35607959,40756929,43531439,4264285,2313801,44511133,3170433,43527908,4178148,2617370,43527995,2001505,43533248,44511271,4238755,2000064,43533353,2313810,2001506,3186322)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% pci) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: PCI"

# Index procedure: CABG
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: CABG"

# Index procedure: Stent #DES 45772824 #BMS 45771110
stent <- c(45771110, 45772824)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% stent) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: stent"

# Index procedure: Stent #DES 45772824 #BMS 45771110
bms <- c(45771110)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% bms) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: BMS"

# Index procedure: Stent #DES 45772824 #BMS 45771110
des <- c(45772824)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% des) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: DES"

# Index drug: heparin 1367571 21600973
heparin <- c(1367571, 21600973)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% heparin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: heparin"

# Index drug: LMWH 35807344 1301025
lmwh <- c(35807344, 1301025)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% lmwh) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: LMWH"

# Index drug: Bivalirudin 19084670
bivalirudin <- c(19084670)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% bivalirudin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Index procedure: bivalirudin"

# ACEi or ARB
aceiarb <- c(21601782)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aceiarb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ACEi or ARB"

# bb
bb <- c(21601664)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% bb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "betablocker"

# statin
statins <- c(21601853) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins"

# CCB
ccb <- c(21601744)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ccb) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ccb"

#aspirin
aspirin <- c(1112807)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"


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
