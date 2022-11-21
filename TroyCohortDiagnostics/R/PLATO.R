# Load RCT information
tcList <- read.csv("C:/git/Troy/TroyCohortDiagnostics/inst/settings/TargetComparatorList.csv") #Should be change path
trials <- "PLATO"
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

covariateSettingsFinalDiagnosis <- createCovariateSettings(useConditionOccurrenceShortTerm = TRUE,
                                                           shortTermStartDays = -7)

covariatesFinalDiagnosis <- getDbCovariateData(connectionDetails = connectionDetails,
                                               cdmDatabaseSchema = cdmDatabaseSchema,
                                               cohortDatabaseSchema = cohortDatabaseSchema,
                                               cohortTable = cohortTable,
                                               cohortId = c(targetNumber, comparatorNumber),
                                               covariateSettings = covariateSettingsFinalDiagnosis)

covariatesFinalDiagnosis$covariates <- covariatesFinalDiagnosis$covariates %>% left_join(covariatesFinalDiagnosis$covariateRef) 

# Result for PLATO trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# age>=75yr
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(.data$covariateId == 1002, .data$covariateValue >= 75) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "age >= 75yr"

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# median bodyweight
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3025315)))

# bodyweight <60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3025315, covariateValue < 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "body weight < 60"

# median BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Hypertension
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# dyslipidemia
dyl <- c(432867,437521,437827,438720,440360,4029258,4029259,4029260,4029261,4029263,4029305,4029890,4029891,4029892,4030586,4030618,4030619,4031945,4031947,4079876,4079885,4079886,4079887,4096215,4104485,4120314,4134862,4142496,4143177,4144326,4144529,4159131,4220010,4223495,4270878,4291436,4292672,4294296,4294297,4295608,4295609,4298010,4298723,4298733,4299409,4300461,4301409,35608140,36674388,36676683,36715325,37016144,37016353,40482885,43530660,43531564,43531651,45757265,45757280,45757432,45757500,45770880)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dyl) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: dyslipidemia"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Medical history: MI

mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: MI"

# Medical history: PCI
pci <- c(43531440,2617369,4225903,4264286,43527999,4265293,45770795,4214516,43527998,4020653,44511273,43527997,4171077,44511532,3190422,44784573,44789455,4175997,2313811,43527994,43531438,37111313,4283892,44512256,43527996,43527909,3180695,44511269,3173558,44511272,4216356,3181809,2313802,43533247,44511270,43533352,4328103,4181025,4006788,44511131,4337738,4329263,44511130,2313803,2313804,4139198,44511268,35607959,40756929,43531439,4264285,2313801,44511133,3170433,43527908,4178148,2617370,43527995,2001505,43533248,44511271,4238755,2000064,43533353,2313810,2001506,3186322)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% pci) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: PCI"

# Medical history: CABG
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: CABG"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# Medical history: Nonhemorrhagic stroke 
ich <- c(372435,377254,379778,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4077086,4108356,4110189,4110190,4110192,4111714,4119140,4129534,4131383,4138327,4141405,4142739,4145897,4146185,4298750,4319146,35610084,35610085,36717605,37110678,37110679,37116473,37395562,40479572,42535465,42535466,42535523,42535524,43530683,43531607,44782773,45767658,45772786,46270031,46270380,46270381,46273649)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% ich) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Nonhemorrhagic stroke "

# Medical history: peripheral arterial disease
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: peripheral arterial disease"

# Medical history: chronic renal disease
ckd <- c(193782,198185,443597,443601,443611,443612,443614,443961,760850,762000,762001,762033,762034,762973,764011,765535,765536,3169303,3177485,3178768,3185897,4030520,4125970,4128067,4128200,4239804,4322556,36679002,36716184,36716455,36716947,36717349,36717534,37017104,37017813,37018761,37018886,37019193,42599648,42599649,42599650,42599651,42872405,43020437,43020455,43020456,43020457,43021835,43021836,43021852,43021853,43021854,43021864,43531559,43531562,43531566,43531577,43531578,43531653,44782429,44782689,44782690,44782691,44782692,44782703,44782717,44784621,44784637,44784638,44784639,44784640,44792226,44792227,44792228,44792229,44792230,44792231,44792232,44792249,44792250,44792251,44792252,44792253,44792254,44792255,45757137,45757139,45757356,45757392,45757393,45757444,45757445,45757446,45757447,45763854,45763855,45768813,45769901,45769902,45769903,45769904,45769906,45771064,45771067,45771075,45772751,45773576,45773688,46270353,46270354,46270355,46270356,46271022,46273164,46273514,46273636,46284566,46284567,46284570,46284572,46284575,46284587,46284588,46284591,46284592,46284593,46284597,46284598,46284599,46284600,46284602,46284603,46286992,46287169)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% ckd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: chronic renal disease"

# Medical history: history of dyspnea
dyspnea <- c(312437,315361,3655952,3656669,4010972,4047610,4059021,4059022,4060052,4060439,4060709,4065915,4066850,4091788,4094132,4097311,4144682,4148800,4178416,4190875,4192279,4193263,4206307,4212233,4217021,4219335,4219740,4220640,4220641,4221237,4222444,4222446,4222908,4223033,4223034,4223035,4223036,4223037,4223038,4228754,4244276,4248284,4253185,4263848,4307188,4308377,4310059,4310172,35610139,35610140,35610141,35610142,35610143,35610144,36685567,36685568,36685569,36685570,36685571,44809070)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dyspnea) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of dyspnea"

# Medical history: COPD
copd <- c(255573,257004,259043,261325,261895,440748,3185549,3654571,3654836,3654837,4046986,4050732,4050733,4050734,4050961,4056405,4083395,4110048,4110056,4110635,4112828,4112836,4115044,4136683,4145496,4148124,4166508,4166517,4177944,4193588,4196712,4209097,4246105,4281815,4286497,4315386,4321845,42574216,42598711,43530693,44791725,44807895,45769389,46269701,46274062)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% copd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: COPD"
# Medical history: asthma
asthma <- c(252658,252942,256448,257581,257583,312950,313236,317009,443801,761844,764677,764949,2108536,3661412,4022592,4051466,4057952,4075237,4080516,4081994,4110051,4119298,4119300,4120261,4120262,4123253,4125022,4138760,4141978,4142738,4143474,4143828,4145356,4145497,4146581,4152913,4155468,4155469,4155470,4191479,4206340,4211530,4212099,4217558,4225553,4225554,4232595,4233784,4235703,4244339,4245292,4245676,4250128,4271333,4277596,4279553,4301938,4309833,4312524,35609846,35609847,36684328,36684335,37108580,37108581,37109103,37116845,37206717,37208352,37310241,40481763,40483397,42535716,42536207,42536208,42536649,42538744,42539549,43530693,43530745,44783610,44810117,45757063,45766727,45766728,45768910,45768911,45768912,45768963,45768964,45768965,45769350,45769351,45769352,45769438,45769441,45769442,45769443,45771045,45772073,45772937,45773005,46269767,46269770,46269771,46269776,46269777,46269784,46269785,46269801,46269802,46270028,46270029,46270030,46270082,46270322,46273452,46273454,46273462,46273487,46273635,46274059,46274062,46274124)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% asthma) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: asthma"
# Medical history: gout
gout <- c(74892,376955,440674,762219,4035436,4035437,4035438,4035751,4035752,4039287,4079745,4083344,4084229,4099013,4101425,4104957,4104958,4104959,4114911,4114912,4116007,4116008,4116009,4153103,4154437,4179065,4182267,4188214,4247334,4281849,4285308,4297662,4299408,35615129,35615130,36674930,36684526,36684527,36684528,36684529,36684530,36684531,36684532,36684533,36684534,36684535,36714445,36715572,37108592,37108594,37209406,37209407,37209408,37209463,37394509,42534842,42534843,42535119,42535120,42535121,42535122,42535127,42535133,42535203,42535204,42535205,42535206,42539209,42539261,42539266,42539288,44783578,45770178,46270412,46270413,46270416,46270417,46270418,46270419,46270420,46270421,46270422,46270423,46270424,46270425,46270426,46270427,46270428,46270429,46270430,46270431,46270432,46270435,46270436,46270437,46270438,46270439,46270440,46270441,46270442,46270443,46270444,46270445,46270446,46270447,46270448,46270449,46270450,46270451,46270452,46270453,46270454,46270455,46270456,46270457,46270458,46270459,46270460,46270461,46270464,46270465,46270466,46270467,46270468,46270470,46270471,46270472,46270473,46270474,46270475,46272380,46272381,46273507,46273508,46273523,46273541,46273638,46273639,46273641,46273648,46273653,46273655)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% gout) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: gout"

# Final diagnosis of STEMI
stemi <- c(312327,319039,434376,436706,438170,438438,438447,441579,444406,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4178129,4243372,4267568,4275436,4296653,4303359,4324413,35610091,35610093,35611570,35611571,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% stemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: ST-elevation MI"

# Final diagnosis of NSTEMI
nstemi <- c(4270024)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% nstemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Non-ST-elevation MI"

# Final diagnosis of Unstable angina
unstableangina <- c(315296)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% unstableangina) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Unstable angina"

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

covariateSettingsFinalDiagnosis <- createCovariateSettings(useConditionOccurrenceShortTerm = TRUE,
                                                           shortTermStartDays = -7)

covariatesFinalDiagnosis <- getDbCovariateData(connectionDetails = connectionDetails,
                                               cdmDatabaseSchema = cdmDatabaseSchema,
                                               cohortDatabaseSchema = cohortDatabaseSchema,
                                               cohortTable = cohortTable,
                                               cohortId = c(targetNumber, comparatorNumber),
                                               covariateSettings = covariateSettingsFinalDiagnosis)

covariatesFinalDiagnosis$covariates <- covariatesFinalDiagnosis$covariates %>% left_join(covariatesFinalDiagnosis$covariateRef) 

# Result for PLATO trial
resultTablePooled <- data.frame()

# Median age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# age>=75yr
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(.data$covariateId == 1002, .data$covariateValue >= 75) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "age >= 75yr"

# Female
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 8532001)))

# median bodyweight
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3025315)))

# bodyweight <60
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3025315, covariateValue < 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "body weight < 60"

# median BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Hypertension
hpt <- c(132685,134414,135601,136743,136760,137940,141084,141639,192684,197930,200157,201912,312648,314090,314103,314423,314958,316866,317895,317898,318437,319826,320128,320456,321074,321080,321638,433536,438490,439077,439393,441922,443771,762994,3169253,3191244,3656115,4006325,4023318,4028741,4028951,4032952,4034031,4034094,4034095,4035655,4048212,4049389,4057976,4057978,4057979,4058987,4061667,4062550,4062811,4062906,4071202,4080325,4083723,4094374,4108213,4110947,4110948,4118910,4146627,4146816,4148205,4151903,4159755,4162306,4167358,4167493,4174979,4178312,4179379,4180283,4199306,4209293,4212496,4215640,4217486,4218088,4219323,4221991,4227607,4242878,4249016,4253928,4262182,4263067,4269358,4276511,4277110,4279525,4283352,4289142,4289933,4291933,4302591,4304837,4305599,4311246,4316372,4321603,4322735,35622939,35624277,36713024,36715087,37016726,37208172,37208293,40481896,42538697,42538946,42873163,43020424,43021830,44783643,44783644,44784483,44784484,44809026,44809027,44809548,44809569,44811110,44811932,44811933,45757119,45757137,45757356,45757444,45757445,45757446,45757447,45757787,45757788,45768449,45771064,45771067)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% hpt) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Hypertension"

# dyslipidemia
dyl <- c(432867,437521,437827,438720,440360,4029258,4029259,4029260,4029261,4029263,4029305,4029890,4029891,4029892,4030586,4030618,4030619,4031945,4031947,4079876,4079885,4079886,4079887,4096215,4104485,4120314,4134862,4142496,4143177,4144326,4144529,4159131,4220010,4223495,4270878,4291436,4292672,4294296,4294297,4295608,4295609,4298010,4298723,4298733,4299409,4300461,4301409,35608140,36674388,36676683,36715325,37016144,37016353,40482885,43530660,43531564,43531651,45757265,45757280,45757432,45757500,45770880)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dyl) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: dyslipidemia"

# Diabets mellitus
dm <- c(201254,201826,443412,3191208,3192767,3193274,3194082,3194119,3194332,4047906,4063042,4063043,4099214,4099215,4099651,4102018,4129519,4130162,4145827,4193704,4230254,4304377,36684827,36717215,42535539,43531008,43531009,43531010,45757474,45757674,45766051,45766052)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dm) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: diabets mellitus"

# Medical history: MI
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
mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: MI"

# Medical history: PCI
pci <- c(43531440,2617369,4225903,4264286,43527999,4265293,45770795,4214516,43527998,4020653,44511273,43527997,4171077,44511532,3190422,44784573,44789455,4175997,2313811,43527994,43531438,37111313,4283892,44512256,43527996,43527909,3180695,44511269,3173558,44511272,4216356,3181809,2313802,43533247,44511270,43533352,4328103,4181025,4006788,44511131,4337738,4329263,44511130,2313803,2313804,4139198,44511268,35607959,40756929,43531439,4264285,2313801,44511133,3170433,43527908,4178148,2617370,43527995,2001505,43533248,44511271,4238755,2000064,43533353,2313810,2001506,3186322)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% pci) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: PCI"

# Medical history: CABG
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: CABG"

# Medical history: CHF
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# Medical history: Nonhemorrhagic stroke 
ich <- c(372435,377254,379778,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4077086,4108356,4110189,4110190,4110192,4111714,4119140,4129534,4131383,4138327,4141405,4142739,4145897,4146185,4298750,4319146,35610084,35610085,36717605,37110678,37110679,37116473,37395562,40479572,42535465,42535466,42535523,42535524,43530683,43531607,44782773,45767658,45772786,46270031,46270380,46270381,46273649)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% ich) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Nonhemorrhagic stroke "

# Medical history: peripheral arterial disease
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: peripheral arterial disease"

# Medical history: chronic renal disease
ckd <- c(193782,198185,443597,443601,443611,443612,443614,443961,760850,762000,762001,762033,762034,762973,764011,765535,765536,3169303,3177485,3178768,3185897,4030520,4125970,4128067,4128200,4239804,4322556,36679002,36716184,36716455,36716947,36717349,36717534,37017104,37017813,37018761,37018886,37019193,42599648,42599649,42599650,42599651,42872405,43020437,43020455,43020456,43020457,43021835,43021836,43021852,43021853,43021854,43021864,43531559,43531562,43531566,43531577,43531578,43531653,44782429,44782689,44782690,44782691,44782692,44782703,44782717,44784621,44784637,44784638,44784639,44784640,44792226,44792227,44792228,44792229,44792230,44792231,44792232,44792249,44792250,44792251,44792252,44792253,44792254,44792255,45757137,45757139,45757356,45757392,45757393,45757444,45757445,45757446,45757447,45763854,45763855,45768813,45769901,45769902,45769903,45769904,45769906,45771064,45771067,45771075,45772751,45773576,45773688,46270353,46270354,46270355,46270356,46271022,46273164,46273514,46273636,46284566,46284567,46284570,46284572,46284575,46284587,46284588,46284591,46284592,46284593,46284597,46284598,46284599,46284600,46284602,46284603,46286992,46287169)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% ckd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: chronic renal disease"

# Medical history: history of dyspnea
dyspnea <- c(312437,315361,3655952,3656669,4010972,4047610,4059021,4059022,4060052,4060439,4060709,4065915,4066850,4091788,4094132,4097311,4144682,4148800,4178416,4190875,4192279,4193263,4206307,4212233,4217021,4219335,4219740,4220640,4220641,4221237,4222444,4222446,4222908,4223033,4223034,4223035,4223036,4223037,4223038,4228754,4244276,4248284,4253185,4263848,4307188,4308377,4310059,4310172,35610139,35610140,35610141,35610142,35610143,35610144,36685567,36685568,36685569,36685570,36685571,44809070)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% dyspnea) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of dyspnea"

# Medical history: COPD
copd <- c(255573,257004,259043,261325,261895,440748,3185549,3654571,3654836,3654837,4046986,4050732,4050733,4050734,4050961,4056405,4083395,4110048,4110056,4110635,4112828,4112836,4115044,4136683,4145496,4148124,4166508,4166517,4177944,4193588,4196712,4209097,4246105,4281815,4286497,4315386,4321845,42574216,42598711,43530693,44791725,44807895,45769389,46269701,46274062)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% copd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: COPD"
# Medical history: asthma
asthma <- c(252658,252942,256448,257581,257583,312950,313236,317009,443801,761844,764677,764949,2108536,3661412,4022592,4051466,4057952,4075237,4080516,4081994,4110051,4119298,4119300,4120261,4120262,4123253,4125022,4138760,4141978,4142738,4143474,4143828,4145356,4145497,4146581,4152913,4155468,4155469,4155470,4191479,4206340,4211530,4212099,4217558,4225553,4225554,4232595,4233784,4235703,4244339,4245292,4245676,4250128,4271333,4277596,4279553,4301938,4309833,4312524,35609846,35609847,36684328,36684335,37108580,37108581,37109103,37116845,37206717,37208352,37310241,40481763,40483397,42535716,42536207,42536208,42536649,42538744,42539549,43530693,43530745,44783610,44810117,45757063,45766727,45766728,45768910,45768911,45768912,45768963,45768964,45768965,45769350,45769351,45769352,45769438,45769441,45769442,45769443,45771045,45772073,45772937,45773005,46269767,46269770,46269771,46269776,46269777,46269784,46269785,46269801,46269802,46270028,46270029,46270030,46270082,46270322,46273452,46273454,46273462,46273487,46273635,46274059,46274062,46274124)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% asthma) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: asthma"
# Medical history: gout
gout <- c(74892,376955,440674,762219,4035436,4035437,4035438,4035751,4035752,4039287,4079745,4083344,4084229,4099013,4101425,4104957,4104958,4104959,4114911,4114912,4116007,4116008,4116009,4153103,4154437,4179065,4182267,4188214,4247334,4281849,4285308,4297662,4299408,35615129,35615130,36674930,36684526,36684527,36684528,36684529,36684530,36684531,36684532,36684533,36684534,36684535,36714445,36715572,37108592,37108594,37209406,37209407,37209408,37209463,37394509,42534842,42534843,42535119,42535120,42535121,42535122,42535127,42535133,42535203,42535204,42535205,42535206,42539209,42539261,42539266,42539288,44783578,45770178,46270412,46270413,46270416,46270417,46270418,46270419,46270420,46270421,46270422,46270423,46270424,46270425,46270426,46270427,46270428,46270429,46270430,46270431,46270432,46270435,46270436,46270437,46270438,46270439,46270440,46270441,46270442,46270443,46270444,46270445,46270446,46270447,46270448,46270449,46270450,46270451,46270452,46270453,46270454,46270455,46270456,46270457,46270458,46270459,46270460,46270461,46270464,46270465,46270466,46270467,46270468,46270470,46270471,46270472,46270473,46270474,46270475,46272380,46272381,46273507,46273508,46273523,46273541,46273638,46273639,46273641,46273648,46273653,46273655)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesHistory$covariates %>% filter(conceptId %in% gout) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: gout"

covariateSettingsFinalDiagnosis <- createCovariateSettings(useConditionOccurrenceShortTerm = TRUE,
                                                           shortTermStartDays = -7)

covariatesFinalDiagnosis <- getDbCovariateData(connectionDetails = connectionDetails,
                                               cdmDatabaseSchema = cdmDatabaseSchema,
                                               cohortDatabaseSchema = cohortDatabaseSchema,
                                               cohortTable = cohortTable,
                                               cohortId = c(targetNumber, comparatorNumber),
                                               covariateSettings = covariateSettingsFinalDiagnosis)

covariatesFinalDiagnosis$covariates <- covariatesFinalDiagnosis$covariates %>% left_join(covariatesFinalDiagnosis$covariateRef) 

# Final diagnosis of STEMI
stemi <- c(312327,319039,434376,436706,438170,438438,438447,441579,444406,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4178129,4243372,4267568,4275436,4296653,4303359,4324413,35610091,35610093,35611570,35611571,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% stemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: ST-elevation MI"

# Final diagnosis of NSTEMI
nstemi <- c(4270024)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% nstemi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Non-ST-elevation MI"

# Final diagnosis of Unstable angina
unstableangina <- c(315296)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariatesFinalDiagnosis$covariates %>% filter(conceptId %in% unstableangina) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: Unstable angina"

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
