# Load RCT information
tcList <- read.csv("C:/git/Troy/TroyCohortDiagnostics/inst/settings/TargetComparatorList.csv") #Should be change path
trials <- "EMPA-REG_OUTCOME"
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

# Result for EMPA-REG_OUTCOME trial
resultTablePooled <- data.frame()

# Mean age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Mean bodyweight
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3025315)))

# Mean BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Medical history: coronary artery disease
acsangina <- c(312327,315296,315830,315832,319039,321318,434376,436706,438170,438438,438447,441579,444406,761735,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4068938,4078531,4108669,4108670,4116486,4119455,4119456,4119457,4119942,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4155008,4155009,4155963,4161456,4161457,4161973,4161974,4178129,4184827,4198141,4201629,4231426,4243372,4262446,4264145,4267568,4270024,4275436,4296653,4303359,4310270,4324413,4324893,35610091,35610093,35611570,35611571,35615052,35615053,36712982,36712983,36712984,37209632,37309713,43020460,43531588,44782712,44782769,45766075,45766076,45766115,45766116,45766150,45766151,45771322,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% acsangina) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "coronary artery disease"

# CV risk factors
multivessel <- c(4124682, 4108673)
mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
stroke <- c(372435,372924,375557,376713,377254,379778,432923,439847,441874,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4077086,4108356,4110189,4110190,4110192,4111714,4119140,4129534,4131383,4138327,4141405,4142739,4145897,4146185,4148906,4298750,4319146,35610084,35610085,36717605,37110678,37110679,37116473,37395562,40479572,42535465,42535466,42535523,42535524,43530683,43530727,43531607,44782773,45767658,45772786,46270031,46270380,46270381,46273649)
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
singlevessel <- c(4111393)
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442)
cvriskfactors <-c(multivessel, mi, cabg, stroke, pvd, singlevessel, chf)

# CV risk factors
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% cvriskfactors) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "cvriskfactors"

# Medical history: multi-vessel artery disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% multivessel) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "multivessel disease"

# Medical history: MI
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "MI"

# Medical history: CABG
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: CABG"

# Medical history: stroke
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: stroke"

# Medical history: peripheral arterial disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: peripheral"

# Medical history: singlevessel artery disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% singlevessel) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: singlevessel artery disease"

# Medical history: CHF
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# glycated hemoglobin
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004410)))

# metformin 
metformin <- c(21600747, 1503297)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "metformin"

# insulin 
insulin  <- c(21600713,35605670,1531601,1567198,35602717,1516976,1502905,1544838,1588986,1550023,1513876,19078608,1590165,1596977,1586346,1513843,1513849,1562586,1586369)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% insulin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "insulin "

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

# Thiazolidinediones
thiazolidinediones <- c(21600779)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% thiazolidinediones) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "thiazolidinediones"

# Hypertension requiring treatment
hptdrugs <- c(1319998,1332418,1314002,1335471,1322081,1338005,1351557,1340128,950370,1346823,19049145,1395058,19050216,19089969,1328165,1341927,1346686,19063575,1353776,1363749,974166,19122327,978555,1347384,1326012,1386957,19004539,19015802,1308216,1367500,19071995,19102106,907013,1307046,1310756,1314577,1318137,1318853,19113063,1319133,1319880,19020061,40226742,19024904,1327978,1373225,1345858,1353766,1331235,1334456,1317640,1342439,1308842,1307863,19010493,19102107)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hptdrugs) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "1>= antihypertension drugs"

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

#	Mineralocorticoids
mineralocorticoids <- c(21602724)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% mineralocorticoids) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "mineralocorticoids"

##  Renin inhibitor 
renininhibitor <- c(1332418, 974166, 40256434)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% renininhibitor) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "renininhibitor"

##  Lipid lowering therapy
lipiddrug <- c(21601853)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% lipiddrug) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "lipiddrug"

# statin
statins <- c(21601853) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins"

# fibrates 
fibrates <- c(21601864)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% fibrates) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "fibrates"

# ezetimibe
ezetimibe <- c(1526475) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ezetimibe) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ezetimibe"

# niacin
niacin <- c(19018419)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% niacin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "niacin"


#anticoagulants
aspirin <- c(1112807)
clopidogrel <- c(1322184)
vka <- c(21600962, 1310149)
anticoagulants <- c(aspirin, clopidogrel, vka)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% anticoagulants) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "anticoagulants"

# aspirin
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

# clopidogrel
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% clopidogrel) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "clopidogrel"

# Vitamin K antagonist
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% vka) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Vitamin K antagonist"

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

# Diastolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3012888)))

# 	Cholesterol [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027114)))

# 	Cholesterol in LDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3028437)))

# 	Cholesterol in HDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3007070)))

#	Triglyceride [Mass/volume] in Serum or Plasma
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3022192)))

#EGFR
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3029859))[1,])

# Renal function, EGFR
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue > 90) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- ">90ml/min"

# Renal function, EGFR
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue < 90, covariateValue > 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "60 to 90<l/min"

# Renal function, EGFR
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue > 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "60</min"

# Normoalbuminuria
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue < 30) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Normoalbuminuria"

# Microalbuminuria
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue >= 30, covariateValue <=300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Microalbuminuria"

# Renal function, creatinine clearance
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3027108, covariateValue > 300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Macroalbuminuria"



#########

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

# Result for EMPA-REG_OUTCOME trial
resultTablePooled <- data.frame()

# Mean age
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter (covariateId == 1002)))

# Mean bodyweight
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3025315)))

# Mean BMI
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 40762636)))

# Medical history: coronary artery disease
acsangina <- c(312327,315296,315830,315832,319039,321318,434376,436706,438170,438438,438447,441579,444406,761735,761736,761737,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4051874,4068938,4078531,4108669,4108670,4116486,4119455,4119456,4119457,4119942,4119943,4119944,4119945,4119946,4119947,4119948,4121464,4121465,4121466,4124684,4124685,4126801,4145721,4151046,4155008,4155009,4155963,4161456,4161457,4161973,4161974,4178129,4184827,4198141,4201629,4231426,4243372,4262446,4264145,4267568,4270024,4275436,4296653,4303359,4310270,4324413,4324893,35610091,35610093,35611570,35611571,35615052,35615053,36712982,36712983,36712984,37209632,37309713,43020460,43531588,44782712,44782769,45766075,45766076,45766115,45766116,45766150,45766151,45771322,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% acsangina) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "coronary artery disease"

# CV risk factors
multivessel <- c(4124682, 4108673)
mi <- c(312327,314666,319039,434376,436706,438170,438438,438447,439693,441579,444406,761736,761737,765132,3189643,3654465,3654466,3654467,3655133,3661502,3661503,3661504,3661520,3661524,3661547,3661641,3661642,3661643,3661644,3661645,3661646,4030582,4051874,4108217,4108218,4108669,4108677,4119456,4119457,4119943,4119944,4119945,4119946,4119947,4119948,4119949,4119950,4121464,4121465,4121466,4121467,4121468,4124684,4124685,4124686,4126801,4138833,4145721,4151046,4170094,4173632,4178129,4200113,4206867,4207921,4209541,4215259,4243372,4267568,4270024,4275436,4296653,4303359,4323202,4324413,4329847,35610087,35610089,35610091,35610093,35611570,35611571,37309626,43020460,44782712,44782769,45766075,45766076,45766113,45766114,45766115,45766116,45766150,45766151,45766241,45771322,45773170,46270158,46270159,46270160,46270161,46270162,46270163,46270164,46273495,46274044)
cabg <- c(2001509,2100872,2100873,2107216,2107217,2107218,2107219,2107220,2107221,2107222,2107223,2107224,2107226,2107227,2107228,2107231,2107242,2107243,2107244,2107250,2108631,2617584,2721131,2721132,2721133,2721134,2724714,2724715,2724716,2724717,2724718,2724719,2724720,2724721,2724722,2724723,2724724,2724725,2724726,2724727,2724728,2724729,2724730,2724731,2724732,2724733,2724734,2724735,2724736,2724737,2724738,2724739,2724740,2724741,2724742,2724743,2724744,2724745,2724746,2724747,2724748,2724749,2724750,2724751,2724752,2724753,2724754,2724755,2724756,2724757,2724758,2724759,2724760,2724761,2724762,2724763,2724764,2724765,2724766,2725042,2725063,2725084,2725105,2725213,2725214,2725215,2725216,2725217,2725218,2725219,2725220,2725221,2725222,2725223,2725224,2725225,2725226,2725227,2725228,2725229,2725230,2725231,2725232,2725233,2725234,2725235,2725236,2725237,2725238,2725239,2725240,2725241,2725242,2725243,2725244,2725245,2725246,2725247,2725248,2725249,2725250,2725251,2725252,2725253,2725254,2725255,2725256,2725257,2725258,2725259,2725260,2725261,2725262,2725263,2725264,2725265,2725266,2725267,2725268,2725269,2725270,2725271,2725272,2725273,2725274,2725275,2725276,2725277,2725278,2725279,2725280,2725281,2725282,2725283,2725284,2725285,2725286,2725287,2725288,2725289,2725290,2725291,2725292,2725293,2725294,2725295,2725296,2725297,2725298,2725299,2725300,2725301,2725302,2725303,2725304,2725305,2725306,2725307,2725308,2725309,2725310,2725311,2725312,2725313,2725314,2725315,2725316,2725317,2725318,2725319,2725320,2725321,2725322,2725323,2725324,2725325,2725326,2725327,2725328,2725329,2725330,2725331,2725332,2725333,2725334,2725335,2725336,2725337,2725338,2725339,2725340,2725341,2725342,2725343,2725344,2725345,2725346,2725347,2725348,2725349,2725350,2725351,2725352,2725353,2725354,2725355,2725356,2725357,2725358,2725359,2725360,2725361,2725362,2725363,2725364,2725365,2725366,2725367,2725368,2725369,2725370,2725371,2725372,2725373,2725374,2725375,2725376,2725377,2725378,2725379,2725380,2725381,2725382,2725383,2725384,2725385,2725386,2725387,2725388,2725389,2725390,2725391,2725392,2725393,2725394,2725395,2725396,2725397,2725398,2725399,2725400,2725401,2725402,2725403,2725404,2725405,2725406,2725407,3655761,3655762,3655763,3655764,3656025,3656026,3656037,3656038,3663267,4000732,4000733,4008625,4011931,4018579,4018692,4018693,4018762,4020213,4031996,4063237,4106548,4140107,4146972,4148779,4168141,4169964,4173645,4189169,4219321,4228304,4228305,4229433,4233420,4233421,4234990,4253805,4284104,4302815,4305509,4309432,4336464,4336465,4336466,4336467,4337056,4337737,4339629,37111313,42537524,42537525,42537526,42537527,42537528,42537529,42537530,42537531,42537532,42537533,42537534,42539671,42894223,42894224,42894225,42894226,42894227,42894228,42894229,42894230,42894231,42894232,42894233,42894234,42894235,42894236,42894237,42894238,42894239,42894240,42894241,42894242,42894243,42894244,42894245,42894246,42894247,42894248,42894249,42894250,42894251,42894252,42894253,42894254,42894255,42894256,42894257,42894258,42894398,42894399,42894400,42894401,42894402,42894403,42894404,42894405,42894406,42894407,42894408,42894409,43528000,43528001,43528002,43528003,43528004,44511074,44511075,44511077,44511079,44511082,44511083,44511085,44511086,44511087,44511088,44511089,44511091,44511092,44511093,44511094,44511095,44511099,44511107,44511110,44511114,44806690)
stroke <- c(372435,372924,375557,376713,377254,379778,432923,439847,441874,443454,443790,443864,444091,761110,762933,762934,762935,762937,762951,763015,765515,3179950,3184775,3185378,3188164,4031045,4043731,4043732,4045735,4045737,4045738,4045740,4045741,4046089,4046090,4046237,4046358,4046359,4046360,4046361,4046362,4048784,4077086,4108356,4110189,4110190,4110192,4111714,4119140,4129534,4131383,4138327,4141405,4142739,4145897,4146185,4148906,4298750,4319146,35610084,35610085,36717605,37110678,37110679,37116473,37395562,40479572,42535465,42535466,42535523,42535524,43530683,43530727,43531607,44782773,45767658,45772786,46270031,46270380,46270381,46273649)
pvd <- c(134380,192657,195834,200132,314965,315558,317309,318712,321052,321822,434961,439836,441044,443407,443408,443729,444264,760869,760870,760871,761206,761428,761429,761430,761431,761432,761433,761434,761463,761464,761465,761466,761467,761740,761741,761742,761743,761744,761748,761749,761804,761805,761812,761813,761814,761815,761816,761822,761823,761824,761825,761826,761827,761828,761829,761853,761854,761858,762053,764026,765079,765080,765226,765284,765300,765301,3171253,3173537,3173780,3176677,3176693,3181089,3184873,3192767,3198350,3654996,4025854,4030663,4033350,4055025,4055720,4108375,4108376,4108377,4108380,4108381,4110331,4110333,4111848,4111850,4112161,4112162,4114011,4115229,4116977,4117679,4117933,4118795,4119141,4119612,4119614,4120098,4121625,4121626,4121627,4124836,4124837,4124838,4124839,4127997,4131908,4136335,4137550,4137551,4139534,4139577,4139578,4139579,4140616,4141106,4141975,4143588,4146012,4151848,4153293,4168060,4173167,4174012,4174013,4174014,4174015,4175569,4175570,4178399,4178400,4178609,4179907,4179908,4179909,4179910,4182191,4193056,4193057,4193302,4194886,4195971,4195972,4195973,4199183,4199887,4200875,4202511,4208071,4208072,4208809,4208944,4226026,4231816,4231826,4235591,4236905,4246643,4253510,4263089,4263648,4264734,4289307,4291464,4294425,4299117,4299127,4301022,4301023,4301024,4316367,4318843,4329498,4337367,4341646,4348039,4348319,35611566,35615018,35615019,35615020,35615021,35615028,35615036,35615054,35615061,35615062,35615063,35615071,35615073,35615075,35615076,35615080,35615081,35615082,35615083,35615086,35615087,35615088,35615089,35615090,35615091,35615092,35615093,35615094,35615095,35615096,35615097,35615098,35615099,35615100,35615106,35615107,35615116,36712805,36712807,36712955,36712987,36713012,36713013,36713014,36713015,36713094,36717256,37016147,37017533,37109921,37109923,37110250,37110251,37116421,37116422,37209668,37312530,40483538,40484541,40484551,40484912,42535143,42535335,42535829,42536634,42536636,42539410,42572961,42597028,42597030,42599607,42599802,42599894,43020443,43021846,43022064,44782775,44782819,44808745,44808746,44808747,44808832,44808833,44813802,44813823,46271459,46271460,46271462)
singlevessel <- c(4111393)
chf <- c(314378,316994,319835,439694,439696,439698,762002,762003,764874,4023479,4139864,4142561,4206009,4215446,4229440,4242669,4284562,4327205,36712927,36712928,36713488,37309625,43021825,43021826,43022068,44782428,44782655,44782713,44782728,44784345,44784442)
cvriskfactors <-c(multivessel, mi, cabg, stroke, pvd, singlevessel, chf)

# CV risk factors
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% cvriskfactors) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "cvriskfactors"

# Medical history: multi-vessel artery disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% multivessel) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "multivessel disease"

# Medical history: MI
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% mi) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "MI"

# Medical history: CABG
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% cabg) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: CABG"

# Medical history: stroke
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: stroke"

# Medical history: peripheral arterial disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% pvd) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "History of atherosclerotic vascular disease: peripheral"

# Medical history: singlevessel artery disease
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% singlevessel) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical history: singlevessel artery disease"

# Medical history: CHF
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% chf) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Medical histroy: Congestive Heart Failure"

# glycated hemoglobin
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004410)))

# metformin 
metformin <- c(21600747, 1503297)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% metformin ) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "metformin"

# insulin 
insulin  <- c(21600713,35605670,1531601,1567198,35602717,1516976,1502905,1544838,1588986,1550023,1513876,19078608,1590165,1596977,1586346,1513843,1513849,1562586,1586369)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% insulin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "insulin "

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

# Thiazolidinediones
thiazolidinediones <- c(21600779)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% thiazolidinediones) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "thiazolidinediones"

# Hypertension requiring treatment
hptdrugs <- c(1319998,1332418,1314002,1335471,1322081,1338005,1351557,1340128,950370,1346823,19049145,1395058,19050216,19089969,1328165,1341927,1346686,19063575,1353776,1363749,974166,19122327,978555,1347384,1326012,1386957,19004539,19015802,1308216,1367500,19071995,19102106,907013,1307046,1310756,1314577,1318137,1318853,19113063,1319133,1319880,19020061,40226742,19024904,1327978,1373225,1345858,1353766,1331235,1334456,1317640,1342439,1308842,1307863,19010493,19102107)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% hptdrugs) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "1>= antihypertension drugs"

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

#	Mineralocorticoids
mineralocorticoids <- c(21602724)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% mineralocorticoids) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "mineralocorticoids"

##  Renin inhibitor 
renininhibitor <- c(1332418, 974166, 40256434)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% renininhibitor) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "renininhibitor"

##  Lipid lowering therapy
lipiddrug <- c(21601853)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% lipiddrug) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "lipiddrug"

# statin
statins <- c(21601853) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% statins) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "statins"

# fibrates 
fibrates <- c(21601864)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% fibrates) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "fibrates"

# ezetimibe
ezetimibe <- c(1526475) 
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% ezetimibe) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "ezetimibe"

# niacin
niacin <- c(19018419)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% niacin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "niacin"


#anticoagulants
aspirin <- c(1112807)
clopidogrel <- c(1322184)
vka <- c(21600962, 1310149)
anticoagulants <- c(aspirin, clopidogrel, vka)
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% anticoagulants) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "anticoagulants"

# aspirin
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% aspirin) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "aspirin"

# clopidogrel
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% clopidogrel) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "clopidogrel"

# Vitamin K antagonist
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %in% vka) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Vitamin K antagonist"

# Systolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3004249)))

# Diastolic blood pressure — mm Hg
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3012888)))

# 	Cholesterol [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3027114)))

# 	Cholesterol in LDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3028437)))

# 	Cholesterol in HDL [Mass/volume] in Serum or Plasma	
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3007070)))

#	Triglyceride [Mass/volume] in Serum or Plasma
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3022192)))

#EGFR
resultTablePooled <- rbind(resultTablePooled, data.frame(statisticsPooled %>% filter(conceptId %like% 3029859))[1,])

# Renal function, EGFR
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue > 90) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- ">90ml/min"

# Renal function, EGFR
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue < 90, covariateValue > 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "60 to 90<l/min"

# Renal function, EGFR
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3029859, covariateValue > 60) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "60</min"

# Normoalbuminuria
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue < 30) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Normoalbuminuria"

# Microalbuminuria
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue >= 30, covariateValue <=300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Microalbuminuria"

# Renal function, creatinine clearance
resultTablePooled[nrow(resultTablePooled)+1, c("n")] <- data.frame(covariates$covariates %>% filter(conceptId %like% 3020682, covariateValue > 300) %>% summarise(n = n_distinct(.data$rowId)))
resultTablePooled[nrow(resultTablePooled), c("covariateName")] <- "Macroalbuminuria"

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
