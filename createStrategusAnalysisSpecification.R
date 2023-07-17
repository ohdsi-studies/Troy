########## Please populate the information below #####################
# Name: "Trials Replication through Observational study by Yonsei - an OHDSI network study"
# Packagename: Troy
# Studyleader: Kyoung Won Kim<kwkim@yuhs.ac>, Seng Chan You<seng@ohdsi.org>
# Createdby: Jaehyeong Cho <jhcho03@yuhs.ac>
# Createddate: 2023-07-14
# Modifiedby: Jaehyeong Cho <jhcho03@yuhs.ac>
# Modifieddate: 
# Organizationname: YUHS
# Description: "Strategus on TROY."

# Install correct versions of HADES packages
# options(install.packages.compile.from.source = "never")
# remotes::install_github("ohdsi/CohortGenerator", ref = "v0.8.0")
# remotes::install_github("ohdsi/CohortDiagnostics", ref = "v3.1.2")
# remotes::install_github("ohdsi/Characterization", ref = "v0.1.1")
# remotes::install_github("ohdsi/CohortIncidence", ref = "v3.1.2")
# remotes::install_github("ohdsi/CohortMethod", ref = "74f017107e0cc1b740a2badc82879ab6ad291b23")
# remotes::install_github("ohdsi/SelfControlledCaseSeries", ref = "15918616814b88137f82bf2aa9986e1dcdf39e74")
# remotes::install_github("ohdsi/PatientLevelPrediction", ref = "v6.3.1")
# remotes::install_github("ohdsi/ROhdsiWebApi")
#
# Interrogate the installed packages to confirm versions above
# installedPackages <- as.data.frame(installed.packages())[,c("Package", "Version")]
# installedPackages[installedPackages$Package %in% c("CohortGenerator", "CohortDiagnostics", "Characterization", "CohortIncidence", "CohortMethod", "SelfControlledCaseSeries", "PatientLevelPrediction", "ROhdsiWebApi"), ]

#"~/Troy/inst"
Sys.setlocale(category="LC_CTYPE", locale="C")

library(CohortGenerator)
cohortDefinitionSet <- getCohortDefinitionSet(
  settingsFileName = "settings/CohortsToCreate.csv",
  jsonFolder = "cohorts",
  sqlFolder = "sql/sql_server",
  packageName = "TroyCohortDiagnostics"
)
ncoCohortSet <- readCsv(file = system.file("testdata/negative_controls_concept_set.csv",
                                           package = "Strategus"
))

kable(cohortDefinitionSet[, c("cohortId", "cohortName")])
kable(ncoCohortSet)

source("https://raw.githubusercontent.com/OHDSI/CohortGeneratorModule/v0.1.0/SettingsFunctions.R")

# Create the cohort definition shared resource element for the analysis specification
cohortDefinitionSharedResource <- createCohortSharedResourceSpecifications(
  cohortDefinitionSet = cohortDefinitionSet
)

# Create the negative control outcome shared resource element for the analysis specification
ncoSharedResource <- createNegativeControlOutcomeCohortSharedResourceSpecifications(
  negativeControlOutcomeCohortSet = ncoCohortSet,
  occurrenceType = "all",
  detectOnDescendants = TRUE
)

studyStartDate <- '19700101' 
studyEndDate <- '20991231'   

# Probably don't change below this line ----------------------------------------

useCleanWindowForPriorOutcomeLookback <- FALSE # If FALSE, lookback window is all time prior, i.e., including only first events
psMatchMaxRatio <- 1 # If bigger than 1, the outcome model will be conditioned on the matched set

# Create the module specification
cohortGeneratorModuleSpecifications <- createCohortGeneratorModuleSpecifications(
  incremental = TRUE,
  generateStats = TRUE
)

source("https://raw.githubusercontent.com/OHDSI/CohortGeneratorModule/v0.1.0/SettingsFunctions.R")

# Create the cohort definition shared resource element for the analysis specification
cohortDefinitionSharedResource <- createCohortSharedResourceSpecifications(
  cohortDefinitionSet = cohortDefinitionSet
)

# Create the negative control outcome shared resource element for the analysis specification
ncoSharedResource <- createNegativeControlOutcomeCohortSharedResourceSpecifications(
  negativeControlOutcomeCohortSet = ncoCohortSet,
  occurrenceType = "all",
  detectOnDescendants = TRUE
)

# Create the module specification
cohortGeneratorModuleSpecifications <- createCohortGeneratorModuleSpecifications(
  incremental = TRUE,
  generateStats = TRUE
)

source("https://raw.githubusercontent.com/OHDSI/CohortDiagnosticsModule/v0.0.7/SettingsFunctions.R")
cohortDiagnosticsModuleSpecifications <- createCohortDiagnosticsModuleSpecifications(
  runInclusionStatistics = TRUE,
  runIncludedSourceConcepts = TRUE,
  runOrphanConcepts = TRUE,
  runTimeSeries = FALSE,
  runVisitContext = TRUE,
  runBreakdownIndexEvents = TRUE,
  runIncidenceRate = TRUE,
  runCohortRelationship = TRUE,
  runTemporalCohortCharacterization = TRUE,
  incremental = FALSE
)

source("https://raw.githubusercontent.com/OHDSI/CharacterizationModule/v0.3.2/SettingsFunctions.R")
characterizationModuleSpecifications <- createCharacterizationModuleSpecifications(
  targetIds = cohortDefinitionSet$cohortId[1:4],
  outcomeIds = cohortDefinitionSet$cohortId[4:length(cohortDefinitionSet$cohortId)],
  covariateSettings = FeatureExtraction::createDefaultCovariateSettings(),
  dechallengeStopInterval = 30,
  dechallengeEvaluationWindow = 30,
  timeAtRisk = data.frame(
    riskWindowStart = c(1, 1),
    startAnchor = c("cohort start", "cohort start"),
    riskWindowEnd = c(0, 365),
    endAnchor = c("cohort end", "cohort end")
  )
)

# CohortMethodModule -----------------------------------------------------------
source("https://raw.githubusercontent.com/OHDSI/CohortMethodModule/v0.1.0/SettingsFunctions.R")
covariateSettings <- FeatureExtraction::createDefaultCovariateSettings(
  addDescendantsToExclude = TRUE # Keep TRUE because you're excluding concepts
)
########### TCO
tcis <- list(
  # Study population: replication study 1 (LEADER trial)
  list(
    targetId = 122, # Liraglutide (LEADER trial)
    comparatorId = 123, # DPP-4 (LEADER trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      40170911, # Liraglutide
      1580747,19122137,40166035,40239216,43009051,43009070,43009089,43013884 # DPP-4
    ) 
  ),
  # Study population: replication study 2 (DECLARE-TIMI 58 trial)
  list(
    targetId = 125, # Dapagliflozin (DECLARE-TIMI 58 trial)
    comparatorId = 126, # DPP-4 (DECLARE-TIMI 58 trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      44785829, # Dapagliflozin
      1580747,19122137,40166035,40239216,43009051,43009070,43009089,43013884 # DPP-4
    ) 
  ),
  # Study population: replication study 3 (EMPA-REG OUTCOME trial)
  list(
    targetId = 127, # Empagliflozin (EMPA-REG OUTCOME trial)
    comparatorId = 128, # DPP-4 (EMPA-REG OUTCOME trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      45774751, # Empagliflozin
      1580747,19122137,40166035,40239216,43009051,43009070,43009089,43013884 # DPP-4
    ) 
  ),
  # Study population: replication study 4 (CANVAS trial)
  list(
    targetId = 129, # Canagliflozin (CANVAS trial)
    comparatorId = 130, # DPP-4 (CANVAS trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      43526465, # Canagliflozin
      1580747,19122137,40166035,40239216,43009051,43009070,43009089,43013884 # DPP-4
    ) 
  ),
  # Study population: replication study 5 (CARMELINA trial)
  list(
    targetId = 135, # Linagliptin (CARMELINA trial)
    comparatorId = 136, # Sulfonylureas (CARMELINA trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      40239216, # Linagliptin
      1502809,1502855,1530014,1559684,1560171,1594973,1597756,19001409,19033498,19059796,19097821,40798860 # SU
    ) 
  ),
  # Study population: replication study 6 (TECOS trial)
  list(
    targetId = 137, # Sitagliptin (TECOS trial)
    comparatorId = 138, # Sulfonylureas (TECOS trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      1580747, # Sitagliptin
      1502809,1502855,1530014,1559684,1560171,1594973,1597756,19001409,19033498,19059796,19097821,40798860 # SU
    ) 
  ),
  # Study population: replication study 7 (SAVOR-TIMI 53 trial)
  list(
    targetId = 139, # Saxagliptin (SAVOR-TIMI 53 trial)
    comparatorId = 140, # Sulfonylureas (SAVOR-TIMI 53 trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      40166035, # Saxagliptin
      1502809,1502855,1530014,1559684,1560171,1594973,1597756,19001409,19033498,19059796,19097821,40798860 # SU
    ) 
  ),
  # Study population: replication study 8 (CAROLINA trial)
  list(
    targetId = 141, # Linagliptin (CAROLINA trial)
    comparatorId = 142, # Sulfonylureas (CAROLINA trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      40239216, # Linagliptin
      1502809,1502855,1530014,1559684,1560171,1594973,1597756,19001409,19033498,19059796,19097821,40798860 # SU
    ) 
  ),
  # Study population: replication study 9 (TRITON-TIMI 38 trial)
  list(
    targetId = 117, # Prasugrel (TRITON-TIMI 38 trial)
    comparatorId = 118, # Clopidogrel (TRITON-TIMI 38 trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      40163718, # Prasugrel
      1322184 # Clopidogrel
    ) 
  ),
  # Study population: replication study 10 (PLATO trial)
  list(
    targetId = 115, # Ticagrelor (PLATO trial)
    comparatorId = 116, # Clopidogrel (PLATO trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      40241186, # Ticagrelor
      1322184 # Clopidogrel
    ) 
  ),
  # Study population: replication study 11 (ROCKET-AF trial)
  list(
    targetId = 107, # Rivaroxaban (ROCKET-AF trial)
    comparatorId = 108, # Warfarin (ROCKET-AF trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      # intervention drugs
      40241331,# Rivaroxaban
      1310149 # Warfarin
    ) 
  ),
  # Study population: replication study 12 (ARISTOTLE trial)
  list(
    targetId = 105, # Apixaban (ARISTOTLE trial)
    comparatorId = 106, # Warfarin (ARISTOTLE trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      43013024,# Apixaban
      1310149 # Warfarin
    ) 
  ),
  # Study population: replication study 13 (ENGAGE AF-TIMI 48 trial)
  list(
    targetId = 109, # Edoxaban (ENGAGE AF-TIMI 48 trial)
    comparatorId = 110, # Warfarin (ENGAGE AF-TIMI 48 trial)
    #indicationId = , # 
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      45892847,# Edoxaban
      1310149 # Warfarin
    ) 
  )
)

outcomes <- tibble(
  cohortId = c(9, 10, 11, 12), # 3P MACE, 4P MACE, HHF + CV death, Storke + systemic embolism 
  cleanWindow = c(0, 0, 0, 0)
)

outcomeList <- append(
  lapply(seq_len(nrow(outcomes)), function(i) {
    if (useCleanWindowForPriorOutcomeLookback)
      priorOutcomeLookback <- outcomes$cleanWindow[i]
    else
      priorOutcomeLookback <- 99999
    CohortMethod::createOutcome(
      outcomeId = outcomes$cohortId[i],
      outcomeOfInterest = TRUE,
      trueEffectSize = NA,
      priorOutcomeLookback = priorOutcomeLookback
    )
  }),
  lapply(ncoCohortSet$cohortId, function(i) {
    CohortMethod::createOutcome(
      outcomeId = i,
      outcomeOfInterest = FALSE,
      trueEffectSize = 1
    )
  })
)

timeAtRisks <- tibble(
  label = c("ITT", "P-P"),
  riskWindowStart  = c(1, 1),
  startAnchor = c("cohort start", "cohort start"),
  riskWindowEnd  = c(0, 9999),
  endAnchor = c("cohort end", "cohort end"),
)

tars <- list()
for (i in seq_len(nrow(timeAtRisks))) {
  tars[[i]] <- CohortIncidence::createTimeAtRiskDef(
    id = i, 
    startWith = gsub("cohort ", "", timeAtRisks$startAnchor[i]), 
    endWith = gsub("cohort ", "", timeAtRisks$endAnchor[i]), 
    startOffset = timeAtRisks$riskWindowStart[i],
    endOffset = timeAtRisks$riskWindowEnd[i]
  )
}

targetComparatorOutcomesList <- list()
for (i in seq_along(tcis)) {
  tci <- tcis[[i]]
  targetId <- tci$targetId
  comparatorId <- tci$comparatorId
  excludedCovariateConceptIds <- tci$excludedCovariateConceptIds
  targetComparatorOutcomesList[[i]] <- CohortMethod::createTargetComparatorOutcomes(
    targetId = targetId,
    comparatorId = comparatorId,
    outcomes = outcomeList,
    excludedCovariateConceptIds = tci$excludedCovariateConceptIds
  )
}


getDbCohortMethodDataArgs <- CohortMethod::createGetDbCohortMethodDataArgs(
  restrictToCommonPeriod = TRUE,
  studyStartDate = studyStartDate,
  studyEndDate = studyEndDate,
  maxCohortSize = 0,
  covariateSettings = covariateSettings
)
createPsArgs = CohortMethod::createCreatePsArgs(
  maxCohortSizeForFitting = 250000,
  errorOnHighCorrelation = TRUE,
  stopOnError = FALSE, # Setting to FALSE to allow Strategus complete all CM operations; when we cannot fit a model, the equipoise diagnostic should fail
  estimator = "att",
  prior = createPrior(
    priorType = "laplace", 
    exclude = c(0), 
    useCrossValidation = TRUE
  ),
  control = createControl(
    noiseLevel = "silent", 
    cvType = "auto", 
    seed = 1, 
    resetCoefficients = TRUE, 
    tolerance = 2e-07, 
    cvRepetitions = 10, 
    startingVariance = 0.01
  )
)
matchOnPsArgs = CohortMethod::createMatchOnPsArgs(
  maxRatio = psMatchMaxRatio,
  caliper = 0.2,
  caliperScale = "standardized logit",
  allowReverseMatch = FALSE,
  stratificationColumns = c()
)
# stratifyByPsArgs <- CohortMethod::createStratifyByPsArgs(
#   numberOfStrata = 5,
#   stratificationColumns = c(),
#   baseSelection = "all"
# )
computeSharedCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(
  maxCohortSize = 250000,
  covariateFilter = NULL
)
computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(
  maxCohortSize = 250000,
  #covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
  covariateFilter = NULL
)
fitOutcomeModelArgs = CohortMethod::createFitOutcomeModelArgs(
  modelType = "cox",
  stratified = psMatchMaxRatio != 1,
  useCovariates = FALSE,
  inversePtWeighting = FALSE,
  prior = createPrior(
    priorType = "laplace", 
    useCrossValidation = TRUE
  ),
  control = createControl(
    cvType = "auto", 
    seed = 1, 
    resetCoefficients = TRUE,
    startingVariance = 0.01, 
    tolerance = 2e-07, 
    cvRepetitions = 10, 
    noiseLevel = "quiet"
  )
)
cmAnalysisList <- list()
for (i in seq_len(nrow(timeAtRisks))) {
  createStudyPopArgs <- CohortMethod::createCreateStudyPopulationArgs(
    firstExposureOnly = FALSE,
    washoutPeriod = 0,
    removeDuplicateSubjects = "keep first",
    censorAtNewRiskWindow = TRUE,
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback = 99999,
    riskWindowStart = timeAtRisks$riskWindowStart[[i]],
    startAnchor = timeAtRisks$startAnchor[[i]],
    riskWindowEnd = timeAtRisks$riskWindowEnd[[i]],
    endAnchor = timeAtRisks$endAnchor[[i]],
    minDaysAtRisk = 1,
    maxDaysAtRisk = 99999
  )
  cmAnalysisList[[i]] <- CohortMethod::createCmAnalysis(
    analysisId = i,
    description = sprintf(
      "Cohort method, %s",
      timeAtRisks$label[i]
    ),
    getDbCohortMethodDataArgs = getDbCohortMethodDataArgs,
    createStudyPopArgs = createStudyPopArgs,
    createPsArgs = createPsArgs,
    matchOnPsArgs = matchOnPsArgs,
    # stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovariateBalanceArgs,
    #computeCovariateBalanceArgs = computeCovariateBalanceArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgs
  )
}
cohortMethodModuleSpecifications <- createCohortMethodModuleSpecifications(
  cmAnalysisList = cmAnalysisList,
  targetComparatorOutcomesList = targetComparatorOutcomesList,
  analysesToExclude = NULL,
  refitPsForEveryOutcome = FALSE,
  refitPsForEveryStudyPopulation = FALSE
)
# Combine across modules -------------------------------------------------------
analysisSpecifications <- Strategus::createEmptyAnalysisSpecificiations() %>%
  Strategus::addSharedResources(cohortDefinitionSharedResource) %>%
  Strategus::addSharedResources(ncoSharedResource) %>%
  Strategus::addModuleSpecifications(cohortGeneratorModuleSpecifications) %>%
  Strategus::addModuleSpecifications(cohortDiagnosticsModuleSpecifications) %>%
  Strategus::addModuleSpecifications(characterizationModuleSpecifications) %>%
  Strategus::addModuleSpecifications(cohortMethodModuleSpecifications)

ParallelLogger::saveSettingsToJson(analysisSpecifications, file.path(getwd(), "analysisSpecification.json"))
