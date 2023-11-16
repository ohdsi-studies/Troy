########## Please populate the information below #####################
# Name: "Trials Replication through Observational study by Yonsei - an OHDSI network study"
# Packagename: Troy
# Studyleader: Kyoung Won Kim<kwkim@yuhs.ac>, Seng Chan You<seng@ohdsi.org>
# Createdby: Jaehyeong Cho <jhcho03@yuhs.ac>
# Createddate: 2023-07-14
# Modifiedby: Jaehyeong Cho <jhcho03@yuhs.ac>
# Modifieddate: 2023-11-1
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
library(dplyr)

outputFolder <- "~/output/TROY"
studyStartDate <- '19700101'
studyEndDate <- '20991231'

########### Cohorts: T and C ##############
tcis <- list(
  # Study population: replication study 1 (LEADER trial)
  list(
    targetId = 324,#1787712, # Liraglutide (LEADER trial)
    comparatorId = 325,#, # DPP-4 (LEADER trial)
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
    targetId = 326,#1787714, # Dapagliflozin (DECLARE-TIMI 58 trial)
    comparatorId = 327,#1787715, # DPP-4 (DECLARE-TIMI 58 trial)
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
    targetId = 328,#1787716, # Empagliflozin (EMPA-REG OUTCOME trial)
    comparatorId = 329,#1787717, # DPP-4 (EMPA-REG OUTCOME trial)
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
    targetId = 330,#1787718, # Canagliflozin (CANVAS trial)
    comparatorId = 331,#1787719, # DPP-4 (CANVAS trial)
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
    targetId = 332,#1787726, # Linagliptin (CARMELINA trial)
    comparatorId = 333,#1787727, # Sulfonylureas (CARMELINA trial)
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
    targetId = 334,#1787722, # Sitagliptin (TECOS trial)
    comparatorId = 335,#1787723, # Sulfonylureas (TECOS trial)
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
    targetId = 336,#1787728, # Saxagliptin (SAVOR-TIMI 53 trial)
    comparatorId = 337,#1787729, # Sulfonylureas (SAVOR-TIMI 53 trial)
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
    targetId = 338,#1787733, # Linagliptin (CAROLINA trial)
    comparatorId = 339,#1787734, # Sulfonylureas (CAROLINA trial)
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
    targetId = 340,#1787737, # Prasugrel (TRITON-TIMI 38 trial)
    comparatorId = 341,#1787738, # Clopidogrel (TRITON-TIMI 38 trial)
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
    targetId = 342,#1787741, # Ticagrelor (PLATO trial)
    comparatorId = 343,#1787742, # Clopidogrel (PLATO trial)
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
    targetId = 344,#1787745, # Rivaroxaban (ROCKET-AF trial)
    comparatorId = 345,#1787746, # Warfarin (ROCKET-AF trial)
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
    targetId = 346,#1787750, # Apixaban (ARISTOTLE trial)
    comparatorId = 347,#1787751, # Warfarin (ARISTOTLE trial)
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
    targetId = 348,#1787754, # Edoxaban (ENGAGE AF-TIMI 48 trial)
    comparatorId = 349,#1787755, # Warfarin (ENGAGE AF-TIMI 48 trial)
    #indicationId = , #
    #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
    #minAge = 35, # All ages In years. Can be NULL
    #maxAge = NULL, # All ages In years. Can be NULL
    excludedCovariateConceptIds = c(
      45892847,# Edoxaban
      1310149 # Warfarin
    )
  ),
  list(
    targetId = 351,#1787763, # Liraglutide (LEADER indication)
    comparatorId = 352,#1787762, # DPP-4 (LEADER indication)
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
    targetId = 353,#1787761, # Dapagliflozin (DECLARE-TIMI 58 indication)
    comparatorId = 354,#1787760, # DPP-4 (DECLARE-TIMI 58 indication)
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
    targetId = 355,#1787758, # Empagliflozin (EMPA-REG OUTCOME trial indication)
    comparatorId = 356,#1787759, # DPP-4 (EMPA-REG OUTCOME trial indication)
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
    targetId = 357,#1787720, # Canagliflozin (CANVAS trial indication)
    comparatorId = 358,#1787721, # DPP-4 (CANVAS trial indication)
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
    targetId = 359,#1787726, # Linagliptin (CARMELINA trial indication)
    comparatorId = 360,#1787727, # Sulfonylureas (CARMELINA trial indication)
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
    targetId = 361,#1787722, # Sitagliptin (TECOS trial indication)
    comparatorId = 362,#1787723, # Sulfonylureas (TECOS trial indication)
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
    targetId = 363,#1787728, # Saxagliptin (SAVOR-TIMI 53 trial indication)
    comparatorId = 364,#1787729, # Sulfonylureas (SAVOR-TIMI 53 trial indication)
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
    targetId = 365,#1787735, # Linagliptin (CAROLINA trial indication)
    comparatorId = 366,#1787736, # Sulfonylureas (CAROLINA trial indication)
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
    targetId = 367,#1787739, # Prasugrel (TRITON-TIMI 38 trial indication)
    comparatorId = 368,#1787740, # Clopidogrel (TRITON-TIMI 38 trial indication)
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
    targetId = 369,#1787743, # Ticagrelor (PLATO trial indication)
    comparatorId = 370,#1787744, # Clopidogrel (PLATO trial indication)
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
    targetId = 371,#1787748, # Rivaroxaban (ROCKET-AF trial indication)
    comparatorId = 372,#1787749, # Warfarin (ROCKET-AF trial indication)
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
    targetId = 373,#1787752, # Apixaban (ARISTOTLE trial indication)
    comparatorId = 374,#1787753, # Warfarin (ARISTOTLE trial indication)
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
    targetId = 375,#1787756, # Edoxaban (ENGAGE AF-TIMI 48 trial indication)
    comparatorId = 376,#1787757, # Warfarin (ENGAGE AF-TIMI 48 trial indication)
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
#
# # TEST ONLY
 # tcis <- list(
 #   # Study population: replication study 12 (ARISTOTLE trial)
 #   list(
 #     targetId = 346,#1787750, # Apixaban (ARISTOTLE trial)
 #     comparatorId = 347,#1787751, # Warfarin (ARISTOTLE trial)
 #     #indicationId = , #
 #     #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
 #     #minAge = 35, # All ages In years. Can be NULL
 #     #maxAge = NULL, # All ages In years. Can be NULL
 #     excludedCovariateConceptIds = c(
 #       43013024,# Apixaban
 #       1310149 # Warfarin
 #     )
 #   ),
 #   list(
 #     targetId = 373,#1787752, # Apixaban (ARISTOTLE trial indication)
 #     comparatorId = 374,#1787753, # Warfarin (ARISTOTLE trial indication)
 #     #indicationId = , #
 #     #genderConceptIds = c(8507, 8532), # use valid genders (remove unknown)
 #     #minAge = 35, # All ages In years. Can be NULL
 #     #maxAge = NULL, # All ages In years. Can be NULL
 #     excludedCovariateConceptIds = c(
 #       43013024,# Apixaban
 #       1310149 # Warfarin
 #     )
 #   )
 # )
########### Cohort: O ##############
outcomes <- tibble(
  cohortId = c(9, 10, 11, 12), # 3P MACE, 4P MACE, HHF + CV death, Storke + systemic embolism
  cleanWindow = c(0, 0, 0, 0)
)

# Probably don't change below this line ----------------------------------------
useCleanWindowForPriorOutcomeLookback <- FALSE # If FALSE, lookback window is all time prior, i.e., including only first events
psMatchMaxRatio <- 4 # If bigger than 1, the outcome model will be conditioned on the matched set

source("https://raw.githubusercontent.com/OHDSI/CohortGeneratorModule/v0.1.0/SettingsFunctions.R")

# Create the module specification
cohortGeneratorModuleSpecifications <- createCohortGeneratorModuleSpecifications(
  incremental = TRUE,
  generateStats = TRUE
)


# Shared Resources -------------------------------------------------------------
baseUrl <- 'http://34.148.35.102:80/WebAPI'  # Sys.getenv("baseUrl")

cohortDefinitionSet <- ROhdsiWebApi::exportCohortDefinitionSet(
  cohortIds =  unique(
    c(
      outcomes$cohortId,
      unlist(sapply(tcis, function(x) c(x$targetId, x$comparatorId, x$indicationId)))
    )
  ),
  generateStats = TRUE,
  baseUrl = baseUrl
)

negativeConceptSetId <- 78  #candidate controls for each trials

negativeControlOutcomeCohortSet <- ROhdsiWebApi::getConceptSetDefinition(
  conceptSetId = negativeConceptSetId,
  baseUrl = baseUrl
) %>%
  ROhdsiWebApi::resolveConceptSet(
    baseUrl = baseUrl
  ) %>%
  ROhdsiWebApi::getConcepts(
    baseUrl = baseUrl
  ) %>%
  rename(outcomeConceptId = "conceptId",
         cohortName = "conceptName") %>%
  mutate(cohortId = row_number() + 1000)

if (any(duplicated(cohortDefinitionSet$cohortId, negativeControlOutcomeCohortSet$cohortId))) {
  stop("*** Error: duplicate cohort IDs found ***")
  rstudioapi::showDialog("Error", "Duplicate cohort IDs found")
}

cohortDefinitionShared <- createCohortSharedResourceSpecifications(cohortDefinitionSet)
negativeControlsShared <- createNegativeControlOutcomeCohortSharedResourceSpecifications(
  negativeControlOutcomeCohortSet = negativeControlOutcomeCohortSet,
  occurrenceType = "first",
  detectOnDescendants = TRUE
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
  lapply(negativeControlOutcomeCohortSet$cohortId, function(i) {
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

# CohortDiagnosticsModule ------------------------------------------------------
source("https://raw.githubusercontent.com/OHDSI/CohortDiagnosticsModule/v0.0.7/SettingsFunctions.R")
library(CohortDiagnostics)
cohortDiagnosticsModuleSpecifications <- createCohortDiagnosticsModuleSpecifications(
  runInclusionStatistics = TRUE,
  runIncludedSourceConcepts = TRUE,
  runOrphanConcepts = TRUE,
  runTimeSeries = FALSE,
  runVisitContext = TRUE,
  runBreakdownIndexEvents = TRUE,
  runIncidenceRate = TRUE,
  runCohortRelationship = TRUE,
  runTemporalCohortCharacterization = FALSE,
  minCharacterizationMean = 0.0001,
  temporalCovariateSettings = getDefaultCovariateSettings(),
  incremental = FALSE,
  cohortIds = cohortDefinitionSet$cohortId)


# CharacterizationModule Settings ---------------------------------------------
source("https://raw.githubusercontent.com/OHDSI/CharacterizationModule/v0.3.2/SettingsFunctions.R")
allCohortIdsExceptOutcomes <- cohortDefinitionSet %>%
  filter(!cohortId %in% outcomes$cohortId) %>%
  pull(cohortId)
characterizationModuleSpecifications <- createCharacterizationModuleSpecifications(
  targetIds = allCohortIdsExceptOutcomes,
  outcomeIds = outcomes$cohortId,
  dechallengeStopInterval = 30,
  dechallengeEvaluationWindow = 30,
  timeAtRisk = timeAtRisks,
  minPriorObservation = 365,
  covariateSettings = FeatureExtraction::createDefaultCovariateSettings()
)

# CohortMethodModule -----------------------------------------------------------
source("https://raw.githubusercontent.com/OHDSI/CohortMethodModule/v0.1.0/SettingsFunctions.R")
covariateSettings <- FeatureExtraction::createDefaultCovariateSettings(
  addDescendantsToExclude = TRUE # Keep TRUE because you're excluding concepts
)

getDbCohortMethodDataArgs <- CohortMethod::createGetDbCohortMethodDataArgs(
  restrictToCommonPeriod = TRUE,
  studyStartDate = studyStartDate,
  studyEndDate = studyEndDate,
  maxCohortSize = 0,
  covariateSettings = covariateSettings
)

createPsArgs = CohortMethod::createCreatePsArgs(
  maxCohortSizeForFitting = 250000,
  errorOnHighCorrelation = FALSE,
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
  covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
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
  Strategus::addSharedResources(cohortDefinitionShared) %>%
  Strategus::addSharedResources(negativeControlsShared) %>%
  Strategus::addModuleSpecifications(cohortGeneratorModuleSpecifications) %>%
  Strategus::addModuleSpecifications(cohortDiagnosticsModuleSpecifications) %>%
  Strategus::addModuleSpecifications(characterizationModuleSpecifications) %>%
  Strategus::addModuleSpecifications(cohortMethodModuleSpecifications)

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

ParallelLogger::saveSettingsToJson(analysisSpecifications, file.path(outputFolder, "analysisSpecification.json"))
