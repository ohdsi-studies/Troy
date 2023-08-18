# install the network package
install.packages('remotes')
remotes::install_github("OHDSI/Strategus", ref="results-upload")
library(Strategus)

##=========== START OF INPUTS ==========
connectionDetailsReference <- "Jmdc"
workDatabaseSchema <- 'scratch_asena5'
cdmDatabaseSchema <- 'cdm_jmdc_v2325'
outputLocation <- '~/output/troy'
minCellCount <- 5
cohortTableName <- "troy"

# the keyring entry should correspond to what you selected in KeyringSetup.R
connectionDetails = DatabaseConnector::createConnectionDetails(
  dbms = keyring::key_get("dbms", keyring = "troy"),
  connectionString = keyring::key_get("connectionString", keyring = "troy"),
  user = keyring::key_get("username", keyring = "troy"),
  password = keyring::key_get("password", keyring = "troy")
)

##=========== END OF INPUTS ==========
##################################
# DO NOT MODIFY BELOW THIS POINT
##################################
analysisSpecifications <- ParallelLogger::loadSettingsFromJson(
  fileName = "inst/analysisSpecification.json"
)

storeConnectionDetails(
  connectionDetails = connectionDetails,
  connectionDetailsReference = connectionDetailsReference,
  keyringName = "troy"
)

executionSettings <- createCdmExecutionSettings(
  connectionDetailsReference = connectionDetailsReference,
  workDatabaseSchema = workDatabaseSchema,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortTableNames = CohortGenerator::getCohortTableNames(cohortTable = cohortTableName),
  workFolder = file.path(outputLocation, connectionDetailsReference, "strategusWork"),
  resultsFolder = file.path(outputLocation, connectionDetailsReference, "strategusOutput"),
  minCellCount = minCellCount
)

# Note: this environmental variable should be set once for each compute node
Sys.setenv("INSTANTIATED_MODULES_FOLDER" = file.path(outputLocation, "StrategusInstantiatedModules"))

execute(
  analysisSpecifications = analysisSpecifications,
  executionSettings = executionSettings,
  executionScriptFolder = file.path(outputLocation, connectionDetailsReference, "strategusExecution"),
  keyringName = "troy"
)
