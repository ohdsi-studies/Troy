# install the network package
#install.packages('remotes')
#remotes::install_github("OHDSI/Strategus", ref="results-upload")
library(Strategus)

##=========== START OF INPUTS ==========
connectionDetailsReference <- "YUHS"
workDatabaseSchema <- 'jaehyeongcho'
cdmDatabaseSchema <- 'YUHS_CDM'
outputLocation <- '~/output' #DO NOT TOUCH
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
  fileName = "~/git/Troy/inst/analysisSpecification.json"
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

source('~/git/Troy/R/troyFunction.R')

troyFunction(cohortDatabaseSchema=executionSettings$workDatabaseSchema,
            cohortTable=executionSettings$cohortTableNames$cohortTable,
            connectionDetails=connectionDetails,
            outputFolder=file.path(outputLocation, connectionDetailsReference, "strategusOutput"))
