# Copyright 2023 Observational Health Data Sciences and Informatics
#
# This file is part of TROY
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
baseUrl <- 'http://34.148.35.102:80/WebAPI'  # Sys.getenv("baseUrl")
referenceRct <- read.csv(file.path('inst/settings', "TargetComparatorList.csv"))
colnames(referenceRct) = c("targetId", "comparatorId", "design", "description", "trialName")
trialNames <- unique(referenceRct$trialName)

# for loop
for (i in 1:length(trialNames)) {
  trial <- unique(trialNames)[i]
  idx <- which(referenceRct[,"trialName"] == trial)
  message(trial)

  # specificationRct <- read.csv(file.path(getwd(), "inst", "csv", paste0(trial, ".csv")))
  specificationRct <- read.csv(file.path('inst/csv/', paste0(trial, ".csv")))
  colnames(specificationRct) <- c("characteristics", "target", "targetSd", "comparator", "comparatorSd", "tag", "isNa", "targetSize", "comparatorSize", "statistics", "analysisId", "conceptIds", "hasValue", "summary", "conceptSetId")
  specificationRct[is.na(specificationRct)]<-""
  specificationRct$conceptSetId <- as.numeric(specificationRct$conceptSetId)

  for (j in 1:nrow(specificationRct)) {
    if(is.na(specificationRct$conceptSetId[j])){
      message(paste0("no concept set for ", specificationRct$characteristics[j]))
      next
    }

  conceptSetId <- ROhdsiWebApi::getConceptSetDefinition(
    conceptSetId = specificationRct$conceptSetId[j],
    baseUrl = baseUrl
  ) %>%
    ROhdsiWebApi::resolveConceptSet(
      baseUrl = baseUrl
    )
  specificationRct$conceptIds[j] <- paste0(conceptSetId, collapse = ";")
  write.csv(specificationRct, paste0("inst/", trial, ".csv"), row.names=FALSE)
  }
}

