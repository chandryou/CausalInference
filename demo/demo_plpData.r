library(DatabaseConnector)
library(FeatureExtraction)

dataFolder<-"/Users/chan/data/test"
cdmDatabaseSchema<-"NHIS_NSC.dbo"
resultsDatabaseSchema<-"NHIS_NSC.dbo"
htn_med<-readRDS("/Users/chan/git/ohdsi/StudyProtocolSandbox/HypertensionCombination/inst/settings/htn_med_list")

connectionDetails<-readRDS(file.path(dataFolder,"connectionDetails.rds"))
# covariateSettings<-createCovariateSettings(useDemographicsGender = FALSE,
#                                   useDemographicsAge = FALSE, useDemographicsAgeGroup = FALSE,
#                                   useDemographicsRace = FALSE, useDemographicsEthnicity = FALSE,
#                                   useDemographicsIndexYear = FALSE, useDemographicsIndexMonth = FALSE,
#                                   useDemographicsPriorObservationTime = FALSE,
#                                   useDemographicsPostObservationTime = FALSE,
#                                   useDemographicsTimeInCohort = FALSE,
#                                   useDemographicsIndexYearMonth = FALSE,
#                                   useConditionOccurrenceAnyTimePrior = FALSE,
#                                   useConditionOccurrenceLongTerm = FALSE,
#                                   useConditionOccurrenceMediumTerm = FALSE,
#                                   useConditionOccurrenceShortTerm = FALSE,
#                                   useConditionOccurrenceInpatientAnyTimePrior = FALSE,
#                                   useConditionOccurrenceInpatientLongTerm = FALSE,
#                                   useConditionOccurrenceInpatientMediumTerm = FALSE,
#                                   useConditionOccurrenceInpatientShortTerm = FALSE,
#                                   useConditionEraAnyTimePrior = FALSE, useConditionEraLongTerm = FALSE,
#                                   useConditionEraMediumTerm = FALSE, useConditionEraShortTerm = FALSE,
#                                   useConditionEraOverlapping = FALSE, useConditionEraStartLongTerm = FALSE,
#                                   useConditionEraStartMediumTerm = FALSE,
#                                   useConditionEraStartShortTerm = FALSE,
#                                   useConditionGroupEraAnyTimePrior = FALSE,
#                                   useConditionGroupEraLongTerm = FALSE,
#                                   useConditionGroupEraMediumTerm = FALSE,
#                                   useConditionGroupEraShortTerm = FALSE,
#                                   useConditionGroupEraOverlapping = FALSE,
#                                   useConditionGroupEraStartLongTerm = FALSE,
#                                   useConditionGroupEraStartMediumTerm = FALSE,
#                                   useConditionGroupEraStartShortTerm = FALSE,
#                                   useDrugExposureAnyTimePrior = FALSE, useDrugExposureLongTerm = FALSE,
#                                   useDrugExposureMediumTerm = FALSE, useDrugExposureShortTerm = FALSE,
#                                   useDrugEraAnyTimePrior = FALSE, useDrugEraLongTerm = FALSE,
#                                   useDrugEraMediumTerm = FALSE, useDrugEraShortTerm = FALSE,
#                                   useDrugEraOverlapping = FALSE, useDrugEraStartLongTerm = FALSE,
#                                   useDrugEraStartMediumTerm = FALSE, useDrugEraStartShortTerm = FALSE,
#                                   useDrugGroupEraAnyTimePrior = FALSE, useDrugGroupEraLongTerm = FALSE,
#                                   useDrugGroupEraMediumTerm = FALSE, useDrugGroupEraShortTerm = FALSE,
#                                   useDrugGroupEraOverlapping = FALSE, useDrugGroupEraStartLongTerm = FALSE,
#                                   useDrugGroupEraStartMediumTerm = FALSE,
#                                   useDrugGroupEraStartShortTerm = FALSE,
#                                   useProcedureOccurrenceAnyTimePrior = FALSE,
#                                   useProcedureOccurrenceLongTerm = FALSE,
#                                   useProcedureOccurrenceMediumTerm = FALSE,
#                                   useProcedureOccurrenceShortTerm = FALSE,
#                                   useDeviceExposureAnyTimePrior = FALSE, useDeviceExposureLongTerm = FALSE,
#                                   useDeviceExposureMediumTerm = FALSE, useDeviceExposureShortTerm = FALSE,
#                                   useMeasurementAnyTimePrior = FALSE, useMeasurementLongTerm = FALSE,
#                                   useMeasurementMediumTerm = FALSE, useMeasurementShortTerm = FALSE,
#                                   useMeasurementValueAnyTimePrior = FALSE,
#                                   useMeasurementValueLongTerm = FALSE,
#                                   useMeasurementValueMediumTerm = FALSE,
#                                   useMeasurementValueShortTerm = FALSE,
#                                   useMeasurementRangeGroupAnyTimePrior = FALSE,
#                                   useMeasurementRangeGroupLongTerm = FALSE,
#                                   useMeasurementRangeGroupMediumTerm = FALSE,
#                                   useMeasurementRangeGroupShortTerm = FALSE,
#                                   useObservationAnyTimePrior = FALSE, useObservationLongTerm = FALSE,
#                                   useObservationMediumTerm = FALSE, useObservationShortTerm = FALSE,
#                                   useCharlsonIndex = FALSE, useDcsi = FALSE, useChads2 = FALSE,
#                                   useChads2Vasc = FALSE, useDistinctConditionCountLongTerm = FALSE,
#                                   useDistinctConditionCountMediumTerm = FALSE,
#                                   useDistinctConditionCountShortTerm = FALSE,
#                                   useDistinctIngredientCountLongTerm = FALSE,
#                                   useDistinctIngredientCountMediumTerm = FALSE,
#                                   useDistinctIngredientCountShortTerm = FALSE,
#                                   useDistinctProcedureCountLongTerm = FALSE,
#                                   useDistinctProcedureCountMediumTerm = FALSE,
#                                   useDistinctProcedureCountShortTerm = FALSE,
#                                   useDistinctMeasurementCountLongTerm = FALSE,
#                                   useDistinctMeasurementCountMediumTerm = FALSE,
#                                   useDistinctMeasurementCountShortTerm = FALSE,
#                                   useVisitCountLongTerm = FALSE, useVisitCountMediumTerm = FALSE,
#                                   useVisitCountShortTerm = FALSE, longTermStartDays = -365,
#                                   mediumTermStartDays = -180, shortTermStartDays = -30, endDays = 0,
#                                   includedCovariateConceptIds = c(), addDescendantsToInclude = FALSE,
#                                   excludedCovariateConceptIds = c(), addDescendantsToExclude = FALSE,
#                                   includedCovariateIds = c())

covariateSettings<-createDefaultCovariateSettings(excludedCovariateConceptIds=htn_med)

plpAC<-PatientLevelPrediction::getPlpData(connectionDetails, 
                                            cdmDatabaseSchema,
                                            oracleTempSchema = NULL, 
                                            cohortId=13180, 
                                            outcomeIds=4320,#list(3,1430),
                                            studyStartDate = "19900101", 
                                            studyEndDate = "20171231",
                                            cohortDatabaseSchema = cdmDatabaseSchema, 
                                            cohortTable = "cohort",
                                            outcomeDatabaseSchema = cdmDatabaseSchema, 
                                            outcomeTable = "cohort",
                                            cdmVersion = "5", excludeDrugsFromCovariates = F,
                                            firstExposureOnly = FALSE, 
                                            washoutPeriod = 0, 
                                            sampleSize = NULL,
                                            covariateSettings)

plpAD<-PatientLevelPrediction::getPlpData(connectionDetails, 
                                            cdmDatabaseSchema,
                                            oracleTempSchema = NULL, 
                                            cohortId=14180, 
                                            outcomeIds=4320,#list(3,1430),
                                            studyStartDate = "19900101", 
                                            studyEndDate = "20171231",
                                            cohortDatabaseSchema = cdmDatabaseSchema, 
                                            cohortTable = "cohort",
                                            outcomeDatabaseSchema = cdmDatabaseSchema, 
                                            outcomeTable = "cohort",
                                            cdmVersion = "5", excludeDrugsFromCovariates = F,
                                            firstExposureOnly = FALSE, 
                                            washoutPeriod = 0, 
                                            sampleSize = NULL,
                                            covariateSettings)

popAC<-PatientLevelPrediction::createStudyPopulation(plpAC, population = NULL, binary = TRUE,outcomeId=4320,
                                                          includeAllOutcomes = T, firstExposureOnly = FALSE, washoutPeriod = 0,
                                                          removeSubjectsWithPriorOutcome = FALSE, priorOutcomeLookback = 99999,
                                                          requireTimeAtRisk = T, minTimeAtRisk = 1, riskWindowStart = 1,
                                                          addExposureDaysToStart = FALSE, riskWindowEnd = 5000,
                                                          addExposureDaysToEnd = F)

popAD<-PatientLevelPrediction::createStudyPopulation(plpAD, population = NULL, binary = TRUE,outcomeId=4320,
                                                          includeAllOutcomes = T, firstExposureOnly = FALSE, washoutPeriod = 0,
                                                          removeSubjectsWithPriorOutcome = FALSE, priorOutcomeLookback = 99999,
                                                          requireTimeAtRisk = T, minTimeAtRisk = 1, riskWindowStart = 1,
                                                          addExposureDaysToStart = FALSE, riskWindowEnd = 5000,
                                                          addExposureDaysToEnd = F)


summary(popAC)
saveRDS(popAC,file.path(dataFolder,"popAC.rds"))
saveCovariateData(covAC, file.path(DataFolder,"covAC"))

summary(popAD)
saveRDS(popAD,file.path(dataFolder,"popAD.rds"))
saveCovariateData(popAD, file.path(DataFolder,"popAD"))

# covACtidy <- tidyCovariateData(covariateData=covAC, 
#                                minFraction = 0.001,
#                                normalize = TRUE,
#                                removeRedundancy = TRUE)
# 
# covADtidy <- tidyCovariateData(covariateData=covAC, 
#                                minFraction = 0.001,
#                                normalize = TRUE,
#                                removeRedundancy = TRUE)
# 
# summary(covACtidy)
# summary(covADtidy)

covariates<-data.frame(covAC$covariates)
covariateref.df<-data.frame(covAC$covariateRef)
subject_id_factor<-as.factor(unique(covariates$rowId))
covariate_id_factor<-as.factor(unique(covariates$covariateId))

#multi-hot vectorization into sparse matrix
covACsparse<-Matrix::sparseMatrix(i=match(covariates$rowId, levels(subject_id_factor)),
                             j=match(covariates$covariateId,levels(covariate_id_factor)),
                             dims=c(max(length(levels(subject_id_factor))),max(length(levels(covariate_id_factor)))),
                             dimnames=list(levels(subject_id_factor),levels(covariate_id_factor))
)
saveRDS(covACsparse,file.path(dataFolder,"covACsparse.rds"))

covariates<-data.frame(covAD$covariates)
covariateref.df<-data.frame(covAD$covariateRef)
subject_id_factor<-as.factor(unique(covariates$rowId))
covariate_id_factor<-as.factor(unique(covariates$covariateId))

covADsparse<-Matrix::sparseMatrix(i=match(covariates$rowId, levels(subject_id_factor)),
                                  j=match(covariates$covariateId,levels(covariate_id_factor)),
                                  dims=c(max(length(levels(subject_id_factor))),max(length(levels(covariate_id_factor)))),
                                  dimnames=list(levels(subject_id_factor),levels(covariate_id_factor))
)
saveRDS(covADsparse,file.path(dataFolder,"covADsparse.rds"))

