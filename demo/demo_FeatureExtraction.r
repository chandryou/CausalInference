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

covAC <- getDbCovariateData(connectionDetails = connectionDetails,
                                    cdmDatabaseSchema = cdmDatabaseSchema,
                                    cohortDatabaseSchema = resultsDatabaseSchema,
                                    cohortTable = "cohort",
                                    cohortId = 13180,
                                    rowIdField = "subject_id",
                                    covariateSettings = covariateSettings)

covAD <- getDbCovariateData(connectionDetails = connectionDetails,
                                      cdmDatabaseSchema = cdmDatabaseSchema,
                                      cohortDatabaseSchema = resultsDatabaseSchema,
                                      cohortTable = "cohort",
                                      cohortId = 14180,
                                      rowIdField = "subject_id",
                                      covariateSettings = covariateSettings)

summary(covAC)
saveRDS(covAC,file.path(dataFolder,"covAC.rds"))
saveCovariateData(covAC, file.path(DataFolder,"covAC"))

summary(covAD)
saveRDS(covAD,file.path(dataFolder,"covAD.rds"))
saveCovariateData(covAD, file.path(DataFolder,"covAD"))

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

