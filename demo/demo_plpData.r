library(DatabaseConnector)
library(FeatureExtraction)

limitCovariatesToPopulation <- function(covariates, rowIds) {
    idx <- !is.na(ffbase::ffmatch(covariates$rowId, rowIds))
    covariates <- covariates[ffbase::ffwhich(idx, idx == TRUE), ]
    return(covariates)
}

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


#################################################################################
###########################TARGET PLP DATA#######################################

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
saveRDS(plpAC,file.path(dataFolder,"plpAC.rds"))

popAC<-PatientLevelPrediction::createStudyPopulation(plpAC, population = NULL, binary = TRUE,outcomeId=4320,
                                                     includeAllOutcomes = T, firstExposureOnly = FALSE, washoutPeriod = 0,
                                                     removeSubjectsWithPriorOutcome = FALSE, priorOutcomeLookback = 99999,
                                                     requireTimeAtRisk = T, minTimeAtRisk = 1, riskWindowStart = 1,
                                                     addExposureDaysToStart = FALSE, riskWindowEnd = 5000,
                                                     addExposureDaysToEnd = F)
saveRDS(popAC,file.path(dataFolder,"popAC.rds"))
#popAC<-readRDS(file.path(dataFolder,"popAC.rds"))

# tidyCovariates<-FeatureExtraction::tidyCovariateData(covariates=plpData$covariates, 
#                                                      covariateRef=plpData$covariateRef,
#                                                      populationSize = nrow(plpData$cohorts),
#                                                      minFraction = 0.001,
#                                                      normalize = TRUE,
#                                                      removeRedundancy = TRUE)

#clone data to prevent accidentally deleting plpData
covariates <- ff::clone(plpAC$covariates)
covariates <- limitCovariatesToPopulation(covariates, ff::as.ff(popAC$rowId))
covariates<-data.frame(covariates)

covariateref.df<-data.frame(plpAC$covariateRef)
subject_id_factor<-as.factor(unique(covariates$rowId))
covariate_id_factor<-as.factor(unique(covariates$covariateId))

#multi-hot vectorization into sparse matrix
covACsparse<-Matrix::sparseMatrix(i=match(covariates$rowId, subject_id_factor), ##remove levels!!!
                                  j=match(covariates$covariateId,covariate_id_factor),
                                  dims=c(max(length(levels(subject_id_factor))),max(length(levels(covariate_id_factor)))),
                                  dimnames=list(subject_id_factor,covariate_id_factor)
)
saveRDS(covACsparse,file.path(dataFolder,"covACsparse.rds"))

covAC<-as.data.frame(as.matrix(covACsparse))

covAC$outcome<-as.logical(popAC$outcomeCount[as.integer(subject_id_factor)])
covAC$treatment <- FALSE

saveRDS(covAC,file.path(dataFolder,"covAC.rds"))

##########################################################################################
###########################COMPARATOR PLP DATA############################################

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
saveRDS(plpAD,file.path(dataFolder,"plpAD.rds"))

popAD<-PatientLevelPrediction::createStudyPopulation(plpAD, population = NULL, binary = TRUE,outcomeId=4320,
                                                          includeAllOutcomes = T, firstExposureOnly = FALSE, washoutPeriod = 0,
                                                          removeSubjectsWithPriorOutcome = FALSE, priorOutcomeLookback = 99999,
                                                          requireTimeAtRisk = T, minTimeAtRisk = 1, riskWindowStart = 1,
                                                          addExposureDaysToStart = FALSE, riskWindowEnd = 5000,
                                                          addExposureDaysToEnd = F)

saveRDS(popAD,file.path(dataFolder,"popAD.rds"))
#saveCovariateData(popAD, file.path(DataFolder,"popAD"))

#clone data to prevent accidentally deleting plpData
covariates <- ff::clone(plpAD$covariates)
covariates <- limitCovariatesToPopulation(covariates, ff::as.ff(popAD$rowId))
covariates<-data.frame(covariates)

covariateref.df<-data.frame(plpAD$covariateRef)
subject_id_factor<-as.factor(unique(covariates$rowId))
covariate_id_factor<-as.factor(unique(covariates$covariateId))

#multi-hot vectorization into sparse matrix
covADsparse<-Matrix::sparseMatrix(i=match(covariates$rowId, subject_id_factor), ##remove levels!!!
                                  j=match(covariates$covariateId,covariate_id_factor),
                                  dims=c(max(length(levels(subject_id_factor))),max(length(levels(covariate_id_factor)))),
                                  dimnames=list(subject_id_factor,covariate_id_factor)
)
saveRDS(covADsparse,file.path(dataFolder,"covADsparse.rds"))

covAD<-as.data.frame(as.matrix(covADsparse))
covAD$outcome<-as.logical(popAD$outcomeCount[as.integer(subject_id_factor)])
covAD$treatment <- FALSE
saveRDS(covAD,file.path(dataFolder,"covAD.rds"))

###################################################
####binding treatment and comparator group ########

#remvove variants which is not in both cohort 
covAC<-covAC[,colnames(covAC) %in% colnames(covAD)]
covAD<-covAD[,colnames(covAD) %in% colnames(covAC)]

#reorder the columns of comparator group
covAD<-covAD[,match(colnames(covAC),colnames(covAD))]

cov.df<-rbind(covAC,covAD)
saveRDS(cov.df,file.path(dataFolder,"cov_df.rds"))
