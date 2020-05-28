setwd("C:/Diana/Desktop/Bayes_project")
#due to the legacy reasons, we can't provide the csv files
#containing the clinical data about the patients.
#the nlst_path_prsn_20180705 file contains variables, 
#that can be subdivided into 10 blocks that contain information about patients
#some the variables, containing in each block are listed below
#for more detailed information about variables, please, consult the 
#open access data dictionaries that you can find in the github or
#by link URL  
#1. STUDY 
#----- patient personal ID 
#----- information about the study branch person was enrolled on
#----- info if patient is elegible of for the study
#2. DEMOGRAPHIC
#----- age, education, gender, height
#----- race, marital status, weight
#3. SMOKING
#----- if patient is a current smoker
#----- years of smoking
#----- measure "pack years" (Total Years Smoked  x Cigarettes Per Day / 20).
#4. SCREENING
#----- information related to computer tomography (CT) screening for lung cancer diagnostics
#5. FOLLOW-UP/PROCEDURES
#----- if patient had invasive procedures (biopsy) related to positive CT screen
#----- if patient had other procedures related to lung cancer
#6. LUNG CANCER
#----- days form beginning of the study until first diagnosis of lung cancer
#----- information about morphology, topology and stage of cancer
#----- if patient received treatment
#7. LAST CONTACT
#----- Days from randomization to date last known alive
#----- status of participant at the end of the study (active/deceased/withdrawn or lost contact)
#8. DEATH
#----- vital status at the end of the study
#----- if the death was related to lung cancer or not
#9. ENDPOINT VERIFICATION PROCESS
#----- technical info about death verification process
#10. WORK HISTORY
#----- was patient doing work, potentially affecting lungs
#----- (asbestos,coal,baking,cotton,chemical,painting)
#----- was patient wearing protection mask/equipment?
#11. DISEASE HISTORY
#----- if person ever head lung-related deseases
#----- (astma,chronic bronchis,tyberculosis,etc...)
#12. PERSONAL CANCER HISTORY
#----- if the person ever had any type of cancer (not only lung)
#13. FAMILY HISTORY OF LUNG CANCER 
#----- if patient's family members had lung cancer
#14. ALCOHOL
#----- info about patient's alcohol consumption (Missing for Adenocarcinoma patients)
#15. CANCERS OF ANY SITE
#----- information about cancer(lung and not lung) of a patient 
#16. ABNORMALITY SUMMARY 
#----- Was a non-calcified nodule/mass >=4mm found on any screen?
#17. PROGRESSION
#----- information of lung cancer progresson since the beginning of the study
#----- information about lung cancer  meastasis
#18. IDENTIFIERS 
#----- tumor ID
#19. TUMOR INFORMATION
#----- diameter,side, topography, WHO class
#----- is likely to have metastasis


data<-read.csv("nlst_path_prsn_20180705.csv")
#data_mod<-read.csv("data_modified.csv")
#taking only the patients that have biopsy tissue image
data <- data[ which(data$has_tissue_image==1), ]


#taking only the patients which have desease type "Adenocarcinoma", because the CNN  for processing the images
#works only with Adenocarcinoma tissue images.
#8140 -- Lung cancer type according to ICD-O-3 morphology (corresonds to "Adenocarcinoma")
data <- data[ which(data$de_type==8140), ] 

#removing patients that are ineligible for the study
#2="Non-smoker or quit > 15 years"
#7="Previous Lung Cancer"
#8="Portion of lung removed"
#9="Cancer within past 5 years"
data<-data[is.na(data$ineligible),]

data$finaldeathlc[is.na(data$finaldeathlc)]<-2
#removing those who died not because of lung cancer
data<-data[which(data$finaldeathlc!=0),]



#write.csv(data,"data_cleaned_diseases.csv")

#--------------------------------------------------------------------------
#Cleaning----1 Study---section done
#--------------------------------------------------------------------------
data[c("cen","dataset_version","elig","has_tissue_image","ineligible","path_selected","study")]<-list(NULL)
#most of the variables are not useful for further analysis after previously
#made filtering
#the variables that are left:
#pid -- person identifier 
#rndgroup -- study arm, Computer Tomography of X-ray

#--------------------------------------------------------------------------
#Cleaning----2 Demographic ---section done
#--------------------------------------------------------------------------
data[c("ethnic")]<-list(NULL)#no variability

data$weight<-data$weight*0.454#pounds to kg
data$height<-data$height*2.54#inches to cm

#--------------------------------------------------------------------------
#Cleaning----3 Smoking ---section done
#--------------------------------------------------------------------------
#pkyr -- Pack years, calculated as: (Total Years Smoked  x Cigarettes Per Day / 20)
#contains information about smoking years and sigarets per days
#so, removing corresponding variables
data[c("smokeyr","smokeday")]<-list(NULL)#no variability

#--------------------------------------------------------------------------
#Cleaning----4 Screening---section 
#--------------------------------------------------------------------------
data[c("scr_lat0","scr_lat1","scr_lat2")]<-list(NULL)
data[c("scr_group","sct_image_has","sct_image_years")]<-list(NULL)
#--------------------------------------------------------------------------
#Cleaning----5 Follow-up/Procedures---section done
#--------------------------------------------------------------------------
#all the considered patients had biopsy and invasive procedures
#because we have chosen the patients that have tissue image(images)
#we are not interest in medical complications during procedures
#all the patients had diagnostic procedures after a posotive screening

data[c("medcomp0","medcomp1","medcomp2","medcomplc","mra_stat0","mra_stat1","mra_stat2")]<-list(NULL)
data[c("no_proc_reas0","no_proc_reas1","no_proc_reas2","proc0","proc1","proc2","proclc")]<-list(NULL)

#--------------------------------------------------------------------------
#Cleaning---- 6 Lung Cancer ---section
#--------------------------------------------------------------------------
data[c("canc_free_days","canc_rpt_link","canc_rpt_source","cancyr")]<-list(NULL)
#canc_free_days is not recommended to use in survival analysis
#because it is not precise

#removing the staging information except of the aggregation parameter "stage"
data[c("clinical_m","clinical_m_7thed","clinical_n","clinical_n_7thed","clinical_t","clinical_stag","clinical_t_7thed","conflc")]<-list(NULL)
data[c("de_grade","de_type","de_stag","path_m","path_m_7thed","path_n","path_n_7thed","path_stag","path_t","path_t_7thed","treatlc")]<-list(NULL)
#we decided to not include "cancer location" boolean variables into analysis
data[c("loccar","loclhil","loclin","locllow","loclmsb","loclup","locmed","locoth","locrhil","locrlow","locrmid","locrmsb","locrup","locunk")]<-list(NULL)

#remaining variables:
#de_stag_7thed -- staging (AJCC7 -- more modern) combining information obtained by non-invasive methods and invasive
#can_scr -- Indicates whether the cancer followed a positive, negative, or missed screen, or whether it occurred after the screening years.
#candx_days -- Days from randomization to first diagnosis of lung cancer 

#--------------------------------------------------------------------------
#FOR CALCULATING countdown(diagnosis) date FOR SURVIVAL ANALYSIS
#--------------------------------------------------------------------------
#subsetting those whose diagnosis don't depend on the screeing of 3 first years of study 
#and defining their countdown date is equal to diagnosis date
#because we don't have better more precise information
subset_data4 <- data[which(data$can_scr!=1),]
subset_data4$diag_date <- subset_data4$candx_days


#for the others the countdown date will equal to the date of first positive
#screening (CT or X-ray)


subset_data1 <- data[which(data$can_scr==1),]

subset_data11<- subset_data1[which(subset_data1$scr_res0==4),]
subset_data11$diag_date <- subset_data11$scr_days0

subset_data12<- subset_data1[which(subset_data1$scr_res0!=4),]
subset_data211<-subset_data12[which(subset_data12$scr_res1==5),]
subset_data212<-subset_data12[which(subset_data12$scr_res1==6),]

subset_data21<-rbind(subset_data211,subset_data212)
subset_data21$diag_date<-subset_data21$scr_days1
rm(subset_data211,subset_data212)

subset_data211<-subset_data12[which(subset_data12$scr_res1==1),]
subset_data212<-subset_data12[which(subset_data12$scr_res1==2),]
subset_data213<-subset_data12[which(subset_data12$scr_res1==3),]
subset_data31<-rbind(subset_data211,subset_data212,subset_data213)
rm(subset_data211,subset_data212,subset_data213)
subset_data31$diag_date<-subset_data31$scr_days2



data<-rbind(subset_data4,subset_data11,subset_data21,subset_data31)
rm(subset_data4,subset_data11,subset_data21,subset_data31,subset_data1,subset_data12)

data[c("biop0","biop1","biop2","bioplc","invas0","invas1","invas2","invaslc")]<-list(NULL)
data[c("scr_res0","scr_res1","scr_res2","scr_iso0","scr_iso1","scr_iso2","can_scr","scr_days0","scr_days1","scr_days2","candx_days")]<-list(NULL)

#--------------------------------------------------------------------------
#Cleaning----7 Last Contact ---section done
#--------------------------------------------------------------------------
data[c("contactstatus","wdlost")]<-list(NULL)
#the remaining variable is:
#fup_days -- follow up days (Days from randomization to date last known alive)

#--------------------------------------------------------------------------
#Cleaning----8 Death---section done
#--------------------------------------------------------------------------
data[c("dcfdeathlc","dcficd","death_days","deathcutoff","deathstat","finaldeathlc","hasdcf","ndicd")]<-list(NULL)
#none of the variables is useful for the further analysis after previously
#made filtering

#--------------------------------------------------------------------------
#Cleaning----9 Endpoint Verification Process ---section done
#--------------------------------------------------------------------------
data[c("evp_revr","evpcert","evpdeath","evpdirect","evpincomplete","evpsel","evpsent")]<-list(NULL)
#none of the variables is useful for the further analysis after previously
#made filtering

#--------------------------------------------------------------------------
#Cleaning----10 Work history ---section 
#--------------------------------------------------------------------------
#we decided to discard the information about wearing a respirator,
#because even from the persons who did a work, dangerous for lungs
#only very few wore respirators
data[c("resasbe","resbaki","resbutc","reschem","rescoal","rescott","resfarm","resfire","resflou","resfoun","reshard","respain","ressand","resweld")]<-list(NULL)
data[c("wrkasbe","wrkbaki","wrkbutc","wrkchem","wrkcoal","wrkcott","wrkfarm","wrkfire","wrkflou","wrkfoun","wrkhard","wrkpain","wrksand","wrkweld")]<-list(NULL)

#creating new variable "dangerous_for_lung_work_years"
subset_data <- data[c("yrsasbe","yrsbaki","yrsbutc","yrschem","yrscoal","yrscott","yrsfarm","yrsfire","yrsflou","yrsfoun","yrshard","yrspain","yrssand","yrsweld")]
data[c("yrsasbe","yrsbaki","yrsbutc","yrschem","yrscoal","yrscott","yrsfarm","yrsfire","yrsflou","yrsfoun","yrshard","yrspain","yrssand","yrsweld")]<-list(NULL)

subset_data[is.na(subset_data)] <- 0
subset_data$lung_dang_work <- rowSums(subset_data)
#manual fixes, because of obvious overlapping occupation years
subset_data$lung_dang_work[7]<-34
subset_data$lung_dang_work[40]<-45
subset_data$lung_dang_work[133]<-30
subset_data$lung_dang_work[154]<-20
lung_dang_work<-subset_data$lung_dang_work

data<-cbind(data,lung_dang_work)
rm(subset_data,lung_dang_work)
#--------------------------------------------------------------------------
#Cleaning----12 Desease history ---section 
#--------------------------------------------------------------------------
data[c("ageasbe","ageadas","agechas","agebron","agechro","agecopd","agediab","ageemph","agefibr","agehear","agepneu","agesarc","agesili","agetube","agehype","agestro")]<-list(NULL)
data[c("diagfibr","diagsili")]<-list(NULL)#no cases of illness reported
data[c("diagchas","diagbron","diagchro","diagcopd","diagdiab","diagemph","diagfibr","diaghear","diagpneu","diagsarc","diagsili","diagtube","diaghype","diagstro")]<-list(NULL)
#--------------------------------------------------------------------------
#Cleaning----12 Personal cancer history ---section done
#--------------------------------------------------------------------------
#we are not interested in other types of cancer of the patient,
#because we concider adenocarcinoma, not metastasis cancer
#and none of the patient had lung cancer diagnosis before the trial(canclung,agelung)
data[c("ageblad","agebrea","agecerv","agecolo","ageesop","agekidn","agelary","agelung","agenasa","ageoral")]<-list(NULL)
data[c("agepanc","agephar","agestom","agethyr","agetran","cancblad","cancbrea","canccerv","canccolo")]<-list(NULL)
data[c("cancesop","canckidn","canclary","canclung","cancnasa","cancoral","cancpanc","cancphar","cancstom","cancthyr","canctran")]<-list(NULL)

#--------------------------------------------------------------------------
#Cleaning----13 Family history of lung cancer done
#--------------------------------------------------------------------------
#we decided to don't include those variables to analysis
data[c("famfather","fammother","fambrother","famsister","famchild")]<-list(NULL)

#--------------------------------------------------------------------
#Cleaning----14 Alcohol---section done
#--------------------------------------------------------------------
#8 variables contain NA, they are not available for these patients
data[c("acrin_alc_ever","acrin_alc_curr","acrin_lastdrink","acrin_drinkyrs_form","acrin_drinknum_form","acrin_drinkyrs_curr","acrin_drinknum_curr","acrin_drink24hr")]<-list(NULL)

#2 variables are left:
#lss_alcohol_freq -- How often do you have a drink containing alcohol?
#lss_alcohol_num -- Number of alcoholic drinks on typical day when drinking

#--------------------------------------------------------------------
#Cleaning----15 Cancers of any site ---section done
#--------------------------------------------------------------------
#we are not including other cancers in analysis, and the info
#about lung cancer we have in 6th section
data[c("confirmed_candxdays1","confirmed_candxdays2","confirmed_candxdays3","confirmed_candxdays4")]<-list(NULL)
data[c("confirmed_conforder1","confirmed_conforder2","confirmed_conforder3","confirmed_conforder4")]<-list(NULL)
data[c("confirmed_icd_behav1","confirmed_icd_behav2","confirmed_icd_behav3","confirmed_icd_behav4")]<-list(NULL)
data[c("confirmed_icd_grade1","confirmed_icd_grade2","confirmed_icd_grade3","confirmed_icd_grade4")]<-list(NULL)
data[c("confirmed_icd_morph1","confirmed_icd_morph2","confirmed_icd_morph3","confirmed_icd_morph4")]<-list(NULL)
data[c("confirmed_icd_topog1","confirmed_icd_topog2","confirmed_icd_topog3","confirmed_icd_topog4")]<-list(NULL)
data[c("confirmed_seer1","confirmed_seer2","confirmed_seer3","confirmed_seer4")]<-list(NULL)
data[c("confirmed_seercat1","confirmed_seercat2","confirmed_seercat3","confirmed_seercat4","num_confirmed")]<-list(NULL)
#--------------------------------------------------------------------
#Cleaning----16 Abnormality summary ---section done
#--------------------------------------------------------------------
data[c("anyscr_has_nodule")]<-list(NULL)

#--------------------------------------------------------------------
#Cleaning---- 17 Progression ---section done
#--------------------------------------------------------------------

data[c("last_progfree_days","prog_days_1st","prog_days_2nd","prog_days_3rd","prog_days_4th","prog_days_5th")]<-list(NULL)
data[c("progsite_adrenal_1st","progsite_adrenal_days","progsite_adrenal_ever","progsite_adrenal_num")]<-list(NULL)
data[c("progsite_bone_1st","progsite_bone_days","progsite_bone_ever","progsite_bone_num","progsite_brain_1st","progsite_brain_days","progsite_brain_ever","progsite_brain_num")]<-list(NULL)
data[c("progsite_liver_1st","progsite_liver_days","progsite_liver_ever","progsite_liver_num")]<-list(NULL)
data[c("progsite_lymph_n1_1st","progsite_lymph_n1_days","progsite_lymph_n1_ever","progsite_lymph_n1_num")]<-list(NULL)
data[c("progsite_lymph_n2_1st","progsite_lymph_n2_days","progsite_lymph_n2_ever","progsite_lymph_n2_num")]<-list(NULL)
data[c("progsite_lymph_n3_1st","progsite_lymph_n3_days","progsite_lymph_n3_ever","progsite_lymph_n3_num")]<-list(NULL)
data[c("progsite_mediastinum_1st","progsite_mediastinum_days","progsite_mediastinum_ever","progsite_mediastinum_num")]<-list(NULL)
data[c("progsite_orig_lung_1st","progsite_orig_lung_days","progsite_orig_lung_ever","progsite_orig_lung_num")]<-list(NULL)
data[c("progsite_other_1st","progsite_other_days","progsite_other_ever","progsite_other_lung_1st","progsite_other_lung_days","progsite_other_lung_ever","progsite_other_lung_num","progsite_other_num")]<-list(NULL)
data[c("progsite_pleura_1st","progsite_pleura_days","progsite_pleura_ever","progsite_pleura_num","progsite_skin_1st","progsite_skin_days","progsite_skin_ever","progsite_skin_num")]<-list(NULL)
data[c("progsite_unk_1st","progsite_unk_days","progsite_unk_ever","progsite_unk_num")]<-list(NULL)

data[c("progressed_ever")]<-list(NULL)
#1 variable is left:
#progression_num, replacing NAs with 0
data$progression_num[is.na(data$progression_num)]<-0

#--------------------------------------------------------------------------
#Cleaning----18 Identifiers ---section done
#--------------------------------------------------------------------------
data[c("tumor_id")]<-list(NULL)


#--------------------------------------------------------------------------
#Cleaning----19 Tumor information ---section done 
#--------------------------------------------------------------------------
data[c("case_highest_inv_grade","case_likely_metastases","case_longest_inv_diam","case_reported_inv_diam","case_side","case_topography","case_topography_descrip","case_WHO_class")]<-list(NULL)


write.csv(data,"data_cleaned.csv")
rm(data)
