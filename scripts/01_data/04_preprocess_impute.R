# =============================================================================
# scripts/01_data/04_preprocess_impute.R
# Purpose : Build the analysis sample.
#   1. Recode WVS non-response codes to NA.
#   2. Drop questions answered by < 30% of respondents.
#   3. Drop respondents answering < 70% of the remaining questions.
#   4. Export the raw geographic variables of the kept respondents.
#   5. Impute residual missing values with a random forest trained on the
#      complete respondents (ranger, 500 trees, seed 123).
#   6. Label the imputed data and split the value block (clustering input)
#      from the political and sociodemographic blocks.
# Thesis  : Section 2.2, Appendix A.2
# Input   : data/raw/WVS_Wave_7_United_States_Csv_v5.1.csv (via 01_import_data.R)
# Output  : data/processed/geo_raw.rds
#           data/processed/data_imputed.rds            (all kept variables)
#           data/processed/data_imputed_selected.rds   (value block, clustering input)
#           data/processed/imputed_labeled.rds         (all variables, labelled)
#           data/processed/imputed_selected_labeled.rds (value block, labelled)
#           data/processed/data_class.rds              (held-out blocks, labelled)
# =============================================================================

source(here::here("R", "00_setup.R"))
source(here("scripts", "01_data", "01_import_data.R"))

# 1 & 2 answer rate ####
# Removing any question that is answered by less than 30% of the respondents
# Removing respondents who answered less than 70% of questions
temp_data <- data %>% 
  mutate( across( everything(), ~ replace(.x,.x %in% missing_codes, NA )))

threshold_q <- 0.3 
threshold_i <- 0.7

q_response_rate <- colMeans(!is.na(temp_data))
i_response_rate <- rowMeans((!is.na(temp_data)))


vars_to_keep <- names(q_response_rate[q_response_rate >= threshold_q])
i_to_keep <-  which(i_response_rate >= threshold_i)

# Filtering 
w_data <- data %>% 
  select(any_of(vars_to_keep)) %>% 
  slice(i_to_keep)
nrow(w_data)

# ── GEO EXPORT ────────────────────────────────────────────────────────────────
# Save geographic variables separately, RAW (untouched by missing-code recoding
# or RF imputation), restricted to the same 2579 respondents kept above.
# All NA/sentinel handling (-5, -999) and labeling is done later in the geo script.

geo_vars <- c("D_INTERVIEW", "G_TOWNSIZE", "G_TOWNSIZE2", "H_SETTLEMENT",
              "H_URBRURAL", "N_REGION_ISO", "N_TOWN",
              "O1_LONGITUDE", "O2_LATITUDE")

stopifnot(all(geo_vars %in% names(data)))   # fail fast: column names must match the codebook

geo_raw <- data %>%
  slice(i_to_keep) %>%
  select(any_of(geo_vars))

stopifnot(nrow(geo_raw) == nrow(w_data))                       # same rows as w_data
stopifnot(nrow(geo_raw) == length(unique(geo_raw$D_INTERVIEW)))# no duplicated IDs

cat("geo_raw: rows saved =", nrow(geo_raw), "\n")

saveRDS(geo_raw, file.path(DIR_DATA, "geo_raw.rds"))

# 5. Imputation ####
# Imputing missing values
# Random forests impute the missing values, using the complete rows as the training sample

#1. defining the train sample 
temp_data2 <- w_data %>%
  mutate(across(everything(), ~ replace(.x, .x %in% missing_codes, NA)))

# changing variables into factor
temp_data2 <- temp_data2 %>%
 mutate(across(everything(), as.factor))


complete_rows <- complete.cases(temp_data2)
T_data <- temp_data2[complete_rows, , drop = FALSE]

cat("Number of complete respondents in T:", nrow(T_data), "\n")
cat("Number of incomplete respondents:", sum(!complete_rows), "\n")


#1. Initialization of Missing Values
# Needed so predictors are available when making predictions
# We use the mode of each variable as initial fill



V_imp <- temp_data2

for (j in names(V_imp)) {
 miss_j <- is.na(V_imp[[j]])
if (any(miss_j)) {
  V_imp[[j]][miss_j] <- get_mode(V_imp[[j]])
}
}


# STEP 3 - ML imputation using T as training sample
# For each variable with missing values:
#   - fit RF on complete respondents T
#   - use all other variables as predictors
#   - predict only where that variable was missing


vars_with_missing <- names(temp_data2)[colSums(is.na(temp_data2)) > 0]

for (target in vars_with_missing) {

 cat("Imputing:", target, "\n")

 #Rows where target was originally missing
miss_idx <- which(is.na(temp_data2[[target]]))

# Predictors = all other variables
predictors <- setdiff(names(V_imp), target)

# Training set: ONLY T (complete respondents from the original matrix)
train_df <- T_data[, c(target, predictors), drop = FALSE]

# New data: rows needing imputation, taken from initialized matrix
new_df <- V_imp[miss_idx, predictors, drop = FALSE]

# Random forest classifier
rf_fit <- ranger(
 formula = as.formula(paste(target, "~ .")),
data = train_df,
classification = TRUE,
num.trees = 500,
mtry = max(1, floor(sqrt(length(predictors)))),
min.node.size = 10,
importance = "none",
probability = FALSE,
seed = 123
)

# Predict missing values
pred <- predict(rf_fit, data = new_df)$predictions

# Replace missing values in the imputed dataset
V_imp[miss_idx, target] <- pred
}

str(V_imp)
# LABELED DATASET ####

# Some levels of the demographic answers are also collapsed because the number
# of respondents choosing them is 1 or 2, which causes problems in cross-validation
data_ord_fct <- V_imp %>% mutate(
  
  # ── SOCIAL VALUES ────────────────────────────────────────────────────────────
  across(Q1:Q6,
         ~ fct_ord(.x, 1:4, c("Very important","Rather important",
                              "Not very important","Not at all important"))),
  
  across(Q7:Q17,
         ~ fct(.x, c(1,2), c("Mentioned","Not mentioned"))),
  
  across(Q18:Q26,
         ~ fct(.x, c(1,2), c("Mentioned","Not mentioned"))),
  
  across(Q27:Q32,
         ~ fct_ord(.x, 1:4, c("Strongly agree","Agree","Disagree","Strongly disagree"))),
  
  across(Q33:Q41,
         ~ fct_ord(.x, 1:5, c("Agree strongly","Agree","Neither agree nor disagree",
                              "Disagree","Disagree strongly"))),
  
  Q42 = fct(Q42, 1:3, c("Revolutionary change","Gradual reform","Defend present society")),
  
  across(Q43:Q45,
         ~ fct(.x, 1:3, c("Good","Don't mind","Bad"))),
  
  # ── HAPPINESS & WELLBEING ────────────────────────────────────────────────────
  Q46 = fct_ord(Q46, 1:4, c("Very happy","Rather happy","Not very happy","Not at all happy")),
  Q47 = fct_ord(Q47, 1:5, c("Very good","Good","Fair","Poor","Very poor")),
  Q48 = fct_10(Q48),
  Q49 = fct_10(Q49),
  Q50 = fct_10(Q50),
  
  across(Q51:Q55,
         ~ fct_ord(.x, 1:4, c("Often","Sometimes","Rarely","Never"))),
  
  Q56 = fct(Q56, 1:3, c("Better off","Worse off","About the same")),
  
  # ── SOCIAL CAPITAL & TRUST ───────────────────────────────────────────────────
  Q57 = fct(Q57, c(1,2), c("Most people can be trusted","Need to be very careful")),
  
  across(Q58:Q63,
         ~ fct_ord(.x, 1:4, c("Trust completely","Trust somewhat",
                              "Do not trust very much","Do not trust at all"))),
  
  across(Q64:Q89,
         ~ fct_ord(.x, 1:4, c("A great deal","Quite a lot","Not very much","None at all"))),
  
  Q90 = fct_10(Q90),
  
  across(Q91:Q93,
         ~ fct(.x, 1:3, c("A","B","C"))),
  
  across(Q94:Q105,
         ~ fct_ord(.x, c(0,1,2), c("Don't belong","Inactive member","Active member"))),
  
  # ── ECONOMIC VALUES ──────────────────────────────────────────────────────────
  across(Q106:Q110, ~ fct_10(.x)),
  
  Q111 = fct(Q111, c(1,2,3), c("Protect environment","Economic growth","Other")),
  
  # ── CORRUPTION ───────────────────────────────────────────────────────────────
  Q112 = fct_10(Q112),
  
  across(Q113:Q117,
         ~ fct_ord(.x, 1:4, c("None of them","Few of them","Most of them","All of them"))),
  
  Q118 = fct_ord(Q118, 1:4, c("Never","Rarely","Frequently","Always")),
  Q119 = fct(Q119, c(1,2,3,4,0),
             c("Strongly agree","Agree","Disagree","Strongly disagree","Hard to say")),
  Q120 = fct_10(Q120),
  
  # ── MIGRATION ────────────────────────────────────────────────────────────────
  Q121 = fct_ord(Q121, 1:5, c("Very bad","Quite bad","Neither good nor bad",
                              "Quite good","Very good")),
  
  across(Q122:Q129,
         ~ fct(.x, c(0,1,2), c("Disagree","Hard to say","Agree"))),
  
  Q130 = fct_ord(Q130, 1:4, c("Let anyone come","Let if jobs available",
                              "Strict limits","Prohibit")),
  
  # ── SECURITY ─────────────────────────────────────────────────────────────────
  Q131 = fct_ord(Q131, 1:4, c("Very secure","Quite secure",
                              "Not very secure","Not at all secure")),
  
  across(Q132:Q138,
         ~ fct_ord(.x, 1:4, c("Very frequently","Quite frequently",
                              "Not frequently","Not at all frequently"))),
  
  across(Q139:Q141,
         ~ fct(.x, c(1,2), c("Yes","No"))),
  
  across(Q142:Q143,
         ~ fct_ord(.x, 1:4, c("Very much","A good deal","Not much","Not at all"))),
  
  across(Q144:Q145,
         ~ fct(.x, c(1,2), c("Yes","No"))),
  
  across(Q146:Q148,
         ~ fct_ord(.x, 1:4, c("Very much","A good deal","Not much","Not at all"))),
  
  Q149 = fct(Q149, c(1,2), c("Freedom","Equality")),
  Q150 = fct(Q150, c(1,2), c("Freedom","Security")),
  Q151 = fct(Q151, c(1,2), c("Yes","No")),
  
  # ── POSTMATERIALIST INDEX ────────────────────────────────────────────────────
  across(Q152:Q153,
         ~ fct(.x, 1:4, c("High economic growth","Strong defense forces",
                          "More say at work/community","Beautiful cities/countryside"))),
  
  across(Q154:Q155,
         ~ fct(.x, 1:4, c("Maintaining order","More say in government",
                          "Fighting rising prices","Protecting freedom of speech"))),
  
  across(Q156:Q157,
         ~ fct(.x, 1:4, c("Stable economy","Less impersonal society",
                          "Ideas count more than money","Fight against crime"))),
  
  # ── SCIENCE & TECHNOLOGY ─────────────────────────────────────────────────────
  across(Q158:Q163, ~ fct_10(.x)),
  
  # ── RELIGIOUS VALUES ─────────────────────────────────────────────────────────
  Q164 = fct_10(Q164),
  
  across(Q165:Q168,
         ~ fct(.x, c(1,2), c("Yes","No"))),
  
  across(Q169:Q170,
         ~ fct_ord(.x, 1:4, c("Strongly agree","Agree","Disagree","Strongly disagree"))),
  
  Q171 = fct_ord(Q171, 1:7, c("More than once a week","Once a week","Once a month",
                              "Only on special holy days","Once a year",
                              "Less often","Never, practically never")),
  Q172 = fct_ord(Q172, 1:8, c("Several times a day","Once a day",
                              "Several times each week",
                              "Only when attending religious services",
                              "Only on special holy days","Once a year",
                              "Less often","Never, practically never")),
  Q173 = fct(Q173, 1:3, c("A religious person","Not a religious person","An atheist")),
  Q174 = fct(Q174, c(1,2), c("Follow religious norms and ceremonies",
                             "Do good to other people")),
  Q175 = fct(Q175, c(1,2), c("Make sense of life after death",
                             "Make sense of life in this world")),
  
  # ── ETHICAL VALUES & NORMS ───────────────────────────────────────────────────
  Q176 = fct_10(Q176),
  
  across(Q177:Q195, ~ fct_10(.x)),
  
  across(Q196:Q198,
         ~ fct_ord(.x, 1:4, c("Definitely should have the right",
                              "Probably should have the right",
                              "Probably should not have the right",
                              "Definitely should not have the right"))),
  
  # ── POLITICAL INTEREST & PARTICIPATION ───────────────────────────────────────
  Q199 = fct_ord(Q199, 1:4, c("Very interested","Somewhat interested",
                              "Not very interested","Not at all interested")),
  Q200 = fct_ord(Q200, 1:3, c("Frequently","Occasionally","Never")),
  
  across(Q201:Q208,
         ~ fct_ord(.x, 1:5, c("Daily","Weekly","Monthly","Less than monthly","Never"))),
  
  across(Q209:Q212,
         ~ fct(.x, 1:3, c("Have done","Might do","Would never do"))),
  
  across(Q213:Q216,
         ~ fct(.x, 1:3, c("Have done","Might do","Would never do"))),
  
  across(Q217:Q220,
         ~ fct(.x, 1:3, c("Have done","Might do","Would never do"))),
  
  across(Q221:Q222,
         ~ fct(.x, c(1,2,3,4), c("Always","Usually","Never","Not allowed to vote"))),
  
  Q223 = fct(Q223, sort(unique(data_clean$Q223[data_clean$Q223 > 0])),
             paste0("Party ", sort(unique(data_clean$Q223[data_clean$Q223 > 0])))),
  
  across(Q224:Q233,
         ~ fct_ord(.x, 1:4, c("Very often","Fairly often","Not often","Not at all often"))),
  
  Q234 = fct_ord(Q234, 1:4, c("Very important","Rather important",
                              "Not very important","Not at all important")),
  
  # ── POLITICAL CULTURE & REGIMES ──────────────────────────────────────────────
  across(Q235:Q239,
         ~ fct_ord(.x, 1:4, c("Very good","Fairly good","Fairly bad","Very bad"))),
  
  Q240 = fct_10(Q240),
  
  across(Q241:Q249, ~ fct_10(.x)),
  
  Q250 = fct_10(Q250),
  Q251 = fct_10(Q251),
  Q252 = fct_10(Q252),
  
  Q253 = fct_ord(Q253, 1:4, c("A great deal of respect","Fairly much respect",
                              "Not much respect","No respect at all")),
  Q254 = fct_ord(Q254, c(1,2,3,4,5),
                 c("Very proud","Quite proud","Not very proud",
                   "Not at all proud","Not this nationality")),
  
  across(Q255:Q259,
         ~ fct_ord(.x, 1:4, c("Very close","Close","Not very close","Not close at all"))),
  
  # ── DEMOGRAPHICS ─────────────────────────────────────────────────────────────
  Q260 = fct(Q260, c(1,2), c("Male","Female")),
  
  across(Q263:Q265,
         ~ fct(.x, c(1,2), c("Born in this country","Immigrant"))),
  
  Q269 = fct(Q269, c(1,2), c("Yes, citizen","Not a citizen")),
  Q270 = fct_ord(Q270, 1:20, as.character(1:20)),
  
  Q271 = fct(Q271, c(1,2,3,4), c("No","Yes, own parents",
                                 "Yes, parents in law","Yes, both")),
  Q273 = fct(Q273, 1:6, c("Married","Living together as married","Divorced",
                          "Separated","Widowed","Single")),
  
  Q274 = fct_ord(Q274, 0:15, c("No children","1","2","3","4","5","6",
                               ">6",">6",">6",">6",">6",">6",">6",">6",">6")),
  
  Q275 = fct_ord(Q275, 0:8, c("No education","Primary","Lower secondary",
                              "Upper secondary","Post-secondary non-tertiary",
                              "Short-cycle tertiary","Bachelor","Master","Doctoral")),
  Q277 = fct_ord(Q277, 0:8, c("No education","Primary","Lower secondary",
                              "Upper secondary","Post-secondary non-tertiary",
                              "Short-cycle tertiary","Bachelor","Master","Doctoral")),
  Q278 = fct_ord(Q278, 0:8, c("No education","Primary","Lower secondary",
                              "Upper secondary","Post-secondary non-tertiary",
                              "Short-cycle tertiary","Bachelor","Master","Doctoral")),
  
  Q279 = fct(Q279, 1:8, c("Full time employee","Part time employee","Self employed",
                          "Retired/pensioned","Housewife not otherwise employed",
                          "Student","Unemployed","Other")),
  
  Q281 = fct(Q281, c(0,1,2,3,4,5,6,7,8,9,10,11),
             c("Never had a job",
               "Professional/technical","Higher administrative","Clerical","Sales",
               "Service","Skilled worker","Semi-skilled worker","Unskilled worker",
               "Farm","Farm","Others")), # maybe check that better 
  
  
  Q282 = fct(Q282, c(0,1,2,3,4,5,6,7,8,9,10),
             c("Never had a job",
               "Professional/technical","Higher administrative","Clerical","Sales",
               "Service","Skilled worker","Semi-skilled worker","Unskilled worker",
               "Farm","Farm")),
  
  Q283 = fct(Q283, c(1,2,3,4,5,6,7,8,9,10,0),
             c("Professional/technical","Higher administrative","Clerical","Sales",
               "Service","Skilled worker","Semi-skilled worker","Unskilled worker",
               "Farm","Farm","Never had a job")),
  
  Q284 = fct(Q284, 1:3, c("Government/public institution",
                          "Private business/industry",
                          "Private non-profit organization")),
  Q285 = fct(Q285, c(1,2), c("Yes","No")),
  Q286 = fct_ord(Q286, 1:4, c("Saved money","Just get by",
                              "Spent some savings","Spent savings and borrowed")),
  Q287 = fct_ord(Q287, 1:5, c("Upper class","Upper middle class",
                              "Lower middle class","Working class","Lower class")),
  Q288 = fct_10(Q288),
  
  Q289 = fct(Q289, c(0,1,2,3,4,5,6,7,8,9),
             c("No denomination","Roman Catholic","Protestant","Orthodox",
               "Jew","Muslim","Hindu","Buddhist","Other","11")),
  Q290 = fct(Q290, c(840001,840002,840003,840004,840005),
             c("White","Black","South Asian","East Asian","Arabic/Central Asian")
))

sum(is.na(data_ord_fct$Q281))
sum(is.na(data_ord_fct))



data_ord_fct <- data_ord_fct %>%
  mutate(across(where(is.factor), fct_drop),
         Q261 = as.numeric(as.character(Q261)),
         Q262 = as.numeric(as.character(Q262)))

table(data_ord_fct$Q282)

# 4 & 5 Political and demo  ####
# Removing Political and Demographics 
vars_to_remove <- c(Political, Demographics)


imputed_selected_labelled <- data_ord_fct %>%
  select(-any_of((vars_to_remove)))
ncol(imputed_selected_labelled)
imputed_selected <- V_imp %>%
  select(-any_of((vars_to_remove)))
ncol(imputed_selected)

# -- Save ---------------------------------------------------------------------
saveRDS(V_imp,                     file.path(DIR_DATA, "data_imputed.rds"))
saveRDS(imputed_selected,          file.path(DIR_DATA, "data_imputed_selected.rds"))
saveRDS(data_ord_fct,              file.path(DIR_DATA, "imputed_labeled.rds"))
saveRDS(imputed_selected_labelled, file.path(DIR_DATA, "imputed_selected_labeled.rds"))


# Classification dataset 
imputed_class <- data_ord_fct %>%
  select(any_of((vars_to_remove))) %>% 
  dplyr::filter(Q276 != "0") %>%
  mutate(Q276 = fct_drop(Q276))

saveRDS(imputed_class, file.path(DIR_DATA, "data_class.rds"))
which(data_ord_fct$Q276 == 0)

sum(is.na(imputed_class))

no_lab <- V_imp %>%
  select(any_of((vars_to_remove))) %>% 
  dplyr::filter(Q276 != "0") %>%
  mutate(Q276 = fct_drop(Q276),
         across(where(is.factor), fct_drop),
         Q261 = as.numeric(as.character(Q261)),
         Q262 = as.numeric(as.character(Q262)),
         Q274 = as.numeric(as.character(Q274)))

