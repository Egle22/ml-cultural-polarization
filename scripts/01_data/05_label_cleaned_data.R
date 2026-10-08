# =============================================================================
# scripts/01_data/05_label_cleaned_data.R
# Purpose : Label the cleaned (non-imputed) data as ordered / nominal factors,
#           keeping the non-response levels, for the codebook. The labels
#           differ slightly from those in 04_preprocess_impute.R (finer
#           occupational and ethnic categories), so the two blocks are kept
#           separate.
# Thesis  : Appendix A.1 (codebook)
# Input   : data/processed/data_clean.rds
# Output  : data/processed/data_ord_fct.rds
# =============================================================================

source(here::here("R", "00_setup.R"))

data_clean <- readRDS(file.path(DIR_DATA, "data_clean.rds"))

data_ord_fct <- data_clean %>% mutate(
  
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
  Q274 = fct_ord(Q274, 0:15, as.character(0:15)),
  
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
  
  Q283 = fct(Q283, c(1,2,3,4,5,6,7,8,9,10,0),
             c("Professional/technical","Higher administrative","Clerical","Sales",
               "Service","Skilled worker","Semi-skilled worker","Unskilled worker",
               "Farm worker","Farm proprietor/manager","Never had a job")),
  
  Q284 = fct(Q284, 1:3, c("Government/public institution",
                          "Private business/industry",
                          "Private non-profit organization")),
  Q285 = fct(Q285, c(1,2), c("Yes","No")),
  Q286 = fct_ord(Q286, 1:4, c("Saved money","Just get by",
                              "Spent some savings","Spent savings and borrowed")),
  Q287 = fct_ord(Q287, 1:5, c("Upper class","Upper middle class",
                              "Lower middle class","Working class","Lower class")),
  Q288 = fct_10(Q288),
  
  Q289 = fct(Q289, c(0,1,2,3,4,5,6,7,8),
             c("No denomination","Roman Catholic","Protestant","Orthodox",
               "Jew","Muslim","Hindu","Buddhist","Other")),
  Q290 = fct(Q290, c(1,2,3,4,5,6),
             c("White","Black","South Asian","East Asian","Arabic/Central Asian","Other"))
)

saveRDS(data_ord_fct, file.path(DIR_DATA, "data_ord_fct.rds"))
