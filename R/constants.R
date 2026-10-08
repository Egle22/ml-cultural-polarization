# =============================================================================
# R/constants.R
# Project-wide constants: cluster colours, WVS missing-value codes, the
# questionnaire sections and the variable blocks (values / political /
# sociodemographic). Sourced by R/00_setup.R.
# =============================================================================

# -- Cluster colours (shared by every figure) ---------------------------------


CLUSTER_COLORS <- c(
  "1"  = "#4C7896", "2"  = "#B5705B", "3"  = "#5A9E6F", "4"  = "#8B6BAE",
  "5"  = "#D4A843", "6"  = "#5E9DB8", "7"  = "#C2614A", "8"  = "#7BAF7A",
  "9"  = "#9B59B6", "10" = "#E8A030", "11" = "#3B7FB5", "12" = "#D46A52",
  "13" = "#4CAF82", "14" = "#7B5EA7", "15" = "#C8922A", "16" = "#4A8FA8"
)

cluster_colors <- CLUSTER_COLORS

# -- WVS missing-value codes ---------------------------------------------------
missing_codes <- c(-1, -2, -3, -4, -5)

missing_labels <- c(
  "-1" = "Don't know",
  "-2" = "No answer / Refused",
  "-3" = "Not applicable (filter)",
  "-4" = "Not asked in survey",
  "-5" = "Missing / Not applicable other reasons"
)

miss_lvl <- c(-5, -3, -2, -1)
miss_lbl <- c("Missing/NA", "Not applicable", "No answer/Refused", "Don't know")

# -- Questionnaire sections and variable blocks --------------------------------
section_labels <- c(
  "Social Values", "Happiness & Wellbeing", "Social Capital & Trust",
  "Economic Values", "Corruption", "Migration", "Security",
  "Postmaterialist Index", "Science & Technology", "Religious Values",
  "Ethical Values & Norms", "Political Interest", "Political Culture", "Demographics"
)

section_colors <- c(
  "Social Values"          = "#4E79A7",
  "Happiness & Wellbeing"  = "#F28E2B",
  "Social Capital & Trust" = "#E15759",
  "Economic Values"        = "#76B7B2",
  "Corruption"             = "#59A14F",
  "Migration"              = "#EDC948",
  "Security"               = "#B07AA1",
  "Postmaterialist Index"  = "#FF9DA7",
  "Science & Technology"   = "#9C755F",
  "Religious Values"       = "#BAB0AC",
  "Ethical Values & Norms" = "#499894",
  "Political Culture"  = "#4E79A7"  , # blue
  "Political Interest" = "#59A14F"   , # green
  "Demographics"       = "#B07AA1"   # purple
)
section_labels_named <- c(
  social_values_attitudes_stereotypes  = "Social Values",
  happiness_wellbeing                  = "Happiness & Wellbeing",
  social_capital_trust_membership      = "Social Capital & Trust",
  economic_values                      = "Economic Values",
  corruption                           = "Corruption",
  migration                            = "Migration",
  security                             = "Security",
  postmaterialist_index                = "Postmaterialist Index",
  science_technology                   = "Science & Technology",
  religious_values                     = "Religious Values",
  ethical_values_norms                 = "Ethical Values & Norms",
  political_interest_participation     = "Political Interest",
  political_culture_regimes            = "Political Culture",
  demographics                         = "Demographics"
)

# QUESTIONS ----
questions <- list(
  
  social_values_attitudes_stereotypes = c(
    "Q1","Q2","Q3","Q4","Q5","Q6",           # importance of life domains
    "Q7","Q8","Q9","Q10","Q11","Q12","Q13",   # child qualities
    "Q14","Q15","Q16","Q17",
    "Q18","Q19","Q20","Q21","Q22","Q23",      # unwanted neighbours
    "Q24","Q25","Q26",
    "Q27","Q28","Q29","Q30","Q31","Q32",      # gender roles (4-point)
    "Q33","Q34","Q35","Q36","Q37","Q38",      # gender roles (5-point)
    "Q39","Q40","Q41",
    "Q42",                                    # attitude to social change
    "Q43","Q44","Q45"                         # evaluation of future changes
  ),
  
  happiness_wellbeing = c(
    "Q46",                                    # subjective happiness
    "Q47",                                    # self-rated health
    "Q48",                                    # freedom of choice
    "Q49",                                    # life satisfaction
    "Q50",                                    # financial satisfaction
    "Q51","Q52","Q53","Q54","Q55",            # material deprivation
    "Q56"                                     # standard of living vs parents
  ),
  
  social_capital_trust_membership = c(
    "Q57",                                    # generalised trust
    "Q58","Q59","Q60","Q61","Q62","Q63",      # trust in specific groups
    "Q64","Q65","Q66","Q67","Q68","Q69",      # confidence in institutions
    "Q70","Q71","Q72","Q73","Q74","Q75",
    "Q76","Q77","Q78","Q79","Q80","Q81",
    "Q82","Q83","Q84","Q85","Q86","Q87",
    "Q88","Q89",
    "Q90",                                    # effectiveness vs democracy of intl. orgs
    "Q91","Q92","Q93",                        # knowledge of international organisations
    "Q94","Q95","Q96","Q97","Q98","Q99",      # organisation membership
    "Q100","Q101","Q102","Q103","Q104","Q105"
  ),
  
  economic_values = c(
    "Q106",                                   # income equality vs incentives
    "Q107",                                   # private vs state ownership
    "Q108",                                   # government vs individual responsibility
    "Q109",                                   # competition: good vs harmful
    "Q110",                                   # hard work vs luck
    "Q111"                                    # environment vs economic growth
  ),
  
  corruption = c(
    "Q112",                                   # perceived corruption level
    "Q113","Q114","Q115","Q116","Q117",       # corruption by group
    "Q118",                                   # personal experience of corruption
    "Q119",                                   # women less corrupt than men
    "Q120"                                    # risk of being held accountable for corruption
  ),
  
  migration = c(
    "Q121",                                   # evaluation of the impact of immigration
    "Q122","Q123","Q124","Q125","Q126",       # specific effects of immigration
    "Q127","Q128","Q129",
    "Q130"                                    # preferred immigration policy
  ),
  
  security = c(
    "Q131",                                   # feeling of personal security
    "Q132","Q133","Q134","Q135","Q136",       # neighbourhood disorder
    "Q137","Q138",
    "Q139","Q140","Q141",                     # security behaviours adopted
    "Q142","Q143",                            # worries: job, education
    "Q144","Q145",                            # crime victimisation
    "Q146","Q147","Q148",                     # geopolitical worries
    "Q149",                                   # freedom vs equality
    "Q150",                                    # freedom vs security
    "Q151"                                   # willingness to fight for the country
  ),
  
  postmaterialist_index = c(
    "Q152","Q153",                            # societal priorities set 1
    "Q154","Q155",                            # societal priorities set 2
    "Q156","Q157"                             # societal priorities set 3
  ),
  
  science_technology = c(
    "Q158",                                   # S&T make life better
    "Q159",                                   # S&T: future opportunities
    "Q160",                                   # too much science, not enough faith
    "Q161",                                   # science erodes morality
    "Q162",                                   # science not important in daily life
    "Q163"                                    # world better/worse off because of S&T
  ),
  
  religious_values = c(
    "Q164",                                   # importance of God
    "Q165","Q166","Q167","Q168",              # religious beliefs
    "Q169","Q170",                            # science vs religion attitudes
    "Q171",                                   # attendance at religious services
    "Q172",                                   # frequency of prayer
    "Q173",                                   # religious self-identification
    "Q174","Q175"                             # meaning of religion
  ),
  
  ethical_values_norms = c(
    "Q176",                                   # difficulty of moral rules
    "Q177","Q178","Q179","Q180","Q181",       # justifiability of behaviours
    "Q182","Q183","Q184","Q185","Q186",
    "Q187","Q188","Q189","Q190","Q191",
    "Q192","Q193","Q194","Q195",
    "Q196","Q197","Q198"                      # government surveillance rights
  ),
  
  political_interest_participation = c(
    "Q199",                                   # interest in politics
    "Q200",                                   # political discussion with friends
    "Q201","Q202","Q203","Q204","Q205",       # information sources
    "Q206","Q207","Q208",
    "Q209","Q210","Q211","Q212",              # conventional political actions
    "Q213","Q214",
    #"Q215",
    "Q216",              # social activism
    "Q217","Q218","Q219","Q220",              # online political actions
    "Q221","Q222",                            # voting frequency
    "Q223",                                   # vote intention
    "Q224","Q225","Q226","Q227","Q228",       # electoral integrity
    "Q229","Q230","Q231","Q232","Q233",
    "Q234"                                  # importance of honest elections
    #"SDG1"                                   # degree of political influence 
  ),
  
  political_culture_regimes = c(
    "Q235","Q236","Q237","Q238","Q239",       # evaluation of regime types
    "Q240",                                   # left-right self-placement
    "Q241","Q242","Q243","Q244","Q245",       # essential characteristics of democracy
    "Q246","Q247","Q248","Q249",
    "Q250",                                   # importance of living in a democracy
    "Q251",                                   # current level of democracy
    "Q252",                                   # satisfaction with the political system
    "Q253",                                   # respect for human rights
    "Q254",                                   # national pride
    "Q255","Q256","Q257","Q258","Q259"        # closeness to local ... world community
  ),
  
  demographics = c(
    "Q260",                                   # sex
    "Q261","Q262",                            # year of birth / age
    "Q263","Q264","Q265",                     # born in country: respondent, mother, father
    "Q269",                                   # citizenship
    "Q270",                                   # household size
    "Q271",                                   # living with parents
    "Q273",                                   # marital status
    "Q274",                                   # number of children
    "Q275","Q276","Q277","Q278",              # education: respondent, spouse, mother, father
    "Q279","Q280",                            # employment: respondent, spouse
    "Q281","Q282","Q283",                     # occupational group
    "Q284",                                   # sector of employment
    "Q285",                                   # chief wage earner
    "Q286",                                   # family savings behaviour
    "Q287",                                   # subjective social class
    "Q288",                                   # income decile
    "Q289",                                   # religious denomination
    "Q290"                                    # ethnic group
  )
  
)

# ORDERED QUESTIONS
# Variables with orderable answers
w_binary_vars <- c(
  paste0("Q", 7:17),
  paste0("Q", 18:26),
  paste0("Q", 58:87),
  paste0("Q", 139:145),
  "Q151",
  "Q57"
)

w_ordinal_scale_vars <- c(
  paste0("Q", 1:6),
  paste0("Q", 27:56),
  paste0("Q", 88:148),
  paste0("Q", 154:208),
  paste0("Q", 209:214),
  paste0("Q", 216:222),
  paste0("Q", 224:259)
)

Demographics <- questions$demographics
Political <- c(questions$political_interest_participation, questions$political_culture_regimes)
Economic_values <- c(questions$economic_values)
Social_trust <- c(questions$social_capital_trust_membership)

demo_cols <- Demographics
pol_cols <- Political

section_lookup <- do.call(rbind, lapply(names(questions), function(sec) {
  data.frame(
    variable      = questions[[sec]],
    section_label = section_labels_named[sec],
    stringsAsFactors = FALSE
  )
})) %>% distinct(variable, .keep_all = TRUE)
