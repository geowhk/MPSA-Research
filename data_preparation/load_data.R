# --- 환경 설정 ---
source("R/data_preparation/setup.R")

# --- 데이터 불러오기 ---
census_api_key("f7a2c390c371e265b4c53f0ccae1f6e900616b74", install = FALSE)

vars <- c(
  # 연속형 변수
  med_income     = "B19013_001",  # Median household income
  med_homeval    = "B25077_001",  # Median home value
  med_age        = "B01002_001",  # Median age
  mean_commute   = "B08119_001",  # Mean commute time
  avg_hh_size    = "B25010_001",  # Average household size
  pop_total      = "B01003_001",  # Total population
  
  # 비율 계산용 count 변수
  renter_count     = "B25003_003",
  owner_count      = "B25003_002",
  housing_total    = "B25003_001",
  
  hs_grad_count    = "B15003_017",
  bachelors_count  = "B15003_022",
  edu_total        = "B15003_001",
  
  white_count      = "B02001_002",
  black_count      = "B02001_003",
  asian_count      = "B02001_005",
  race_total       = "B02001_001",
  
  hispanic_count   = "B03002_012",
  hispanic_total   = "B03002_001",
  
  vacant_count     = "B25002_003",
  housing_all      = "B25002_001",
  
  unemployment     = "B23025_005",
  labor_force      = "B23025_003",
  
  # 산업군 종사자 수 (범주형 변수 생성용)
  total_workers    = "C24050_001",
  ag_forest        = "C24050_003",
  manufacturing    = "C24050_005",
  retail           = "C24050_007",
  transport        = "C24050_008",
  education        = "C24050_014",
  information      = "C24050_017",
  finance          = "C24050_018"
)


franklin_tracts <- get_acs(
  geography = "tract",
  state = "OH",
  county = "Franklin",
  variables = vars,
  year = 2020,
  survey = "acs5",
  geometry = TRUE
)

franklin <- franklin_tracts |> 
  select(GEOID, NAME, variable, estimate, geometry) |> 
  pivot_wider(
    names_from = variable,
    values_from = estimate
  )

franklin <- franklin |> 
  mutate(
    # ✅ 비율 변수
    pct_renter        = renter_count / housing_total * 100,
    pct_owner         = owner_count / housing_total * 100,
    pct_vacant        = vacant_count / housing_all * 100,
    pct_white         = white_count / race_total * 100,
    pct_black         = black_count / race_total * 100,
    pct_asian         = asian_count / race_total * 100,
    pct_hispanic      = hispanic_count / hispanic_total * 100,
    hs_grad_rate      = hs_grad_count / edu_total * 100,
    bachelors_rate    = bachelors_count / edu_total * 100,
    unemployment_rate = unemployment / labor_force * 100,
    
    # ✅ 산업군 비율
    p_ag_forest     = ag_forest / total_workers,
    p_manufacturing = manufacturing / total_workers,
    p_retail        = retail / total_workers,
    p_transport     = transport / total_workers,
    p_education     = education / total_workers,
    p_information   = information / total_workers,
    p_finance       = finance / total_workers,
    
    # ✅ 주요 산업군(factor): 가장 비율이 높은 산업군을 범주형으로 지정
    main_industry = case_when(
      p_ag_forest     >= pmax(p_ag_forest, p_manufacturing, p_retail, p_transport, p_education, p_information, p_finance, na.rm = TRUE) ~ "agriculture",
      p_manufacturing >= pmax(p_ag_forest, p_manufacturing, p_retail, p_transport, p_education, p_information, p_finance, na.rm = TRUE) ~ "manufacturing",
      p_retail        >= pmax(p_ag_forest, p_manufacturing, p_retail, p_transport, p_education, p_information, p_finance, na.rm = TRUE) ~ "retail",
      p_transport     >= pmax(p_ag_forest, p_manufacturing, p_retail, p_transport, p_education, p_information, p_finance, na.rm = TRUE) ~ "transport",
      p_education     >= pmax(p_ag_forest, p_manufacturing, p_retail, p_transport, p_education, p_information, p_finance, na.rm = TRUE) ~ "education",
      p_information   >= pmax(p_ag_forest, p_manufacturing, p_retail, p_transport, p_education, p_information, p_finance, na.rm = TRUE) ~ "information",
      p_finance       >= pmax(p_ag_forest, p_manufacturing, p_retail, p_transport, p_education, p_information, p_finance, na.rm = TRUE) ~ "finance",
      TRUE ~ NA_character_
    ),
    main_industry = factor(main_industry, levels = c("agriculture", "manufacturing", "retail", "transport", "education", "information", "finance"))
  )

vars_to_drop <- c(
  # 기존 count 및 분모 변수
  "renter_count", "owner_count", "housing_total",
  "hs_grad_count", "bachelors_count", "edu_total",
  "white_count", "black_count", "asian_count", "race_total",
  "hispanic_count", "hispanic_total",
  "vacant_count", "housing_all",
  "unemployment", "labor_force",
  
  # 산업군 count
  "ag_forest", "manufacturing", "retail", "transport",
  "education", "information", "finance", "total_workers",
  
  # 산업군 비율 (중간 계산)
  "p_ag_forest", "p_manufacturing", "p_retail", "p_transport",
  "p_education", "p_information", "p_finance"
)


franklin <- franklin |> 
  select(-all_of(vars_to_drop))

# 공항(모든 컬럼 NA) 제거
franklin <- franklin |> 
  filter(GEOID != "39049980000")

# 결측치 대체
df_numeric <- franklin |> 
  st_drop_geometry() |> 
  select(-GEOID, -NAME)

df_imputed <- df_numeric |> 
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

franklin <- franklin |> 
  select(GEOID, NAME, geometry) |> 
  bind_cols(df_imputed)

qtm(franklin)

write_rds(franklin, "data/franklin.rds")
