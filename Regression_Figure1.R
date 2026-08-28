library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(georefdatar)
library(lubridate)
library(corrplot)
library(ggpubr)
library(purrr)
library(broom)
library(ggthemes)
library(patchwork)
library(ggrepel)
library(writexl)
library(scales)

################################### Data import ####################################

#DL Sediment data from Heard et al. 2024 DOI:10.1016/j.gca.2024.07.037
#And Zenodo Data Repository DOI:10.5281/zenodo.21953894
DL_sediments <- read.csv("2022_Deming_Sediment.csv")

# Total carbon and nitrogen data for sediment transect
#Data from Ostrander et al 2025 DOI:https://doi.org/10.1016/j.gca.2024.12.008
#And Zenodo Data Repository DOI: 10.5281/zenodo.21953894
DLtoc_tn <- read.csv("TOC_TN_data.csv", stringsAsFactors = FALSE)

#Brownie Lake Sediment from Swanner et al
# DOI: 10.6073/pasta/68b50baa0a767ab33f2b7dd91948036e

BL_sediment <- read_csv("BL_sediment_data_v2.csv")

# Archean and early Proterozoic Shale 
# Compilation is on Zenodo Data Repository DOI:10.5281/zenodo.21953894

Archean_shale <- read_xlsx("Compiled_Shale.xlsx")

# Modern hydrothermal vent and BIF data
#Compilation is on Zenodo Data Repository DOI:10.5281/zenodo.21953894

IFs <- read_csv("Compiled_IF.csv")

# UCC from Rudnick and Gao 2014 https://doi.org/10.1016/j.gca.2013.11.006 for some elements
data(CC_Upper__Rudnick_Gao__2014)

# UCC through history from Ptáček et al 2020 https://doi.org/10.1016/j.epsl.2020.116090
ArcheanUCC <- read_xlsx("Ptacek_TableS5_crustthroughtime.xlsx")

# Chad Wittkop's MN surficial till https://doi.org/10.1016/j.chemgeo.2019.119390
CW_norm <- read_csv("CW_norm.csv")

################################### Data wrangling ####################################

###Deming Lake Sediments###

# atomic weights for molar conversions
  AW_Al = 26.982
  AW_P = 30.974
  AW_Ti = 47.867
  AW_Fe = 55.845

# Convert from ug/g to mol/kg

DL_sediments <- DL_sediments %>%
  mutate(
    Al_molkg = pAl..ug.g. / (AW_Al * 1000),
    Fe_molkg = pFe..ug.g. / (AW_Fe * 1000),
    P_molkg  = pP..ug.g.  / (AW_P  * 1000),
    Ti_molkg = pTi..ug.g. / (AW_Ti * 1000)
  )

# Add column to classify measurement as sediment transect or sediment core

DL_sediments <- DL_sediments %>%
  mutate(
    sample_type = case_when(
      str_detect(Sample.Type, "-Transect$") ~ "Transect",
      str_detect(Sample.Type, "-Core$") ~ "Core",
      TRUE ~ NA_character_
    )
  )

# Pull local till values
  
  Fe2O3_till <- CW_norm %>%
    filter(Component == "Fe2O3 %") %>%
    pull(`Average carbonate-free till (<63 um)`)
  
  P2O5_till <- CW_norm %>%
    filter(Component == "P2O5 %") %>%
    pull(`Average carbonate-free till (<63 um)`)
  
  TiO2_till <- CW_norm %>%
    filter(Component == "TiO2 %") %>%
    pull(`Average carbonate-free till (<63 um)`)
  
  Al2O3_till <- CW_norm %>%
    filter(Component == "Al2O3 %") %>%
    pull(`Average carbonate-free till (<63 um)`)

# Molar conversions for oxides
  MW_Al2O3 = 101.96
  MW_P2O5 = 141.94
  MW_TiO2 = 79.87
  MW_Fe2O3 = 159.69
  
  Fe_till <- (Fe2O3_till / MW_Fe2O3) * 2 # Multiplied by molar ratio
  P_till  <- (P2O5_till  / MW_P2O5)  * 2
  Ti_till <- (TiO2_till  / MW_TiO2)
  Al_till <- (Al2O3_till / MW_Al2O3) * 2

# Ratios for authigenic calculations
  
  FeTi_till <- Fe_till / Ti_till
  PTi_till  <- P_till / Ti_till
  FeAl_till <- Fe_till / Al_till
  PAl_till  <- P_till / Al_till

# Save CW local till ratios for reference
  CW_till_endmember <- tibble(
    FeTi_till = FeTi_till,
    PTi_till = PTi_till
  )
  
  CW_till_endmember

# Calculate detrital and authigenic fractions in Deming using local till ratios
  
  DL_sediments <- DL_sediments %>%
    mutate(
      Fe_det = Ti_molkg * FeTi_till,
      P_det  = Ti_molkg * PTi_till,
      
      Fe_auth = Fe_molkg - Fe_det,
      P_auth  = P_molkg  - P_det,
      
      Fe_auth_pct = 100 * Fe_auth / Fe_molkg,
      P_auth_pct  = 100 * P_auth  / P_molkg
    )  

# Convert UCC data to long
CC_long <- CC_Upper__Rudnick_Gao__2014 %>%
  pivot_longer(
    cols = everything(),
    names_to = "Component",
    values_to = "value"
  )

###Brownie Lake Sediments###

#Pull out relevant data from dataset
BL_FP <- BL_sediment %>%
  filter(
    `Sample Type` == "FC",
    Analyte %in% c(
      "Al2O3",
      "Fe2O3",
      "P2O5",
      "TiO2"
    )
  )


#Average replicate measurements in BL dataset
BL_FP <- BL_FP %>%
  group_by(`Sample Type`, Date, Depth, Analyte, Unit) %>%
  summarise(
    Value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

#convert to wide format 

BL_FP <- BL_FP %>%
  pivot_wider(
    names_from = Analyte,
    values_from = Value
  )

#Convert oxides to elemental molar

BL_FP <- BL_FP %>%
  mutate(
    
    Al_molkg =
      ((Al2O3 * 10) /  MW_Al2O3) * 2,
    
    Fe_molkg =
      ((Fe2O3 * 10) / MW_Fe2O3) * 2,
    
    P_molkg =
      ((P2O5 * 10) / MW_P2O5) * 2,
    
    Ti_molkg =
      (TiO2 * 10) / MW_TiO2
    
  )

#Calculate Brownie lake authigenic and detrital fractions 

BL_FP <- BL_FP %>%
  mutate(
    
    Fe_det = Ti_molkg * FeTi_till,
    P_det  = Ti_molkg * PTi_till,
    
    Fe_auth = Fe_molkg - Fe_det,
    P_auth  = P_molkg - P_det,
    
    Fe_auth_pct = 100 * Fe_auth / Fe_molkg,
    P_auth_pct  = 100 * P_auth / P_molkg
    
  )

### Brownie and Deming Lake Regression ###

#Stitch lake data together
DL_reg <- DL_sediments %>%
  select(Fe_auth, P_auth) %>%
  mutate(Lake = "Deming")

BL_reg <- BL_FP %>%
  select(Fe_auth, P_auth) %>%
  mutate(Lake = "Brownie")

FP_compare <- bind_rows(DL_reg, BL_reg)

# Linear regression and table of Brownie and Deming data

Lake_regression_table <- FP_compare %>%
  group_by(Lake) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    model = map(data, ~ lm(P_auth ~ Fe_auth, data = .x)),
    glance = map(model, broom::glance),
    tidy = map(model, broom::tidy)
  ) %>%
  mutate(
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    # 95% CI for slope
    slope_ci_low = map_dbl(model, ~ confint(.x)["Fe_auth", 1]),
    slope_ci_high = map_dbl(model, ~ confint(.x)["Fe_auth", 2]),
    slope_95CI = paste0(
      round(slope_ci_low, 3),
      " to ",
      round(slope_ci_high, 3)
    ),
    residual_se = map_dbl(glance, "sigma"),
    r_squared = map_dbl(glance, "r.squared"),
    # Pearson correlation coefficient
    pearson_r = sign(slope) * sqrt(r_squared),
    p_value = map_dbl(
      tidy,
      ~ .x %>%
        filter(term == "Fe_auth") %>%
        pull(p.value)
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    linetype = if_else(
      significance == "ns",
      "longdash",
      "solid"
    )
  ) %>%
  select(
    Lake,
    n,
    slope,
    slope_95CI,
    intercept,
    r_squared,
    p_value,
    residual_se,
    pearson_r,
    significance,
    linetype
  ) %>%
  arrange(desc(r_squared))

Lake_regression_table

# DL regression excluding two data points from lake depths 14.0 and 16.5 m

DL_FP_no_outliers <- DL_sediments %>%
  filter(!Lake.Depth..m. %in% c(14.0, 16.5))

DL_no_outliers_regression_table<- DL_FP_no_outliers %>%
  filter(
    !is.na(Fe_auth),
    !is.na(P_auth)
  ) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    model = map(data, ~ lm(P_auth ~ Fe_auth, data = .x)),
    glance = map(model, broom::glance),
    tidy = map(model, broom::tidy)
  ) %>%
  mutate(
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    # 95% CI for slope
    slope_ci_low = map_dbl(model, ~ confint(.x)["Fe_auth", 1]),
    slope_ci_high = map_dbl(model, ~ confint(.x)["Fe_auth", 2]),
    slope_95CI = paste0(
      round(slope_ci_low, 3),
      " to ",
      round(slope_ci_high, 3)
    ),
    residual_se = map_dbl(glance, "sigma"),
    r_squared = map_dbl(glance, "r.squared"),
    # Pearson correlation coefficient
    pearson_r = sign(slope) * sqrt(r_squared),
    p_value = map_dbl(
      tidy,
      ~ .x %>%
        filter(term == "Fe_auth") %>%
        pull(p.value)
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    linetype = if_else(
      significance == "ns",
      "longdash",
      "solid"
    )
  ) %>%
  select(
    n,
    slope,
    slope_95CI,
    intercept,
    r_squared,
    p_value,
    residual_se,
    pearson_r,
    significance,
    linetype
  ) %>%
  arrange(desc(r_squared))

DL_no_outliers_regression_table


tol_colors <- c(
  "#332288",  # dark blue
  "#88CCEE",  # light blue
  "#44AA99",  # teal
  "#117733",  # green
  "#999933",  # olive
  "#DDCC77",  # sand
  "#CC6677",  # rose
  "#882255",  # wine
  "#AA4499"   # purple
)

lakes_label_data <- tibble(
  `Lake` = unique(FP_compare$`Lake`)
) %>%
  mutate(
    x = max(FP_compare$Fe_auth, na.rm = TRUE) * -0.1,
    y = seq(
      max(FP_compare$P_auth, na.rm = TRUE),
      max(FP_compare$P_auth, na.rm = TRUE) * 0.95,
      length.out = n()
    )
  )

lake_linetypes <- setNames(
  Lake_regression_table$linetype,
  Lake_regression_table$Lake
)

lakes_plot <- ggplot(
  FP_compare,
  aes(
    x = Fe_auth,
    y = P_auth,
    color = Lake
  )
) +
  geom_point(size = 2.5, alpha = 0.4) +
  
  geom_smooth(aes(linetype = Lake),
    method = "lm",
    se = FALSE,
    linewidth = 1.75,
    lineend = "round"
  ) +
  geom_abline(
    intercept = coef(DL_lm_no_outliers)[1],
    slope = coef(DL_lm_no_outliers)[2],
    linetype = "longdash",
    lineend = "round",
    color = "#88CCEE",
    linewidth = 1.75
  ) +
  scale_color_manual(values = tol_colors) +
  scale_linetype_manual(values = lake_linetypes) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.05)) +
  scale_x_continuous(breaks = seq(0, 0.6, by = 0.2)) +
  geom_text(
    data = lakes_label_data,
    aes(
      x = x,
      y = y,
      label = Lake,
      color = Lake
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +

  labs(
    x = expression("Authigenic Fe (mol kg"^{-1}*")"),
    y = expression("Authigenic P (mol kg"^{-1}*")"),
    color = "Lake",
    title = "Ferruginous Lakes",
    tag = "A"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5),
        axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1))

print(lakes_plot)

### Archean Shale ###

Archean_FP <- Archean_shale %>%
  mutate(
    Al_molkg =
      ((`Al2O3 [wt.%]` * 10) / MW_Al2O3) * 2,
    Fe_molkg =
      ((`Fe2O3 [wt.%]` * 10) / MW_Fe2O3) * 2,
    P_molkg =
      ((`P2O5 [wt.%]` * 10) / MW_P2O5) * 2,
    Ti_molkg =
      (`TiO2 [wt.%]` * 10) / MW_TiO2
  )

# Make lookup table for Archean UCC values

ArcheanUCC_lookup <- ArcheanUCC %>%
  transmute(
    ucc_age = `TIME (MA)`,
    
    UCC_FeAl =
      (FeOt / 71.844) /
      ((Al2O3 / MW_Al2O3) * 2),
    
    UCC_PAl =
      ((P2O5 / MW_P2O5) * 2) /
      ((Al2O3 / MW_Al2O3) * 2)
  )

#Group Archean shale record to nearest age in Archean UCC
Archean_FP <- Archean_FP %>%
  rowwise() %>%
  mutate(
    ucc_age =
      ArcheanUCC_lookup$ucc_age[
        which.min(
          abs(ArcheanUCC_lookup$ucc_age - `Age (MA)`)
        )
      ]
  ) %>%
  ungroup()

#Join Archean shale with Archean UCC based on nearest age

Archean_FP <- Archean_FP %>%
  left_join(
    ArcheanUCC_lookup,
    by = "ucc_age"
  )

# Calculate authigenic and detrital fractions of Archean shale
Archean_FP <- Archean_FP %>%
  mutate(
    Fe_det = Al_molkg * UCC_FeAl,
    P_det  = Al_molkg * UCC_PAl,
    
    Fe_auth = Fe_molkg - Fe_det,
    P_auth  = P_molkg - P_det,
    
    Fe_auth_pct = 100 * Fe_auth / Fe_molkg,
    P_auth_pct  = 100 * P_auth / P_molkg
  )

# Linear regression of auth_Fe vs auth_P Archean shale record

formation_models <- Archean_FP %>%
  filter(
    !is.na(Fe_auth),
    !is.na(P_auth)
  ) %>%
  group_by(`Formation/ODP Site`) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    model = map(data, ~ lm(P_auth ~ Fe_auth, data = .x)),
    glance = map(model, broom::glance),
    tidy = map(model, broom::tidy)
  )

ArcheanShale_regression_table <- formation_models %>%
  mutate(
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    
    # 95% CI for slope
    slope_ci_low = map_dbl(model, ~ confint(.x)["Fe_auth", 1]),
    slope_ci_high = map_dbl(model, ~ confint(.x)["Fe_auth", 2]),
    slope_95CI = paste0(
      round(slope_ci_low, 3),
      " to ",
      round(slope_ci_high, 3)
    ),
    
    # Model statistics
    r_squared = map_dbl(glance, "r.squared"),
    residual_se = map_dbl(glance, "sigma"),
    
    # p-value for slope
    p_value = map_dbl(
      tidy,
      ~ .x %>%
        filter(term == "Fe_auth") %>%
        pull(p.value)
    ),
    
    # Pearson correlation coefficient
    pearson_r = sign(slope) * sqrt(r_squared)
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    linetype = if_else(
      significance == "ns",
      "longdash",
      "solid"
    )
  ) %>%
  select(
    Formation = `Formation/ODP Site`,
    n,
    slope,
    slope_95CI,
    intercept,
    r_squared,
    p_value,
    residual_se,
    pearson_r,
    significance,
    linetype
  ) %>%
  arrange(desc(r_squared))

ArcheanShale_regression_table

# Filter the six formations with most sample data points

top6_shaleformations <- Archean_FP %>%
  count(`Formation/ODP Site`, sort = TRUE) %>%
  slice_head(n = 6)

top6_shaleformations

top6_shalenames <- top6_shaleformations$`Formation/ODP Site`

Archean_top6 <- Archean_FP %>%
  filter(`Formation/ODP Site` %in% top6_shalenames)

#Edit names
Archean_top6<- Archean_top6 %>%
  mutate(
    Formation_label = `Formation/ODP Site` %>%
      stringr::str_remove_all("Formation") %>%
      stringr::str_squish()
  )

#Filter top 6 from  regression table
regression_table_top6 <- ArcheanShale_regression_table %>%
  filter(`Formation/ODP Site` %in% top6_shalenames)

Archean_top6_label_data <- tibble(
  `Formation_label` = unique(Archean_top6$`Formation_label`)
) %>%
  mutate(
    x = max(Archean_top6$Fe_auth, na.rm = TRUE) * -0.25,
    y = seq(
      max(Archean_top6$P_auth, na.rm = TRUE),
      max(Archean_top6$P_auth, na.rm = TRUE) * 0.65,
      length.out = n()
    )
  )

shale_linetypes <- regression_table_top6 %>%
  mutate(
    Formation_label = str_remove(
      `Formation/ODP Site`,
      " Formation$"
    )
  ) %>%
  select(Formation_label, linetype)

shale_linetypes <- setNames(
  shale_linetypes$linetype,
  shale_linetypes$Formation_label
)


Top6_shaleplot <- ggplot(
  Archean_top6,
  aes(
    x = Fe_auth,
    y = P_auth,
    color = `Formation_label`
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.4
  ) +
  scale_color_manual(values = tol_colors) +
  scale_linetype_manual(values = shale_linetypes) +
  scale_y_continuous(breaks = seq(0, 0.08, by = 0.04 )) +
  scale_x_continuous(breaks = seq(0, 6, by = 2)) +
  geom_text(
    data = Archean_top6_label_data,
    aes(
      x = x,
      y = y,
      label = `Formation_label`,
      color = `Formation_label`
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +

  geom_smooth(aes(linetype = Formation_label),
    method = "lm",
    se = FALSE,
    linewidth = 1.75,
    lineend = "round"
  ) +
  
  labs(
    x = expression("Authigenic Fe (mol kg"^{-1}*")"),
    y = expression("Authigenic P (mol kg"^{-1}*")"),
    color = "Formation",
    title = "Archean Shales",
    tag = "C"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 11, hjust = 0.5),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1))

print(Top6_shaleplot)


### IF Data ###

# Remove ca from ages, choose midpoint of ranges, and convert from GA to MA

IFs <- IFs %>%
  mutate(
    Age_clean = str_remove(Age_Ga, "^ca\\.\\s*"),
    
    Age_Ga_num = if_else(
      str_detect(Age_clean, "-"),
      (
        parse_number(str_split_fixed(Age_clean, "-", 2)[,1]) +
          parse_number(str_split_fixed(Age_clean, "-", 2)[,2])
      ) / 2,
      parse_number(Age_clean)
    ),
    
    Age_Ma = Age_Ga_num * 1000
  )

# Classify as modern hydrothermal vent or Archean

IFs <- IFs %>%
  filter(
    !str_detect(Sample, regex("average", ignore_case = TRUE))
  ) %>%
  mutate(
    Age_group = case_when(
      Age_Ma > 2400 ~ "Archean IF",
      Age_Ma > 15 & Age_Ma <= 2400 ~ "Intermediate",
      Age_Ma <= 15 ~ "Hydrothermal",
      TRUE ~ NA_character_
    ),
    Age_group = factor(
      Age_group,
      levels = c(
        "Archean IF",
        "Intermediate",
        "Hydrothermal"
      )
    )
  )

# Convert to mol/kg

IFs <- IFs %>%
  mutate(
    Fe_molkg = (`Fe (%)` * 10) / AW_Fe,
    P_molkg  = (`P (%)` * 10) / AW_P,
    Al_molkg = (`Al (%)` * 10) / AW_Al
  )

#Group to nearest age in Archean UCC dataset

#Group Archean shale record to nearest age in Archean UCC
IFs<- IFs %>%
  rowwise() %>%
  mutate(
    ucc_age =
      ArcheanUCC_lookup$ucc_age[
        which.min(
          abs(ArcheanUCC_lookup$ucc_age - `Age_Ma`)
        )
      ]
  ) %>%
  ungroup()

#Join IF with Archean UCC based on nearest age

IFs <- IFs %>%
  left_join(
    ArcheanUCC_lookup,
    by = "ucc_age"
  )

# Calculate authigenic and detrital fractions of IFs
IFs <- IFs %>%
  mutate(
    Fe_det = Al_molkg * UCC_FeAl,
    P_det  = Al_molkg * UCC_PAl,
    
    Fe_auth = Fe_molkg - Fe_det,
    P_auth  = P_molkg - P_det,
    
    Fe_auth_pct = 100 * Fe_auth / Fe_molkg,
    P_auth_pct  = 100 * P_auth / P_molkg
  )

# Linear regressions of IFs

Hydrothermal_regression_table <- IFs %>%
  filter(Age_group == "Hydrothermal") %>%
  group_by(Formation_Location) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    model = map(data, ~ lm(P_molkg ~ Fe_molkg, data = .x)),
    glance = map(model, broom::glance),
    tidy = map(model, broom::tidy)
  ) %>%
  mutate(
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    # 95% CI for slope
    slope_ci_low = map_dbl(model, ~ confint(.x)["Fe_molkg", 1]),
    slope_ci_high = map_dbl(model, ~ confint(.x)["Fe_molkg", 2]),
    slope_95CI = paste0(
      round(slope_ci_low, 3),
      " to ",
      round(slope_ci_high, 3)
    ),
    residual_se = map_dbl(glance, "sigma"),
    r_squared = map_dbl(glance, "r.squared"),
    pearson_r = sign(slope) * sqrt(r_squared),
    p_value = map_dbl(
      tidy,
      ~ .x %>%
        filter(term == "Fe_molkg") %>%
        pull(p.value)
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    linetype = if_else(
      significance == "ns",
      "longdash",
      "solid"
    )
  ) %>%
  select(
    Formation_Location,
    n,
    slope,
    slope_95CI,
    intercept,
    r_squared,
    p_value,
    residual_se,
    pearson_r,
    significance,
    linetype
  ) %>%
  arrange(desc(n))

Hydrothermal_regression_table

ArcheanIF_regression_table <- IFs %>%
  filter(Age_group == "Archean IF") %>%
  group_by(`Formation_Location`) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow)
  ) %>%
  filter(n >= 3) %>%
  mutate(
    model = map(data, ~ lm(P_molkg ~ Fe_molkg, data = .x)),
    glance = map(model, broom::glance),
    tidy = map(model, broom::tidy)
  ) %>%
  mutate(
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    # 95% CI for slope
    slope_ci_low = map_dbl(model, ~ confint(.x)["Fe_molkg", 1]),
    slope_ci_high = map_dbl(model, ~ confint(.x)["Fe_molkg", 2]),
    slope_95CI = paste0(
      round(slope_ci_low, 3),
      " to ",
      round(slope_ci_high, 3)
    ),
    residual_se = map_dbl(glance, "sigma"),
    r_squared = map_dbl(glance, "r.squared"),
    pearson_r = sign(slope) * sqrt(r_squared),
    p_value = map_dbl(
      tidy,
      ~ filter(.x, term == "Fe_molkg") %>%
        pull(p.value)
    ),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    linetype = if_else(
      significance == "ns",
      "longdash",
      "solid"
    )
  ) %>%
  select(
    Formation_Location,
    n,
    slope,
    slope_95CI,
    intercept,
    r_squared,
    p_value,
    residual_se,
    pearson_r,
    significance,
    linetype
  ) %>%
  arrange(desc(n))

ArcheanIF_regression_table

# Filter the six Archean formations and top two hydrothermal with most sample data points

top6_archean <- IFs %>%
  filter(Age_group == "Archean IF") %>%
  count(`Formation_Location`, sort = TRUE) %>%
  slice_head(n = 6)

top6_archean

top6_names <- top6_archean$`Formation_Location`

IFs_top6_archean <- IFs %>%
  filter(`Formation_Location` %in% top6_names)

#Remove iron formation from location names
IFs_top6_archean <- IFs_top6_archean %>%
  mutate(
    Formation_label = `Formation_Location` %>%
      stringr::str_remove_all("[Ii]ron [Ff]ormation") %>%
      stringr::str_remove_all("[Bb]elt") %>%
      stringr::str_squish()
  )

top2_hydro <- IFs %>%
  filter(Age_group == "Hydrothermal") %>%
  count(`Formation_Location`, sort = TRUE) %>%
  slice_head(n = 2)

top2_hydro

top2_names <- top2_hydro$`Formation_Location`

Hydro_top2<- IFs %>%
  filter(`Formation_Location` %in% top2_names)

selected_formations <- c(
  top6_archean$`Formation_Location`,
  top2_hydro$`Formation_Location`
)

IFs_top <- IFs %>%
  filter(`Formation_Location` %in% selected_formations)

#label for plotting

label_data <- tibble(
  `Formation_label` = unique(IFs_top6_archean$`Formation_label`)
) %>%
  mutate(
    x = max(IFs_top6_archean$Fe_molkg, na.rm = TRUE) * 1.02,
    y = seq(
      max(IFs_top6_archean$P_molkg, na.rm = TRUE),
      max(IFs_top6_archean$P_molkg, na.rm = TRUE) * 0.725,
      length.out = n()
    )
  )

#Make linestypes for plot

IF_linetypes <- ArcheanIF_regression_table %>%
  select(`Formation_Location`, linetype)

IF_linetypes <- setNames(
  IF_linetypes$linetype,
  IF_linetypes$`Formation_Location`
)

Archean_IF_plot <- ggplot(
IFs_top6_archean,
  aes(
    x = Fe_molkg,
    y = P_molkg,
    color = `Formation_label`
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.4
  ) +
  scale_color_manual(values = tol_colors) +
  scale_linetype_manual(values = IF_linetypes) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.04)) +
  geom_text(
    data = label_data,
    aes(
      x = x,
      y = y,
      label = `Formation_label`,
      color = `Formation_label`
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 3
  ) +

  geom_smooth(aes(linetype = `Formation_Location`),
               method = "lm",
               se = FALSE,
               linewidth = 1.75, lineend = "round"
  ) +

  labs(
    x = expression("Fe (mol kg"^{-1}*")"),
    y = expression("P (mol kg"^{-1}*")"),
    color = "Formation",
    title = "Archean IFs",
    tag = "D"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1))

print(Archean_IF_plot)

hydro_label_data <- tibble(
  `Formation_Location` = unique(Hydro_top2$`Formation_Location`)
) %>%
  mutate(
    x = max(Hydro_top2$Fe_molkg, na.rm = TRUE) * 0,
    y = seq(
      max(Hydro_top2$P_molkg, na.rm = TRUE),
      max(Hydro_top2$P_molkg, na.rm = TRUE) * 0.94,
      length.out = n()
    )
  )

#Hydrothermal vent IF plot

HF_linetypes <- Hydrothermal_regression_table %>%
  select(`Formation_Location`, linetype)

HF_linetypes <- setNames(
  HF_linetypes$linetype,
  HF_linetypes$`Formation_Location`
)

Hydro_IF_plot <- ggplot(
  Hydro_top2,
  aes(
    x = Fe_molkg,
    y = P_molkg,
    color = `Formation_Location`
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.4
  ) +
  geom_text(
    data = hydro_label_data,
    aes(
      x = x,
      y = y,
      label = `Formation_Location`,
      color = `Formation_Location`
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +
  scale_color_manual(values = tol_colors) +
  scale_x_continuous(breaks = seq(0, 2.5, by =1)) +
  scale_linetype_manual(values = HF_linetypes) +
  geom_smooth( aes(linetype = `Formation_Location`),
    method = "lm",
    se = FALSE,
    linewidth = 1.75,
    lineend = "round"
  ) +
  
  labs(
    x = expression("Fe (mol kg"^{-1}*")"),
    y = expression("P (mol kg"^{-1}*")"),
    color = "Formation",
    title = "Hydrothermal IFs",
    tag = "B"
  ) +
  
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(size = 8),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1))

print(Hydro_IF_plot)

###Combined figure###

combined_plot_fig1 <-
  (lakes_plot | Hydro_IF_plot) /
  (Top6_shaleplot | Archean_IF_plot) 

combined_plot_fig1

ggsave(
  "Fe_P_comparison_figure.pdf",
  combined_plot_fig1,
  width = 180,
  height = 180,
  units = "mm",
  dpi = 300
)

#Add significance labels to regression tables
ArcheanShale_regression_table <- ArcheanShale_regression_table %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

ArcheanIF_regression_table <- ArcheanIF_regression_table %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

Hydrothermal_regression_table <- Hydrothermal_regression_table %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

Lake_regression_table <- Lake_regression_table %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

DL_no_outliers_regression_table <- DL_no_outliers_regression_table %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )


#Save regression tables

write.csv(
  ArcheanShale_regression_table,
  "ArcheanShale_regression_table.csv",
  row.names = FALSE
)

write.csv(
  ArcheanIF_regression_table,
  "ArcheanIF_regression_table.csv",
  row.names = FALSE
)

write.csv(
  Hydrothermal_regression_table,
  "Hydrothermal_regression_table.csv",
  row.names = FALSE
)

write.csv(
  Lake_regression_table,
  "Lake_regression_table.csv",
  row.names = FALSE
)

write.csv(
  DL_no_outliers_regression_table,
  "DL_no_outliers_regression_table.csv",
  row.names = FALSE
)


write_xlsx(
  list(
    Archean_Shale = ArcheanShale_regression_table,
    Archean_IF = ArcheanIF_regression_table,
    Hydrothermal = Hydrothermal_regression_table,
    Lake = Lake_regression_table,
    Deming_NoOutliers =  DL_no_outliers_regression_table
    
  ),
  "Regression_Tables.xlsx"
)
