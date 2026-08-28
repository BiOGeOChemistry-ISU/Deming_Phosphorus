# Deming Metals Data source: Rico et al. (2024)
# DOI: https://doi.org/10.1594/PANGAEA.970399
#Chlorophyll data source: Swanner et al. (2026)
#DOI: https://doi.org/10.6073/pasta/05b42443cdb707bef311e3da2a7892e9

# Load Packages
library(tidyverse)
library(readxl)
library(broom)
library(grid)
library(stringr)
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(scico)
library(fuzzyjoin)
library(patchwork)


#DL water data from Heard et al. 2024 DOI:10.1016/j.gca.2024.07.037
#And Zenodo Data Repository DOI:10.5281/zenodo.21953894
DL_water <- read.csv("082022_Deming_Water.csv")

# Chad Wittkop's MN surficial till https://doi.org/10.1016/j.chemgeo.2019.119390
CW_norm <- read_csv("CW_norm.csv")

itasca <- read.csv("Itasca_database_2026.csv")


################################### Data wrangling ####################################

# atomic weights for molar conversions
AW_Al = 26.982
AW_P = 30.974
AW_Ti = 47.867
AW_Fe = 55.845

###Deming Lake Water Column###

#Convert dissolved to molar

DL_water <- DL_water %>%
  mutate(
    dFe_umol = dFe..ng.ml. / (AW_Fe),
    pFe_umol = pFe..ng.ml. / (AW_Fe),
    dP_umol  = dP..ng.ml.  / (AW_P),
    pP_umol  = pP..ng.ml.  / (AW_P),
    pTi_umol  = pTi..ng.ml.  / (AW_Ti),
  ) %>%
  #Rename columns
  rename(
    Depth = Depth..m.)


chl <- itasca %>%
  filter(
    Lake == "Deming",
    Date == "2022-08-09",
    Device == "YSI ProDSS",
    Measurement == "Chlorophyll"
  ) %>%
  select(
    chl_depth = Depth,
    Chlorophyll = Value
  )

### Calculate authigenic portions of particulates ###
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

DL_water <- DL_water %>%
  mutate(
    pFe_det = pTi_umol * FeTi_till,
    pP_det  = pTi_umol * PTi_till,
    
    pFe_auth = pFe_umol - pFe_det,
    pP_auth  = pP_umol  - pP_det,
    
    pFe_auth_pct = 100 * pFe_auth / pFe_umol,
    pP_auth_pct  = 100 * pP_auth  / pP_umol
  ) 


### Panel A ###
col_dP  <- "#332288"  
col_dFe <- "#882255"

#Define layer depths and labels 
layer_labels <- data.frame(
  Depth = c(2, 7.5, 11.5),
  label = c("L1", "L2", "L3"),
  col = c("#44AA99", "#CC6677", "#88CCEE")  # replace with your colors
)

#arrange by depth for flipped y-axis (Depth 0 at top)
DL_water <- DL_water %>%
  arrange(Depth)

#Figure plot

Fig2a <- ggplot() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = 5.75,
    fill = "#44AA99",
    alpha = 0.08
  )+
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 5.75,
    ymax = 10,
    fill = "#CC6677",
    alpha = 0.08
  )+
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 10,
    ymax = Inf,
    fill = "#88CCEE",
    alpha = 0.08
  )+
  # dP
  geom_path(
    data = DL_water %>% filter(!is.na(dP_umol)),
    aes(x = dP_umol, y = Depth, color = "dP_umol"),
    linewidth = 0.8,
    linetype = "solid"
  ) + 
  geom_point(
    data = DL_water %>% filter(!is.na(dP_umol)),
    aes(x = dP_umol, y = Depth, color = "dP_umol", shape = "dP_umol"),
    size = 2.8,  alpha = 0.9
  ) +
  
  # dFe 
  geom_path(
    data = DL_water %>% filter(!is.na(dFe_umol)),
    aes(x = dFe_umol * scale_factor, y = Depth, color = "dFe_umol"),
    linetype = "solid",
    linewidth = 0.8
  ) +
  geom_point(
    data = DL_water %>% filter(!is.na(dFe_umol)),
    aes(x = dFe_umol * scale_factor, y = Depth,
        color = "dFe_umol", shape = "dFe_umol"),
    size = 2.8,  alpha = 0.9
  ) +
  geom_hline(yintercept = c(5.75, 10),
             linetype = "dashed",
             color = "black",
             linewidth = 0.2) +
  
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = Depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  # Reverse y axis
  scale_y_reverse(
    breaks = seq(0, 16, by = 2),
    name = "Depth (m)"
  ) +
  #Double x axis (dFe top and dP bottom)
  scale_x_continuous(
    name = expression(paste("Dissolved P (", mu, "M)")),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = expression(paste("Dissolved Fe (", mu, "M)"))
    )
  ) +

  # dP and dFe Colors
  scale_color_manual(
    values = c("dP_umol" = col_dP, "dFe_umol" = col_dFe),
    labels = c("dP_umol" = expression("P"), "dFe_umol" = expression("Fe")),
    name = NULL
  ) +
  #Shapes
  scale_shape_manual(
    values = c("dP_umol" = 16,   # filled circle
               "dFe_umol" = 17), # filled triangle
    labels = c("dP_umol" = expression("P"), "dFe_umol" = expression("Fe")),
    name = NULL
  ) +
  labs(tag = "A") +
  theme_bw() +
  theme(
    legend.position = c(0.3, 0.8),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1)
  )

print(Fig2a)

#save figure
ggsave("Fig2a.pdf", Fig2a, width = 3.35, height = 5)

### Panel B ###

#grouping chl and particulate data by depth using a nearest join depth
#of 0.15 m. since depths do not perfectly match between datasets
DL_dat_joined <- difference_left_join(
  DL_water,
  chl,
  by = c("Depth" = "chl_depth"),
  max_dist = 0.15
)

# Define primary (bottom) axis using pP & pFe, and compute Chl transform

prim_min <- min(c(DL_dat_joined$pP_auth, DL_dat_joined$pFe_auth), na.rm = TRUE)
prim_max <- max(c(DL_dat_joined$pP_auth, DL_dat_joined$pFe_auth), na.rm = TRUE)

chl_min  <- min(DL_dat_joined$Chlorophyll, na.rm = TRUE)
chl_max  <- max(DL_dat_joined$Chlorophyll, na.rm = TRUE)

if (!is.finite(prim_min) || !is.finite(prim_max) || prim_min == prim_max) {
  stop("pP_auth and pFe_auth must have a non-zero range.")
}
if (!is.finite(chl_min) || !is.finite(chl_max) || chl_min == chl_max) {
  stop("Chlorophyll must have a non-zero range.")
}

# Affine mapping: Chlorophyll -> primary scale (µM)
scale_chl  <- (prim_max - prim_min) / (chl_max - chl_min)
offset_chl <- prim_min - scale_chl * chl_min

# Limits to show everything
# (We’ll recompute later using the interpolated series, but set a provisional version now)
x_lims <- c(prim_min, prim_max)

# Interpolate Chl to fill gaps

dat_interp <- DL_dat_joined %>%
  arrange(Depth)

ok <- complete.cases(dat_interp$Depth, dat_interp$Chlorophyll)
if (sum(ok) >= 2) {
  dat_interp$Chlorophyll_filled <- approx(
    x = dat_interp$Depth[ok],
    y = dat_interp$Chlorophyll[ok],
    xout = dat_interp$Depth,
    method = "linear",
    ties = "ordered"
  )$y
} else {
  # Fallback if too few valid points to interpolate
  dat_interp$Chlorophyll_filled <- dat_interp$Chlorophyll
}

# Transform the FILLED Chl for plotting on the primary axis
dat_interp <- dat_interp %>%
  mutate(Chl_trans = offset_chl + scale_chl * Chlorophyll_filled)

# Now recompute x-limits to ensure the ribbon fits fully
x_lims <- c(
  min(prim_min, min(dat_interp$Chl_trans, na.rm = TRUE)),
  max(prim_max, max(dat_interp$Chl_trans, na.rm = TRUE))
)

# Ensure ordered by depth right before plotting
dat <- DL_dat_joined %>% arrange(Depth)

# Colors
col_pP  <- "#332288"  
col_pFe <- "#882255"
col_Chl <- "#117733"  

layer_labels <- data.frame(
  Depth = c(3, 8, 12),
  label = c("L1", "L2", "L3"),
  col = c("#44AA99", "#CC6677", "#88CCEE")  # replace with your colors
)

Fig2b <- ggplot() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = 5.75,
    fill = "#44AA99",
    alpha = 0.08
  )+
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 5.75,
    ymax = 10,
    fill = "#CC6677",
    alpha = 0.08
  )+
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 10,
    ymax = Inf,
    fill = "#88CCEE",
    alpha = 0.08
  )+
  geom_ribbon(
    data = dat_interp,
    aes(y = Depth, xmin = ribbon_baseline, xmax = Chl_trans, fill = "Chl"),
    alpha = 0.3,
    color = NA,
    na.rm = TRUE,
    key_glyph = "rect"
  ) +
  geom_path(
    data = dat %>% dplyr::filter(!is.na(pP_auth)),
    aes(x = pP_auth, y = Depth, color = "pP", group = 1),
    linewidth = 0.8, na.rm = TRUE
  ) + 
  geom_point(
    data = dat %>% dplyr::filter(!is.na(pP_auth)),
    aes(x = pP_auth, y = Depth, color = "pP", shape = "pP"),
    size = 2.8, na.rm = TRUE,
    alpha = 0.9
  ) +
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = Depth, label = label),
    hjust = 1.5,              # pushes slightly inside panel
    size = 3.25,
    color = layer_labels$col
  ) +
  
  geom_point(
    data = dat %>% dplyr::filter(!is.na(pFe_auth)),
    aes(x = pFe_auth, y = Depth, color = "pFe", shape = "pFe"),
    size = 2.8, na.rm = TRUE,
    alpha = 0.9
  ) +
  geom_path(
    data = dat %>% dplyr::filter(!is.na(pFe_auth)),
    aes(x = pFe_auth, y = Depth, color = "pFe", group = 1),
    linewidth = 0.8, na.rm = TRUE
  ) +
  geom_hline(yintercept = c(5.75, 10),
             linetype = "dashed",
             color = "black",
             linewidth = 0.2) +
  # Axes
  scale_y_reverse(
    breaks = seq(0, 16, by = 2),
    name = "Depth (m)"
  ) +
  scale_x_continuous(
    name = expression(paste("Authigenic pP & pFe (", mu, "M)")),
    sec.axis = sec_axis(~ (. - offset_chl) / scale_chl, 
                        name = expression(paste("Chlorophyll (", mu, "g ", L^{-1}, ")"))),
    breaks = pretty(x_lims, n = 5)
  ) +
  
  scale_color_manual(
    values = c("pP" = col_pP, "pFe" = col_pFe),
    labels = c("pP" = expression("P"), "pFe" = expression("Fe")),
    name = NULL        # <-- same title as fill
  ) +
  scale_fill_manual(values = c("Chl" = col_Chl)) +
  scale_shape_manual(
    values = c("pP" = 16, "pFe" = 17),  
    labels = c("pP" = expression("P"), "pFe" = expression("Fe")),
    name = NULL
  ) +
  guides(
    fill  = guide_legend(order = 1, override.aes = list(alpha = 0.5)),
    color = guide_legend(order = 2, override.aes = list(fill = NA)),
    shape = guide_legend(order = 2),
    guide = guide_legend(
      override.aes = list(
        alpha = 0.3,
        colour = NA
      )
    )
  ) +
  labs(tag = "B") +
  theme_bw() +
  theme(
    legend.position     = "right",
    legend.box.spacing  = grid::unit(0, "pt"),   # remove extra space between boxes
    legend.spacing.y    = grid::unit(0, "pt"),   # tighten vertical spacing between keys
    legend.key.height   = grid::unit(8, "pt"),
    legend.key.width    = grid::unit(14, "pt"),
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1),
    axis.text = element_text(size = 8)
     )

print(Fig2b)

### Panel C ###

#Add and define layers

DL_water  <- DL_water %>%
  mutate(
    L = case_when(
      Depth <= 5.75 ~ 1,
      Depth <= 10 ~ 2,
      TRUE ~ 3
    ),
    L = factor(L, levels = c(1, 2, 3))
  )

palette3 <- c("#44AA99", "#CC6677", "#88CCEE")

# --- Linear regression and format table ---

Water_column_regression_table <-DL_water %>%
  group_by(L) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    model = map(data, ~ lm(pP_auth ~ pFe_auth, data = .x)),
    glance = map(model, broom::glance),
    tidy = map(model, broom::tidy)
  ) %>%
  mutate(
    intercept = map_dbl(model, ~ coef(.x)[1]),
    slope = map_dbl(model, ~ coef(.x)[2]),
    # 95% CI for slope
    slope_ci_low = map_dbl(model, ~ confint(.x)["pFe_auth", 1]),
    slope_ci_high = map_dbl(model, ~ confint(.x)["pFe_auth", 2]),
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
        filter(term == "pFe_auth") %>%
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
    L,
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
  arrange(desc(L))

Water_column_regression_table 

r2_dat_DL_Aug2022_Water <- DL_water %>%
  group_by(L) %>%
  summarise(
    R2   = summary(lm(pP_auth ~ pFe_auth))$r.squared,
    n    = n(),
    pval = summary(lm(pP_auth ~ pFe_auth))$coefficients[2,4]
  ) %>%
  ungroup() %>%
  arrange(L) %>%
  mutate(
    line_type = if_else(pval <= 0.05, "solid", "longdash"),
    lab = paste0(
      "L", L,
      "\nR² = ", formatC(R2, digits = 4, format = "f"),
      "\np = ", ifelse(pval < 0.05, "< 0.05",
                       formatC(pval, digits = 3, format = "f"))
    ),
    x_label = c(1.6, 3.4, 5.2),
    y_label = 3.9
  )

DL_Aug2022_Water_plot <- DL_Aug2022_Water %>%
  left_join(
    r2_dat_DL_Aug2022_Water %>% select(L, line_type),
    by = "L"
  )


#plot figure
Fig2C <- ggplot(DL_Aug2022_Water_plot, aes(pFe_auth, pP_auth, color = L)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_smooth(aes(linetype = line_type), method = "lm", se = FALSE, linewidth = 1.25, lineend = "round") +
  geom_label(
    data = r2_dat_DL_Aug2022_Water,
    aes(x = x_label, y = y_label,
      label = lab,
      color = L,
    ),
    fill = alpha("white", 0.85),
    label.size = 0.5,
    size = 3,
    fontface = "bold",
    label.padding = unit(0.4, "lines"),
    show.legend = FALSE,
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 1.05
  ) +
  
  scale_color_manual(values = palette3, name = "L") +
  scale_y_continuous(breaks = seq(0,4, by = 1)) +
  guides(color = "none") +
  labs(
    x = expression(paste("Authigenic pFe (", mu, "M)")),
    y = expression(paste("Authigenic pP (", mu, "M)"))
  ) +
  scale_linetype_identity() +
  guides(
    color = "none",
    linetype = "none"
  )+
  labs(tag = "C") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0,1),
        strip.text = element_text(face = "bold"))

print(Fig2C)

#save figure
ggsave("Figure3.pdf",
       plot = Figure3, width = 3.35, height = 3.35, device = "pdf")




#combined figure

###Combined figure###
Fig2a <- Fig2a +
  theme(legend.position = "none")

combined_plot_fig2ab <-
  (Fig2a| Fig2b) + plot_layout(guides = "collect") 

combined_plot_fig2 <-
  combined_plot_fig2ab /
  ((plot_spacer() | Fig2C | plot_spacer()) +
  plot_layout(widths = c(1, 6, 1)))


combined_plot_fig2

#save figure
ggsave("combined_plot_fig2.pdf", combined_plot_fig2,   width = 180,
       height = 210,
       units = "mm",
       dpi = 300)

write.csv(
  Water_column_regression_table,
  "Water_column_regression_table.csv",
  row.names = FALSE
)
