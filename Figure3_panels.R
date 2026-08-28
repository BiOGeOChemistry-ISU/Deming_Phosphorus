# Deming Metals Data source: Rico et al. (2024)
# DOI: https://doi.org/10.1594/PANGAEA.970399

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

################################### Data import ####################################

#DL Sediment data from Heard et al. 2024 DOI:10.1016/j.gca.2024.07.037
#And Zenodo Data Repository DOI:10.5281/zenodo.21953894
DL_sediments <- read.csv("2022_Deming_Sediment.csv")


# Chad Wittkop's MN surficial till https://doi.org/10.1016/j.chemgeo.2019.119390
CW_norm <- read_csv("CW_norm.csv")


# Total carbon and nitrogen data for sediment transect from Ostrander et al 2025
#https://doi.org/10.1016/j.gca.2024.12.008

DLtoc_tn <- read.csv("TOC_TN_data.csv", stringsAsFactors = FALSE)

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

#Assign layers

DL_sediments <- DL_sediments %>%
  mutate(
    layer = case_when(
      Lake.Depth..m. <= 5 ~ "L1",
      Lake.Depth..m. > 5 &  Lake.Depth..m. <= 10 ~ "L2",
      Lake.Depth..m. > 10 &  Lake.Depth..m. <= 17 ~ "L3"
    )
  )

#Regression
DL_model <- lm(P_auth ~ Fe_auth, data = DL_sediments)

DL_stats <- tibble(
  R2 = summary(DL_model)$r.squared,
  n = nrow(DL_FP),
  pval = summary(DL_model)$coefficients[2,4]
)

#Annotation Label

DL_stats <- DL_stats %>%
  mutate(
    lab = paste0(
      "R² = ", formatC(R2, digits = 4, format = "f"),
      "\np = ",
      ifelse(
        pval < 0.001,
        "< 0.001",
        formatC(pval, digits = 3, format = "f")
      )
    )
  )

palette3 <- c("#44AA99", "#CC6677", "#88CCEE")

# Regression plot
Fig3D <- ggplot(
  DL_sediments,
  aes(x = Fe_auth, y = P_auth, color = layer)
) +
  geom_point(size = 3.25, alpha = 0.8) +
  
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    linewidth = 1.25,
    color = "#544c4a",
    lineend = "round",
    alpha = 0.8
  ) +
  
  geom_label(
    data = DL_stats,
    aes(
      x = 0.2,
      y = 0.045,
      label = lab
    ),
    fill = alpha("white", 0.85),
    label.size = 0.5,
    size = 3,
    fontface = "bold",
    label.padding = unit(0.6, "lines"),
    show.legend = FALSE,
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 1.0,
    color = "#544c4a"
  ) +
  
  scale_color_manual(
    values = palette3,
    breaks = c("L1", "L2", "L3"),
    name = "Layer"
  ) +
  
  labs(tag = "D",
    x = expression(paste("Authigenic Fe (mol ", kg^{-1}, ")")),
    y = expression(paste("Authigenic P (mol ", kg^{-1}, ")"))
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = c(0.15, 0.8),
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      color = "#544c4a"
    ),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(color = "black", size = 8),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1)
  )

Fig3D

#DL Transect

DL_transect <- DL_sediments %>%
  filter(Sample.Type == "Transect", !is.na(Sample.Type))

### Panel A ###
col_pP  <- "#332288"  
col_pFe <- "#882255"

#Define layer depths and labels 
layer_labels <- data.frame(
  depth = c(2, 7.5, 11.5),
  label = c("L1", "L2", "L3"),
  col = c("#44AA99", "#CC6677", "#88CCEE")
)

#arrange by depth for flipped y-axis (Depth 0 at top)
DL_transect <- DL_transect %>%
  arrange( Lake.Depth..m.)

#Figure plot

Fig3a <- ggplot() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = 5,
    fill = "#44AA99",
    alpha = 0.08
  )+
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 5,
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
  # pP
  geom_path(
    data = DL_transect %>% filter(!is.na(P_auth)),
    aes(x = P_auth, y =  Lake.Depth..m.),
    color = "#332288",
    linewidth = 0.8,
    linetype = "solid"
  ) + 
  geom_point(
    data = DL_transect %>% filter(!is.na(P_auth)),
    aes(x = P_auth, y =  Lake.Depth..m.),
   color = "#332288", shape = 16,
    size = 2.8,  alpha = 0.9
  ) +
  
  geom_hline(yintercept = c(5, 10),
             linetype = "dashed",
             color = "black",
             linewidth = 0.2) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
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
    name = expression(paste("Authigenic P (mol ", kg^{-1}, ")"))
  ) +
  labs(tag = "A") +
  theme_bw() +
  theme(
    legend.text = element_blank(),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1)
  )


print(Fig3a)
  

Fig3b <- ggplot() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = 5,
    fill = "#44AA99",
    alpha = 0.08
  )+
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 5,
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
  # pP
  geom_path(
    data = DL_transect %>% filter(!is.na(Fe_auth)),
    aes(x = Fe_auth, y =  Lake.Depth..m.),
    color = "#882255",
    linewidth = 0.8,
    linetype = "solid"
  ) + 
  geom_point(
    data = DL_transect %>% filter(!is.na(Fe_auth)),
    aes(x = Fe_auth, y =  Lake.Depth..m.),
    color = "#882255", shape = 17,
    size = 2.8,  alpha = 0.9
  ) +
  
  geom_hline(yintercept = c(5, 10),
             linetype = "dashed",
             color = "black",
             linewidth = 0.2) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  # Reverse y axis
  scale_y_reverse(
    breaks = seq(0, 16, by = 2)
  ) +
  scale_x_continuous(
    name = expression(paste("Authigenic Fe (mol ", kg^{-1}, ")"))
  ) +
  labs(tag = "B") +
  theme_bw() +
  theme(
    legend.text = element_blank(),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(color = NA),
    axis.title.x  = element_text(size = 10),
    axis.title.y = element_blank(),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1)
  )


Fig3b

##TOC profile

DLtoc_tn_clean <- DLtoc_tn %>%
  rename(
    Depth_m = Depth,
    TOC_pct = TOC
  ) %>%
  filter(X != "Sample") %>%   # remove header row that was imported as data
  mutate(
    Depth_m = as.numeric(Depth_m),
    TOC_pct = as.numeric(TOC_pct)
  ) %>%
  filter(!is.na(TOC_pct))

Fig3c <- ggplot() +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = 5,
    fill = "#44AA99",
    alpha = 0.08
  )+
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = 5,
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
  # pP
  geom_path(
    data = DLtoc_tn_clean,
    aes(x = TOC_pct, y = Depth_m),
    color = "#117733",
    linewidth = 0.8,
    linetype = "solid"
  ) + 
  geom_point(
    data = DLtoc_tn_clean,
    aes(x = TOC_pct, y = Depth_m),
    color = "#117733", shape = 18,
    size = 2.8,  alpha = 0.9
  ) +
  
  geom_hline(yintercept = c(5, 10),
             linetype = "dashed",
             color = "black",
             linewidth = 0.2) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  # Reverse y axis
  scale_y_reverse(
    breaks = seq(0, 16, by = 2)
  ) +
  scale_x_continuous(
    name = "TOC %") +
  labs(tag = "C") +
  theme_bw() +
  theme(
    legend.text = element_blank(),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(color = NA),
    axis.title.x  = element_text(size = 10),
    axis.title.y = element_blank(),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1)
  )

Fig3c


combined_plot_fig3abc <-
  (Fig3a| Fig3b | Fig3c)

combined_plot_fig3 <-
  combined_plot_fig3abc /
  ((plot_spacer() | Fig3D | plot_spacer()) +
     plot_layout(widths = c(1, 6, 1)))


combined_plot_fig3

#save figure
ggsave("combined_plot_fig3.pdf", combined_plot_fig3,   width = 180,
       height = 210,
       units = "mm",
       dpi = 300)

