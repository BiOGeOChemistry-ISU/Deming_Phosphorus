library(readxl)
library(tidyverse)
library(patchwork)

# Data extracted from EDI repository Swanner 2026 DOI:https://doi.org/10.6073/pasta/05b42443cdb707bef311e3da2a7892e9

#Load data
df <- read_excel("Aug_YSI.xlsx") %>%
  rename(
    DO = `ODO (mg/L)`,
    Temp = `Temp (C)`,
    SpC = `Sp Cond (µS/cm)`
  )


POC_data <- read_excel("Aug_2023_POC.xlsx")

# --- Individual plots ---

layer_labels <- data.frame(
  depth = c(2, 7.5, 11.5),
  label = c("L1", "L2", "L3"),
  col = c("#44AA99", "#CC6677", "#88CCEE")
)


p_pH <- ggplot(df, aes(x = pH, y = Depth)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = 5.75,
    fill = "#44AA99", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 5.75, ymax = 10,
    fill = "#CC6677", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = Inf,
    fill = "#88CCEE", alpha = 0.08
  ) +
  geom_path( color = "#363636",linewidth = 0.8,
             linetype = "solid") +
  geom_point(color = "#363636", shape = 16,
                size = 2.8,  alpha = 0.9) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  scale_y_reverse(limits = c(14,0),breaks = seq(0, 16, by = 2)) +
  geom_hline(
    yintercept = c(5.75, 10),
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  labs(x = "pH", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
  )

p_pH

p_DO <- ggplot(df, aes(x = DO, y = Depth)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = 5.75,
    fill = "#44AA99", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 5.75, ymax = 10,
    fill = "#CC6677", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = Inf,
    fill = "#88CCEE", alpha = 0.08
  ) +
  geom_path( color = "#363636",linewidth = 0.8,
             linetype = "solid") +
  geom_point(color = "#363636", shape = 16,
             size = 2.8,  alpha = 0.9) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  scale_y_reverse(limits = c(14,0),breaks = seq(0, 16, by = 2)) +
  geom_hline(
    yintercept = c(5.75, 10),
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  labs(
    x = expression(paste("DO (mg ", L^{-1}, ")")),
    y = "Depth (m)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
  )

p_DO 

p_Temp <- ggplot(df, aes(x = Temp, y = Depth)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = 5.75,
    fill = "#44AA99", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 5.75, ymax = 10,
    fill = "#CC6677", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = Inf,
    fill = "#88CCEE", alpha = 0.08
  ) +
  geom_path( color = "#363636",linewidth = 0.8,
             linetype = "solid") +
  geom_point(color = "#363636", shape = 16,
             size = 2.8,  alpha = 0.9) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 3,
    size = 3.25,
    color = layer_labels$col
  ) +
  scale_y_reverse(limits = c(14,0),breaks = seq(0, 16, by = 2)) +
  geom_hline(
    yintercept = c(5.75, 10),
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  labs(x = "Temperature (°C)", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
  )

p_Temp

p_SpC <- ggplot(df, aes(x = SpC, y = Depth)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = 5.75,
    fill = "#44AA99", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 5.75, ymax = 10,
    fill = "#CC6677", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = Inf,
    fill = "#88CCEE", alpha = 0.08
  ) +
  geom_path( color = "#363636",linewidth = 0.8,
             linetype = "solid") +
  geom_point(color = "#363636", shape = 16,
             size = 2.8,  alpha = 0.9) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  scale_y_reverse(limits = c(14,0),breaks = seq(0, 16, by = 2)) +
  geom_hline(
    yintercept = c(5.75, 10),
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  labs(
    x = expression(paste("Sp Conductance (", mu, "S/cm)")),
    y = "Depth (m)"
  )  +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
  )

p_SpC


p_PAR <- ggplot(df, aes(x = PAR, y = Depth)) +
  annotate(
    "rect",
    xmin = 0, xmax = Inf,
    ymin = -Inf, ymax = 5.75,
    fill = "#44AA99", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = 0, xmax = Inf,
    ymin = 5.75, ymax = 10,
    fill = "#CC6677", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = 0, xmax = Inf,
    ymin = 10, ymax = Inf,
    fill = "#88CCEE", alpha = 0.08
  ) +
  geom_path( color = "#363636",linewidth = 0.8,
             linetype = "solid") +
  geom_point(color = "#363636", shape = 16,
             size = 2.8,  alpha = 0.9) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  scale_x_log10(
    breaks = scales::breaks_log(n = 6),
    labels = scales::label_number()
  )+
  scale_y_reverse(limits = c(14,0),breaks = seq(0, 16, by = 2)) +
  geom_hline(
    yintercept = c(5.75, 10),
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  labs(
    x = expression(paste("PAR (", mu, "mol photons ", m^{-2}, " ", s^{-1}, ")")),
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
  )

p_PAR 
  
p_POC <- ggplot(POC_data, aes(x = POC, y = Depth)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -Inf, ymax = 5.75,
    fill = "#44AA99", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 5.75, ymax = 10,
    fill = "#CC6677", alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = 10, ymax = Inf,
    fill = "#88CCEE", alpha = 0.08
  ) +
  geom_path( color = "#363636",linewidth = 0.8,
             linetype = "solid") +
  geom_point(color = "#363636", shape = 16,
             size = 2.8,  alpha = 0.9) +
  # Layer labels
  geom_text(
    data = layer_labels,
    aes(x = Inf, y = depth, label = label),
    inherit.aes = FALSE,
    hjust = 1.5,
    size = 3.25,
    color = layer_labels$col
  ) +
  scale_y_reverse(limits = c(14,0),
    breaks = seq(0, 16, by = 2)) +
  geom_hline(
    yintercept = c(5.75, 10),
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  labs(
    x = expression(paste("POC (", mu, "M)")),
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10),
  )

p_POC

# --- Combine horizontally ---
SupFig1a<-  p_SpC + p_Temp  + p_pH

SupFig1B <- p_DO + p_PAR + p_POC

SupFig1 <- SupFig1a/SupFig1B

SupFig1

ggsave("SupFig1.pdf", SupFig1,   width = 180,
       height = 180,
       units = "mm",
       dpi = 300)
