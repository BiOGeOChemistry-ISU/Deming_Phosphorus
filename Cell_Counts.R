library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(dbplyr)
library(scales)
require(scales)


# Data from Chamberlain and Swanner 2025. DOI: https://doi.org/10.25380/iastate.27629178
cellcounts <- read.csv("CellCounts.csv")

ps_cellcounts <- cellcounts %>%
  filter(Variable %in% c("Photosynthetic", "Non-Photosynthetic"))

ps_cellcounts$Depth <- factor(
  ps_cellcounts$Depth,
  levels = rev(sort(unique(ps_cellcounts$Depth)))
)

cellcount_plot <- ggplot(
  ps_cellcounts,
  aes(x = Depth, y = Value, color = Variable)
) +
  
  geom_point(
    size = 4, alpha = 0.8
  ) +
  
  coord_flip() +
  
  scale_color_manual(
    values = c(
      "Photosynthetic" = "#117733",
      "Non-Photosynthetic" = "#88CCEE"
    )
  ) +
  
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(
    x = "Depth (m)",
    y = expression(Cells~mL^-1)
  ) +
  
  theme_bw() +
  theme(
    legend.position = c(0.32, 0.9),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.margin = margin(2, 2, 2, 2),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

cellcount_plot

ggsave(
  "cellcount_plot.pdf",
  plot = cellcount_plot,
  width = 90,
  height = 90,
  units = "mm", dpi = 300
)


total_cellcounts <- filter(cellcounts, Variable == "Total")
total_cellcounts$Depth <- factor(total_cellcounts$Depth, levels = rev(sort(unique(total_cellcounts$Depth))))

total_cellcount_plot <- ggplot(total_cellcounts, aes(x = Depth, y = Value)) + 
  geom_point(color = "black", size = 3) +
  geom_errorbar(aes(ymin = Value - stderror, ymax = Value + stderror), width = .02) +
  theme(
    plot.title = element_blank(), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 14), 
    panel.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white"),
    axis.title = element_text(size = 12),
    axis.line = element_line(color = "black")
  ) +
  labs(x = "Depth (m)", y = expression("Cells mL"^-1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_flip()

# Display the plot
print(total_cellcount_plot)

ggsave("total_cellcount_plot.tiff", width = 3, height = 3, device = "tiff")
total_cellcount_plot
dev.off()
dev.off()
