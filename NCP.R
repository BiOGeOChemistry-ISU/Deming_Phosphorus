library(tidyverse)

# Data from Iowa State Figshare Chamberlain and Swanner 2025 DOI: 10.25380/iastate.29109704.v1

dat <- read_csv("C_Fixatoin_long.csv")

nep <- dat %>%
  filter(Variable == "NEP") %>%
  mutate(
    Depth = factor(as.character(Depth),
                   levels = c("7", "5.75", "3.5"))
  )

NCP_plot <- ggplot(nep, aes(x = Value, y = Depth)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(height = 0.01, width = 0, size = 3, alpha = 0.8) +
  labs(
    x = expression(NCP ~ (mu*mol~C~L^{-1}~d^{-1})),
    y = "Depth (m)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text = element_text(size = 8),
    axis.title  = element_text(size = 10)
  )
  

NCP_plot

#save figure
ggsave("NCP_plot.pdf", NCP_plot,   width = 90,
       height = 90,
       units = "mm",
       dpi = 300)

