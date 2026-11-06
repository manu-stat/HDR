# Import ------------------------------------------------------------------
source("1 HDR measures.r")
source("2 Main HDR.r")
source("4. HDR error.r")

measures = c("KDE", "KNN", "NP_Cop", "P_Cop")
sample_sizes = c(50, 100, 250, 500, 1000)
mc = 100

# Error 2D ----------------------------------------------------------------
Error_2D = do.call(rbind, c(setNames(lapply(sample_sizes, function(n) {
  MC_error(MC = mc, size = n, scenarios_file = "3a. Scenarios 2D.r", measures = measures)$Summary
  }), sample_sizes), make.row.names = FALSE))

# --- Error Plot ---
library(dplyr); library(ggplot2)
lapply(unique(Error_2D$Scenario), function(nm){
  fig = Error_2D %>% filter(Scenario == nm & Metric != "time") %>% ggplot(aes(x = N, y = Mean, colour = Measure, group = Measure)) + 
    geom_line(aes(linetype = Measure), position = position_dodge(width = 20)) + 
    geom_point(position = position_dodge(width = 20)) + 
    labs(title = nm) + 
    theme_minimal() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), legend.position = "bottom") +
    facet_wrap(~Metric, scale = "free_y")
  plot(fig)
  readline("Next plot: ENTER")
})


# Error 3D ----------------------------------------------------------------
Error_3D = do.call(rbind, c(setNames(lapply(sample_sizes, function(n) {
  MC_error(MC = mc, size = n, scenarios_file = "3b. Scenarios 3D.r", measures = measures)$Summary
}), sample_sizes), make.row.names = FALSE))

# --- Error Plot ---
library(dplyr); library(ggplot2)
lapply(unique(Error_2D$Scenario), function(nm){
  fig = Error_3D %>% filter(Scenario == nm & Metric != "time") %>% ggplot(aes(x = N, y = Mean, colour = Measure, group = Measure)) + 
    geom_line(aes(linetype = Measure), position = position_dodge(width = 20)) + 
    geom_point(position = position_dodge(width = 20)) + 
    labs(title = nm) + 
    theme_minimal() + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), legend.position = "bottom") +
    facet_wrap(~Metric, scale = "free_y")
  plot(fig)
  readline("Next plot: ENTER")
})



# Save -----
save("Error_2D", "Error_3D", file = "Error_mc100.RData")

