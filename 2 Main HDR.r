HDR = function(data, measure = "KDE", coverage_prob = 0.95, margin = "norm", k = NULL, H = NULL, eps = NULL, CDF = NULL, plot = TRUE, ...){
  
  library(ggplot2); library(plotly)
  
  n = dim(data)[1]; d = dim(data)[2]
  data = data.frame(data); colnames(data) = paste0("X", 1:d)
  
  # density estimation 
  data$f_joint = switch(measure,
                        "KDE"        = KDE(data, H),
                        "KNN"        = KNN(data, k),
                        "CDF_KNN"    = CDF_KNN(data, k, CDF),
                        "CDF_e"      = CDF_e(data, eps),
                        "NP_Cop"     = NP_Cop(data),
                        "P_Cop"      = P_Cop(data, margin, ...),
                        "CDF_P_Cop"  = CDF_P_Cop(data, margin, eps, ...),
                        
                        stop("Measure not specified or invalid. \n Choose: KDE, KNN, CDF_KNN, CDF_e, NP_Cop, P_Cop, CDF_P_Cop")
  )
  
  # quantile-approach
  alpha = quantile(data$f_joint, 1-coverage_prob, na.rm = TRUE)
  data$IN = data$f_joint > alpha
  
  # figure
  if(plot & d == 2){
    fig = ggplot(data = data, aes(x = X1, y = X2, color = IN)) + 
      geom_point(size = 2, alpha = 0.8) + 
      scale_color_manual(values = c("#42314A", "#7FC6A4"), labels = c("Outside", "Inside")) +
      labs(x = "X", y = "Y", title = bquote("Estimated "*.(coverage_prob*100)*"% HDR with "*.(measure)), color = "HDR") + 
      theme_minimal() + theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), legend.position = "bottom") 
  } else if(plot & d == 3){
    fig = plot_ly(data = data, x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers", color = ~IN, colors = c("#42314A", "#7FC6A4"), 
                    marker = list(size = 4, sizemode = "diameter")) %>% 
      layout(title = list(text = paste0("Estimated ", coverage_prob*100, "% HDR with ", measure), font = list(size = 16), x = 0.5, y = 0.95, xanchor = "center"),
             scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")),
             legend = list(title = list(text = "HDR", xanchor = "center"), orientation = "v", font = list(size = 13), yanchor = "center", y = 0, 
             itemsizing = "constant"), margin = list(r = 100))
  } else{ fig = NULL }
  
  
  results = list(df = data, IN = data[data$IN == TRUE, ], OUT = data[data$IN == FALSE, ], alpha = alpha)
  if(plot) results$plot = fig
  
  return(results)
}
