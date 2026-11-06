HDR_error = function(data, measure = "KDE", f_true, coverage_prob = 0.95, margin = "norm", plot = TRUE, ...){
  
  n = dim(data)[1]; d = dim(data)[2]
  hdr = HDR(data = data, measure = measure, coverage_prob = coverage_prob, margin = margin, plot = plot, ...)
  
  alpha = quantile(f_true, 1-coverage_prob)
  IN_true = f_true > alpha
  
  TP = sum(IN_true & hdr$df$IN)
  TN = sum(!IN_true & !hdr$df$IN)
  FP = sum(!IN_true & hdr$df$IN)
  FN = sum(IN_true & !hdr$df$IN)
  ER = FP+FN
  
  error = round(data.frame(
    ERR = ER/n,
    FPR = FP/(FP+TN),
    FNR = FN/(FN+TP),
    Accuracy = 1-ER/n,
    F1_score = 2*TP/(2*TP+FP+FN) + 2*TN/(2*TN+FP+FN),
    MCC = (TP*TN - FP*FN)/sqrt(round((TP+FP)*(TP+FN),5)*round((TN+FP)*(TN+FN),5))
  ), 2)
  
  results = list(Error = error)
  if(plot) results$plot = hdr$plot
    
  return(results)
}





# Monte Carlo Error
MC_error = function(MC = 10, worker = parallelly::availableCores()-1, size = 500, scenarios_file, measures){
  
  # MC Simulation
  library(future.apply); library(future)
  
  old_plan = plan()
  on.exit(plan(old_plan), add = TRUE)
  
  plan(multisession, workers = worker); start_time = Sys.time()
  Results = future_lapply(seq_len(MC), function(mc){
    
    seed = mc; N = size
    source(file.path(scenarios_file), local = TRUE)
    
    setNames(lapply(measures, function(measure) {
      setNames(lapply(names(scenarios), function(nm) {
        start_time = Sys.time()
        data.frame(
          HDR_error(data = scenarios[[nm]][[1]], measure = measure, f_true = scenarios[[nm]][[2]], coverage_prob = 0.95, margin = scenarios[[nm]][[3]], plot = FALSE)$Error,
          time = as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      }), names(scenarios))
    }), measures)
    
  }, future.seed = NULL); end_time = Sys.time()
  
  # Elapsed time: hh:mm:ss
  elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
  time = sprintf("%02d:%02d:%02d mins", 
                 floor(elapsed / 3600),
                 floor((elapsed %% 3600) / 60),
                 round(elapsed %% 60))
  
  # MC mean error
  env = new.env(); source(scenarios_file, local = env)
  scenario_names = names(env$scenarios)
  
  
  Summary_list = setNames(lapply(scenario_names, function(nm) {
    setNames(lapply(measures, function(measure) {
      error_mat = do.call(rbind, lapply(seq_len(MC), function(mc) Results[[mc]][[measure]][[nm]]))
      
      data.frame(
        N = size,
        Scenario =  nm,
        Measure = measure,
        Metric = names(error_mat),
        Mean = colMeans(error_mat, na.rm = TRUE),
        SD = apply(error_mat, 2, sd, na.rm = TRUE),
        row.names = NULL
        
      )}), measures)
  }), scenario_names)
  
  Summary = do.call(rbind, c(unlist(
    lapply(scenario_names, function(nm) {
      lapply(measures, function(measure) {
        error_mat = do.call(rbind, lapply(seq_len(MC), function(mc) Results[[mc]][[measure]][[nm]]))
        
        data.frame(
          N = size,
          Scenario =  nm,
          Measure = measure,
          Metric = names(error_mat),
          Mean = round(colMeans(error_mat, na.rm = TRUE), 3),
          SD = round(apply(error_mat, 2, sd, na.rm = TRUE), 3),
          row.names = NULL
        )
      })
    }), recursive = FALSE), make.row.names = FALSE)
  )
  
  message("MC simulation completed in ", time, " (hh:mm:ss)")
  return(list(Summary = Summary, time = time, Results = Results, Results_list = Summary_list))
}
