##Fit new models with uniform prior
load("output/resSMunif.RData")

##Fit models and plot model output
for (i in 1:length(joys)) {
  #Set JOY
  joy <- joys[i]
  print(joy)
  
  #Identify countries
  ctrys_unif <- resSMunif[[i]]$ctry$Code
  print(ctrys_unif)
  
  #Identify uniform prior bounds
  cfr_lu <- resSMunif[[i]]$bounds[1:2,"sd"] %>% unlist() %>% round(4)
  ppr_lu <- resSMunif[[i]]$bounds[-c(1:2),"sd"] %>% unlist() %>% round(4)
  print(cfr_lu)
  print(ppr_lu)
  
  #Fit models for each country
  for (ctry in ctrys_unif) {
    print(ctry)
    joy_sel <- joy_sel_lst[[ctry]]
    red_sel <- red_sel_lst[[ctry]]
    red <- red_sel[which(joy_sel==joy)]
    #Set cohort range (if possible, 40 cohorts are considered: 11 complete cohorts and age 15 only 
    #observed in the last cohort; otherwise, 35 cohorts are considered, with only 6 complete cohorts)
    if (red==F) cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
    if (red==T) cohrange <- (joy-max(agerange)-5):(joy-min(agerange))
    
    ##Load stan data
    if (red==F) load(paste0("data/standata_",ctry,"_",joy,".RData"))
    if (red==T) load(paste0("data/standata_",ctry,"_",joy,"_red.RData"))
    
    #Add bounds to standata
    standata$cfr_lu <- cfr_lu
    standata$ppr_lu <- ppr_lu
    
    #Fit aggregate model
    source("r/6_fit_mod_a_SMunif.r")
    
    #Fit parity-specific model
    source("r/6_fit_mod_p_SMunif.r")
    
    #Plot model output
    source("r/6_plot_model_output_ctry_joy_SMunif.r")
  }
}
