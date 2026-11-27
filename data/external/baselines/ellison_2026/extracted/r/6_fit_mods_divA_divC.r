###Fit models with different basis function knot spacing
load("output/resSMunif.RData")

##Set number of knots for the age and cohort basis functions
divA <- 7.5 #Number of age knots = NkA = floor(Na/divA)
divC <- 7.5 #Number of cohort knots = NkC = floor(Nc/divC)
#Values we consider (number of basis functions NbA = NkA + 2, NbC = NkC + 2)
#divA = 5, divC = 5 (NkA = 6, NkC = 8; NbA = 8, NbC = 10) - original, already fitted
#divA = 3, divC = 3 (NkA = 10, NkC = 13; NbA = 12, NbC = 15)
#divA = 5, divC = 3 (NkA = 6, NkC = 13; NbA = 8, NbC = 15)
#divA = 3, divC = 5 (NkA = 10, NkC = 8; NbA = 12, NbC = 10)
#divA = 7.5, divC = 7.5 (NkA = 4, NkC = 5; NbA = 6, NbC = 7)


##Prepare data for fitting models, fit models, plot model output
for (i in 1:1) {
  #Set JOY
  joy <- joys[i]
  print(joy)
  
  #Identify countries
  ctrys_unif <- setdiff(unique(resSMunif[[i]]$bounds$Code),
                        resSMunif[[i]]$ctry$Code)
  
  print(ctrys_unif)
  
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
    
    #Prepare stan data
    source("r/6_prepare_stan_data_divA_divC.r")

    #Fit aggregate model
    source("r/6_fit_mod_a_divA_divC.r")
    
    #Fit parity-specific model
    source("r/6_fit_mod_p_divA_divC.r")
    
    #Plot model output
    source("r/6_plot_model_output_ctry_joy_divA_divC.r")
  }
}
