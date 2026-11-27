#############################################################################################
#                                                                                           #
# Code to replicate results in "Can Incorporating Parity Information Improve the            #
# Reliability of Completed Cohort Fertility Projections? Insights from a Bayesian           #
# Generalized Additive Model Approach"                                                      #
# Joanne Ellison, Jakub Bijak, and Erengul Dodd*                                            #
# *E.Dodd@soton.ac.uk (corresponding author)                                                #
#                                                                                           #    
#############################################################################################

##0. Set-up
#To run this code you will need R version 3.4.0 or above and be able to install
#packages, in particular rstan, which is required to implement the Hamiltonian
#Monte Carlo methodology to fit the models. This will probably require you to also 
#install the package Rtools, explained in the rstan installation instructions that
#can be found at "https://github.com/stan-dev/rstan/wiki".

#You will also need to register with the Human Fertility Database (HFD) for free at
#https://www.humanfertility.org/cgi-bin/registration.php, and then download 'births.zip',
#'exposure.zip', 'asfr.zip' and 'parityexp.zip' from https://www.humanfertility.org/Data/ZippedDataFiles
#(rows 1-3 and 10 of 'Data by type' table; 'births', 'exposure', 'asfr' and 'parityexp' hyperlinks) and save
#the zip files in the 'data' subfolder.

#Please note that the steps are designed to be followed sequentially.

#Install required packages (if not already installed)
packages <- c("tidyverse",
              "curl",
              "readxl",
              "viridisLite",
              "magrittr",
              "directlabels",
              "rstudioapi",
              "HMDHFDplus",
              "rstan",
              "flextable",
              "officer",
              "ggh4x",
              "scales",
              "countrycode",
              "devtools",
              "bayesTFR",
              "bayesPop",
              "viridis",
              "gridExtra",
              "cowplot")
ind <- setdiff(1:length(packages), which(packages %in% rownames(installed.packages())))

if (length(ind)>0) {
  install.packages(packages[ind],
                   repos="https://cloud.r-project.org",
                   dependencies = T)
}

#Load packages
lapply(packages, require, character.only = TRUE)
windowsFonts(Calibri=windowsFont("Calibri"))

#Install wpp2024 package (if not already installed)
if (!("wpp2024" %in% installed.packages())) {
  options(timeout = 600)
  install_github("PPgp/wpp2024")
}
library(wpp2024)

#Set working directory to the location of the source file
setwd(dirname(getActiveDocumentContext()$path)) # set to location of 'source_file.r'
getwd() # check working directory is correct

#Source functions file
source("r/0_functions.r")


##1. Process data (only need to run once)
#Countries chosen based on information in 'HFD_data_availability.xlsx' in the 'data' subfolder
#Save country names (abbreviations, full names and numeric codes)
ctrys <- c("BLR","CAN","CZE","DNK","EST","HUN","IRL","JPN","LTU",
           "NLD","NOR","RUS","SVK","SWE","TWN","UKR","USA")
ctrys_long <- countrycode(ctrys, origin = "genc3c", destination = "country.name")
names(ctrys_long) <- ctrys
ctrys_code <- countrycode(ctrys, origin = "genc3c", destination = "genc3n")
names(ctrys_code) <- ctrys

#Cohorts with quality issues identified by HFD and therefore removed from analysis
#Details provided in 'HFD_data_availability.xlsx' in the 'data' subfolder - cohorts
#removed are summarised in the 'Decision' column of this spreadsheet
coh_iss <- list()
coh_iss[["BLR"]] <- 1949:1962
coh_iss[["EST"]] <- 1949:1962
coh_iss[["LTU"]] <- 1957:1960
coh_iss[["NLD"]] <- 1935:1944
coh_iss[["RUS"]] <- 1944:1954 
coh_iss[["UKR"]] <- c(1944:1954,1958:1961)

source("r/1_process_WPP_data.r")
# Running this code processes the World Population Prospects (WPP) population estimates and projections
# for the countries. It saves 'WPP_pop.RData' in the 'data' subfolder.

for (ctry in ctrys) source("r/1_process_HFD_data.r")
#Running this code extracts 'birthsVH.txt' and 'birthsVHbo.txt' from 'births.zip' in the 'data' subfolder,
#'exposVH.txt' from 'exposure.zip' in the 'data' subfolder, and saves '<ctry>_rates.RData' in the 'data' subfolder.
 
#Plot rates for all countries
for (ctry in ctrys) source("r/1_plot_rates.r")
#Running this code saves '<ctry>_rates_1946_2000.png' in the 'plots' subfolder.
#If there are cohorts within the 1946-2000 range that have had previous quality issues identified,
#the code also saves '<ctry>_rates_1946_2000_iss.png' in the 'plots' subfolder, where these cohorts
#are plotted with a dashed line.
 
#Comparison of approximate parity-specific exposures with census/register-based estimates
for (ctry in ctrys) source("r/1_comp_expos.r")
#Running this code saves '<ctry>_comp_expos_age.png' and '<ctry>_comp_expos_yr.png' in the
#'plots' subfolder, which plot the exposure estimates against each other with either the ages
#or years indicated by the colour of the points, and open circles indicating that at least one
#of the cohorts involved in the estimate have a previous quality issue identified.
#The code also saves 'ALL_comp_expos_cor_age.png' and 'ALL_comp_expos_cor_yr.png' in the
#'plots' subfolder, which plot the correlation between the exposure estimates for each country-age
#or country-year combination respectively.
 
#Comparison of within-cohort exposures for consecutive ages to assess the effect of scaling
source("r/1_comp_expos_scaling.r")
#Running this code saves four plots in the 'plots' subfolder:
#'ALL_comp_expos_scaling_dens.png': a density plot of the ratios of the exposures across all countries,
#ages and cohorts
#'ALL_comp_expos_scaling_age.png': the consecutive exposures plotted against each other, with the ages
#indicated by the colour of the points
#'ALL_comp_expos_scaling_yr.png': the consecutive exposures plotted against each other, with the years
#indicated by the colour of the points
#'ALL_comp_expos_scaling_ctry.png': the consecutive exposures plotted against each other for each country

#Identify available JOYs for each country
agerange <- 15:44 #Set age range
joys <- seq(2000,2015,5) #Create vector of jump-off years to be considered
source("r/1_identify_joys.r")
#Running this code identifies for each country whether there is sufficient data for each JOY
#to be considered. It outputs the list 'joy_sel_lst', with country-specific elements giving
#the available JOYs, and the list 'red_sel_lst', with country-specific elements indicating
#whether the reduced setup (35 instead of 40 cohorts, see below) is required. 
ctrys_fit <- lapply(joy_sel_lst,length) %>% unlist() %>% as.data.frame() %>%
  filter(.>0) %>% rownames() #Extract list of countries with at least one available JOY

#Plot Figure 1
source("r/1_plot_fig1.r")

##Steps 2-4: Prepare data for fitting models, fit models, plot model output
for (ctry in ctrys_fit) { #For each country
  print(ctry)
  joy_sel <- joy_sel_lst[[ctry]] #Extract available JOYs
  red_sel <- red_sel_lst[[ctry]] #Extract corresponding indicator of reduced setup
  for (i in 1:length(joy_sel)) { #For each JOY
    #Set JOY
    joy <- joy_sel[i]; print(joy)
    red <- red_sel[i]
    #Set cohort range (if possible, 40 cohorts are considered: 11 complete cohorts and age 15 only 
    #observed in the last cohort; otherwise, 35 cohorts are considered, with only 6 complete cohorts)
    if (red==F) cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
    if (red==T) cohrange <- (joy-max(agerange)-5):(joy-min(agerange))
    
    ##2. Prepare data for fitting models in stan
    source("r/2_prepare_stan_data.r")
    #Running this code extracts 'asfrVHbo.txt' from 'asfr.zip' in the 'data' subfolder, saves
    #'standata_<ctry>_<joy>.RData' in the 'data' subfolder for the given country and JOY
    #(with "_red" added to the filename if red=T), and saves 'HFD_sum_original_<joy>.png',
    #'HFD_sum_original_2row_<joy>.png' and 'HFD_sum_differenced_<joy>.png' in the 'plots' subfolder.
    #It also saves 'ppr01_prior_<joy>.png' in the 'plots' subfolder.
    
    
    ##3. Fit models
    #Fit aggregate model
    source("r/3_fit_mod_a.r")
    #Running this code fits the aggregate model 'mod_a.stan' to the data for the given country and JOY,
    #and saves the model output as 'mod_a_<ctry>_<joy>.RData' in the 'output' subfolder.
    #It also saves a plot of the rhat coefficients in 'mod_a_<ctry>_rhat_<joy>.png' and
    #trace plots of selected parameters in 'mod_a_<ctry>_trace1_<joy>.png' in the 'plots' subfolder.
    #"_red" is added to the filenames if red=T. 
    
    #Fit parity-specific model
    source("r/3_fit_mod_p.r")
    #Running this code fits the parity-specific model 'mod_p.stan' to the data for the given country and JOY,
    #and saves the model output as 'mod_p_<ctry>_<joy>.RData' in the 'output' subfolder.
    #It also saves a plot of the rhat coefficients in 'mod_p_<ctry>_rhat_<joy>.png' and
    #trace plots of selected parameters in 'mod_p_<ctry>_trace1_<joy>.png' and 'mod_p_<ctry>_trace2_<joy>.png'
    #in the 'plots' subfolder. "_red" is added to the filenames if red=T.
    
    
    ##4. Plot model output
    source("r/4_plot_model_output_ctry_joy.r")
    #Running this code generates and saves various plots of the model output in the 'plots' subfolder
    #for the given country and JOY:
    #'mod_ap_<ctry>_beta_med_<joy>.png': a plot of the posterior median basis function coefficients
    #for both model types
    #'mod_ap_<ctry>_rates_med_<joy>.png': a plot of the posterior median age-specific fertility rates
    #(at the aggregate and parity-specific levels) for both model types
    #'mod_ap_cfr_post_<joy>.png': a plot of the posterior distribution of the CFR for both model types
    #'mod_ap_ppr_post_<joy>.png': a plot of the posterior distribution of the PPRs for the parity-specific model
    #'mod_ap_fsd_post_<joy>.png': a plot of the posterior distributions of the FSD for the parity-specific model
  } 
}

#Plot model output for multiple countries/JOYs
source("r/4_plot_model_output_ctry_joy_mult.r")

#Plot computation times and parameter estimates 
source("r/4_plot_comp_time_param_est.r")
#Running this code generates and saves the following plots, for all countries and JOYs and both model types:
#'mod_ap_comp_time.png': a plot of the computation times
#'mod_ap_sigma_par_est.png': a plot of the age and cohort smoothing parameter estimates and 95% intervals
#'mod_ap_phi_par_est.png': a plot of the negative binomial dispersion parameter estimates and 95% intervals


##5. Compare with existing models - update
#Fit comparator models
source("r/5a_fit_comparator_models.r")
#Running this code fits the comparator models, i.e. the freeze rates model, the Myrskylä et al. model and
#the Schmertmann et al. model for all jump-off years and countries. For the Schmertmann et al. model, it is
#necessary for users to first download the code from the project website and make some minor amendments
#which are described in lines 84-104. The code also processes the output for the proposed models (aggregate
#and parity-specific) generated in Step 3. The code also computes two scoring rules: the continuous ranked
#probability score (CRPS) and the logarithmic score (LogS, only for the probabilistic models). The fit and
#scores are computed on the ASFR and CFR levels.
#We note that for the CFR, in order for the measures of predictive accuracy for the probabilistic models
#to be comparable with those from the deterministic models (which fit to the observed ASFRs perfectly),
#we compute the point estimates and uncertainty about the CFR in two different ways. The first way simply
#uses the posterior distribution of the CFR. The second way treats the ASFRs observed up to and including
#the jump-off year as known, hence only using the posterior distribution of the ASFRs for age-cohort comb-
#inations after the jump-off year to compute the CFR point estimates and uncertainty. This ensures a fairer
#comparison with the deterministic models. When plotting the posterior distribution of the CFR, we use the
#unadjusted version.

#Compare models
source("r/5b_compare_models.r")
#Running this code combines the ASFR and CFR fitted values and forecasts for all models. The error is then
#computed, allowing summary statistics to be calculated. These are saved as 'ASFR_CFR.docx' in the 'tables'
#subfolder. Various comparative plots are then generated nd saved in the 'plots' subfolder:
#'CRPS_ACtcj_bar.png': a plot of the average CRPS values for all model-country-jump-off-year combinations,
#for the ASFRs and CFRs
#'CRPS_ACtcj_bar_diff.png': the above plot expressed in terms of differences from the average CRPS for the
#proposed (parity-specific) model
#'vis_CFR_<joy1>_<joy2>.png': a plot of the CFR forecast distributions for all models and countries,
#and a pair of jump-off years chosen by the user (joy1 and joy2, chosen to be 2000 and 2005)
#'vis_ASFR_<coh2>.png': a plot of the ASFR forecast distributions for all models and countries, for a single
#cohort chosen by the user (coh2, chosen to be 1975), across three jump-off years chosen by the user (chosen
#to be 2000, 2005 and 2010).
#'vis_CFR_joys_<CFR_joys>_ASFR_<coh2>_joys_<ASFR_joys>.png': a plot of the CFR forecast distributions for all
#models and countries and jump-off years chosen by the user (CFR_joys, separated by underscores in the file name),
#combined with a plot of the ASFR forecast distributions for all models and countries, for a single cohort
#chosen by the user (coh2), across jump-off years chosen by the user (ASFR_joys,separated by underscores in the
#file name).


##6. Sensitivity analyses
#Summary fertility measure prior

source("r/6_process_HFD_sum_data_SMunif.r")
#Done

source("r/6_fit_mods_SMunif.r")
#Done

source("r/6_plot_sd_dens_SMunif.r")
#Done

source("r/6_comp_mods_SMunif.r")
#Done

#One density plot
#Remove 50% coverage plot
#Have ASFR and CFR plots for CRPS, 50% and 90% coverage
#Add to appendix and make conclusions

#Basis function dimension
source("r/6_fit_mods_divA_divC.r")
#Done

source("r/6_comp_mods_divA_divC.r")
#Done

#Using observed exposures
source("r/6_use_expos_pac.r")
#Done

##CFR tables by JOY/ctry?
##Get ASFR tables and plots, make sure we are happy that the functions
##are working correctly and we aren't missing anything. Shouldn't be
##as we only need to generate from aggregate data for these summary
##statistics so don't need to read in correct parity-specific exposures,
##just generating from aggregate model.
##Check consistent with plots and our visual comparisons
##Also can present some of the comparison plots, do one for 
##all 3 models? Or choose one of the variants to present?

#Change models 2 and 3 - pick one?


