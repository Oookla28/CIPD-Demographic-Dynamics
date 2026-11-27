###Identify available JOYs for each country
##Create lists to store results
joy_sel_lst <- red_sel_lst <- list()

for (ctry in ctrys) {
  ##Read in rates
  load(paste0("data/",ctry,"_rates.RData"))
  
  
  ##Process data
  #Identify cohorts observed from age 15
  coh_obs <- alldata %>% filter(par==1) %>%
    group_by(Cohort) %>%
    mutate(mina=min(Age),maxa=max(Age),len=n(),len2=maxa-mina+1) %>%
    mutate(diff=len2-len) %>%
    mutate(Age2=mina[1]+0:(len[1]-1)) %>%
    filter(mina==15) %>%
    filter(Age==Age2) %>%
    filter(Age==15) %>% pull(Cohort)
  
  #Remove cohorts with quality issues
  if (!is.null(coh_iss[[ctry]])) coh_obs <- setdiff(coh_obs,coh_iss[[ctry]])
  
  #Identify minimum and maximum JOYs
  min.yr <- min(coh_obs)+15
  max.yr <- max(coh_obs)+15
  
  #Check suitability for JOYs w/ 11 complete cohorts (40 cohort setup)
  joy_poss <- function(joy) {
    out <- ifelse(
      max.yr>=joy &
        all((joy-max(agerange)-10):(joy-min(agerange)) %in% coh_obs), T, F)
    out
  }
  joy_ind <- sapply(joys,joy_poss)
  
  #Check suitability for JOYs w/ 6 complete cohorts (35 cohort setup)
  joy_poss_red <- function(joy) {
    out <- ifelse(
      max.yr>=joy &
        all((joy-max(agerange)-5):(joy-min(agerange)) %in% coh_obs), T, F)
    out
  }
  joy_ind_red <- sapply(joys,joy_poss_red)
  
  
  ##Construct selection variables
  joy_sel <- joys[which(as.numeric(joy_ind+joy_ind_red)>0)]
  red_func <- function(joy) {
    i <- which(joys==joy)
    out <- ifelse(joy_ind[i]==T, F, T)
    out
  }
  red_sel <- sapply(joy_sel,red_func)
  
  
  ##Output joy_sel and red_sel
  joy_sel_lst[[ctry]] <- joy_sel
  red_sel_lst[[ctry]] <- red_sel
}
