###Process WPP data

##Extract WPP female population estimates and projections for selected countries
data(popAge1dt); data(popprojAge1dt)
ACest <- popAge1dt %>% filter(country_code %in% ctrys_code) %>%
  select(country_code,name,year,age,popF)
ACproj <- popprojAge1dt %>% filter(country_code %in% ctrys_code) %>%
  select(country_code,name,year,age,popF)


##Process WPP data
#Combine into a single data frame
ACall <- rbind(ACest, ACproj)

#Filter to reproductive age range
ACall %<>% filter(age %in% 15:44)

#Convert year to numeric variable
ACall$year <- as.numeric(ACall$year)

#Convert counts from "thousands" to raw values, renaming counts as "N"
ACall %<>% rename(N=popF) %>% mutate(N=1000*N)

#Convert to midyear population estimates 
ACall %<>% group_by(name,age) %>% mutate(N2=lag(N)) %>%
  mutate(Navg=0.5*(N+N2)) %>%
  select(country_code,name,year,age,Navg) %>%
  rename(N=Navg) %>% filter(!is.na(N))

#For a given year y and age a, we assume that half the exposure is from
#the cohort y-a (turning age a in year y), and half is from the cohort y-(a+1)
#(turning age a+1 in year y)

#Add cohort column (y-a, the younger cohort)
ACall$Cohort <- ACall$year-ACall$age

#Duplicate ACall, with the duplication having the older cohort (y-(a+1))
ACall <- rbind(ACall,ACall %>% mutate(Cohort=Cohort-1))

#Halve the exposures
ACall %<>% mutate(N=0.5*N)

#Add a column to count the number of copies of each age-cohort combination
ACall %<>% group_by(name,age,Cohort) %>% mutate(ind=n())

#Remove age-cohort combinations that are only partially observed
ACall %<>% filter(ind==2)

#Sum exposures across country, age and cohort
ACall %<>% group_by(country_code,name,age,Cohort) %>% summarise(N=sum(N))

#Arrange by cohort
ACall %<>% arrange(name,Cohort,age)

#Add abbreviated country names
ACall$Country <- ctrys[match(ACall$country_code,ctrys_code)]


##Save WPP data
save(ACall,file="data/WPP_pop.RData")
