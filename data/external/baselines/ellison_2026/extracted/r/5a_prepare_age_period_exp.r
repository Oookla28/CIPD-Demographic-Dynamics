###Prepare age-period exposures
##Prepare HFD data
expos  <-  read.table("data/exposRR.txt",skip=2,as.is=T, header=T) %>% as_tibble()
expos %<>% filter(Code %in% ctrys_fit, Age %in% 15:44)


##Prepare WPP data
#Extract WPP female population estimates and projections for selected countries
data(popAge1dt); data(popprojAge1dt)
APest <- popAge1dt %>% filter(country_code %in% ctrys_code[ctrys_fit]) %>%
  select(country_code,name,year,age,popF)
APproj <- popprojAge1dt %>% filter(country_code %in% ctrys_code[ctrys_fit]) %>%
  select(country_code,name,year,age,popF)

#Process WPP data
#Combine into a single data frame
APall <- rbind(APest, APproj)

#Filter to reproductive age range
APall %<>% filter(age %in% 15:44)

#Convert year to numeric variable
APall$year <- as.numeric(APall$year)

#Convert counts from "thousands" to raw values, renaming counts as "N"
APall %<>% rename(N=popF) %>% mutate(N=1000*N)

#Convert to midyear population estimates 
APall %<>% group_by(name,age) %>% mutate(N2=lag(N)) %>%
  mutate(Navg=0.5*(N+N2)) %>%
  select(country_code,name,year,age,Navg) %>%
  rename(N=Navg) %>% filter(!is.na(N))

#Add abbreviated country names
APall$Country <- ctrys[match(APall$country_code,ctrys_code)]


##Combine exposures from HFD and WPP data
APcomb <- data.frame()
ctrynames <- unique(tdat$Code)
ctrycodes <- unique(tdat$country_code)
for (ctry in ctrys_fit) {
  tmp <- expos %>% filter(Code==ctry)
  tmp2 <- APall %>% filter(Country==ctry,year>max(tmp$Year))
  tmp3 <- rbind(tmp %>% left_join(data.frame(Code=ctrynames,country_code=ctrycodes)) %>%
                  select(country_code,Code,Year,Age,Exposure) %>% rename(N=Exposure),
                tmp2 %>% ungroup() %>% select(country_code,Country,year,age,N) %>%
                  rename(Code=Country,Year=year,Age=age))
  tmp3 %<>% filter(Year %in% (min(joys)-30):(max(joys)+30)) #restrict to years of interest
  APcomb <- rbind(APcomb,tmp3)
}
