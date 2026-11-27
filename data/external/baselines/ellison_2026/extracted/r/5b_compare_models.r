###Compare models
##Combine model predictions
#ASFRs
combadat <- rbind(data.frame(Model="Freeze rates",fradat %>% select(-mincoh,-red)),
                  data.frame(Model="Myrskylä et al.",myrskadat %>% select(-mincoh,-red)),
                  data.frame(Model="Schmertmann et al.",schmertadat %>% select(-mincoh,-red)),
                  data.frame(Model="Ševčíková et al.",unadat %>% select(-mincoh,-mean,-red)),
                  data.frame(Model="Aggregate",
                             propadat %>% filter(Model=="Aggregate") %>% select(-Model)),
                  data.frame(Model="Parity-specific",
                             propadat %>% filter(Model=="Parity-specific") %>% select(-Model))) %>%
  mutate(Model=factor(Model,levels=c("Parity-specific","Aggregate","Ševčíková et al.","Schmertmann et al.","Myrskylä et al.","Freeze rates"))) %>%
  arrange(Model,joy,ctry,Cohort,Age)
combadat2 <- combadat %>% filter(Year>joy,!is.na(x))
combadat2$err <- combadat2$x-combadat2$pred

#CFRs
combcdat <- rbind(data.frame(Model="Freeze rates",frcdat),
                  data.frame(Model="Myrskylä et al.",myrskcdat),
                  data.frame(Model="Schmertmann et al.",schmertjcdat[,c(2,1,3,8,7,9,10,13,11,14,12,15,16,17)]),
                  data.frame(Model="Ševčíková et al.",unjcdat[,c(2,1,3,8,7,10,11:18)]),
                  data.frame(Model="Aggregate",
                             propjcdat %>% filter(Model=="Aggregate") %>% select(-Model)) %>%
                               select(Model,ctry,Cohort,joy,pred,x,sd,lwr95,upr95,lwr90,upr90,lwr50,upr50,crps,logs),
                  data.frame(Model="Parity-specific",
                             propjcdat %>% filter(Model=="Parity-specific") %>% select(-Model)) %>%
                               select(Model,ctry,Cohort,joy,pred,x,sd,lwr95,upr95,lwr90,upr90,lwr50,upr50,crps,logs)) %>%
  mutate(Model=factor(Model,levels=c("Parity-specific","Aggregate","Ševčíková et al.","Schmertmann et al.","Myrskylä et al.","Freeze rates"))) %>%
  arrange(Model,joy,ctry,Cohort)
combcdat2 <- combcdat %>% filter(Cohort+max(agerange)>joy,!is.na(x))
combcdat2$err <- combcdat2$x-combcdat2$pred


##Create tables of overall summary statistics
#Set-up
files <- c('tables/ASFR_CFR.docx')
sect_properties <- prop_section(
  page_size = page_size(
    orient = "portrait",
    width = 8.3, height = 11.7
  ),
  type = "continuous",
  page_margins = page_mar()
)

#ASFRs
tmp <- combadat2 %>% group_by(Model) %>%
  summarise(
    MAE=mean(abs(err)),
    MAPE=100*mean(abs((err)/x)),
    RMSE=sqrt(mean(err^2)),
    Bias=mean(-err),
    Rel_bias=100*mean(-err/x),
    cov90=100*mean(x>=lwr90 & x<=upr90),
    cov50=100*mean(x>=lwr50 & x<=upr50),
    CRPS=mean(crps)) %>%
  flextable() %>%
  set_formatter(
    MAE = function(x) formatC(x, digits = 4, format = 'f'),
    MAPE = function(x) formatC(x, digits = 1, format = 'f'),
    RMSE = function(x) formatC(x, digits = 3, format = 'f'),
    Bias = function(x) formatC(x, digits = 3, format = 'f'),
    Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
    cov90 = function(x) formatC(x, digits = 0, format = 'f'),
    cov50 = function(x) formatC(x, digits = 0, format = 'f'),
    CRPS = function(x) formatC(x, digits = 4, format = 'f')) %>%
  bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
  bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
  bold(~ min(MAE) == MAE, ~ MAE) %>%
  bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
  bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
  bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
  bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
  bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
  font(part = 'all', fontname = 'Calibri') %>%
  fontsize(part = 'all', size = 10) %>%
  bold(part = 'header') %>%
  align(part = 'header', align = 'center') %>%
  set_caption("ASFR summary statistics by model: all countries, all JOYs")

#CFRs
tmp2 <- combcdat2 %>% group_by(Model) %>%  
  summarise(
    MAE=mean(abs(err)),
    MAPE=100*mean(abs((err)/x)),
    RMSE=sqrt(mean(err^2)),
    Bias=mean(-err),
    Rel_bias=100*mean(-err/x),
    cov90=100*mean(x>=lwr90 & x<=upr90),
    cov50=100*mean(x>=lwr50 & x<=upr50),
    CRPS=mean(crps)) %>%
  flextable() %>%
  set_formatter(
    MAE = function(x) formatC(x, digits = 3, format = 'f'),
    MAPE = function(x) formatC(x, digits = 1, format = 'f'),
    RMSE = function(x) formatC(x, digits = 3, format = 'f'),
    Bias = function(x) formatC(x, digits = 3, format = 'f'),
    Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
    cov90 = function(x) formatC(x, digits = 0, format = 'f'),
    cov50 = function(x) formatC(x, digits = 0, format = 'f'),
    CRPS = function(x) formatC(x, digits = 3, format = 'f')) %>%
  bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
  bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
  bold(~ min(MAE) == MAE, ~ MAE) %>%
  bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
  bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
  bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
  bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
  bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
  font(part = 'all', fontname = 'Calibri') %>%
  fontsize(part = 'all', size = 10) %>%
  bold(part = 'header') %>%
  align(part = 'header', align = 'center') %>%
  set_caption("CFR summary statistics by model: all countries, all JOYs")

#Save tables to file
save_as_docx(tmp, tmp2, pr_section=sect_properties, path = files)


##Create tables of summary statistics by JOY
#ASFRs
afunc <- function(joy2) {
  tmp <- combadat2 %>% filter(joy==joy2) %>%
    group_by(Model) %>%
    summarise(
      MAE=mean(abs(err)),
      MAPE=100*mean(abs((err)/x)),
      RMSE=sqrt(mean(err^2)),
      Bias=mean(-err),
      Rel_bias=100*mean(-err/x),
      cov90=100*mean(x>=lwr90 & x<=upr90),
      cov50=100*mean(x>=lwr50 & x<=upr50),
      CRPS=mean(crps)) %>%
    flextable() %>%
    set_formatter(
      MAE = function(x) formatC(x, digits = 4, format = 'f'),
      MAPE = function(x) formatC(x, digits = 1, format = 'f'),
      RMSE = function(x) formatC(x, digits = 3, format = 'f'),
      Bias = function(x) formatC(x, digits = 3, format = 'f'),
      Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
      cov90 = function(x) formatC(x, digits = 0, format = 'f'),
      cov50 = function(x) formatC(x, digits = 0, format = 'f'),
      CRPS = function(x) formatC(x, digits = 4, format = 'f')) %>%
    bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
    bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
    bold(~ min(MAE) == MAE, ~ MAE) %>%
    bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
    bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
    bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
    bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
    bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
    font(part = 'all', fontname = 'Calibri') %>%
    fontsize(part = 'all', size = 10) %>%
    bold(part = 'header') %>%
    align(part = 'header', align = 'center') %>%
    set_caption(paste0("ASFR summary statistics by model: all countries, JOY = ",joy2))
  tmp
}

tmpa <- list()
for (i in 1:length(joys)) {
  tmpa[[i]] <- afunc(joys[i])
}

#CFRs
cfunc <- function(joy2) {
  tmp2 <- combcdat2 %>% filter(joy==joy2) %>%
    group_by(Model) %>%  
    summarise(
      MAE=mean(abs(err)),
      MAPE=100*mean(abs((err)/x)),
      RMSE=sqrt(mean(err^2)),
      Bias=mean(-err),
      Rel_bias=100*mean(-err/x),
      cov90=100*mean(x>=lwr90 & x<=upr90),
      cov50=100*mean(x>=lwr50 & x<=upr50),
      CRPS=mean(crps)) %>%
    flextable() %>%
    set_formatter(
      MAE = function(x) formatC(x, digits = 3, format = 'f'),
      MAPE = function(x) formatC(x, digits = 1, format = 'f'),
      RMSE = function(x) formatC(x, digits = 3, format = 'f'),
      Bias = function(x) formatC(x, digits = 3, format = 'f'),
      Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
      cov90 = function(x) formatC(x, digits = 0, format = 'f'),
      cov50 = function(x) formatC(x, digits = 0, format = 'f'),
      CRPS = function(x) formatC(x, digits = 3, format = 'f')) %>%
    bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
    bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
    bold(~ min(MAE) == MAE, ~ MAE) %>%
    bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
    bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
    bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
    bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
    bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
    font(part = 'all', fontname = 'Calibri') %>%
    fontsize(part = 'all', size = 10) %>%
    bold(part = 'header') %>%
    align(part = 'header', align = 'center') %>%
    set_caption(paste0("CFR summary statistics by model: all countries, JOY = ",joy2))
  tmp2
}

tmp2a <- list()
for (i in 1:length(joys)) {
  tmp2a[[i]] <- cfunc(joys[i])
}

#Combine tables
tmp3a <- c(tmpa,tmp2a)

#Save tables to file
write_word_table <- function(table, doc){
  doc %>%
    body_add_flextable(table) %>%
    body_add_par("")
}

my_doc <- read_docx()

for (i in 1:length(tmp3a)) {
  my_doc <- write_word_table(tmp3a[[i]], my_doc)
}

print(my_doc, target = "tables/ASFR_CFR_JOY.docx")


##Create tables of summary statistics by country
#ASFRs
afunc <- function(ctry2) {
  tmp <- combadat2 %>% filter(ctry==ctry2) %>%
    group_by(Model) %>%
    summarise(
      MAE=mean(abs(err)),
      MAPE=100*mean(abs((err)/x)),
      RMSE=sqrt(mean(err^2)),
      Bias=mean(-err),
      Rel_bias=100*mean(-err/x),
      cov90=100*mean(x>=lwr90 & x<=upr90),
      cov50=100*mean(x>=lwr50 & x<=upr50),
      CRPS=mean(crps)) %>%
    flextable() %>%
    set_formatter(
      MAE = function(x) formatC(x, digits = 4, format = 'f'),
      MAPE = function(x) formatC(x, digits = 1, format = 'f'),
      RMSE = function(x) formatC(x, digits = 3, format = 'f'),
      Bias = function(x) formatC(x, digits = 3, format = 'f'),
      Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
      cov90 = function(x) formatC(x, digits = 0, format = 'f'),
      cov50 = function(x) formatC(x, digits = 0, format = 'f'),
      CRPS = function(x) formatC(x, digits = 4, format = 'f')) %>%
    bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
    bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
    bold(~ min(MAE) == MAE, ~ MAE) %>%
    bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
    bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
    bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
    bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
    bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
    font(part = 'all', fontname = 'Calibri') %>%
    fontsize(part = 'all', size = 10) %>%
    bold(part = 'header') %>%
    align(part = 'header', align = 'center') %>%
    set_caption(paste0("ASFR summary statistics by model: all JOYs, ",ctry2))
  tmp
}

tmpa <- list()
for (i in 1:length(ctrys_fit)) {
  tmpa[[i]] <- afunc(ctrys_fit[i])
}

#CFRs
cfunc <- function(ctry2) {
  tmp2 <- combcdat2 %>% filter(ctry==ctry2) %>%
    group_by(Model) %>%  
    summarise(
      MAE=mean(abs(err)),
      MAPE=100*mean(abs((err)/x)),
      RMSE=sqrt(mean(err^2)),
      Bias=mean(-err),
      Rel_bias=100*mean(-err/x),
      cov90=100*mean(x>=lwr90 & x<=upr90),
      cov50=100*mean(x>=lwr50 & x<=upr50),
      CRPS=mean(crps)) %>%
    flextable() %>%
    set_formatter(
      MAE = function(x) formatC(x, digits = 3, format = 'f'),
      MAPE = function(x) formatC(x, digits = 1, format = 'f'),
      RMSE = function(x) formatC(x, digits = 3, format = 'f'),
      Bias = function(x) formatC(x, digits = 3, format = 'f'),
      Rel_bias = function(x) formatC(x, digits = 1, format = 'f'),
      cov90 = function(x) formatC(x, digits = 0, format = 'f'),
      cov50 = function(x) formatC(x, digits = 0, format = 'f'),
      CRPS = function(x) formatC(x, digits = 3, format = 'f')) %>%
    bold(~ min(MAPE) == MAPE, ~ MAPE) %>%
    bold(~ min(RMSE) == RMSE, ~ RMSE) %>%
    bold(~ min(MAE) == MAE, ~ MAE) %>%
    bold(~ min(abs(Bias)) == abs(Bias), ~ Bias) %>%
    bold(~ min(abs(Rel_bias)) == abs(Rel_bias), ~ Rel_bias) %>%
    bold(~ min(abs(cov90-90),na.rm=T) == abs(cov90-90), ~ cov90) %>%
    bold(~ min(abs(cov50-50),na.rm=T) == abs(cov50-50), ~ cov50) %>%
    bold(~ min(CRPS) == CRPS, ~ CRPS) %>%
    font(part = 'all', fontname = 'Calibri') %>%
    fontsize(part = 'all', size = 10) %>%
    bold(part = 'header') %>%
    align(part = 'header', align = 'center') %>%
    set_caption(paste0("CFR summary statistics by model: all JOYs, ",ctry2))
  tmp2
}

tmp2a <- list()
for (i in 1:length(ctrys_fit)) {
  tmp2a[[i]] <- cfunc(ctrys_fit[i])
}

#Combine tables
tmp3a <- c(tmpa,tmp2a)

#Save tables to file
write_word_table <- function(table, doc){
  doc %>%
    body_add_flextable(table) %>%
    body_add_par("")
}

my_doc <- read_docx()

for (i in 1:length(tmp3a)) {
  my_doc <- write_word_table(tmp3a[[i]], my_doc)
  }

print(my_doc, target = "tables/ASFR_CFR_ctry.docx")


##Plot CRPS values for CFR forecasts
#Prepare dataset
combacdat2 <- rbind(data.frame(
  type="ASFR",
  combadat2 %>% select(Model,Age,ctry,Cohort,joy,crps) %>%
    group_by(Model,ctry,joy) %>% summarise(crps=mean(crps))),
  data.frame(
    type="CFR",
    combcdat2 %>% select(Model,ctry,Cohort,joy,crps) %>%
      group_by(Model,ctry,joy) %>% summarise(crps=mean(crps))))

#Prepare dataset - JOY-specific
combacdat2 <- data.frame(type="CFR",
                         combcdat2 %>% select(Model,ctry,Cohort,joy,crps) %>%
                           group_by(Model,ctry,joy) %>% summarise(crps=mean(crps)))

#Prepare dataset - all JOYs combined
combacdat2_avg <- data.frame(type="CFR",
                         combcdat2 %>% select(Model,ctry,Cohort,joy,crps) %>%
                           group_by(Model,ctry) %>% summarise(crps=mean(crps)))

#Prepare objects for plotting
combacdat2_avg_order <- combacdat2_avg %>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long)) %>%
  group_by(Model) %>% arrange(desc(crps)) %>% mutate(ind=1:n()) %>%
  arrange(Model) %>%
  group_by(ctry) %>% mutate(ind=ind[1]) %>% mutate(joy=0)

combacdat2_order2 <- combacdat2 %>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long)) %>%
  group_by(Model,joy) %>% arrange(desc(crps)) %>% mutate(ind=1:n()) %>%
  arrange(Model) %>% group_by(joy,ctry) %>% mutate(ind=ind[1]) %>%
  select(type,Model,ctry,crps,ctryl,ind,joy)

combacdat2_avg_order_joy2 <- rbind(combacdat2_avg_order,combacdat2_order2)
combacdat2_avg_order_joy2 %<>% mutate(ind=ind+joy*10)
combacdat2_avg_order_joy2$joy2 <- factor(combacdat2_avg_order_joy2$joy,levels=c(0,joys))
combacdat2_avg_order_joy2$joy3 <- factor(combacdat2_avg_order_joy2$joy)
levels(combacdat2_avg_order_joy2$joy3) <- 
  c("All JOYs",
    expression(paste(italic("JOY")," = 2000")),
    expression(paste(italic("JOY")," = 2005")),
    expression(paste(italic("JOY")," = 2010")),
    expression(paste(italic("JOY")," = 2015")))

#Plot average CRPS values for each JOY and overall
supp.labj <- c("All JOYs",joys)
names(supp.labj) <- c(0,joys)
combacdat2_avg_order_joy3 <- combacdat2_avg_order_joy2 %>% filter(joy!=0) %>%
  left_join(ctry_joy_df %>% select(-mincoh)) %>%
  mutate(ctryl2=ifelse(red==F,ctryl,paste0(ctryl,"*")))

p1 <- combacdat2_avg_order_joy2 %>% filter(joy==0) %>%
  ggplot(aes(x=ind,y=crps,color=Model,shape=Model,linetype=Model)) +
  scale_x_continuous(breaks=combacdat2_avg_order_joy2$ind,minor_breaks=NULL,
                     labels=combacdat2_avg_order_joy2$ctryl) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,0.12,0.02)) +
  geom_line(aes(group=Model),show.legend = F) +
  geom_point(show.legend = F) +
  labs(y="Average CRPS",x="Country",color="Model",fill="Model") +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.08)) +
  scale_shape_manual(values=c(19,1,15,0,17,8)) +
  scale_linetype_manual(values=c(1,2,1,2,2,1)) +
  scale_color_manual(values=c("blue","blue","orange","orange","black","black")) +
  guides(color=guide_legend(nrow=1),linetype=guide_legend(nrow=1),shape=guide_legend(nrow=1)) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(color="black"),
        axis.title=element_text(face="bold"),
        plot.title=element_text(hjust=0.5,size=12),
        plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
  labs(title = expression(bold("a. All")~bolditalic("JOY")~bold("values combined")))

p2 <- combacdat2_avg_order_joy3 %>%
  ggplot(aes(x=ind,y=crps,color=Model,shape=Model,linetype=Model)) +
  scale_x_continuous(breaks=combacdat2_avg_order_joy3$ind,minor_breaks=NULL,
                     labels=combacdat2_avg_order_joy3$ctryl2) +
  geom_line(aes(group=Model)) + geom_point() +
  facet_wrap(~joy3,scales="free",labeller=label_parsed) +
  labs(y="Average CRPS",x="Country",color="Model",fill="Model") +
  theme_bw() +
  scale_shape_manual(values=c(19,1,15,0,17,8)) +
  scale_linetype_manual(values=c(1,2,1,2,2,1)) +
  scale_color_manual(values=c("blue","blue","orange","orange","black","black")) +
  guides(color=guide_legend(nrow=2,byrow=T),
         linetype=guide_legend(nrow=2,byrow=T),
         shape=guide_legend(nrow=2,byrow=T)) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text=element_text(color="black"),
        axis.title=element_text(face="bold"),
        plot.title=element_text(hjust=0.5,size=12),
        plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
  labs(title = expression(bold("b. Each individual"~bolditalic("JOY")~bold("value"))))

png(file=paste0("plots/CRPS_Ctcj_line_all_JOYs.png"),width=15,height=22,units="cm",res=400,family="Calibri")
print({plot_grid(p1, p2, nrow = 2, rel_heights = c(1/3, 2/3)) +
    theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"))})
dev.off()

tiff(file=paste0("plots/CRPS_Ctcj_line_all_JOYs.tif"),
     width=15,height=22,units="cm",res=600,family="Calibri",compression="none")
print({plot_grid(p1, p2, nrow = 2, rel_heights = c(1/3, 2/3)) +
    theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"))})
dev.off()


##Plot CFR forecasts for a given pair of jump-off years and four countries
#Set-up
joydat <- data.frame(joy=joys,xint=joys-max(agerange))
joy2 <- joys[1:2] #Choose two jump-off years to plot
ctry2 <- c("CAN","CZE","NLD","USA") #choose countries

#Prepare dataset - replace proposed models and Schmertmann et al. model with true predictions
tmp <- combcdat %>% filter(joy %in% joy2, ctry %in% ctry2) %>%
  left_join(propcdat %>% select(joy,ctry,Model,coh,med) %>% rename("Cohort"=coh)) %>%
  mutate(pred=ifelse(Model %in% c("Parity-specific","Aggregate"),med,pred)) %>% select(-med) %>%
  left_join(schmertcdat %>% select(joy,ctry,Cohort,pred) %>% rename("med"=pred) %>% mutate(Model="Schmertmann et al.")) %>%
  mutate(pred=ifelse(Model=="Schmertmann et al.",med,pred)) %>% select(-med) %>%
  mutate(pred=ifelse(Model=="Ševčíková et al." & is.na(pred),x,pred))
tmp %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long)) %>%
  mutate(Model=factor(Model,levels=c("Parity-specific","Aggregate","Ševčíková et al.","Schmertmann et al.","Myrskylä et al.","Freeze rates")))

#Save plot
png(file=paste0("plots/vis_CFR_",joy2[1],"_",joy2[2],"_",paste(ctry2,collapse="_"),".png"),width=16,height=10.15,units="cm",res=400,family="Calibri")
print({ggplot(tmp,aes(x=Cohort,color=Model,fill=Model)) +
    geom_hline(aes(yintercept=2),color="grey50") +
    geom_ribbon(aes(ymin=lwr95,ymax=upr95),alpha=0.3,linewidth=0.3,linetype=2) +
    geom_line(aes(y=pred)) +
    facet_grid(paste0("JOY = ",joy)~ctryl,scales="free_x") +
    scale_x_continuous(breaks=seq(1900,2005,5),minor_breaks=1946:2010,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=seq(0,3,0.2)) +
    geom_point(aes(x=Cohort,y=x),color=1,size=0.5,show.legend=F) +
    geom_vline(data=joydat %>% filter(joy %in% joy2),aes(xintercept=xint,group=joy),color=1) +
    labs(x="Cohort",y="CFR",color="Model",fill="Model") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom",
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))})
dev.off()


##Plot ASFR forecasts for a given cohort across multiple jump-off years
#Set-up
coh2 <- 1975 #Choose cohort to plot
joy2 <- joys[1:3] #Choose jump-off years to plot
joydat2 <- data.frame(expand.grid(Cohort=coh2,joy=joy2))
joydat2$age <- joydat2$joy-joydat2$Cohort

#Prepare dataset
tmp <- combadat %>% filter(Cohort %in% coh2, joy %in% joy2, ctry %in% ctry2) %>%
  left_join(data.frame(ctry=ctrys,ctryl=ctrys_long)) %>%
  mutate(pred=ifelse(Model=="Ševčíková et al." & is.na(pred),x,pred))

#Save plot
png(file=paste0("plots/vis_ASFR_",coh2,"_",paste(ctry2,collapse="_"),".png"),width=16,height=10,units="cm",res=400,family="Calibri")
print({ggplot(tmp,aes(x=Age,color=Model,fill=Model)) +
    geom_ribbon(aes(ymin=lwr95,ymax=upr95,color=Model,fill=Model),alpha=0.3,linewidth=0.3,linetype=2) +
    geom_line(aes(y=pred,color=Model)) +
    facet_grid(paste0("JOY = ",joy)~ctryl) +
    scale_x_continuous(breaks=seq(15,45,5),minor_breaks=agerange,expand=c(0.025,0.025)) +
    scale_y_continuous(breaks=seq(0,1,0.05)) +
    coord_cartesian(ylim=c(0,0.2)) +
    geom_point(aes(x=Age,y=x),color=1,size=0.5,show.legend=F) +
    geom_vline(data=joydat2,aes(xintercept=age),color=1) +
    labs(x="Age",y="ASFR",color="Model",fill="Model") +
    theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom") })
dev.off()


##Plot CFR and ASFR posterior distributions together
#Set-up
cfr_joy <- joys[1:2] #Choose CFR jump-off years to plot
asfr_coh <- 1975 #Choose ASFR cohort to plot
asfr_joy <- joys[1:3] #Choose ASFR jump-off years to plot
joydat_cfr <- data.frame(joy=cfr_joy,xint=cfr_joy-max(agerange))
joydat_cfr$joy2 <- factor(joydat_cfr$joy)
levels(joydat_cfr$joy2) <- c(expression(paste(italic("JOY")," = 2000")),
                             expression(paste(italic("JOY")," = 2005"))) # adjust as necessary
intdat_cfr_x <- data.frame(joy=cfr_joy,yint=2)
joydat_asfr <- data.frame(expand.grid(Cohort=asfr_coh,joy=asfr_joy))
joydat_asfr$cov <- joydat_asfr$joy-joydat_asfr$Cohort
joydat_asfr$joy2 <- factor(joydat_asfr$joy)
levels(joydat_asfr$joy2) <- c(expression(paste(italic("JOY")," = 2000")),
                              expression(paste(italic("JOY")," = 2005")),
                              expression(paste(italic("JOY")," = 2010"))) # adjust as necessary
ctry2 <- c("CAN","CZE","NLD","USA") #choose countries

#Prepare datasets
tmp_cfr <- combcdat %>% filter(joy %in% cfr_joy, ctry %in% ctry2) %>%
  left_join(propcdat %>% select(joy,ctry,Model,coh,med) %>% rename("Cohort"=coh)) %>%
  mutate(pred=ifelse(Model %in% c("Parity-specific","Aggregate"),med,pred)) %>% select(-med) %>%
  left_join(schmertcdat %>% select(joy,ctry,Cohort,pred) %>% rename("med"=pred) %>% mutate(Model="Schmertmann et al.")) %>%
  mutate(pred=ifelse(Model=="Schmertmann et al.",med,pred)) %>% select(-med) %>%
  mutate(pred=ifelse(Model=="Ševčíková et al." & is.na(pred),x,pred))
tmp_cfr %<>% left_join(data.frame(ctry=ctrys,ctryl=ctrys_long)) %>%
  mutate(Model=factor(Model,levels=c("Parity-specific","Aggregate","Ševčíková et al.","Schmertmann et al.","Myrskylä et al.","Freeze rates")))
tmp_cfr %<>% rename("cov"="Cohort")

tmp_asfr <- combadat %>% filter(Cohort %in% asfr_coh, joy %in% asfr_joy, ctry %in% ctry2) %>%
  left_join(data.frame(ctry=ctrys,ctryl=ctrys_long))  %>%
  mutate(pred=ifelse(Model=="Ševčíková et al." & is.na(pred),x,pred))
tmp_asfr %<>% rename("cov"="Age")

#Prepare objects for plotting
tmp1 <- data.frame(tmp_cfr[,c("Model","ctryl","cov","joy","pred","lwr95","upr95","x")])
tmp1$joy2 <- factor(tmp1$joy)
levels(tmp1$joy2) <- c(expression(paste(italic("JOY")," = 2000")),
                       expression(paste(italic("JOY")," = 2005")))
tmp1$ctryl2 <- factor(tmp1$ctryl)
levels(tmp1$ctryl2) <- c(expression(paste("Canada")),
                         expression(paste("Czechia")),
                         expression(paste("Netherlands")),
                         expression(paste("United States")))

tmp2 <- data.frame(type="2",tmp_asfr[,c("Model","ctryl","cov","joy","pred","lwr95","upr95","x")])
tmp2$joy2 <- factor(tmp2$joy)
levels(tmp2$joy2) <- c(expression(paste(italic("JOY")," = 2000")),
                       expression(paste(italic("JOY")," = 2005")),
                       expression(paste(italic("JOY")," = 2010")))
tmp2$ctryl2 <- factor(tmp2$ctryl)
levels(tmp2$ctryl2) <- c(expression(paste("Canada")),
                         expression(paste("Czechia")),
                         expression(paste("Netherlands")),
                         expression(paste("United States")))

#Plot CFR and ASFR forecasts
p1 <- ggplot(tmp1,aes(x=cov,color=Model,fill=Model)) +
  geom_hline(data=intdat_cfr_x,aes(yintercept=yint),color="grey50") +
  geom_ribbon(aes(ymin=lwr95,ymax=upr95,color=Model,fill=Model),alpha=0.3,linewidth=0.3,linetype=2) +
  geom_line(aes(y=pred,color=Model)) +
  facet_grid(joy2~ctryl2,labeller=label_parsed,axes="all_x") +
  geom_point(aes(x=cov,y=x),color=1,size=0.5,show.legend=F) +
  geom_vline(data=joydat_cfr,aes(xintercept=xint),color=1) +
  labs(x="Cohort",y="CFR",color="Model",fill="Model") +
  theme_bw() + theme(text = element_text("Calibri"),legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     axis.title=element_text(face="bold"),
                     axis.text=element_text(color="black"),
                     plot.title=element_text(face="bold",hjust=0.5,size=12),
                     plot.margin=unit(c(5.5, 5.5, 11, 5.5), "pt")) +
  guides(color=guide_legend(byrow=T)) +
  scale_x_continuous(limits=range(tmp_cfr$cov),breaks=seq(1900,2005,5),minor_breaks=1946:2010,expand=c(0.02,0.02)) +
  scale_y_continuous(breaks=seq(0,3,0.4),limits=c(0.7,2.71)) +
  labs(title = "a. CFR forecasts")

p2 <- ggplot(tmp2,aes(x=cov,color=Model,fill=Model)) +
  geom_ribbon(aes(ymin=lwr95,ymax=upr95,color=Model,fill=Model),alpha=0.3,linewidth=0.3,linetype=2) +
  geom_line(aes(y=pred,color=Model)) +
  facet_grid(joy2~ctryl2,labeller=label_parsed,axes="all_x") +
  geom_point(aes(x=cov,y=x),color=1,size=0.5,show.legend=F) +
  geom_vline(data=joydat_asfr,aes(xintercept=cov),color=1) +
  labs(x="Age",y="ASFR",color="Model",fill="Model") +
  theme_bw() + theme(text = element_text("Calibri"),legend.position = "bottom",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                     axis.title=element_text(face="bold"),
                     axis.text=element_text(color="black"),
                     plot.title=element_text(face="bold",hjust=0.5,size=12),
                     plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt")) +
  guides(color=guide_legend(byrow=T)) +
  scale_x_continuous(breaks=seq(15,45,5),minor_breaks=agerange,expand=c(0.025,0.025)) +
  scale_y_continuous(breaks=seq(0,1,0.05),limits=c(0,0.2)) +
  labs(title = "b. ASFR forecasts")

#Save plot
png(file=paste0("plots/vis_CFR_joys_",paste0(cfr_joy,collapse="_"),"_ASFR_",asfr_coh,"_joys_",paste0(asfr_joy,collapse="_"),"_",paste(ctry2,collapse="_"),".png"),
    width=16,height=22.6,units="cm",res=400,family="Calibri")
print({plot_grid(p1, p2, nrow = 2, rel_heights = c(2.25/5,2.75/5)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))})
dev.off()

tiff(file=paste0("plots/vis_CFR_joys_",paste0(cfr_joy,collapse="_"),"_ASFR_",asfr_coh,"_joys_",paste0(asfr_joy,collapse="_"),"_",paste(ctry2,collapse="_"),".tif"),
    width=16,height=22.6,units="cm",res=600,family="Calibri",compression="none")
print({plot_grid(p1, p2, nrow = 2, rel_heights = c(2.25/5,2.75/5)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))})
dev.off()
