###Plot rates for Fig 1

##Set-up
ctry <- "CAN" #choose country
agerange <- 15:44
joy <- joys[4] #choose JOY
cohrange <- (joy-max(agerange)-10):(joy-min(agerange))
load(paste0("data/",ctry,"_rates.RData")) #read in rates


##Lexis surfaces
#Set-up
stepfunc <- function(s1,s2) {
  dat <- data.frame(y=(s1+30):s1+0.5)
  dat$x <- seq(s2,by=1,length=nrow(dat))+0.5
  dat
}

#Aggregate
png(file=paste0("plots/",ctry,"_rates_lexis_",min(cohrange),"_",max(cohrange),".png"),
    width=9,height=6,units="cm",res=400,family="Calibri")
print({ggplot(alldata %>% filter(par==0,Cohort %in% cohrange, Age+Cohort<=joy),
              aes(x=Cohort, y=Age, fill=x)) +
    geom_tile() +
    labs(x = "Cohort", y = "Age", fill="Rate")+
    scale_fill_viridis(option ="D", guide = guide_colorbar(barheight = 7.5, barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(0,0.5,0.02), limits=c(0,0.13)) +
    scale_x_continuous(breaks=seq(1905,2105,10),minor_breaks=setdiff(1900:2100,seq(1905,2105,10)),expand=c(0,0)) +
    scale_y_continuous(limits=c(14.5,44.5), breaks=seq(15,44,5), minor_breaks=setdiff(15:44,seq(15,44,5)), expand=c(0,0)) +
    geom_vline(xintercept=cohrange+0.5,color="lightgray",linewidth=0.25) +
    geom_hline(yintercept=15.5:43.5,color="lightgray",linewidth=0.25) +
    geom_step(data=stepfunc(14,max(cohrange)-30), mapping=aes(x=x,y=y,fill=NULL),direction="hv") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text("Calibri"),
          axis.text=element_text(color="black"), axis.title.y = element_text(margin = margin(r = unit(5, "pt"))))})
dev.off()

#Parity-specific
png(file=paste0("plots/",ctry,"_rates_lexis_p14_",min(cohrange),"_",max(cohrange),".png"),
    width=20,height=6,units="cm",res=400,family="Calibri")
print({ggplot(alldata %>% filter(par>0,Cohort %in% cohrange,Age+Cohort<=joy) %>%
                mutate(par2=paste0("Parity ",c(0,1,2,"3+")[as.numeric(par)])) %>%
                mutate(x=ifelse(x>0.3,0.3,x)),
              aes(x=Cohort, y=Age, fill=x)) +
    geom_tile() +
    facet_wrap(~par2,nrow=1) +
    labs(x = "Cohort", y = "Age", fill="Rate")+
    scale_fill_viridis(option ="D", guide = guide_colorbar(barheight = 7.5, barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(0,0.5,0.05), limits=c(0,0.3001)) +
    scale_x_continuous(breaks=seq(1905,2105,10),minor_breaks=setdiff(1900:2100,seq(1905,2105,10)),expand=c(0,0)) +
    scale_y_continuous(limits=c(14.5,44.5), breaks=seq(15,44,5), minor_breaks=setdiff(15:44,seq(15,44,5)), expand=c(0,0)) +
    geom_vline(xintercept=cohrange+0.5,color="lightgray",size=0.1) +
    geom_hline(yintercept=15.5:43.5,color="lightgray",size=0.1) +
    geom_step(data=stepfunc(14,max(cohrange)-30), mapping=aes(x=x,y=y,fill=NULL),direction="hv") +
    theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text("Calibri"),
                       axis.text=element_text(color="black"), axis.title.y = element_text(margin = margin(r = unit(5, "pt"))))})
dev.off()


##Age curves
#Aggregate
png(file=paste0("plots/",ctry,"_rates_curves_",min(cohrange),"_",max(cohrange),".png"),
    width=9,height=6,units="cm",res=400,family="Calibri")
print({ggplot(alldata %>% filter(par==0,Cohort %in% cohrange, Age+Cohort<=joy),
              aes(x=Age, y=x, color=Cohort)) +
    geom_line(aes(group=Cohort)) +
    labs(x = "Age", y = "Rate", color="Cohort")+
    scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 7.5,barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(1925,2005,10), limits=range(cohrange)) +
    scale_x_continuous(breaks=seq(15,40,5),minor_breaks=agerange,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=seq(0,0.13,0.02),minor_breaks=NULL,limits=c(0,0.13)) +
    theme_bw() + theme(text = element_text("Calibri"),axis.text=element_text(color="black"))})
dev.off()

#Parity-specific
png(file=paste0("plots/",ctry,"_rates_curves_p14_",min(cohrange),"_",max(cohrange),".png"),
    width=20,height=6,units="cm",res=400,family="Calibri")
print({ggplot(alldata %>% filter(par>0,Cohort %in% cohrange,Age+Cohort<=joy) %>%
                mutate(par2=paste0("Parity ",c(0,1,2,"3+")[as.numeric(par)])),
              aes(x=Age, y=x, color=Cohort)) +
    geom_line(aes(group=Cohort),size=0.3) +
    labs(x = "Age", y = "Rate", color="Cohort")+
    facet_wrap(~par2,nrow=1) +
    coord_cartesian(ylim=c(0,0.3)) +
    scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 7.5,barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(1925,2005,10), limits=range(cohrange)) +
    scale_x_continuous(breaks=seq(15,40,5),minor_breaks=agerange,expand=c(0.02,0.02)) +
    scale_y_continuous(breaks=seq(0,0.3,0.1)) +
    theme_bw() + theme(text = element_text("Calibri"),axis.text=element_text(color="black"))})
dev.off()


##Combined plot
p1 <- ggplot(alldata %>% filter(par==0,Cohort %in% cohrange, Age+Cohort<=joy),
             aes(x=Cohort, y=Age, fill=x)) +
  geom_tile() +
  labs(x = "Cohort", y = "Age", fill="Rate")+
  scale_fill_viridis(option ="D", guide = guide_colorbar(barheight = 7.5, barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(0,0.5,0.02), limits=c(0,0.13)) +
  scale_x_continuous(breaks=seq(1905,2105,10),minor_breaks=setdiff(1900:2100,seq(1905,2105,10)),expand=c(0,0)) +
  scale_y_continuous(limits=c(14.5,44.5), breaks=seq(15,44,5), minor_breaks=setdiff(15:44,seq(15,44,5)), expand=c(0,0)) +
  geom_vline(xintercept=cohrange+0.5,color="lightgray",linewidth=0.25) +
  geom_hline(yintercept=15.5:43.5,color="lightgray",linewidth=0.25) +
  geom_step(data=stepfunc(14,max(cohrange)-30), mapping=aes(x=x,y=y,fill=NULL),direction="hv") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text("Calibri"),
                     axis.text=element_text(color="black"),
                     axis.title=element_text(face="bold"),
                     axis.title.y = element_text(margin = margin(r = unit(5, "pt"))),
                     plot.margin=unit(c(5.5, 5.5, 11, 5.5), "pt"))

p2 <- ggplot(alldata %>% filter(par==0,Cohort %in% cohrange, Age+Cohort<=joy),
             aes(x=Age, y=x, color=Cohort)) +
  geom_line(aes(group=Cohort)) +
  labs(x = "Age", y = "Rate", color="Cohort")+
  scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 7.5,barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(1925,2005,10), limits=range(cohrange)) +
  scale_x_continuous(breaks=seq(15,40,5),minor_breaks=agerange,expand=c(0.02,0.02)) +
  scale_y_continuous(breaks=seq(0,0.13,0.02),minor_breaks=NULL,limits=c(0,0.13)) +
  theme_bw() + theme(text = element_text("Calibri"),
                     axis.text=element_text(color="black"),
                     axis.title=element_text(face="bold"),
                     plot.margin=unit(c(5.5, 5.5, 11, 5.5), "pt"))

p1_p2 <- plot_grid(p1, p2) +
  theme(plot.title=element_text(hjust=0.5,size=12),
        text = element_text("Calibri")) +
  labs(title = expression(bold("a. Canada age-specific fertility rates")))

p3 <- ggplot(alldata %>% filter(par>0,Cohort %in% cohrange,Age+Cohort<=joy) %>%
               mutate(par2=paste0("Parity ",c(0,1,2,"3+")[as.numeric(par)])) %>%
               mutate(x=ifelse(x>0.3,0.3,x)),
             aes(x=Cohort, y=Age, fill=x)) +
  geom_tile() +
  facet_wrap(~par2,nrow=1) +
  labs(x = "Cohort", y = "Age", fill="Rate")+
  scale_fill_viridis(option ="D", guide = guide_colorbar(barheight = 7.5, barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(0,0.5,0.05), limits=c(0,0.3001)) +
  scale_x_continuous(breaks=seq(1905,2105,10),minor_breaks=setdiff(1900:2100,seq(1905,2105,10)),expand=c(0,0)) +
  scale_y_continuous(limits=c(14.5,44.5), breaks=seq(15,44,5), minor_breaks=setdiff(15:44,seq(15,44,5)), expand=c(0,0)) +
  geom_vline(xintercept=cohrange+0.5,color="lightgray",size=0.1) +
  geom_hline(yintercept=15.5:43.5,color="lightgray",size=0.1) +
  geom_step(data=stepfunc(14,max(cohrange)-30), mapping=aes(x=x,y=y,fill=NULL),direction="hv") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text("Calibri"),
        axis.text=element_text(color="black"),
        axis.title=element_text(face="bold"),
        axis.title.y = element_text(margin = margin(r = unit(5, "pt"))),
        plot.title=element_text(hjust=0.5,size=12),
        plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
  labs(title = expression(bold("b. Canada conditional parity-specific fertility rates")))

p4 <- ggplot(alldata %>% filter(par>0,Cohort %in% cohrange,Age+Cohort<=joy) %>%
               mutate(par2=paste0("Parity ",c(0,1,2,"3+")[as.numeric(par)])),
             aes(x=Age, y=x, color=Cohort)) +
  geom_line(aes(group=Cohort),size=0.3) +
  labs(x = "Age", y = "Rate", color="Cohort")+
  facet_wrap(~par2,nrow=1) +
  coord_cartesian(ylim=c(0,0.3)) +
  scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 7.5,barwidth=1,frame.colour="black",ticks.colour="black"), breaks = seq(1925,2005,10), limits=range(cohrange)) +
  scale_x_continuous(breaks=seq(15,40,5),minor_breaks=agerange,expand=c(0.02,0.02)) +
  scale_y_continuous(breaks=seq(0,0.3,0.1)) +
  theme_bw() + theme(text = element_text("Calibri"),
                     axis.text=element_text(color="black"),
                     axis.title=element_text(face="bold"),
                     plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"))

png(file=paste0("plots/",ctry,"_rates_curves_p04_",min(cohrange),"_",max(cohrange),"_redo.png"),
    width=20,height=20,units="cm",res=400,family="Calibri")
print({plot_grid(p1_p2, p3, p4, ncol=1, rel_heights = c(0.35, 0.65/2, 0.65/2)) +
    theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))})
dev.off()

tiff(file=paste0("plots/",ctry,"_rates_curves_p04_",min(cohrange),"_",max(cohrange),"_redo.tif"),
    width=20,height=20,units="cm",res=600,family="Calibri",compression="none")
print({plot_grid(p1_p2, p3, p4, ncol=1, rel_heights = c(0.35, 0.65/2, 0.65/2)) +
    theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))})
dev.off()
