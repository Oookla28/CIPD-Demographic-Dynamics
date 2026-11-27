###Plot observed rates
##Read in rates
load(paste0("data/",ctry,"_rates.RData"))


##Plot rates
supp.labs1 <- c("ASFR","Parity 0","Parity 1","Parity 2","Parity 3+")
names(supp.labs1) <- c(0,1,2,3,4)

png(file=paste0("plots/",ctry,"_rates_1946_2000.png"),width=20,height=6.5,units="cm",res=400,family="Calibri")
print(ggplot(alldata %>% filter(Cohort %in% 1946:2000),
             aes(x=Age, y=x, color = Cohort)) +
        geom_line(aes(group=Cohort),lwd=0.3) +
        labs(x = "Age", y = "Rate", color = "Cohort") +
        scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks = seq(1800,2000,10), limits=c(1946,2000)) +
        scale_y_continuous(breaks=seq(0,0.5,0.1)) +
        scale_x_continuous(breaks=seq(15,45,5), minor_breaks=setdiff(15:45,seq(15,45,5)), expand=c(0.02,0.02)) +
        coord_cartesian(ylim=c(0,0.5)) +
        theme_bw() + theme(text = element_text("Calibri")) +
        facet_grid(~par,labeller=labeller(par=supp.labs1)))
dev.off()


##Plot rates, highlighting cohorts with accuracy issues
if (!is.null(coh_iss[[ctry]]) & length(intersect(1946:2000, coh_iss[[ctry]]))>0) {
  tmp <- alldata %>% filter(Cohort %in% 1946:2000)
  tmp %<>% mutate(iss=ifelse(Cohort %in% coh_iss[[ctry]],T,F))
  png(file=paste0("plots/",ctry,"_rates_1946_2000_iss.png"),width=20,height=6.5,units="cm",res=400,family="Calibri")
  print(ggplot(tmp, aes(x=Age, y=x, color = Cohort)) +
          geom_line(data=tmp %>% filter(iss==T), aes(group=Cohort),lwd=0.3, linetype=2) +
          geom_line(data=tmp %>% filter(iss==F), aes(group=Cohort),lwd=0.3) +
          labs(x = "Age", y = "Rate", color = "Cohort") +
          scale_color_gradientn(colours=viridis(100)[100:1], guide = guide_colorbar(barheight = 10,frame.colour="black",ticks.colour="black"), breaks = seq(1800,2000,10), limits=c(1946,2000)) +
          scale_y_continuous(breaks=seq(0,0.5,0.1)) +
          scale_x_continuous(breaks=seq(15,45,5), minor_breaks=setdiff(15:45,seq(15,45,5)), expand=c(0.02,0.02)) +
          coord_cartesian(ylim=c(0,0.5)) +
          theme_bw() + theme(text = element_text("Calibri")) +
          facet_grid(~par,labeller=labeller(par=supp.labs1)))
  dev.off()
}


