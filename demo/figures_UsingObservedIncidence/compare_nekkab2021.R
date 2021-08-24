#=============================================================================
# 			VIVAX ODE MODELS
#
#  Comparing the outputs of the simple ODE model
#  with the effect sizes presented in Nekkab et al. 2020
#
#
#  created by : Clara Champagne (clara.champagne@swisstph.ch)
#  originally created : 2020
#
#=============================================================================
rm(list=ls())
# charge libraries

library(VivaxModelR)
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)


mydata0=data.frame(incidence=c(23,112,267)) %>%
  mutate(h=incidence/1000/365,
         lambda=-1, I=NA, prop_import=0,
         rho=c(0.18,0.13,0.08) ,
         PQ_effect=c(0.431,0.429,0.422),
         TQ_effect=c(0.59,0.605,0.619),
         alpha=0.95*rho,
         scenario="S1")

mydata=rbind(mydata0,
             mydata0 %>% mutate(scenario="S2",
                                PQ_effect=c(0.431*0.9/0.667,0.429*0.9/0.667,0.422*0.9/0.667)),
             mydata0 %>%  mutate(scenario="S3",
                                 PQ_effect=c(0.431*0.3/0.667,0.429*0.3/0.667,0.422*0.3/0.667))) %>%
  mutate(beta=PQ_effect, id=paste(scenario, incidence, sep="_"))

r=1/60
gamma=1/383
f=1/69

mydata_lambda=calculate_r0_rc_fromdata(df=mydata, f=f, gamma =gamma, r=r,
                                   return.all = T)



simul.PQ=simulate_from_equilibrium_fromdata(mydata_lambda %>% mutate(alpha.old=alpha,
                                                                     beta.old=beta),
                                        f=f, gamma =gamma, r=r,
                                        maxtime=2000,year=T)
names(simul.PQ)=c("time",paste0(names(simul.PQ)[c(-1,-10)], "_PQ"),"id")

simul.TQ=simulate_from_equilibrium_fromdata(mydata_lambda %>% mutate(alpha.old=alpha,
                                                                     beta.old=beta,
                                                                     beta.new=TQ_effect),
                                        f=f, gamma =gamma, r=r,
                                        maxtime=2000,year=T)
names(simul.TQ)=c("time",paste0(names(simul.TQ)[c(-1,-10)], "_TQ"), "id")

db_compare=simul.PQ %>% left_join(simul.TQ) %>%
  mutate(scenario=strtrim(id,2), incidence=str_sub(id,start=4))
db_compare2 = db_compare%>% mutate(effect_size=-(h_TQ-h_PQ)/h_PQ,
                                   year=time/365+2020)


db_compare2%>% filter(time == 365*5) %>%
  mutate(scenario=factor(scenario, levels=c("S3", "S1", "S2")),
         incidence_lab=factor(paste0("Annual incidence = ",incidence), levels=c("Annual incidence = 23", "Annual incidence = 112", "Annual incidence = 267"))
  ) %>%
  ggplot() +
  geom_bar(aes(x=scenario, y=effect_size, fill=scenario),stat='identity') +
  coord_flip() +
  scale_fill_manual(values=c("darkolivegreen3", "goldenrod1","forestgreen"))+
  geom_text(aes(x=scenario,y=effect_size,label=paste0(round(effect_size,digits = 2)*100,"%")),vjust=0, hjust=-0.1) +
  facet_wrap( incidence_lab ~., ncol=1)
