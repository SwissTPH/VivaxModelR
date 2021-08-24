#=============================================================================
# 			VIVAX ODE MODELS
#
#  Validity conditions
#
#
#  created by : Clara Champagne (clara.champagne@swisstph.ch)
#  originally created : 2020
#
#
#
#=============================================================================
rm(list=ls())
# charge libraries
library(VivaxModelR)
library(ggplot2)
library(dplyr)
library(viridis)

my_directory="./"
dirPlots=my_directory

##################################################################
# VALIDITY CONDITIONS
##################################################################

full=list()
full$myalpha=c(0,0.4,0.8)
full$mybeta=c(0,0.4,0.8)
full$myincidence=seq(1,500,1)/1000/365
full$myp=seq(0,99,0.1)/100
simul_db = expand.grid( full )
simul_db$ID=rownames(simul_db)

my.solve.lamba=function(incidence, myalpha, myp, mybeta){
  my.lambda=solve_lambda_vivax(h=incidence, r=1/60,  gamma=1/223 ,  f=1/72,
                               alpha=myalpha , beta=mybeta,  rho=1, p=myp,omega = 1)
  my.lambda=ifelse( (length(my.lambda)==0), NA, my.lambda)

  return(my.lambda)
}


for(i in 1:nrow(simul_db)){
  print(i)
  simul_db$lambda[i]=my.solve.lamba(incidence= simul_db$myincidence[i],
                                    myalpha=simul_db$myalpha[i],
                                    myp=simul_db$myp[i],
                                    mybeta=simul_db$myb[i])
}


simul_db %>% filter(myalpha %in% c(0,0.4,0.8), mybeta %in% c(0,0.4,0.8)) %>%
  mutate(alpha=paste0("alpha:",myalpha),beta=paste0("beta:",mybeta)) %>%
  ggplot()+
  geom_tile(aes(x=myp, y=myincidence*365*1000, fill=lambda))+
  scale_fill_viridis(option="viridis", na.value="lightgrey")+
  labs(y="Incidence (per 1000 person-year)", x="Proportion of imported cases (p)", fill=expression(lambda))+
  facet_grid(beta ~ alpha, labeller=label_parsed )+
  theme(legend.title  = element_text(size=25))+
  theme(strip.background =element_rect(fill="white", color="black"),
        strip.text = element_text(size=20))+
  theme(axis.title  = element_text(size=20)  ,axis.text = element_text(size=15))+
  theme(panel.background = element_blank())
ggsave(file.path(dirPlots, "validity_conditions.pdf"), width=14, height = 7)
ggsave(file.path(dirPlots, "validity_conditions.png"), width=14, height = 7)
