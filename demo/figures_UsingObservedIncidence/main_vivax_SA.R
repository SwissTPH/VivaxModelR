#=============================================================================
# 			VIVAX ODE MODELS
#
#  Sensitivity analysis on model parameters
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
library(tidyr)
library(sensitivity)
library(lhs)
library(dplyr)
library(cowplot)

my_directory="./"
source(file.path(my_directory,"demo/figures_UsingObservedIncidence//functions_for_SA.R"))
dirPlots=my_directory



# LHS design
mysize=10000
X1_u = data.frame(as.matrix( maximinLHS( n = mysize , k =  6) ))
X2_u = data.frame(as.matrix( maximinLHS( n = mysize , k = 6 ) ))
colnames(X1_u) = c("r", "gamma", "f","CM" ,"prop_radCure", "rho")
colnames(X2_u) = c("r", "gamma", "f","CM", "prop_radCure", "rho")
X1 = X1_u %>% mutate(r=qunif(r, min = 1/85, max =1/35),
                     gamma=qunif(gamma, min = 1/500, max =1/200),
                     f=qunif(f, min = 1/125, max =1/40),
                     CM=qunif(CM, min = 0.1, max =0.8), # parameterise in terms of alphatilde, alphafinal, proplivercleared
                     prop_radCure=qunif(prop_radCure, min = 0, max =1),
                     rho=qunif(rho, min = 0.1, max = 0.8))

X2 = X2_u %>% mutate(r=qunif(r, min = 1/85, max =1/35),
                     gamma=qunif(gamma, min = 1/500, max =1/200),
                     f=qunif(f, min = 1/125, max =1/40),
                     CM=qunif(CM, min = 0.1, max =0.8), # parameterise in terms of alphatilde, alphafinal, proplivercleared
                     prop_radCure=qunif(prop_radCure, min = 0, max =1),
                     rho=qunif(rho, min = 0.1, max = 0.8))





h.values=c(5,50,100)/365/1000
p.values=c(0,0.1)


#=============================================================================
#        SA with CM
#=============================================================================

myvars=c("r", "gamma", "f", "prop_radCure", "rho", 'CM')


################
# UNCERTAINTY AND PRCC

all_uncertainty=data.frame()

for(i in 1:3){
  my.h.fixed=h.values[i]
  print(my.h.fixed)
  for(j in 1:2){
    my.p.fixed=p.values[j]
    print(my.p.fixed)

    # calculate uncertainty
    uncertainty_cm_lambda=sensi_fit_cm_import_lambda(X1[myvars])
    uncertainty_cm_r0=sensi_fit_cm_import(X1[myvars])
    uncertainty_cm_rc=sensi_fit_cm_import_RC(X1[myvars])
    uncertainty_cm_relapse=sensi_fit_cm_import_relapses(X1[myvars])

    this_uncertainty=data.frame(lambda=uncertainty_cm_lambda,R0=uncertainty_cm_r0,
                                Rc=uncertainty_cm_rc, prop_relapse=uncertainty_cm_relapse,
                                myh=my.h.fixed, myp=my.p.fixed)

    all_uncertainty=rbind(all_uncertainty, this_uncertainty)
  }
}

all_uncertainty_plot=all_uncertainty %>%
  pivot_longer(cols=c(lambda, R0, Rc, prop_relapse)) %>%
  filter(name !="lambda")%>%
  mutate(name=factor(name, levels=c("lambda","R0", "Rc", "prop_relapse"), labels = c("lambda","R0", "Rc","Prop. of relapses")),
         incidence=factor(myh*365*1000, levels = c("5", "50", "100"))) %>%
  ggplot( )+
  geom_boxplot(aes(y=value, x=incidence, group=interaction(name, myp,incidence), color=paste0(myp*100,"%")), size=1)+theme_bw()+
  facet_grid(name ~ ., scales="free", switch="y")+
  scale_color_manual(values=c( "cyan3","darkblue"))+
  theme( panel.border = element_rect(colour = "black", fill=NA), legend.position = "none")+
  labs(y="", x="Incidence per 1000 person-year", color="Prop. of imports")+
  theme(strip.background =element_rect(fill="white", color="white"), strip.text = element_text(size=12))+
  theme(strip.placement = "outside", legend.position="bottom")


#################
# SOBOL
#################

all_sobol=data.frame()
all_sobol_format=data.frame()

for(i in 1:3){
  my.h.fixed=h.values[i]
  print(my.h.fixed)
  for(j in 1:2){
    my.p.fixed=p.values[j]
    print(my.p.fixed)

    # calculate sobol indices for R0
    sobol_cm=soboljansen(sensi_fit_cm_import, X1[myvars], X2[myvars], nboot = 100 )
    for_plot_effects=cbind(sobol_cm$X, y=sobol_cm$y)  %>% mutate(outcome="R0", incidence=my.h.fixed, p=my.p.fixed)
    all_sobol=rbind(all_sobol,for_plot_effects)
    sobol_ref_cm=reformat_sobol_output(sobol_cm) %>% mutate(outcome="R0", incidence=my.h.fixed, p=my.p.fixed)
    all_sobol_format=rbind(all_sobol_format,sobol_ref_cm)

    # calculate sobol indices for Rc
    sobol_cm_Rc=soboljansen(sensi_fit_cm_import_RC, X1[myvars], X2[myvars], nboot = 100 )
    for_plot_effects=cbind(sobol_cm_Rc$X, y=sobol_cm_Rc$y)  %>% mutate(outcome="Rc", incidence=my.h.fixed, p=my.p.fixed)
    all_sobol=rbind(all_sobol,for_plot_effects)
    sobol_ref_cm=reformat_sobol_output(sobol_cm_Rc) %>% mutate(outcome="Rc", incidence=my.h.fixed, p=my.p.fixed)
    all_sobol_format=rbind(all_sobol_format,sobol_ref_cm)


    # calculate sobol indices for Relapse
    sobol_cm_relapse=soboljansen(sensi_fit_cm_import_relapses, X1[myvars], X2[myvars], nboot = 100 )
    for_plot_effects=cbind(sobol_cm_relapse$X, y=sobol_cm_relapse$y)  %>% mutate(outcome="Prop. of relapses", incidence=my.h.fixed, p=my.p.fixed)
    all_sobol=rbind(all_sobol,for_plot_effects)
    sobol_ref_cm=reformat_sobol_output(sobol_cm_relapse) %>% mutate(outcome="Prop. of relapses", incidence=my.h.fixed, p=my.p.fixed)
    all_sobol_format=rbind(all_sobol_format,sobol_ref_cm)

  }
}



#########
# SOBOL PLOTS

all_sobol_format$outcome[all_sobol_format$outcome=="Prop. of relapses"]="Prop. of\nrelapses"

sobol5=plot_sobol(all_sobol_format%>%filter(incidence==5/365/1000))
sobol50=plot_sobol(all_sobol_format%>%filter(incidence==50/365/1000))
sobol100=plot_sobol(all_sobol_format%>%filter(incidence==100/365/1000))


all_sobol$outcome[all_sobol$outcome=="Prop. of relapses"]="Prop. of\nrelapses"

dplot_5_minmax= plot_effect_minmax(all_sobol %>% filter(incidence==5/365/1000))
dplot_50_minmax= plot_effect_minmax(all_sobol %>% filter(incidence==50/365/1000))
dplot_100_minmax= plot_effect_minmax(all_sobol %>% filter(incidence==100/365/1000))

plot5_minmax=plot_grid(sobol5, dplot_5_minmax, labels = c("B.", "C."), rel_widths = c(1,2))

plot50_minmax=plot_grid(sobol50, dplot_50_minmax, labels = c("B.", "C."), rel_widths = c(1,2))

plot100_minmax=plot_grid(sobol100, dplot_100_minmax, labels = c("B.", "C."), rel_widths = c(1,2))


plot5_full_minmax=plot_grid(all_uncertainty_plot, plot5_minmax, labels = c("A.", ""), ncol=1, rel_heights = c(1.8,2))
ggsave(plot=plot5_full_minmax,file.path(dirPlots, "full_full_sobol_inc5_minmax.png"), width=12, height = 12)

plot50_full_minmax=plot_grid(all_uncertainty_plot, plot50_minmax, labels = c("A.", ""), ncol=1, rel_heights = c(1.8,2))
ggsave(plot=plot50_full_minmax,file.path(dirPlots, "full_full_sobol_inc50_minmax.png"), width=12, height = 12)

plot100_full_minmax=plot_grid(all_uncertainty_plot, plot100_minmax, labels = c("A.", ""), ncol=1, rel_heights = c(1.8,2))
ggsave(plot=plot100_full_minmax,file.path(dirPlots, "full_full_sobol_inc100_minmax.png"), width=12, height = 12)

