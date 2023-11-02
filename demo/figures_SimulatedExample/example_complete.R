################################################################
#       Simulation of intervention effects with VivaxModelR
#
#         Clara Champagne, June 2022
#
################################################################

rm(list=ls())
library(VivaxModelR)
library(AnophelesModel)  #devtools::install_github("SwissTPH/AnophelesModel")
library(ggplot2)
library(dplyr)
library(tidyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions_IRS_parameters.R")
source("helper_functions.R")
output_dir=(".")


####################################################################################################
# Step 1: Specify baseline parameters for all three areas and calibrate the model
####################################################################################################

PQ_success_proba=0.76
# create imaginary settings
mydata=data.frame(id=c("Area1","Area2", "Area3"),
                  cases=c(6,95, 540),
                  cases_import=c(1,0, 27),
                  population=c(5000, 5000, 5000),
                  alpha=0.5*0.95*0.7,
                  beta=c(0.86*PQ_success_proba,0.86*PQ_success_proba,0.86*PQ_success_proba),
                  rho=0.5*0.95*0.7,
                  omega=1,
                  sigma=1/10, N=5000,
                  iota=50/7/5000,nu=5, eta=0.95, tau=5)  # treatment delays last on average 15 days
mydata$cases_local=mydata$cases-mydata$cases_import

nsim=1000
mydata_uncertainty=sample_uncertainty_incidence_import(mydata, ndraw = nsim)%>%
  dplyr::mutate(iter=rep(1:nsim, each=nrow(mydata)),
                id0=id, id=paste0(id,"_",iter))

mydata$incidence=1000*mydata$cases/mydata$population
mydata$prop_import=(mydata$cases-mydata$cases_local)/mydata$cases
mydata$h=incidence_year2day(mydata$incidence)

# calibrate the model
mydata_withR0RC=calibrate_vivax_equilibrium(df=mydata, f=1/72, gamma=1/223, r=1/60, h.cutoff = 1e-08, return.all = TRUE,
                                            delay = TRUE, rcd=T)


mydata_withR0RC_uncertainty=calibrate_vivax_equilibrium(df=mydata_uncertainty, f=1/72, gamma=1/223, r=1/60, h.cutoff = 1e-08, return.all = TRUE,
                                                        delay = TRUE, rcd=T)


mydata_withR0RC_uncertainty %>%filter(lambda>0) %>%
  group_by(id0) %>%
  summarise(R0_median=median(R0), R0_q025=q025(R0), R0_q975=q975(R0),
            Rc_median=median(Rc), Rc_q025=q025(Rc), Rc_q975=q975(Rc))
####################################################################################################
# Step 2: Specific functions for interventions
####################################################################################################

### RCD
# Function for time varying tau from Chitnis et al. 2019
my_tau=function(pr){
  return(varying_tau(nu=10, pr=pr, N=5000))
}


### IRS
mosquito_species="Anopheles albimanus"

# Bionomics from Briet 2019
bionomics_albimanus=data.frame(species_name=mosquito_species,
                               M=0.484,Chi=.054,A0=.405,zeta.3=1,td=.33,tau=3,ts=10,to=5,endophily=.34,endophagy=.4)
# Rhythms  from Briet 2019
combination_activities=list(HBI=c(0.22,0.21,0.22,0.10,0.13,0.17,
                                  0.08,0.12,0.03,0.12,0.17,0.18,0.25),
                            HBO=c(0.25,0.35,0.37,0.33,0.36,0.32,
                                  0.13,0.14,0.09,0.15,0.36,0.23,0.25),
                            humans_in_bed=c(5.628637, 14.084496, 28.558507, 47.761203, 67.502964, 83.202297, 92.736838, 96.734455, 96.637604,
                                            92.983136, 85.264515, 73.021653, 57.079495)/100) # Briet et al 2019, Knutson et al. 2014)
combination_activities$humans_indoors=combination_activities$humans_in_bed

# Calculate impact with AnophelesModel package
int_list = list(IRS_Bendiocarb=list(id="IRS",description="IRS albimanus Bendiocarb",parameterisation="IRS06"))

vec_params = def_vector_params(mosquito_species=mosquito_species,vector_table = bionomics_albimanus)
hosts_params = def_host_params()
model_params = build_model_obj(vec_params, hosts_params,activity=combination_activities, total_pop =2000)
intervention_vec = def_interventions_effects(intervention_list=int_list,model_p=model_params,num_ip_points =100)
impacts = calculate_impact(intervention_vec, coverage_vec = seq(0,.99, by=.01), model_p = model_params,Nv0 = 10000, num_ip_points = 100)

# Extract database of vector effects
exp_df=data.frame(id=c("IRS 60%"),int=c("IRS albimanus Bendiocarb"),cov=c(.60),rounds=c(20),time_btw=c(1)) %>%
  mutate(f_name=paste0(id,"\n",int,"\n",rounds," rounds, ",cov*100,"% coverage"))
IRS_bendiocarb_60_outdoor=format_VCimpact_omega(impacts=impacts,exp_df = exp_df, max_time=20)


####################################################################################################
# Step 3: Define and simulate interventions (Case management, RCD, IRS)
####################################################################################################


# simulate interventions
intervention_object0=list(intervention_name="Baseline", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "rho.new"=NA, "iota.new"=NA, "nu.new"=NA, "tau.new"=NA, "eta.new"=NA)
intervention_objectCM=list(intervention_name="Improved case management", "alpha.new"=0.8*0.7*0.95, "beta.new"=0.9*0.76, "omega.new"=NA, "sigma.new"=1/5, "rho.new"=0.8*0.7*0.95, "iota.new"=NA, "nu.new"=NA, "tau.new"=NA, "eta.new"=NA )
intervention_objectRCD=list(intervention_name="RCD1", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA,"sigma.new"=NA, "iota.new"=50/7/5000, "nu.new"=10, "tau.new"=5, "eta.new"=0.95 )
intervention_objectRCD2=list(intervention_name="RCD2", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "rho.new"=NA,"sigma.new"=NA, "iota.new"=50/7/5000, "nu.new"=10, "tau.new"=my_tau, "eta.new"=0.95 )
intervention_objectIRS=list(intervention_name="Foutdoors", "alpha.new"=NA, "beta.new"=NA, "rho.new"=NA,"omega.new"=IRS_bendiocarb_60_outdoor, "sigma.new"=NA, "iota.new"=NA, "nu.new"=NA, "tau.new"=NA, "eta.new"=NA )


my_intervention_list=list(intervention_object0,
                          intervention_objectCM, intervention_objectRCD, intervention_objectRCD2, intervention_objectIRS)

simulation_model= simulate_vivax_interventions(df=mydata_withR0RC, my_intervention_list,
                                               delay = TRUE, rcd = TRUE, referral = FALSE,
                                               maxtime = 6*365,  year = FALSE, rcd_at_baseline = TRUE)

simulation_model_uncertainty= simulate_vivax_interventions(df=mydata_withR0RC_uncertainty%>%
                                                             dplyr::filter(lambda>0) , my_intervention_list,
                                                           delay = TRUE, rcd = TRUE, referral = FALSE,
                                                           maxtime = 6*365,  year = TRUE, rcd_at_baseline = T)


# postprocess and visualise the results
simulation_model_uncertainty$area=substring(simulation_model_uncertainty$id,1,5)
simulation_model_uncertainty$simulid=substring(simulation_model_uncertainty$id,7,9)

simulation_model_uncertainty_central=
  compute_centrality_score_all_nothing(data=simulation_model_uncertainty,
                                       n_curves_sample=50,n_samples=100,nsim=500, myseed=1)

simulation_model_uncertainty_reformat=simulation_model_uncertainty%>%
  ungroup()%>%
  group_by(time,intervention, area) %>%
  summarise(incidence_median=median(incidence)) %>%
  left_join(simulation_model_uncertainty_central)%>%
  left_join(mydata %>% mutate(time=0, area=id, incidence_data=incidence) %>% select(time, area, incidence_data) )%>%
  mutate(Intervention=factor(intervention, levels = c("Baseline",
                                                      "Improved case management","RCD1", "RCD2", "Foutdoors" ),
                             labels = c("Baseline",
                                        "Improved case management","Increased RCD", "Increased RCD", "IRS" )),
         Intervention2=factor(intervention, levels = c("Baseline",
                                                       "Improved case management","RCD1", "RCD2", "Foutdoors" ),
                              labels = c("Baseline",
                                         "Improved case management","Increased RCD", "Increased RCD2", "IRS" )),
         tau=ifelse(intervention=="RCD2", "time varying", "fixed"),
         area=gsub("Area", "Area ",area)
  )


plot_cm_uncertainty=simulation_model_uncertainty_reformat %>%
  ggplot()+
  geom_line(aes(x=time/365,y=incidence_median, color=Intervention, linetype=tau), lwd=1)+
  geom_ribbon(aes(x=time/365,ymin=incidence_inf,ymax=incidence_sup, fill=Intervention2), alpha=0.1)+
  geom_point( aes(x=time,y=incidence_data), color="black")+
  facet_wrap(area ~., scales="free")+xlab("Year since baseline")+labs(color="", fill="", linetype="Targeting ratio", y="Reported annual incidence (per 1000 inh.)")+
  scale_color_manual(values=c( "black","red","orange","dodgerblue","darkblue"))+
  scale_fill_manual(values=c( "black","red","orange","orange","dodgerblue","darkblue"))+
  theme_bw()+
  theme(legend.position = "bottom",   axis.text=element_text(size=12), legend.title =element_text(size=18),
        axis.title=element_text(size=18), legend.text =element_text(size=18), strip.text =element_text(size=18) ,
        legend.box="vertical")+ guides(fill="none")+
  xlim(0,5)+ylim(0, NA)

plot_cm_uncertainty
ggsave(plot_cm_uncertainty, filename = file.path(output_dir, "cm_improvements_uncertainty_incidence.jpg"), width=12, height = 7)


####################################################################################################
# Step 4: Define and simulate interventions (MDA)
####################################################################################################

# define MDA interventions
intervention_objectMDA_0=list(intervention_name="baseline_0_0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0, "MDAp_length.new"=30, "MDArad_cure.new"=0, "rho.new"=NA , "iota.new"=NA, "nu.new"=NA, "tau.new"=NA, "eta.new"=NA )
intervention_objectMDA_80_0=list(intervention_name="cov_80_0", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.8, "MDAp_length.new"=30, "MDArad_cure.new"=0, "rho.new"=NA  , "iota.new"=NA, "nu.new"=NA, "tau.new"=NA, "eta.new"=NA)
intervention_objectMDA_80_50=list(intervention_name="cov_80_50", "alpha.new"=NA, "beta.new"=NA, "omega.new"=NA, "sigma.new"=NA, "MDAcov.new"=0.8, "MDAp_length.new"=30, "MDArad_cure.new"=0.5*PQ_success_proba, "rho.new"=NA , "iota.new"=NA, "nu.new"=NA, "tau.new"=NA, "eta.new"=NA )

my_intervention_list_MDA=list(intervention_objectMDA_0, intervention_objectMDA_80_0,
                              intervention_objectMDA_80_50)

# simulate MDA interventions

# before MDA
simulation_model_day= simulate_vivax_interventions(df=mydata_withR0RC, my_intervention_list_MDA, year=T,delay=T, maxtime = 365, sto = T, runs = 5000, rcd_at_baseline = T, rcd = T)
# round 1
simulation_model2_MDA= simulate_vivax_interventions(df=mydata_withR0RC, my_intervention_list_MDA, previous_simulation = simulation_model_day, year=T, mda = T, delay=T, maxtime = 1*365, sto = T, runs =5000, rcd_at_baseline = T, rcd = T)
# round 2
simulation_model3_MDA= simulate_vivax_interventions(df=mydata_withR0RC, my_intervention_list_MDA, previous_simulation = simulation_model2_MDA, year=T, mda = T, delay=T, maxtime = 1*365, sto = T, runs = 5000, rcd_at_baseline = T, rcd = T)
# round 3
simulation_model4_MDA= simulate_vivax_interventions(df=mydata_withR0RC, my_intervention_list_MDA, previous_simulation = simulation_model3_MDA, year=T, mda = T, delay=T, maxtime = 4*365, sto = T, runs = 5000, rcd_at_baseline = T, rcd = T)

# postprocess and visualise the results
simulation_model4_MDA$simulid=simulation_model4_MDA$run

simulation_MDA_central=
  compute_centrality_score_all_nothing(data=simulation_model4_MDA,
                                       n_curves_sample=50,n_samples=100,
                                       nsim=5000, myseed=1)

simulation_MDA=simulation_model4_MDA %>% mutate(area=id) %>%
  group_by(intervention, time, area) %>%
  summarise(incidence_median=median(incidence)) %>%
  left_join(simulation_MDA_central) %>%
  left_join(mydata %>% mutate(time=0, incidence_data=incidence, area=id) %>% select(time, area, incidence_data) )%>%
  separate(intervention, into = c("type", "cov", "radcure"), sep='_') %>%
  mutate(intervention=ifelse(cov==0, "Baseline", ifelse(radcure==0, "MDA, no PQ", "MDA, with PQ")),
         area=gsub("Area", "Area ",area))


plot_mda=simulation_MDA%>%
  ggplot()+
  geom_line(aes(x=time/365,y=incidence_median, color=intervention), lwd=1)+
  geom_ribbon(aes(x=time/365,ymin=incidence_inf, ymax=incidence_sup, fill=intervention), alpha=0.2)+
  geom_point(aes(x=time/365,y=incidence_data), color="black")+
  facet_wrap(area ~., scales="free")+xlab("Year since baseline")+ylab("Reported annual incidence (per 1000 inh.)")+
  scale_color_manual(values=c("black", "lightgreen", "forestgreen", "darkgreen"), name="")+
  scale_fill_manual(values=c("black", "lightgreen", "forestgreen", "darkgreen"), name="")+
  theme_bw()+
  theme(legend.position = "bottom",   axis.text=element_text(size=12), legend.title =element_text(size=18),
        axis.title=element_text(size=18), legend.text =element_text(size=18), strip.text =element_text(size=18) , legend.box="vertical")
plot_mda
ggsave(plot_mda, filename = file.path(output_dir, "mda_sim2.jpg"), width=12, height = 7)


# Calculate elimination probability
simuldat_elim=simulation_model4_MDA %>%
  mutate(year=ceiling(time/365)) %>%
  group_by(intervention, year, id, run) %>%
  summarise(h=sum(h), hh=sum(hh), hhl=sum(hhl), hl=sum(hl))


simuldat_elim7=simuldat_elim %>%
  filter(year==7) %>%
  group_by(id, intervention) %>%
  summarise(elim=sum(hh==0)/sum(!is.na(hh)), tot=sum(!is.na(hh)), tot_elim=sum(hh==0),
            elim_loc=sum(hhl==0)/sum(!is.na(hhl)), tot_loc=sum(!is.na(hhl)), tot_elim_loc=sum(hhl==0))
simuldat_elim7 # elimination probability by year 7

