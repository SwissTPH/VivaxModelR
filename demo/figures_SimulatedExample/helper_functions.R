q025=function(x){return(quantile(x, probs=0.025))}
q975=function(x){return(quantile(x, probs=0.975))}
q90=function(x){return(quantile(x, probs=0.9))}


compute_centrality_score_all_nothing=function(data=simulation_model_uncertainty,
                                              n_curves_sample=10,
                                              n_samples=100,nsim=500, myseed=1){

  set.seed(myseed)

  mydata= data%>%
    select(id, simulid, time, area, intervention, incidence,I)
  final.score=mydata%>%
    mutate(centrality_score_inc=0,centrality_score_pr=0)
  for(i in 1:n_samples){
    #pick curves and create enveloped
    mysample=sample(1:nsim,n_curves_sample)
    df=mydata %>% dplyr::filter(simulid %in% mysample) %>%
      group_by(time, intervention, area)%>%
      mutate(env_sup_inc=max(incidence), env_inf_inc=min(incidence),
             env_sup_pr=max(incidence), env_inf_pr=min(incidence))%>%
      select(time, intervention, area, env_sup_inc, env_inf_inc,env_sup_pr, env_inf_pr)%>%
      ungroup() %>% data.frame() %>%
      unique()

    # check if each simulation is in the envelope
    this.score=mydata %>%
      left_join(df) %>%
      mutate(in_envelope_inc=ifelse(incidence>=env_inf_inc & incidence <= env_sup_inc, 0, 1),
             in_envelope_pr=ifelse(I>=env_inf_pr & I <= env_sup_pr, 0, 1))%>%
      group_by(intervention, area, simulid) %>%
      summarise(centrality_score_inc_update=as.numeric(sum(in_envelope_inc)==0),
                centrality_score_pr_update=as.numeric(sum(in_envelope_pr)==0))

    final.score= this.score %>%
      right_join(final.score) %>%
      mutate(centrality_score_inc=centrality_score_inc+centrality_score_inc_update,
             centrality_score_pr=centrality_score_pr+centrality_score_pr_update) %>%
      select(-centrality_score_inc_update,-centrality_score_pr_update)
    # print(max(final.score$centrality_score))
  }

  envelope90=
    final.score %>%
    group_by(intervention, area) %>%
    mutate(central=q90(centrality_score_inc),
           in_central=(centrality_score_inc>=central)) %>%
    filter(in_central==TRUE)%>%
    ungroup() %>%
    group_by(area, time, intervention) %>%
    dplyr::summarise(incidence_inf=min(incidence),
                     incidence_sup=max(incidence))
  return(envelope90)
}










# simulation_model_uncertainty_reformat=simulation_model_uncertainty%>%
#   mutate(Intervention=factor(intervention, levels = c("Baseline",
#                                                       "Improved case management","RCD1", "RCD2", "Foutdoors" ),
#                              labels = c("Baseline",
#                                         "Improved case management","Increased RCD", "Increased RCD", "IRS" )),
#          Intervention2=factor(intervention, levels = c("Baseline",
#                                                        "Improved case management","RCD1", "RCD2", "Foutdoors" ),
#                               labels = c("Baseline",
#                                          "Improved case management","Increased RCD", "Increased RCD2", "IRS" )),
#          tau=ifelse(intervention=="RCD2", "time varying", "fixed")) %>%
#   ungroup()%>%
#   group_by(time, Intervention,Intervention2, intervention, area, tau) %>%
#   summarise(incidence_mean=mean(incidence), I_mean=mean(I),
#             I_inf=q025(I), I_sup=q975(I),
#             incidence_inf=q025(incidence), incidence_sup=q975(incidence))


# plot_cm_uncertainty=truc2 %>%
#   ggplot()+
#   #geom_line(aes(x=time/365,y=I_mean*100, color=Intervention, linetype=tau), lwd=1)+
#   #geom_ribbon(aes(x=time/365,ymin=I_inf*100,ymax=I_sup*100, fill=Intervention), alpha=0.1)+
#   geom_line(aes(x=time/365,y=incidence_mean, color=Intervention, linetype=tau), lwd=1)+
#   geom_ribbon(aes(x=time/365,ymin=incidence_inf,ymax=incidence_sup, fill=Intervention2), alpha=0.1)+
#   geom_point(data=mydata %>% mutate(time=0, area=id), aes(x=time,y=incidence), color="black")+
#   facet_wrap(area ~., scales="free")+xlab("Year since baseline")+labs(color="", fill="", linetype="Targeting ratio", y="Annual reported incidence (per 1000 inh.)")+
#   scale_color_manual(values=c( "black","red","orange","dodgerblue","darkblue"))+
#   scale_fill_manual(values=c( "black","red","orange","orange","dodgerblue","darkblue"))+
#   theme_bw()+
#   theme(legend.position = "bottom",   axis.text=element_text(size=12), legend.title =element_text(size=18),
#         axis.title=element_text(size=18), legend.text =element_text(size=18), strip.text =element_text(size=18) , legend.box="vertical")+
#   xlim(0,5)+ylim(0, NA)

#plot_cm_uncertainty
#ggsave(plot_cm_uncertainty, filename = file.path(output_dir, "cm_improvements_uncertainty_incidence.jpg"), width=12, height = 7)
#
# plot_cm=simulation_model%>%
#   mutate(Intervention=factor(intervention, levels = c("Baseline",
#                                                       "Improved case management","RCD1", "RCD2", "Foutdoors" ),
#                              labels = c("Baseline",
#                                         "Improved case management","Increased RCD", "Increased RCD", "IRS" )),
#          tau=ifelse(intervention=="RCD2", "time varying", "fixed"),
#          area=gsub("Area", "Area ",id)) %>%
#   ggplot()+
#   geom_line(aes(x=time/365,y=I*100, color=Intervention, linetype=tau), lwd=1)+
#   facet_wrap(area ~., scales="free")+xlab("Year since baseline")+labs(color="", linetype="Targeting ratio", y="Prevalence (%)")+
#   scale_color_manual(values=c( "black","red","orange","dodgerblue","darkblue"))+
#   theme_bw()+
#   theme(legend.position = "bottom",   axis.text=element_text(size=12), legend.title =element_text(size=18),
#         axis.title=element_text(size=18), legend.text =element_text(size=18), strip.text =element_text(size=18) , legend.box="vertical")+
#   xlim(0,5)+ylim(0, NA)
#
# plot_cm
# #ggsave(plot_cm, filename = file.path(output_dir, "cm_improvements.jpg"), width=12, height = 7)
#
# simulation_model%>%
#   mutate(Intervention=factor(intervention, levels = c("Baseline",
#                                                       "Improved case management","RCD1", "RCD2", "Foutdoors" ),
#                              labels = c("Baseline",
#                                         "Improved case management","Increased RCD", "Increased RCD", "IRS" )),
#          tau=ifelse(intervention=="RCD2", "time varying", "fixed"),
#          area=gsub("Area", "Area ",id)) %>%
#   ggplot()+
#   geom_line(aes(x=time/365,y=incidence, color=Intervention, linetype=tau), lwd=1)+
#   facet_wrap(area ~., scales="free")+xlab("Year since baseline")+labs(color="", linetype="Targeting ratio", y="Prevalence (%)")+
#   scale_color_manual(values=c( "black","red","orange","dodgerblue","darkblue"))+
#   theme_bw()+
#   theme(legend.position = "bottom",   axis.text=element_text(size=12), legend.title =element_text(size=18),
#         axis.title=element_text(size=18), legend.text =element_text(size=18), strip.text =element_text(size=18) , legend.box="vertical")+
#   xlim(0,5)+ylim(0, NA)
#

## postprocess and visualise the results
# simuldat=simulation_model4_MDA %>%
#   group_by(intervention, time, id) %>%
#   summarise(Imean=mean(incidence), Imin=quantile(incidence, probs = 0.025), Imax=quantile(incidence, probs = 0.975)) %>%
#   separate(intervention, into = c("type", "cov", "radcure"), sep='_')
#
# simuldat2=simuldat %>% mutate(intervention=ifelse(cov==0, "Baseline", ifelse(radcure==0, "MDA, no PQ", "MDA, with PQ")))
# simuldat2$area=gsub("Area", "Area ",simuldat2$id)

#
# plot_mda=simuldat2%>%
#   ggplot()+
#   geom_line(aes(x=time/365,y=Imean, color=intervention), lwd=1)+
#   geom_ribbon(aes(x=time/365,ymin=Imin, ymax=Imax, fill=intervention), alpha=0.2)+
#   facet_wrap(area ~., scales="free")+xlab("Year since baseline")+ylab("Prevalence (%)")+
#   scale_color_manual(values=c("black", "lightgreen", "forestgreen", "darkgreen"), name="")+
#   scale_fill_manual(values=c("black", "lightgreen", "forestgreen", "darkgreen"), name="")+
#   theme_bw()+
#   theme(legend.position = "bottom",   axis.text=element_text(size=12), legend.title =element_text(size=18),
#         axis.title=element_text(size=18), legend.text =element_text(size=18), strip.text =element_text(size=18) , legend.box="vertical")
# plot_mda
# ggsave(plot_mda, filename = file.path(output_dir, "mda_sim.jpg"), width=12, height = 7)
#
