################################################################
#       Simulation of intervention effects with VivaxModelR
#
#         Jeanne Lemant, 2021-2022
#         Clara Champagne, 2023
#
################################################################

############
# This function converts the outputs of AnophelesModel into inputs for VivaxModelR
# It can be adapted by the user to match other specific uses cases

format_VCimpact_omega = function(impacts, exp_df,
                                 max_time=20){

  max_time=max_time*365


  time_impact_0=as.data.frame(impacts$interventions_vec[[1]]$effects$impact)
  colnames(time_impact_0)=paste0("cov_",impacts$interventions_vec[[1]]$coverages)

  # extract the duration of IRS06
  d= interventions_param$interventions_summary %>%
    filter(Parameterisation==impacts$interventions_vec[[1]]$parameterisation) %>%
    select(Duration) %>% pull()


  time_impact = time_impact_0 %>%
    mutate(time=seq(from=0,to=d,length.out = dim(time_impact_0)[1])) %>%
    pivot_longer(cols=starts_with("cov"),names_to="coverage"
                 ,names_prefix="cov_",values_to="impact") %>%
    mutate(coverage=as.numeric(coverage)
           ,intervention_name = impacts$interventions_vec[[1]]$description
           ,duration=d) %>%
    rowwise() %>%
    mutate(impact=max(impact,0)) %>%
    mutate(impact=min(impact,1)) %>%
    ungroup() %>%
    mutate(value=1-impact
           ,time=365*time
           ,round=1
           ,duration=duration*365
    )


  intervention_df = time_impact %>%
    filter(coverage==exp_df$cov) %>%
    mutate(id=exp_df$id
           ,nb_rounds=exp_df$rounds
           ,time_btw_rounds=exp_df$time_btw*365
           ,f_name=exp_df$f_name)

  if(max(intervention_df$time)<max_time){#need to extend until max time
    time_vector = sort(unique(intervention_df$time))
    time_interval=max(intervention_df$time)-time_vector[length(time_vector)-1] #in days
    additional_time=seq(max(intervention_df$time),max_time,by = time_interval)[-1]
    extended_intervention_df=NULL

    given_cov =intervention_df$coverage
    specific_df = intervention_df
    additional_df=cbind.data.frame(time=additional_time,specific_df[1,] %>% select(-time)) %>%
      mutate(value=1,round=Inf) #no vector control after the end of the duration
    extended_specific_df=rbind(specific_df,additional_df)

    if(unique(specific_df$nb_rounds)>1){#several rounds
      omega_intermediate_round=extended_specific_df %>%
        filter(time<=time_btw_rounds) %>% select(value) %>% pull()
      length_round = length(extended_specific_df %>%
                              filter(time<=time_btw_rounds & time<=duration) %>% select(time) %>% pull())
      length_break=length(additional_df %>% filter(time<=time_btw_rounds) %>% select(time) %>% pull())
      if(length_break>0){
        numeroted_rounds=c()
        for (r in 1:(unique(specific_df$nb_rounds)-1)){
          n=c(rep(r,length_round),rep(Inf,length_break))
          numeroted_rounds=append(numeroted_rounds,n)
        }
      }else{
        numeroted_rounds=rep(1:(unique(specific_df$nb_rounds)-1),each=length(omega_intermediate_round))
      }

      replacing_omega1=data.frame(value=rep(omega_intermediate_round,times=unique(specific_df$nb_rounds)-1)
                                  ,round=numeroted_rounds)
      replacing_omega2=data.frame(value=specific_df$value,
                                  round=unique(specific_df$nb_rounds))#last round
      replacing_omega=rbind(replacing_omega1,replacing_omega2)
      extended_specific_df$value[1:dim(replacing_omega)[1]]<-replacing_omega$value
      extended_specific_df$round[1:dim(replacing_omega)[1]]<-replacing_omega$round
      extended_intervention_df=extended_specific_df
    }

  }else{
    extended_intervention_df=intervention_df
  }
  return(extended_intervention_df)
}


