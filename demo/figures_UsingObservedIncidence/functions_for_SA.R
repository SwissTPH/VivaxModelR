#=============================================================================
# 			VIVAX MODELS: FUNCTIONS FOR SENSITIVITY ANALYSIS
#
#  Functions to encapsulate ODE functions for Senstivity analysis
#  and formatting the outputs
#
#  created by : Clara Champagne (clara.champagne@swisstph.ch)
#  originally created : 2020
#
#=============================================================================
##########################################
# ENCAPSULATED FUNCTIONS FOR SA

#lambda
sensi_fit_cm_import_lambda=function(mydf){
  y=c()
  for(i in 1:nrow(mydf)){
    lambda=solve_lambda_vivax(h=my.h.fixed,r=mydf[i,]$r,  gamma=mydf[i,]$gamma,
                              f= mydf[i,]$f, p=my.p.fixed,
                              alpha=mydf[i,]$CM, beta=mydf[i,]$prop_radCure,
                              rho=mydf[i,]$rho, omega=1)

    y=c(y, lambda)
  }

  return(y)
}

# CALCUALTE R0 WITH CM
sensi_fit_cm_import=function(mydf){
  y=c()
  for(i in 1:nrow(mydf)){
    lambda=solve_lambda_vivax(h=my.h.fixed,
                                 r=mydf[i,]$r,  gamma=mydf[i,]$gamma,
                                 f= mydf[i,]$f, p=my.p.fixed,
                                 alpha=mydf[i,]$CM, beta=mydf[i,]$prop_radCure,
                                 rho=mydf[i,]$rho, omega=1)

    y=c(y, get_r0_vivax(lambda,  f=mydf[i,]$f, r=mydf[i,]$r,gamma= mydf[i,]$gamma))
  }

  return(y)
}

# CALCUALTE Rc WITH CM
sensi_fit_cm_import_RC=function(mydf){
  y=c()
  for(i in 1:nrow(mydf)){
    lambda=solve_lambda_vivax(h=my.h.fixed,
                                 r=mydf[i,]$r,  gamma=mydf[i,]$gamma,
                                 f= mydf[i,]$f, p=my.p.fixed,
                                 alpha=mydf[i,]$CM, beta=mydf[i,]$prop_radCure,
                                 rho=mydf[i,]$rho, omega = 1)

    y=c(y, get_rc_vivax(lambda,  f=mydf[i,]$f, r=mydf[i,]$r,gamma= mydf[i,]$gamma,
                           alpha=mydf[i,]$CM, beta=mydf[i,]$prop_radCure))
  }

  return(y)
}


# CALCULATE proportion of relapses with CM
sensi_fit_cm_import_relapses=function(mydf){
  y=c()
  for(i in 1:nrow(mydf)){
    lambda=solve_lambda_vivax(h=my.h.fixed,
                                 r=mydf[i,]$r,  gamma=mydf[i,]$gamma,
                                 f= mydf[i,]$f,  p=my.p.fixed,
                                 alpha=mydf[i,]$CM, beta=mydf[i,]$prop_radCure,
                                 rho=mydf[i,]$rho, omega=1)

    y=c(y, get_prop_relapse(lambda, h=my.h.fixed, f=mydf[i,]$f, r=mydf[i,]$r,gamma= mydf[i,]$gamma,
                                                    alpha=mydf[i,]$CM, beta=mydf[i,]$prop_radCure,
                                                    p=my.p.fixed, omega=1, rho = mydf[i,]$rho))
  }
  return(y)
}

#######################3
# OTHER HELPER FUNCTIONS

reformat_sobol_output=function(sobol.obj, this.name=""){
  my.sobol.obj=sobol.obj
  my.sobol.df=rbind(data.frame("param"=row.names(my.sobol.obj$S), "myvalue"=(my.sobol.obj$S), "index"="First order"),
                    data.frame("param"=row.names(my.sobol.obj$T), "myvalue"=c(my.sobol.obj$T), "index"="Total")) %>%
    mutate(sobol_name=this.name)
  return(my.sobol.df)
}


plot_sobol = function(df){
  my.p=  df %>% mutate(param=factor(param, levels=c("rho", "prop_radCure", "CM", "f","gamma" , "r")),
                 outcome=factor(outcome, levels=c("R0", "Rc", "Prop. of\nrelapses"))) %>%
    #filter(outcome != "lambda") %>%
    ggplot()+
    geom_point(aes(x=myvalue.original, y=param, col=paste0(p*100, "%"), shape=index), size=3)+
    geom_linerange(aes(xmin=myvalue.min..c.i., xmax=myvalue.max..c.i., y=param, col=paste0(p*100, "%")))+
    labs(y='',x='Sobol indices' , col="Prop. of imports", shape="")+
    scale_color_manual(values=c( "cyan3","darkblue"))+
    scale_y_discrete(labels = c("rho"=expression(rho),"gamma"=expression(gamma), "prop_radCure"=expression(beta), "CM"=expression(alpha)))+
    facet_grid( outcome~., switch="both" )+theme_bw()+
    theme(legend.position = "bottom", legend.box="vertical", legend.text = element_text(size = 12 ), legend.title = element_text(size = 15 ))+
    theme(strip.background =element_rect(fill="white", color="white"))+
    theme(strip.placement = "outside", legend.position = "bottom",
          strip.text = element_text(size=12), axis.text.y = element_text(size=10))+
    theme(axis.text.y=element_text(size=15))
  return(my.p)
}



plot_effect = function(df){
  my.p= df %>% rename("alpha"=CM, "beta"=prop_radCure, "rho"=rho) %>%
    mutate( outcome=factor(outcome, levels=c("R0", "Rc", "Prop. of\nrelapses"))) %>%
    pivot_longer(cols=c("alpha", "beta", "rho"))%>% group_by(name) %>%
    dplyr::mutate(my_cut=cut(value, breaks=50) ) %>%
    group_by(my_cut, name, outcome, p) %>%
    summarise(median_y=median(y),
              q975_y=quantile(y, probs = 0.975), q025_y=quantile(y, probs = 0.025),
              median_x=median(value),
              q975_x=quantile(value, probs = 0.975), q025_x=quantile(value, probs = 0.025)) %>%
    ungroup() %>%
    ggplot()+
    geom_line(aes(x=median_x, y=median_y, color=paste0(p*100, "%")))+
    geom_ribbon(aes(ymin=q025_y, ymax=q975_y, x=median_x, fill=paste0(p*100, "%")), alpha=0.3,  color=NA)+
    scale_color_manual(values=c("cyan3","darkblue"))+
    scale_fill_manual(values=c("cyan3","darkblue"))+
    facet_grid( outcome~name, scales="free", switch="both", labeller=labeller(name=label_parsed))+
    labs(x="", fill="Proportion of imports", color="Proportion of imports", y="")+ theme_bw()+
    theme(strip.background =element_rect(fill="white", color="white"))+
    theme(strip.placement = "outside", legend.position = "bottom", strip.text.x = element_text(size=15),
          strip.text.y = element_text(size=12), axis.text.y = element_text(size=10))
  return(my.p)
}

plot_effect_minmax = function(df){
  my.p= df %>% rename("alpha"=CM, "beta"=prop_radCure, "rho"=rho) %>%
    mutate( outcome=factor(outcome, levels=c("R0", "Rc", "Prop. of\nrelapses"))) %>%
    pivot_longer(cols=c("alpha", "beta", "rho"))%>% group_by(name) %>%
    dplyr::mutate(my_cut=cut(value, breaks=50) ) %>%
    group_by(my_cut, name, outcome, p) %>%
    summarise(median_y=median(y),
              q975_y=max(y), q025_y=min(y),
              median_x=median(value),
              q975_x=max(value), q025_x=min(value)) %>%
    ungroup() %>%
    ggplot()+
    geom_line(aes(x=median_x, y=median_y, color=paste0(p*100, "%")))+
    geom_ribbon(aes(ymin=q025_y, ymax=q975_y, x=median_x, fill=paste0(p*100, "%")), alpha=0.3,  color=NA)+
    scale_color_manual(values=c("cyan3","darkblue"))+
    scale_fill_manual(values=c("cyan3","darkblue"))+
    facet_grid( outcome~name, scales="free", switch="both", labeller=labeller(name=label_parsed))+
    labs(x="", fill="Proportion of imports", color="Proportion of imports", y="")+ theme_bw()+
    theme(strip.background =element_rect(fill="white", color="white"))+
    theme(strip.placement = "outside", legend.position = "bottom", strip.text.x = element_text(size=15),
          strip.text.y = element_text(size=12), axis.text.y = element_text(size=10))
  return(my.p)
}
