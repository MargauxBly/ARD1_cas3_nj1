library(shiny)
library(MASS) 
library(ggplot2)
library(SuppDists)
library(gridExtra)
library(png) 
library(doParallel)


al=24 ;  ale =runif(1)*10^3

 ARD1<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio){
  Xts<-NULL 
  nb_maint=floor((n-1)/(tau_periode+1)) #nombre de maintenances cas 1
  n=n-nb_maint
  if(perio){
    t <- seq(0,(n*p-1),by=p)
  }else{
    set.seed(al)
    t <- c(0,sort(sample(1:(n*5),(n-1))))
    set.seed(NULL)
  }
 
  set.seed(ale)
  inc=rnorm(n-1,mu*diff(t),sqrt(sigma2*diff(t)))
  Xt=c(0,cumsum(inc))
  if((tau_periode+1)<=(length(t)-1)){
    taus<-c(0,t[seq(tau_periode+1,length(t)-1, by = tau_periode)],t[length(t)])
    
  Xt0<-Xt[t<=taus[2]]
  ti<-sort(c(t,t[pmatch(taus[-c(1,length(taus))],t)]))
  for ( i in 1:(length(taus)-2)){
    Xt1<-Xt[t<=taus[i+1] & t>=taus[i]]
    Xt2<-Xt[t<=taus[i+2] & t>=taus[i+1]]
    Xtf=Xt2-rho*Xt1[length(Xt1)]
    Xts<-c(Xts,Xtf)
  }
  Xts<-c(Xt0,Xts)
  }else{
    taus<-0
    ti<-t
    Xts<-Xt
  }
  df_Xt<-data.frame(Temps=ti,Degradation=Xts)
  list(dfXt=df_Xt,t=t,taus=taus,Xt=Xt)
}

 ARD1_bis<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio){
   nb_maint=floor((n-1)/(tau_periode+1)) 
   df=ARD1(al,ale,mu,sigma2,rho,n+1,p,tau_periode,perio)
   list(dfXt=df$dfXt[1:n,],taus=df$taus,Xt=df$Xt[1:(n-nb_maint)])
 }
 
 
 

#ale ne change qu'en cliquant sur le bouton resimuler, il faut que alea soit le meme pour les courbes et les estimations

#### cas 1 ####

dfs<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    d=list();ard1=list()
    #nb_maint=floor((n-1)/(tau_periode+1))
    plott<-ARD1_bis(al,ale,mu,sigma2,rho,n,p,tau_periode,perio)
    alea=sort(sample((floor(ale)+1):(ale+nb_simu)))
    d[[1]]<-plott$dfXt
    ard1[[1]]<-plott
    if(nb_simu>1){
    for(i in 2:nb_simu){
        ard1[[i]]<-ARD1_bis(al,alea[i],mu,sigma2,rho,n,p,tau_periode,perio)
        d[[i]]<-ard1[[i]]$dfXt
    }
    }
    list(d=d,alea=alea,plott=plott,ard1=ard1)
}


plotARD1<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){ 
    dfs1<-dfs(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    ints<--rho*mu*(dfs1$plott$taus[-length(dfs1$plott$taus)])
     if(perio){ 
        sxc<-scale_x_continuous(breaks = dfs1$plott$taus[-c(1,length(dfs1$plott$taus))])
        gab<-geom_abline(intercept = ints,slope = rep(mu,length(dfs1$plott$taus)-1),colour="light blue")
        gg<-ggplot(dfs1$d[[1]],aes(x = Temps,y=Degradation))+sxc+gab
          

    }else{
    
        gab<-geom_abline(intercept = ints,slope = rep(mu,length(dfs1$plott$taus)-1),colour="light blue")
      
        gg<-ggplot(dfs1$d[[1]],aes(x = Temps,y=Degradation))+gab
          
    }
    for(i in 1:nb_simu){
        if(i!=num_simu){
            gg=gg+geom_line(data=dfs1$d[[i]],col="light grey")+geom_point(data=dfs1$d[[i]],col="light grey")
        }
    }
    gg=gg+geom_line(data=dfs1$d[[num_simu]])+geom_point(data=dfs1$d[[num_simu]])+theme_bw()+
      theme(panel.grid.minor.y = element_blank(),panel.grid = element_blank())+labs(x="time")
    
    return(gg) 
}



ARD1_para<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    mus<-NULL ; sigs<-NULL
    dfard<-dfs(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    alea<-dfard$alea
    for (i in 1:nb_simu){
        df<-dfard$d[[i]]
        DT<-c(0,diff(df$Temps))
        DD<-c(0,diff(df$Degradation))
        deltay<-data.frame(DTemps=DT,DDegradation=DD)
        deltay<-deltay[-which(deltay$DTemps==0),] #sans les sauts
        mu0<-sum(deltay$DDegradation)/sum(deltay$DTemps)
        #mu0<-dfard$ard1[[i]]$Xt[length(dfard$ard1[[i]]$Xt)]/df$Temps[length(df$Temps)] 
        sigs<-c(sigs,(sum(((deltay$DDegradation-mu0*deltay$DTemps)^2)/deltay$DTemps)/(dim(deltay)[1]))) # -1 pour debiaiser
        mus<-c(mus,mu0)
    }
    list(dfy=deltay,mu1=mus,sig1=sigs)
}


plot_ARD1_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    franz<-data.frame(mus=ARD1_para(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$mu1)
    gf<-ggplot(franz,aes(y=mus))+geom_boxplot(color="dark blue",fill="light grey")+
        labs(title="Estimation de mu",y="mu",x=paste0("mean=",round(mean(franz$mus),2),", sd=",round(sd(franz$mus),2)))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
        geom_hline(yintercept=mu, color = "red") +
        geom_hline(yintercept =franz$mus[num_simu],linetype="dashed")+
      theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    list(gf=gf,franz=franz)
} 

plot_ARD1_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    sabina<-data.frame(mus=ARD1_para(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$sig1)
    gs<-ggplot(sabina,aes(y=mus))+geom_boxplot(color="dark green",fill="light grey")+
        labs(title="Estimation de sigma^2",y="sigma^2",x=paste0("mean=",round(mean(sabina$mus),2),", sd=",round(sd(sabina$mus),2)))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
        geom_hline(yintercept=sigma2, color = "red") +
        geom_hline(yintercept=sabina$mus[num_simu], color = "dark green",linetype="dashed")+
      theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    list(gs=gs,sabina=sabina)
    } 


para_mu_C1<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    ard1_para_mu=ARD1_para(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$mu1
    return(ard1_para_mu)
},vectorize.args = "n")


para_sig_C1<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    ard1_para_sig=ARD1_para(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$sig1
    return(ard1_para_sig)
},vectorize.args = "n")

biais_C1_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    Mmu=para_mu_C1(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    mean_mus<-apply(Mmu,2,mean)
    var_mus<-apply(Mmu,2,var)
    mean_mus<-mean_mus-mu
    Bmu=data.frame(taille=n,biais=mean_mus)
    gg1=ggplot(Bmu,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de mu",
                                                                          x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")
    list(gg1=gg1,bmu=mean_mus,vmu=var_mus)
}

biais_C1_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    Msig=para_sig_C1(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    mean_sigs<-apply(Msig,2,mean)
    var_sigs<-apply(Msig,2,var) 
    mean_sigs<-mean_sigs-sigma2
    Bsig=data.frame(taille=n,biais=mean_sigs)
    gg2=ggplot(Bsig,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de sigma^2",
                                                                           x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")
    list(gg2=gg2,bsig=mean_sigs,vsig=var_sigs)
}

EQM_C1_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    bc1<-biais_C1_mu(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    eqm_mu<-bc1$bmu^2+bc1$vmu
    DEQM=data.frame(taille=n,EQM1=eqm_mu)
    gg3=ggplot(DEQM,aes(x=taille,y=EQM1))+geom_point()+geom_smooth()+labs(title="EQM de mu",
                                                                          x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")
    list(gg3=gg3,deqm=DEQM)
}

EQM_C1_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    bc1<-biais_C1_sig(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    eqm_var<-bc1$bsig^2+bc1$vsig
    DEQM=data.frame(taille=n,EQM2=eqm_var)
    gg4=ggplot(DEQM,aes(x=taille,y=EQM2))+geom_point()+geom_smooth()+labs(title="EQM de sigma^2",
                                                                          x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")
    list(gg4=gg4,deqm=DEQM)
}

#### cas 2 ####

ARD1_C2<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio){
    nb_maint=floor((n-2)/(tau_periode))
    ard1<-ARD1_bis(al,ale,mu,sigma2,rho,n+nb_maint,p,tau_periode,perio)
    Xt_C2<-ard1$Xt
    dinit0<-ard1$dfXt #donnees completes
    dinit<-ard1$dfXt
    taus0<-ard1$taus
    taus1<-taus0[-c(1,length(taus0))]
    df<-dinit[-(pmatch(taus1,dinit$Temps)+1),] #sans les maintenances
    Degradation1<-dinit$Degradation
    dinit$Degradation[(pmatch(taus1,dinit$Temps)+1)]<-NA
    df0<-data.frame(Temps=dinit$Temps,Degradation=dinit$Degradation) #non obs en NA
    df1<-cbind(df0,Degradation1)
    list(df=df,df0=df0,dinit=dinit,dinit0=dinit0,df1=df1,taus1=taus1,taus0=taus0,Xt2=Xt_C2)
}

dfs_C2<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    d=list();g=list();h=list();ardc2=list()
    #alea=ale
    plott<-ARD1_C2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio)
    d[[1]]<-plott$df0 # non obs en NA
    g[[1]]<-plott$dinit0 # donnees completes
    h[[1]]<-plott$df # sans maintenance , sans non obs
    ardc2[[1]]<-plott
    alea=sort(sample((floor(ale)+1):(ale+nb_simu)))
    if(nb_simu>1){
        for(i in 2:nb_simu){
            #alea<-c(alea,runif(1)*10^3)
            ardc2[[i]]<-ARD1_C2(al,alea[i],mu,sigma2,rho,n,p,tau_periode,perio)
            d[[i]]<-ardc2[[i]]$df0
            g[[i]]<-ardc2[[i]]$dinit0
            h[[i]]<-ardc2[[i]]$df
        }
    }
    list(d=d,g=g,h=h,alea=alea,plott=plott,ardc2=ardc2)
}

plotARD1_C2<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    plott<-dfs_C2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    ints<-c(0,-rho*mu*plott$plott$taus1)

        if(perio){
            sxc<-scale_x_continuous(breaks = plott$plott$taus1)
            gab<-geom_abline(intercept = ints,slope = rep(mu,length(ints)),colour="light blue")
         
        gg<-ggplot(plott$d[[1]],aes(x = Temps,y=Degradation))+sxc+gab
          
          
        }else{
            gg<-ggplot(plott$d[[1]],aes(x = Temps,y=Degradation))
        }
    
    for(i in 1:nb_simu){
        if(i!=num_simu){
            gg=gg+geom_line(data=plott$g[[i]],col="light grey")+geom_point(data=plott$g[[i]],col="light grey")
            }
        }
            gg=gg+geom_point(data = plott$d[[num_simu]])+geom_line(data=plott$d[[num_simu]])+geom_line(data=plott$g[[num_simu]],linetype="dotted")+
              theme_bw()+
              theme(panel.grid.minor.x=element_blank(),panel.grid = element_blank())+labs(x="time")
            
    return(gg)
}


ARD1_paraC2<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo){
    mus<-NULL ; sigs<-NULL ;logvrais=NULL
    df0<-dfs_C2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    for (i in 1:nb_simu){
        df<-df0$h[[i]]
        DT<-c(0,diff(df$Temps))
        DD<-c(0,diff(df$Degradation))
        deltay<-data.frame(DTemps=DT,DDegradation=DD)
        deltay<-deltay[-which(deltay$DTemps==0),]
        
        ytim<-df$Degradation[pmatch(df0$ardc2[[i]]$taus1,df$Temps)] #y (tau_i -)
        ytip<-df$Degradation[(pmatch(df0$ardc2[[i]]$taus1,df$Temps)+1)] #y (tau_i +)
        
        deltay1<-deltay[-pmatch(df0$ardc2[[i]]$taus1,df$Temps),] # sans les sauts
        dtj1<-deltay$DTemps[(pmatch(df0$ardc2[[i]]$taus1,df$Temps))] #delta tj,1
        
        func_rhoo<-Vectorize(function(rhoo){
            su=NULL ;suu=NULL ;su2= NULL ; suu2= NULL
              if((n-1)%%tau_periode !=0){
                 res=floor((n-1)/tau_periode)
             }else{
                 res=floor((n-1)/tau_periode)-1
             }
            
            for (j in 1:res){ 
                su<-NULL
                for(i in 1:j){
                    su<-c(su,rhoo^(j-i)*ytim[i])
                }
                suu<-c(suu,ytip[j]-(1-rhoo)*sum(su))
            }
            suuu<-sum(suu)
            mu0<-(sum(deltay1$DDegradation)+suuu)/sum(deltay$DTemps)
            
            su0<-sum(((deltay1$DDegradation-mu0*deltay1$DTemps)^2)/deltay1$DTemps)
            for (j in 1:res){ 
                su2<-NULL
                for(ii in 1:j){
                    su2<-c(su2,rhoo^(j-ii)*ytim[ii])
                }
                suu2<-c(suu2,((ytip[j]-mu0*dtj1[j]-(1-rhoo)*sum(su2))^2)/dtj1[j])
            }
            suuu2<-sum(suu2)
            sigs0<-(su0+suuu2)/(dim(df)[1]-1) #on enleve la valeur initiale + -1 si on debiaise
            
            #logvrais<-sum(log(1/(sqrt(2*pi*sigs0*deltay1$DTemps)))-((deltay1$DDegradation-mu0*deltay1$DTemps)^2/(2*sigs0*deltay1$DTemps))
            #)+sum(log(1/(sqrt(2*pi*sigs0*dtj1)))-(suu2/(2*sigs0)))
            
            list(mu1=mu0,sig1=sigs0,lgv=logvrais)
        })
        if(length(rhoo)==nb_simu){    # sinon juste pour la courbe de la log vrais
            mus=c(mus,as.numeric(func_rhoo(rhoo)["mu1",i]))
            sigs=c(sigs,as.numeric(func_rhoo(rhoo)["sig1",i]))
            #logvrais=c(logvrais,as.numeric(func_rhoo(rhoo)["lgv",i]))
        }else{
            mus=c(mus,as.numeric(func_rhoo(rhoo)["mu1",]))
            sigs=c(sigs,as.numeric(func_rhoo(rhoo)["sig1",]))
            #logvrais=c(logvrais,as.numeric(func_rhoo(rhoo)["lgv",]))
        }
    } 
    
    list(dfy=deltay,dfy1=deltay1,mu1=mus,sig1=sigs)
}


# plot_estim_rho<-function( al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo){
#     lgv_rho0<-ARD1_paraC2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo)$lgv
#     lgv_rho<-lgv_rho0[(length(rhoo)*(nb_simu-1)+1):(length(rhoo)*nb_simu)]
#     dflg<-data.frame(rho=rhoo,lgv=lgv_rho)
#     ggplot(dflg,aes(x=rho,y=lgv))+geom_line()+
#         xlab("Rho")+
#         ylab("Log Vraisemblance")+
#         ggtitle("Log-vraisemblance du modele en fonction de rho")
# }



optim_test<-function(rhoo,al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  lgvr<-ARD1_paraC2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu=1,rhoo)$sig1
  return(lgvr)
}


bxp_estim_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){ 
  rho_hat=NULL
  alea<-dfs_C2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$alea
  for(i in 1:nb_simu){
    rho_hat<-c(rho_hat,optimize(optim_test,c(0,1.1),maximum=F,tol=0.0001,
                                al,alea[i],mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$minimum)
  }
  rho_hat<-data.frame(rhohat=rho_hat)
  gb<-ggplot(rho_hat,aes(y=rhohat))+geom_boxplot(color="dark red",fill="light grey")+
    labs(title="Estimation de rho",y="rho",x=paste0("mean=",round(mean(rho_hat$rhohat),2),", sd=",round(sd(rho_hat$rhohat),2)))+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
    geom_hline(yintercept=rho, color = "red") +
    geom_hline(yintercept=rho_hat$rhohat[num_simu], linetype="dashed")+
    theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
  
  list(gb=gb,rhohat=rho_hat)
}


plot_ARD1C2_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    rhoh<-c(bxp_estim_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$rhohat) 
    ardpar<-ARD1_paraC2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo=rhoh$rhohat)
    lou<-data.frame(mus=ardpar$mu1)
    glou=ggplot(lou,aes(y=mus))+geom_boxplot(color="dark blue",fill="light grey")+
        labs(title="Estimation de mu",y="mu",x=paste0("mean=",round(mean(lou$mus),2),", sd=",round(sd(lou$mus),2)))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
        geom_hline(yintercept=mu, color = "red") +
        geom_hline(yintercept=lou$mus[num_simu], linetype="dashed")+
      theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    list(glou=glou,lou=lou)
    } 


plot_ARD1C2_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    rhoh<-c(bxp_estim_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$rhohat) 
    ardpar<-ARD1_paraC2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo=rhoh$rhohat)
    myr<-data.frame(sigs=ardpar$sig1)
    gmyr=ggplot(myr,aes(y=sigs))+geom_boxplot(color="dark green",fill="light grey")+
        labs(title="Estimation de sigma^2",y="Sigma^2",x=paste0("mean=",round(mean(myr$sigs),2),", sd=",round(sd(myr$sigs),2)))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
        geom_hline(yintercept=sigma2, color = "red") +
        geom_hline(yintercept=myr$sigs[num_simu],linetype="dashed")+
      theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    list(gmyr=gmyr,myr=myr)
} 


para_mu_C2<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    ard1_para_mu=plot_ARD1C2_mu(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$lou$mus
    return(ard1_para_mu)
},vectorize.args = "n")


para_sig_C2<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    ard1_para_sig=plot_ARD1C2_sig(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$myr$sigs
    return(ard1_para_sig)
},vectorize.args = "n",SIMPLIFY = T)

para_rho_C2<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    parar=bxp_estim_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$rhohat$rhohat
},vectorize.args = "n")


biais_C2_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    Mmu=para_mu_C2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    mean_mus<-apply(Mmu,2,mean)
    var_mus<-apply(Mmu,2,var)
    mean_mus<-mean_mus-mu
    Bmu=data.frame(taille=n,biais=mean_mus)
    gg1=ggplot(Bmu,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de mu",
                                                                          x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")
    list(gg1=gg1,bmu=mean_mus,vmu=var_mus)
}

biais_C2_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    Msig=para_sig_C2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    mean_sigs<-apply(Msig,2,mean)
    var_sigs<-apply(Msig,2,var)
    mean_sigs<-mean_sigs-sigma2
    Bsig=data.frame(taille=n,biais=mean_sigs)
    gg2=ggplot(Bsig,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de sigma^2",
                                                                           x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")
    list(gg2=gg2,bsig=mean_sigs,vsig=var_sigs)
}

biais_C2_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    dfrho=para_rho_C2(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    Mrho=apply(dfrho,2,mean)
    Vrho=apply(dfrho,2,var)
    brho=Mrho-rho
    Brho=data.frame(taille=n,biais=brho)
    ggr=ggplot(Brho,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de rho",
                          x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")
    list(ggr=ggr,brho=brho,vrho=Vrho)
}

EQM_C2_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    bc1<-biais_C2_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    eqm_rho<-bc1$brho^2+bc1$vrho
    DEQM=data.frame(taille=n,EQM1=eqm_rho)
    gg3=ggplot(DEQM,aes(x=taille,y=EQM1))+geom_point()+geom_smooth()+labs(title="EQM de rho",
                                                                          x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")
    list(gg3=gg3,deqm=DEQM)
}


EQM_C2_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    bc1<-biais_C2_mu(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    eqm_mu<-bc1$bmu^2+bc1$vmu
    DEQM=data.frame(taille=n,EQM1=eqm_mu)
    gg3=ggplot(DEQM,aes(x=taille,y=EQM1))+geom_point()+geom_smooth()+labs(title="EQM de mu",
                                                                          x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")

    list(gg3=gg3,deqm=DEQM)
}

EQM_C2_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    bc1<-biais_C2_sig(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    eqm_var<-bc1$bsig^2+bc1$vsig
    DEQM=data.frame(taille=n,EQM2=eqm_var)
    gg4=ggplot(DEQM,aes(x=taille,y=EQM2))+geom_point()+geom_smooth()+labs(title="EQM de sigma^2",
    x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")
    
    list(gg4=gg4,deqm=DEQM)
}


#### cas 3 #### 

ARD1_C3<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio){
  nb_maint=floor((n-1)/tau_periode)
    ard1<-ARD1_bis(al,ale,mu,sigma2,rho,n+nb_maint,p,tau_periode,perio)
    Xt<-ard1$Xt[1:n]
    dinit0<-ard1$dfXt
    dinit<-ard1$dfXt
    taus0<-ard1$taus
    # if(length(taus0)==(nb_maint+2)){
    #   taus1<-taus0[-c(1,length(taus0))]
    # }else{
    #   taus1<-taus0[-1]
    # }
    taus1<-taus0[-c(1,length(taus0))]
    df<-dinit[-(pmatch(taus1,dinit$Temps)),]
    Degradation1<-dinit$Degradation
    dinit$Degradation[(pmatch(taus1,dinit$Temps))]<-NA
    df0<-data.frame(Temps=dinit$Temps,Degradation=dinit$Degradation)
    df1<-cbind(df0,Degradation1)
    list(df1=df1,df=df,taus1=taus1,dinit0=dinit0,df0=df0,ard1=ard1,Xt=Xt)
}

dfs_C3<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    d=list();g=list();h=list();ardc3=list()  ; xt=list()
    plott<-ARD1_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio)
    xt[[1]]<-plott$Xt
    d[[1]]<-plott$df0 # non obs en NA
    g[[1]]<-plott$dinit0 # donnees completes
    h[[1]]<-plott$df # sans les donnees non obs
    ardc3[[1]]<-plott
    alea=sort(sample((floor(ale)+1):(ale+nb_simu)))
    if(nb_simu>1){
        for(i in 2:nb_simu){
            ardc3[[i]]<-ARD1_C3(al,alea[i],mu,sigma2,rho,n,p,tau_periode,perio)
            d[[i]]<-ardc3[[i]]$df0
            g[[i]]<-ardc3[[i]]$dinit0
            h[[i]]<-ardc3[[i]]$df
            xt[[i]]<-ardc3[[i]]$Xt
        }
    }
    list(d=d,g=g,h=h,alea=alea,plott=plott,alea=alea,ardc3=ardc3,xt=xt)
}

plotARD1_C3<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    plotti<-dfs_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    ints<-c(0,-rho*mu*plotti$plott$taus1)
        if(perio){
            sxc<-scale_x_continuous(breaks = plotti$plott$taus1)
            gab<-geom_abline(intercept = ints,slope = rep(mu,length(ints)),colour="light blue")
          
        gg<-ggplot(plotti$d[[1]],aes(x = Temps,y=Degradation))+sxc+gab
         
        }else{
            gg<-ggplot(plotti$d[[1]],aes(x = Temps,y=Degradation))
        }
    
    for(i in 1:nb_simu){
        if(i!=num_simu){
            gg=gg+geom_line(data=plotti$g[[i]],col="light grey")+geom_point(data=plotti$g[[i]],col="light grey")
            }
        }
            gg=gg+geom_point(data = plotti$d[[num_simu]])+geom_line(data=plotti$d[[num_simu]])+geom_line(data=plotti$g[[num_simu]],linetype="dotted")+
              theme_bw()+
              theme(panel.grid.minor.x=element_blank(),panel.grid = element_blank())+labs(x="time")
            
    return(gg)
} 



ARD1_paraC3<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo){
    mus<-NULL ; sigs<-NULL ; logvrais=NULL
    nb_maint=floor((n-1)/tau_periode)
    df0<-dfs_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    for (i in 1:nb_simu){
        dfxt<-df0$xt[[i]]
        df<-df0$h[[i]] #sans les donnees non obs
        dfg<-df0$g[[i]] #donnees completes
        DT<-diff(df$Temps)
        DD<-diff(df$Degradation)
        deltay<-data.frame(DTemps=DT,DDegradation=DD)
        
        jumps<-seq(tau_periode,nb_maint*(tau_periode),by=tau_periode) #indices des sauts en increment
        deltay1<-deltay[-jumps,] # sans les sauts
        
        deltay2<-deltay[jumps,] #juste les sauts

        if(tau_periode!=1){ 
        deltay3<-deltay1[1:(jumps[length(jumps)]-length(jumps)),] #increments precedents les sauts
        
        jumps0<-c(0,jumps)+1 ; vec_d3<-NULL
        
        for (k in 1:length(jumps)){
          vec_d3<-c(vec_d3,sum(deltay$DDegradation[jumps0[k]:(jumps[k]-1)]))
        }
        su=vec_d3
        
        }else{
          deltay3<-deltay1 ; su<-0
        }
        
        func_rhoo<-Vectorize(function(rhoo){
        #mu0<-(sum(deltay1$DDegradation)
             # +(1/(1-rhoo))*(sum(deltay2$DDegradation)
                             #+rhoo*sum(deltay3$DDegradation)))/(sum(deltay$DTemps))
        
        mu0<-dfxt[length(dfxt)]/df$Temps[length(df$Temps)] 
        
        sig0<-(sum(((deltay1$DDegradation-mu0*deltay1$DTemps)^2)/deltay1$DTemps)+
                   sum((((deltay2$DDegradation-mu0*(1-rhoo)*deltay2$DTemps)+
                             rhoo*su)^2)/((1-rhoo)^2*deltay2$DTemps)))/(dim(df)[1]-1) # moins la valeur initiale moins 1 pour debiaiser
        
        lgv0<-sum(log(1/(sqrt(2*pi*sig0*deltay1$DTemps)))-((deltay1$DDegradation-mu0*deltay1$DTemps)^2/(2*sig0*deltay1$DTemps))
        )+sum(log(1/(sqrt(2*pi*sig0*deltay2$DTemps*(1-rhoo)^2)))-(((deltay2$DDegradation-mu0*(1-rhoo)*deltay2$DTemps)+
                                                                   rhoo*su)^2)/((1-rhoo)^2*deltay2$DTemps*2*sig0))
        
        #lgv0<-sum(log(1/(sqrt(2*pi*sig0*deltay1$DTemps))))+sum(log(1/(sqrt(2*pi*sig0*deltay2$DTemps*(1-rhoo)^2))))-((length(jumps)-1)/2) # en debiaisant pareil que lgv0
        
        #logl3_t1<--log(sig0*(1-rhoo)^2)*length(jumps)/2-sum(log(2*pi*deltay2$DTemps))/2-(length(jumps)/2) #lorsque tau_periode=1
        #lgv0<--(dim(df)[1]-1)*(log(sig0)/2)-length(jumps)*log(1-rhoo)-log(2*pi)*(dim(df)[1]-1)-sum(log(deltay1$DTemps))-sum(log(deltay2$DTemps))-((dim(df)[1]-1)/2)

            list(mu1=mu0,sig1=sig0,lgv=lgv0)
        })
       
        
        if(length(rhoo)==nb_simu){    # sinon juste pour la courbe de la log vrais
            mus=c(mus,as.numeric(func_rhoo(rhoo)["mu1",i]))
            sigs=c(sigs,as.numeric(func_rhoo(rhoo)["sig1",i]))
            logvrais=c(logvrais,as.numeric(func_rhoo(rhoo)["lgv",i]))
        }else{
            mus=c(mus,as.numeric(func_rhoo(rhoo)["mu1",]))
            sigs=c(sigs,as.numeric(func_rhoo(rhoo)["sig1",]))
            logvrais=c(logvrais,as.numeric(func_rhoo(rhoo)["lgv",]))
        }
    }
    
    list(dfy=deltay,mu1=mus,sig1=sigs,lgv1=logvrais,df=df,nbmaint=nb_maint)
}


# plot_estim_rho_C3<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo){
#     lgv_rho0<-ARD1_paraC3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo)$lgv
#     lgv_rho<-lgv_rho0[(length(rhoo)*(nb_simu-1)+1):(length(rhoo)*nb_simu)]
#     dflg<-data.frame(rho=rhoo,lgv=lgv_rho)
#     ggplot(dflg,aes(x=rho,y=lgv))+geom_line()+
#         xlab("Rho")+
#         ylab("Log Vraisemblance")+
#         ggtitle("Log-vraisemblance du modele en fonction de rho")
# }

optim_test_C3<-function(rhoo,al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    paraC3<-ARD1_paraC3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu=1,rhoo)
    #sig_r<-paraC3$sig1
    #fsig<-(dim(paraC3$df)[1]-1)*(log(sig_r)/2)+(paraC3$nbmaint*log(1-rhoo))
    #argm<-sig_r*(1-rhoo)^2
    return(paraC3$lgv1)
}


bxp_estim_rho_C3<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    rho_hat=NULL
    alea<-dfs_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$alea
    for(i in 1:nb_simu){
        rho_hat<-c(rho_hat,optimize(optim_test_C3,c(0,1),maximum=T,tol=0.0001,
                                    al,alea[i],mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$maximum)
    }
    rho_hat<-data.frame(rhohat=rho_hat)
    dfg=ggplot(rho_hat,aes(y=rhohat))+geom_boxplot(color="dark red",fill="light grey",outlier.shape = NA)+
        labs(title="Estimation de rho < 1 ",y="rho",x=paste0("mean=",round(mean(rho_hat$rhohat),2),", sd=",round(sd(rho_hat$rhohat),2)))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
        geom_hline(yintercept=rho, color = "red") +
        geom_hline(yintercept=rho_hat$rhohat[num_simu], linetype="dashed")+
      theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    list(dfg=dfg,rhohat=rho_hat)
}

plot_ARD1C3_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    rhoh<-bxp_estim_rho_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$rhohat 
    ardpar<-ARD1_paraC3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo=rhoh$rhohat)
    paul<-data.frame(mus=ardpar$mu1)
    dgp<-ggplot(paul,aes(y=mus))+geom_boxplot(color="dark blue",fill="light grey")+
        labs(title="Estimation de mu",y="mu",x=paste0("mean=",round(mean(paul$mus),2),", sd=",round(sd(paul$mus),2)))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
        geom_hline(yintercept=mu, color = "red") +
        geom_hline(yintercept=paul$mus[num_simu], linetype="dashed", color = "dark green")+
      theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    list(dgp=dgp,paul=paul)
} 



plot_ARD1C3_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    rhoh<-c(bxp_estim_rho_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$rhohat) 
    ardpar<-ARD1_paraC3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo=rhoh$rhohat)
    myr<-data.frame(sigs=ardpar$sig1)
    ggm<-ggplot(myr,aes(y=sigs))+geom_boxplot(color="dark green",fill="light grey")+
        labs(title="Estimation de sigma^2",y="Sigma^2",x=paste0("mean=",round(mean(myr$sigs),2),", sd=",round(sd(myr$sigs),2)))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
        geom_hline(yintercept=sigma2, color = "red") +
        geom_hline(yintercept=myr$sigs[num_simu], linetype="dashed", color = "dark green")+
      theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    list(ggm=ggm,myr=myr)
} 

para_mu_C3<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    ard1_para_mu=plot_ARD1C3_mu(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$paul$mus
    return(ard1_para_mu)
},vectorize.args = "n")


para_sig_C3<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    ard1_para_sig=plot_ARD1C3_sig(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$myr$sigs
    return(ard1_para_sig)
},vectorize.args = "n",SIMPLIFY = T)

para_rho_C3<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    parar=bxp_estim_rho_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)$rhohat$rhohat
},vectorize.args = "n")

biais_C3_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    Mmu=para_mu_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    mean_mus<-apply(Mmu,2,mean)
    var_mus<-apply(Mmu,2,var)
    mean_mus<-mean_mus-mu
    Bmu=data.frame(taille=n,biais=mean_mus)
    gg1=ggplot(Bmu,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de mu",
                                                                          x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")

    list(gg1=gg1,bmu=mean_mus,vmu=var_mus)
}

biais_C3_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    Msig=para_sig_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    mean_sigs<-apply(Msig,2,mean)
    var_sigs<-apply(Msig,2,var)
    mean_sigs<-mean_sigs-sigma2
    Bsig=data.frame(taille=n,biais=mean_sigs)
    gg2=ggplot(Bsig,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de sigma^2",
                                                                           x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")
    list(gg2=gg2,bsig=mean_sigs,vsig=var_sigs)
}

biais_C3_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    dfrho=para_rho_C3(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    Mrho=apply(dfrho,2,mean)
    Vrho=apply(dfrho,2,var)
    brho=Mrho-rho
    Brho=data.frame(taille=n,biais=brho)
    ggr=ggplot(Brho,aes(x=taille,y=biais))+geom_point()+geom_smooth()+labs(title="Biais de rho",
             x ="taille de l'echantillon", y = "Biais")+geom_hline(yintercept = 0,col="red")
    list(ggr=ggr,brho=brho,vrho=Vrho)
}

EQM_C3_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    bc1<-biais_C3_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    eqm_rho<-bc1$brho^2+bc1$vrho
    DEQM=data.frame(taille=n,EQM1=eqm_rho)
    gg3=ggplot(DEQM,aes(x=taille,y=EQM1))+geom_point()+geom_smooth()+labs(title="EQM de rho",
               x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")
    list(gg3=gg3,deqm=DEQM)
}

EQM_C3_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    bc1<-biais_C3_mu(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    eqm_mu<-bc1$bmu^2+bc1$vmu
    DEQM=data.frame(taille=n,EQM1=eqm_mu)
    gg3=ggplot(DEQM,aes(x=taille,y=EQM1))+geom_point()+geom_smooth()+labs(title="EQM de mu",
                                                                          x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")
    list(gg3=gg3,deqm=DEQM)
}

EQM_C3_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    bc1<-biais_C3_sig(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu)
    eqm_var<-bc1$bsig^2+bc1$vsig
    DEQM=data.frame(taille=n,EQM2=eqm_var)
    gg4=ggplot(DEQM,aes(x=taille,y=EQM2))+geom_point()+geom_smooth()+labs(title="EQM de sigma^2",
                                                                          x ="taille de l'echantillon", y = "EQM")+geom_hline(yintercept = 0,col="red")
    list(gg4=gg4,deqm=DEQM)
}
 

#### CAS 4 ####

ARD1_C4<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio){
    nb_maint=floor((n-2)/(tau_periode-1))
    ard1<-ARD1_bis(al,ale,mu,sigma2,rho,n+2*nb_maint,p,tau_periode,perio) # n+2*nb_maint
    dinit0<-ard1$dfXt
    dinit<-ard1$dfXt
    taus0<-ard1$taus
    taus1<-taus0[-c(1,length(taus0))]
    df<-dinit[-(pmatch(taus1,dinit$Temps)),]
    df<-df[-(pmatch(taus1,df$Temps)),]
    Degradation1<-dinit$Degradation
    dinit$Degradation[c(pmatch(taus1,dinit$Temps),(pmatch(taus1,dinit$Temps)+1))]<-NA
    df0<-data.frame(Temps=dinit$Temps,Degradation=dinit$Degradation)
    df1<-cbind(df0,Degradation1)
    list(df1=df1,df=df,taus1=taus1,dinit0=dinit0,df0=df0)
}

dfs_C4<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
    d=list();g=list();h=list();ardc4=list()
    plott<-ARD1_C4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio)
    d[[1]]<-plott$df0 # non obs en NA
    g[[1]]<-plott$dinit0 # donnees completes
    h[[1]]<-plott$df # sans les donnees non obs
    ardc4[[1]]<-plott
    alea=sort(sample((floor(ale)+1):(ale+nb_simu)))
    if(nb_simu>1){
        for(j in 2:nb_simu){
            ardc4[[j]]<-ARD1_C4(al,alea[j],mu,sigma2,rho,n,p,tau_periode,perio)
            d[[j]]<-ardc4[[j]]$df0
            g[[j]]<-ardc4[[j]]$dinit0
            h[[j]]<-ardc4[[j]]$df
        }
      
    }
    list(d=d,g=g,h=h,alea=alea,plott=plott,ardc4=ardc4)

}


plotARD1_C4<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
    plotti<-dfs_C4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
    ints<-c(0,-rho*mu*plotti$plott$taus1)
    if(perio){
        sxc<-scale_x_continuous(breaks = plotti$plott$taus1)
        gab<-geom_abline(intercept = ints,slope = rep(mu,length(ints)),colour="light blue")
      
        gg<-ggplot(plotti$d[[1]],aes(x = Temps,y=Degradation))+sxc+gab
          
            }else{
        gg<-ggplot(plotti$d[[1]],aes(x = Temps,y=Degradation))
    }
    
    for(i in 1:nb_simu){
        if(i!=num_simu){
            gg=gg+geom_line(data=plotti$g[[i]],col="light grey")+geom_point(data=plotti$g[[i]],col="light grey")
        }
    }
    gg=gg+geom_point(data = plotti$d[[num_simu]])+geom_line(data=plotti$d[[num_simu]])+geom_line(data=plotti$g[[num_simu]],linetype="dotted")+
      theme_bw()+
      theme(panel.grid.minor.x=element_blank(),panel.grid = element_blank())+labs(x="time")
    
    return(gg)
} 


ARD1_paraC4<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo){
  nb_maint=floor((n-2)/(tau_periode-1))
  Ardc4<-dfs_C4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
  mus<- NULL ; sigs<-NULL ; logvrais<-NULL
  for ( i in 1:nb_simu){
  dfxt<-Ardc4$g[[i]] #donnees completes
  dfxt2<-Ardc4$h[[i]] #sans les sauts verticaux
  
  deg4<-dfxt2$Degradation
  taus<-Ardc4$ardc4[[i]]$taus1 ; taus0=c(0,taus,dfxt$Temps[length(dfxt$Temps)]) 
  jumps<-seq(tau_periode,((tau_periode-1)*nb_maint)+1,by=tau_periode-1) #indices sauts des increments observes
  jumps0<-c(0,jumps) ; jumps_c<-seq(tau_periode+1,nb_maint*(tau_periode+1),by=tau_periode+1) #indices sauts donnees completes
  #deltay0<-diff(dfxt2$Degradation)[1:(jumps[length(jumps)]-1)][-jumps] #delta y (j-1) jusqu a la derniere maintenance
  dtj1<-diff(dfxt$Temps)[jumps_c+1] ; dtjf <-diff(dfxt$Temps)[jumps_c-1] ; dtj1_bis<-c(0,dtj1[-length(dtj1)])
  dy_noj<-diff(deg4)[-jumps] #sans les increments des sauts
  sum_dynoj<-sum(dy_noj)
  dt_noj<-diff(dfxt2$Temps)[-jumps] # dt sans les increments des sauts
  sum_dtnoj<-sum(dt_noj)
  dys<-diff(deg4)[jumps] # que les incrememnts des sauts
  
  func_rhoo<-Vectorize(function(rhoo){
  vec_d4<-NULL
  
  if(sum(diff(jumps)!=1)==0){
    vec_d4<-c(diff(dfxt2$Degradation)[1],rep(0,length(jumps)-1))
  }else{
    for (k in 1:length(jumps)){
      vec_d4<-c(vec_d4,sum(diff(dfxt2$Degradation)[(jumps0[k]+1):(jumps[k]-1)]))
    }
  }
  v_vec<-rhoo*vec_d4
  u_vec<-dtj1+(1-rhoo)*dtjf-rhoo*dtj1_bis
  s_vec<-dtj1+(1-rhoo)^2*dtjf+(rhoo^2)*dtj1_bis
  SIG<-diag(s_vec)
  SIG[col(SIG)-row(SIG)==1]<--rhoo*dtj1[-length(dtj1)] ; SIG[row(SIG)-col(SIG)==1]<--rhoo*dtj1[-length(dtj1)]
    
    deg4<-dfxt2$Degradation
     taus0<-c(0,taus,n-1) 
    
    muh<-as.numeric((sum_dynoj+(t(u_vec)%*%solve(SIG)%*%dys)+(t(u_vec)%*%solve(SIG)%*%v_vec))/((t(u_vec)%*%solve(SIG)%*%u_vec)+sum_dtnoj))
    sigh<-as.numeric((t(dys-muh*u_vec+v_vec)%*%solve(SIG)%*%(dys-muh*u_vec+v_vec)+sum((dy_noj-muh*dt_noj)^2/dt_noj))/(dim(dfxt2)[1]-1))    # - la 1ere valeur et - 1 si on debiaise                                                                      
    rems<-as.numeric(t(dys-muh*u_vec+v_vec)%*%solve(SIG)%*%(dys-muh*u_vec+v_vec))
    lgv<--(1/2)*(log((2*pi)^length(jumps)*sigh^length(jumps)*det(SIG))+(rems/sigh)+
                   log(2*pi*sigh)*(dim(dfxt2)[1]-length(jumps))+sum(log(dt_noj))+(sum((dy_noj-muh*dt_noj)^2/(dt_noj*sigh))))
    
    
    list(muh=muh,sigh=sigh,lgv=lgv)
  })
  
  mus=c(mus,as.numeric(func_rhoo(rhoo)["muh",i]))
  sigs=c(sigs,as.numeric(func_rhoo(rhoo)["sigh",i]))
  logvrais=c(logvrais,as.numeric(func_rhoo(rhoo)["lgv",i]))
   
 
   }
  
  list(muh=mus,sigh=sigs,lgv=logvrais,ardC4=Ardc4)
}


  #estimation de rho


# logvrais_C42<-function(rhoo,al,ale,mu,sigma2,rho,n,p,tau_periode,perio){
#   lgv<-ARD1_paraC4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu=1,rhoo)$lgv
#   return(lgv)
# }


logvrais_C4<-function(rhoo,al,ale,mu,sigma2,rho,n,p,tau_periode,perio){
  paraC4<-ARD1_paraC4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu=1,rhoo)
  #sig4<-paraC4$sigh
  #argmin<-log(sig4)*((dim(paraC4$ardC4$h[[1]])[1]-1)/2)+log(solve(paraC4$SIG))
  return(paraC4$lgv)
   
}


estim_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  alea=dfs_C4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$alea
  
  Ncpus <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(cl)
  
  iter <- itertools::isplitIndices(n=nb_simu, chunks = Ncpus)
  res <- foreach(i=iter, .combine='c') %dopar% {
    source("server.R")
    rrho<-NULL
    for (j in i){
      rrho<-c(rrho,optimize(logvrais_C4,c(0,1),maximum=T,tol=0.0001,
                            al,ale=alea[j],mu,sigma2,rho,n,p,tau_periode,perio)$maximum
      )
    }
    return(rrho)
  }
  parallel::stopCluster(cl)
  return(res)
}
  


plot_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  lgvr<-logvrais_C4(rhoo=seq(0,1,0.01),al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
  martin<-data.frame(rho=seq(0,1,0.01),lgv=lgvr)
  ggplot(martin,aes(x=rho,y=lgv))+geom_line()+geom_vline(xintercept = rho,color="red")+
    xlab("Rho")+ylab("Log Vraisemblance")
}

bxp_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu){
  rho_hat<-estim_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
  dfrho<-data.frame(rho=rho_hat)
  ggrho<-ggplot(dfrho,aes(y=rho))+geom_boxplot(color="dark red",fill="light grey")+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
    labs(title = "Estimation de rho",x=paste0("mean=",round(mean(dfrho$rho),2),", sd=",round(sd(dfrho$rho),2)))+
    geom_hline(yintercept=rho, color = "red") +
    geom_hline(yintercept=dfrho$rho[num_simu], linetype="dashed") +
    theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
  list(ggrho=ggrho,rhos=rho_hat)
}


bxp_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu,rhoo){
  estp<-ARD1_paraC4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo)
  mus<-estp$muh ; sigs<-estp$sigh
  mu1<-mus[num_simu] 
  michel<-data.frame(mus=mus,sigs=sigs)
  ggm<-ggplot(michel,aes(y=mus))+geom_boxplot(color="dark blue",fill="light grey")+
    labs(title="Estimation de mu",y="mu",x=paste0("mean=",round(mean(michel$mus),2),", sd=",round(sd(michel$mus),2)))+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
    geom_hline(yintercept=mu, color = "red") +
    geom_hline(yintercept=michel$mus[num_simu], linetype="dashed")+
    theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
  list(ggm=ggm,mu1=mu1,michel=michel,mus=mus)
}


bxp_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu,rhoo){
  michel<-bxp_mu(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,num_simu,rhoo)$michel
  sig1<-michel$sigs[num_simu]
  ggs<-ggplot(michel,aes(y=sigs))+
    geom_boxplot(color="dark green",fill="light grey",outlier.shape = NA)+
    labs(title="Estimation de sigma2",y="sigma2",x=paste0("mean=",round(mean(michel$sigs),2),", sd=",round(sd(michel$sigs),2)))+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
    geom_hline(yintercept=sigma2, color = "red") +
    geom_hline(yintercept=michel$sigs[num_simu], linetype="dashed")+
    theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())
    
  list(ggs=ggs,sig1=sig1,sigs=michel$sigs)
}

mus_n<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  mus_sig<-Vectorize(function(n,mu){
  rhoo=estim_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
  estp<-ARD1_paraC4(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu,rhoo=rhoo)
  mus<-estp$muh ; sign<-estp$sigh
  if(mu){
    return(mus)
  }else{
    return(sign)
  }
  })
  mus<-mus_sig(n,mu=T) ; sign<-mus_sig(n,mu=F)
list(mus=mus,sign=sign)
}


biais_C4_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  Mn<-mus_n(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$mus
  colnames(Mn)<-paste0("n=",n) ; rownames(Mn)<-paste0("rho",1:nb_simu)
  biais<-apply(Mn,2,mean)-mu
  varest<-apply(Mn,2,var)
  eqm<-biais^2+varest
  remi<-data.frame(n=n,biais=biais,eqm=eqm)
  gg<-ggplot(remi,aes(x=n,y=biais))+geom_smooth()+labs(title = "Biais de mu")+
    geom_hline(yintercept = 0,color="red")
  list(gg=gg,remi=remi)

}

EQM_C4_mu<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  bmu<-biais_C4_mu(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$remi
  ggplot(bmu,aes(x=n,y=eqm))+geom_smooth()+labs("EQM de mu")+geom_hline(yintercept = 0,color="red")
  }


biais_C4_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  Mn<-mus_n(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$sign
  colnames(Mn)<-paste0("n=",n) ; rownames(Mn)<-paste0("rho=",1:nb_simu)
  biais<-apply(Mn,2,mean)-sigma2
  varest<-apply(Mn,2,var)
  eqm<-biais^2+varest
  remi<-data.frame(n=n,biais=biais,eqm=eqm)
  gg<-ggplot(remi,aes(x=n,y=biais))+geom_smooth()+labs(title = "Biais de sigma2")+
    geom_hline(yintercept = 0,color="red")
  list(gg=gg,remi=remi)
  
}

EQM_C4_sig<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  bmu<-biais_C4_sig(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$remi
  ggplot(bmu,aes(x=n,y=eqm))+geom_smooth()+labs("EQM de sigma2")+geom_hline(yintercept = 0,color="red")
}

rhos_n<-Vectorize(function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  rhos<-estim_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
  return(rhos)
},vectorize.args = "n")


biais_C4_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  Mn<-rhos_n(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)
  biais<-apply(Mn,2,mean)-rho
  varest<-apply(Mn,2,var)
  eqm<-biais^2+varest
  remi<-data.frame(n=n,biais=biais,eqm=eqm)
  gg<-ggplot(remi,aes(x=n,y=biais))+geom_smooth()+labs(title = "Biais de Rho")+
    geom_hline(yintercept = 0,color="red")
   list(gg=gg,remi=remi)
}



EQM_C4_rho<-function(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu){
  bmu<-biais_C4_rho(al,ale,mu,sigma2,rho,n,p,tau_periode,perio,nb_simu)$remi
  ggplot(bmu,aes(x=n,y=eqm))+geom_smooth()+labs("EQM de rho")+geom_hline(yintercept = 0,color="red")
}




shinyServer(function(input, output,session) {
   
        output$Text1 <-renderText({if(input$but=="CAS 1"){
            "
            CAS 1 : Les degradations avant et apres chaque maintenance sont observees 
            \n"
        }else if (input$but=="CAS 2"){
            "CAS 2 : Les degradations apres chaque maintenance ne sont pas observees \n" 
        }else if (input$but=="CAS 3"){
            "CAS 3 : Les degradations avant chaque maintenance ne sont pas observees \n"
        }else{
            "CAS 4 : Les degradations avant et apres la maintenance ne sont pas observees \n"
        }
            
            })
          

        output$Text3 <- renderText({if(input$but=="CAS 1"){
            "
            CAS 1 : Les degradations avant et apres chaque maintenance sont observees 
            \n"
        }else if (input$but2=="CAS 2"){
            "CAS 2 : Les degradations apres chaque maintenance ne sont pas observees \n" 
        }else if (input$but2=="CAS 3"){
            "CAS 3 : Les degradations avant chaque maintenance ne sont pas observees \n"
        }else{
            "CAS 4 : Les degradations avant et apres chaque maintenance ne sont pas observees \n"
        }
            
        })
    
    
    #values <- reactiveValues(df1 = NULL, df2 = NULL, Temps = NULL, Boolean = FALSE)
    val<-reactiveValues(al=24,ale=100)

    observeEvent(input$button, {
        b<-val$ale+input$nbtest+1
        val$ale<-b
        # mu_x <- input$mu_x
        # var1 <- input$var1
        # p <- input$p
        # n <- input$n
        # h <- input$h
        # taille_echatillon = input$taille_echatillon
        # nombre_des_echantillons = input$nombre_des_echantillons
        # updateNumericInput(session, "mu_x", value=0)
        # updateNumericInput(session, "var1", value=0)
        # updateNumericInput(session, "p", value=0)
        # updateNumericInput(session, "n", value=0)
        # updateNumericInput(session, "h", value=0)
        # updateNumericInput(session, "taille_echatillon", value=0)
        # updateNumericInput(session, "nombre_des_echantillons", value=0)
        # 
        # updateNumericInput(session, "mu_x", value=mu_x)
        # updateNumericInput(session, "var1", value=var1)
        # updateNumericInput(session, "p", value=p)
        # updateNumericInput(session, "n", value=n)
        # updateNumericInput(session, "h", value=h)
        # updateNumericInput(session, "taille_echatillon", value=taille_echatillon)
        # updateNumericInput(session, "nombre_des_echantillons", value=nombre_des_echantillons)
        
    })
      
    output$mytable <- renderPlotly({
        validate()
        
        if(input$but == "CAS 1"){
            plotARD1(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)
          }
        else if(input$but == "CAS 2"){
            plotARD1_C2(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)
        }
        else if (input$but=="CAS 3"){
            plotARD1_C3(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)
          
           }else if (input$but=="CAS 4"){
            plotARD1_C4(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)
                  }
    })
    

    
    # output$tabes<-renderTable({
    #     if(input$but=="CAS 1"){
    #         data.frame(mu=
    #    ARD1_para(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest)$mu1[input$num]
    # ,sigma2=ARD1_para(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest)$sig1[input$num]
    #         ) 
    #     }else if (input$but=="CAS 2"){
    #         data.frame(mu=plot_ARD1C2_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)$lou[input$num,1],
    #                    sigma2=plot_ARD1C2_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)$myr[input$num,1],
    #                    rho=bxp_estim_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,
    #                                            input$pm_periode,perio=T,input$nbtest,input$num)$rhohat[input$num,1]
    #                          
    #         )
    #     }else if (input$but=="CAS 3"){ 
    #         data.frame(mu=plot_ARD1C3_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)$paul[input$num,1],
    #         sigma2=plot_ARD1C3_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)$myr[input$num,1],
    #         rho=bxp_estim_rho_C3(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)$rhohat[input$num,1]
    #         )                          
    #     }
      # else if (input$but=="CAS 4"){
      #     rhohat<-estim_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,
      #                       input$pm_periode,perio=T,input$nbtest,num_simu=1:input$nbtest)
      #     data.frame(mu=estim_paras(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,
      #                               input$pm_periode,perio=T,input$num,rhoo=rhohat[input$num])$muh
      #                ,sigma2=estim_paras(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,
      #                                   input$pm_periode,perio=T,input$num,rhoo=rhohat[input$num])$sigh
      #         
      #                ,rho=rhohat[input$num]
      #   ) }
       #  })
    
    # output$tabes2<-renderTable({
    #     if(input$but2=="CAS 1"){
    #         data.frame(mu=
    #                        ARD1_para(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2)$mu1[input$num]
    #                    ,sigma2=ARD1_para(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2)$sig1[input$num]
    #         ) 
    #     }else if (input$but2=="CAS 2"){
    #         data.frame(mu=plot_ARD1C2_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)$lou[input$num2,1],
    #                    sigma2=plot_ARD1C2_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)$myr[input$num2,1],
    #                    rho=bxp_estim_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)$rhohat[input$num2,1]
    #                    
    #         )
    #     }else if (input$but2=="CAS 3"){ 
    #         data.frame(mu=plot_ARD1C3_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)$paul[input$num2,1],
    #                    sigma2=plot_ARD1C3_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)$myr[input$num2,1],
    #                    rho=bxp_estim_rho_C3(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)$rhohat[input$num2,1]
    #         )                          
    #     } 
    # 
    #   
    # })
    
    
    output$bxARD1<-renderPlot({
        validate()
      
      
        if(input$but == "CAS 1"){
            p1<-plot_ARD1_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                             nb_simu = input$nbtest,input$num)[[1]]
            p2<-plot_ARD1_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                              nb_simu = input$nbtest,input$num)[[1]]
            grid.arrange(p1,p2,ncol=2, nrow = 1)
          
        }
        else if(input$but == "CAS 2"){
         
            p3<-plot_ARD1C2_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                               nb_simu = input$nbtest,input$num)[[1]]
            p4<-plot_ARD1C2_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                               nb_simu = input$nbtest,input$num)[[1]]
            #g3<-plot_estim_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,rhoo=seq(0,1,0.01))
            g4<-bxp_estim_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,input$num)[[1]]
            grid.arrange(p3,p4,g4,ncol=2, nrow = 2)
       
        }else if (input$but=="CAS 3"){
            p5<-plot_ARD1C3_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                               nb_simu = input$nbtest,input$num)[[1]]
            p6<-plot_ARD1C3_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                                nb_simu = input$nbtest,input$num)[[1]]
            #g5<-plot_estim_rho_C3(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                                  #nb_simu = input$nbtest,rhoo=seq(0,1,0.01)[-length(seq(0,1,0.01))])
            g6<-bxp_estim_rho_C3(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                                 nb_simu = input$nbtest,input$num)[[1]]
            grid.arrange(p5,p6,g6,ncol=2, nrow = 2)
        
        }
      else if (input$but=="CAS 4"){
        rhohat<-estim_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                          nb_simu = input$nbtest)
        p7<-bxp_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,input$nbtest,
                        input$num,rhoo=rhohat)[[1]]
        p8<-bxp_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                            input$nbtest,input$num,rhoo=rhohat)[[1]]
        #g7<-plot_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                     #nb_simu = input$nbtest)
        g8<-bxp_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,input$n,input$p,input$pm_periode,perio=T,
                   nb_simu = input$nbtest,input$num)[[1]]
        grid.arrange(p7,p8,g8,ncol=2, nrow = 2)
        
      }

    })
    
    
    
    output$mytable2 <- renderPlotly({
        validate()
        
        if(input$but2 == "CAS 1"){
            plotARD1(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)
          
             }
        else if (input$but2 == "CAS 2"){
            plotARD1_C2(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)
         
             }else if (input$but2 == "CAS 3"){
          
                 
            plotARD1_C3(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)
               
                   }else if (input$but2=="CAS 4"){
                         plotARD1_C4(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)
                    
                       }
    })
    
    
    
    output$bxARD2<-renderPlot({
        validate()
        if(input$but2 == "CAS 1"){
          
            p1<-plot_ARD1_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                             nb_simu = input$nbtest2,input$num2)[[1]]
            p2<-plot_ARD1_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                              nb_simu = input$nbtest2,input$num2)[[1]]
            grid.arrange(p1,p2,ncol=2, nrow = 1)
  
        }else if(input$but2 == "CAS 2"){
            
            p3<-plot_ARD1C2_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                             nb_simu = input$nbtest2,input$num2)[[1]]
            p4<-plot_ARD1C2_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                              nb_simu = input$nbtest2,input$num2)[[1]]
            #g3<-plot_estim_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,rhoo=seq(0,1,0.01))
            g4<-bxp_estim_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,input$nbtest2,input$num2)[[1]]
            
            grid.arrange(p3,p4,g4,ncol=2, nrow = 2)
         
        }else if (input$but2=="CAS 3"){
            
            p5<-plot_ARD1C3_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                               nb_simu = input$nbtest2,input$num2)[[1]]
            p6<-plot_ARD1C3_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                                nb_simu = input$nbtest2,input$num2)[[1]]
            #g5<-plot_estim_rho_C3(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                                  #nb_simu = input$nbtest2,rhoo=seq(0,1,0.01)[-length(seq(0,1,0.01))])
            g6<-bxp_estim_rho_C3(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                                 nb_simu = input$nbtest2,input$num2)[[1]]
            grid.arrange(p5,p6,g6,ncol=2, nrow = 2)
            
        }else if(input$but2=="CAS 4"){
          rhohat<-estim_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                  nb_simu = input$nbtest2)
          p7<-bxp_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                     input$num2,rhoo=rhohat)[[1]]
          p8<-bxp_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                      input$num2,rhoo=rhohat)[[1]]
        #g7<-plot_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                     # nb_simu = input$nbtest2)
          g8<-bxp_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,input$n2,input$p,input$pm_periode2,perio=F,
                      nb_simu = input$nbtest2,input$num2)[[1]]
          grid.arrange(p7,p8,g8,ncol=2, nrow = 2)
        
        }
        
    })
    
    
    observeEvent(input$button2, {
        a<-sample(10000,1)
        val$al<-a
        # b<-runif(1)*10^3
        # val$ale<-b
        
        # mu_x2 <- input$mu_x2
        # var2 <- input$var2
        # n2 <- input$n2
        # h2 <- input$h2
        # 
        # taille_echatillon2 = input$taille_echatillon2
        # nombre_des_echantillons2 = input$nombre_des_echantillons2
         # values$Boolean = FALSE
        # updateNumericInput(session, "mu_x2", value=0)
        # updateNumericInput(session, "var2", value=0)
        # updateNumericInput(session, "h2", value=0)
        # updateNumericInput(session, "n2", value=0)
        # updateNumericInput(session, "taille_echatillon2", value=0)
        # updateNumericInput(session, "nombre_des_echantillons2", value=0)
        # 
        # updateNumericInput(session, "mu_x2", value=mu_x2)
        # updateNumericInput(session, "var2", value=var2)
        # updateNumericInput(session, "h2", value=h2)
        # updateNumericInput(session, "n2", value=n2)
        # updateNumericInput(session, "taille_echatillon2", value=taille_echatillon2)
        # updateNumericInput(session, "nombre_des_echantillons2", value=nombre_des_echantillons2)
        # 
          })
    
    
    
    observeEvent(input$bu2, {
        b<-val$ale+input$nbtest2+1
        val$ale<-b
        # mu_x2 <- input$mu_x2
        # var2 <- input$var2
        # n2 <- input$n2
        # h2 <- input$h2
        # taille_echatillon2 = input$taille_echatillon2
        # nombre_des_echantillons2 = input$nombre_des_echantillons2
        # values$Boolean = TRUE
        # updateNumericInput(session, "mu_x2", value=0)
        # updateNumericInput(session, "var2", value=0)
        # updateNumericInput(session, "h2", value=0)
        # updateNumericInput(session, "n2", value=0)
        # updateNumericInput(session, "taille_echatillon2", value=0)
        # updateNumericInput(session, "nombre_des_echantillons2", value=0)
        # 
        # updateNumericInput(session, "mu_x2", value=mu_x2)
        # updateNumericInput(session, "var2", value=var2)
        # updateNumericInput(session, "h2", value=h2)
        # updateNumericInput(session, "n2", value=n2)
        # updateNumericInput(session, "taille_echatillon2", value=taille_echatillon2)
        # updateNumericInput(session, "nombre_des_echantillons2", value=nombre_des_echantillons2)
          })
    
    output$perf<-renderPlot({
        validate()
        if(input$but=="CAS 1"){
            bbc1m<-biais_C1_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=2),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest)
            bbc1s<-biais_C1_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=2),input$p,input$pm_periode,perio=T,
                               nb_simu = input$nbtest)
            beqm1m<-EQM_C1_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=2),input$p,input$pm_periode,perio=T,
                          nb_simu = input$nbtest)
            beqm1s<-EQM_C1_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=2),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest)
            
            gg1<-bbc1m$gg1
            gg2<-bbc1s$gg2
            gg3<-beqm1m$gg3
            gg4<-beqm1s$gg4
            grid.arrange(gg1,gg3,gg2,gg4,ncol=2,nrow=2)
        }else if (input$but=="CAS 2"){
            bbc2m<-biais_C2_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=8),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest,input$num)
            bbc2s<-biais_C2_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=8),input$p,input$pm_periode,perio=T,
                              nb_simu = input$nbtest,input$num)
            bbc2r<-biais_C2_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=8),input$p,input$pm_periode,perio=T,
                                nb_simu = input$nbtest,input$num)
            beqm2m<-EQM_C2_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=8),input$p,input$pm_periode,perio=T,
                          nb_simu = input$nbtest,input$num)
            beqm2s<-EQM_C2_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=8),input$p,input$pm_periode,perio=T,
                             nb_simu = input$nbtest,input$num)
            beqm2r<-EQM_C2_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=8),input$p,input$pm_periode,perio=T,
                               nb_simu = input$nbtest,input$num)
            gg1<-bbc2m$gg1
            gg2<-bbc2s$gg2
            gg3<-beqm2m$gg3
            gg4<-beqm2s$gg4
            gg5<-bbc2r$ggr
            gg6<-beqm2r$gg3
            grid.arrange(gg1,gg3,gg2,gg4,gg5,gg6,ncol=2,nrow=3)
            
        }else if (input$but=="CAS 3"){
            bbc3m<-biais_C3_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=10),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest,input$num)
            bbc3s<-biais_C3_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=10),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest,input$num)
            beqm3m<-EQM_C3_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=10),input$p,input$pm_periode,perio=T,
                          nb_simu = input$nbtest,input$num)
            beqm3s<-EQM_C3_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=10),input$p,input$pm_periode,perio=T,
                          nb_simu = input$nbtest,input$num)
            bbc3r<-biais_C3_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=10),input$p,input$pm_periode,perio=T,
                                nb_simu = input$nbtest,input$num)
            beqm3r<-EQM_C3_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=10),input$p,input$pm_periode,perio=T,
                                nb_simu = input$nbtest,input$num)
            
            gg1<-bbc3m$gg1
            gg2<-bbc3s$gg2
            gg3<-beqm3m$gg3
            gg4<-beqm3s$gg4
            gg5<-bbc3r$ggr
            gg6<-beqm3r$gg3
            grid.arrange(gg1,gg3,gg2,gg4,gg5,gg6,ncol=2,nrow=3) 
        }
      else if (input$but=="CAS 4"){
        
        bbc4m<-biais_C4_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=15),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest)
        bbc4s<-biais_C4_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=15),input$p,input$pm_periode,perio=T,
                            nb_simu = input$nbtest)
        beqm4m<-EQM_C4_mu(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=15),input$p,input$pm_periode,perio=T,
                          nb_simu = input$nbtest)
        beqm4s<-EQM_C4_sig(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=15),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest)
        bbc4r<-biais_C4_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=15),input$p,input$pm_periode,perio=T,
                            nb_simu = input$nbtest)
        beqm4r<-EQM_C4_rho(val$al,val$ale,input$mu_x,input$var1,input$rho,seq(input$pm_periode+2,input$n,by=15),input$p,input$pm_periode,perio=T,
                           nb_simu = input$nbtest)
        
        gg1<-bbc4m$gg
        gg2<-bbc4s$gg
        gg3<-beqm4m
        gg4<-beqm4s
        gg5<-bbc4r$gg
        gg6<-beqm4r
        grid.arrange(gg1,gg3,gg2,gg4,gg5,gg6,ncol=2,nrow=3) 
        
      }
    })
    
    output$perf2<-renderPlot({
        validate()
        if(input$but2=="CAS 1"){
            bbc1m<-biais_C1_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=2),input$p,input$pm_periode2,perio=F,
                           nb_simu = input$nbtest2)
            bbc1s<-biais_C1_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=2),input$p,input$pm_periode2,perio=F,
                               nb_simu = input$nbtest2)
            beqm1m<-EQM_C1_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=2),input$p,input$pm_periode2,perio=F,
                          nb_simu = input$nbtest2)
            beqm1s<-EQM_C1_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=2),input$p,input$pm_periode2,perio=F,
                              nb_simu = input$nbtest2)
            
            gg1<-bbc1m$gg1
            gg2<-bbc1s$gg2
            gg3<-beqm1m$gg3
            gg4<-beqm1s$gg4
            grid.arrange(gg1,gg3,gg2,gg4,ncol=2,nrow=2)
        }else if (input$but2=="CAS 2"){
            bbc2m<-biais_C2_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=8),input$p,input$pm_periode2,perio=F,
                           nb_simu = input$nbtest2,input$num2)
            bbc2s<-biais_C2_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=8),input$p,input$pm_periode2,perio=F,
                               nb_simu = input$nbtest2,input$num2)
            beqm2m<-EQM_C2_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=8),input$p,input$pm_periode2,perio=F,
                          nb_simu = input$nbtest2,input$num2)
            beqm2s<-EQM_C2_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=8),input$p,input$pm_periode2,perio=F,
                          nb_simu = input$nbtest2,input$num2)
            bbc2r<-biais_C2_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=8),input$p,input$pm_periode2,perio=F,
                                nb_simu = input$nbtest2,input$num2)
            beqm2r<-EQM_C2_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=8),input$p,input$pm_periode2,perio=F,
                                nb_simu = input$nbtest2,input$num2)
            
            gg1<-bbc2m$gg1
            gg2<-bbc2s$gg2
            gg3<-beqm2m$gg3
            gg4<-beqm2s$gg4
            gg5<-bbc2r$ggr
            gg6<-beqm2r$gg3
            grid.arrange(gg1,gg3,gg2,gg4,gg5,gg6,ncol=2,nrow=3)
            
        }else if (input$but2=="CAS 3"){
            bbc3m<-biais_C3_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                           nb_simu = input$nbtest2,input$num2)
            bbc3s<-biais_C3_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                           nb_simu = input$nbtest2,input$num2)
            beqm3m<-EQM_C3_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                          nb_simu = input$nbtest2,input$num2)
            beqm3s<-EQM_C3_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                          nb_simu = input$nbtest2,input$num2)
            bbc3r<-biais_C3_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                                nb_simu = input$nbtest2,input$num2)
            beqm3r<-EQM_C3_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                                nb_simu = input$nbtest2,input$num2)
            
            gg1<-bbc3m$gg1
            gg2<-bbc3s$gg2
            gg3<-beqm3m$gg3
            gg4<-beqm3s$gg4
            gg5<-bbc3r$ggr
            gg6<-beqm3r$gg3
            grid.arrange(gg1,gg2,gg3,gg4,gg5,gg6,ncol=2,nrow=3) 
            

        }else if (input$but2=="CAS 4"){
          
          bbc4m<-biais_C4_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                             nb_simu = input$nbtest2)
          bbc4s<-biais_C4_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                              nb_simu = input$nbtest2)
          beqm4m<-EQM_C4_mu(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                            nb_simu = input$nbtest2)
          beqm4s<-EQM_C4_sig(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                             nb_simu = input$nbtest2)
          bbc4r<-biais_C4_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                              nb_simu = input$nbtest2)
          beqm4r<-EQM_C4_rho(val$al,val$ale,input$mu_x2,input$var2,input$rho2,seq(input$pm_periode2+2,input$n2,by=10),input$p,input$pm_periode2,perio=F,
                             nb_simu = input$nbtest2)
          
          gg1<-bbc4m$gg
          gg2<-bbc4s$gg
          gg3<-beqm4m
          gg4<-beqm4s
          gg5<-bbc4r$gg
          gg6<-beqm4r
          grid.arrange(gg1,gg3,gg2,gg4,gg5,gg6,ncol=2,nrow=3) 
          
        }
      
    })
    
    
   
    
}
)













