# Use pred-prey model (2 focal prey sp + 1 "other" nuiscence sp) to predict
# potential equilibrium seaa otter density for Drakes Estero based on
# Model fit to equivalent data in Elkhorn Slough
require(readxl)
require(gtools)
require(parallel)
require(bayesplot)
require(ggplot2)
require(fitdistrplus)
require(stats)
require(cmdstanr)
require(posterior)
require(cowplot)
require(deSolve)
require(dplyr)
require(tidyr)
rstan::rstan_options(javascript=FALSE)
# 
# Load data ----------------------------------------------
# df_Cr1 = read_excel("Crab_historical_data.xlsx")
load("elkhorn_PredPrey_rslt_2021-12-09.rdata")
df_Cr = read_excel("Crab_data_Elkhorn_Tomales__Drakes_2011-2016.xlsx")
df_Cl = read_excel("Benthic Data_Drakes.xlsx")
df_crab_biomass_fxn = read_excel("Crab_biomass_params.xlsx")
# df_Ot = read_excel("Elkhorn Slough time series.xlsx")
# df_Fdat = read_excel("Foraging_Sum_elkhorn.xlsx")
#
Area_eelgr = 6.18
Area_bare = 19.02
Area_all = 25.2
#  
r_adjust = 1 
# Process data ------------------------------------------
# Crab counts (prey X1)
df_Cr_all = df_Cr
df_Cr = df_Cr_all[which(df_Cr_all$Site=="Tomales Bay"),]
df_Cr$Size_mm = as.numeric(df_Cr$`Size (mm)`)
tmp = df_Cr %>%
  left_join(df_crab_biomass_fxn,by="Species") %>%
  select(a_par,b_par)
df_Cr$Biomass = rep(0,nrow(df_Cr))
ii = which(df_Cr$Num>0)  
df_Cr$Biomass[ii] = df_Cr$Num[ii] * 1.2 * (tmp$a_par[ii] * df_Cr$Size_mm[ii]^tmp$b_par[ii])
df_Cr$Num_traps = rep(1,nrow(df_Cr))
df_Cr_sum = df_Cr %>% group_by(Record_ID) %>% 
  summarise(year = mean(Year),
            num_cr = sum(Num),
            BM_tot = sum(Biomass),
            num_trp = mean(Num_traps))
df_Cr_sum$BM_units = round(df_Cr_sum$BM_tot/100)
df_Cr_sum$Record_ID = as.numeric(as.factor(df_Cr_sum$Record_ID))
crab_counts = round(df_Cr_sum$BM_units/df_Cr_sum$num_trp)
# Clam counts (holes that could be clam siphons)
ii = which(df_Cl$Species=="Fat")
df_Cl = df_Cl[-ii,]
ii = which(df_Cl$`Eelgrass or Bare`=="E")
clam_counts_E = df_Cl$Count[ii]
ii = which(df_Cl$`Eelgrass or Bare`=="B")
clam_counts_B = df_Cl$Count[ii]
# Fit negative binomial distributions to prey abundances
invisible(capture.output(ftCr <- fitdist(crab_counts,"nbinom")))
invisible(capture.output(ftClE <- fitdist(clam_counts_E,"nbinom")))
invisible(capture.output(ftClB <- fitdist(clam_counts_B,"nbinom")))
NBpars = data.frame(Prey = c("Crab","Crab","Clam","Clam"),
                    Habitat = c("Eelgrass","Bare","Eelgrass","Bare"),
                    mu = c(ftCr$estimate[2],ftCr$estimate[2]*.2,
                           ftClE$estimate[2],ftClB$estimate[2] ),
                    mu_sd = c(ftCr$sd[2],ftCr$sd[2],
                           ftClE$sd[2],ftClB$sd[2] ),
                    size = c(ftCr$estimate[1],ftCr$estimate[1],
                             ftClE$estimate[1],ftClB$estimate[1] ),
                    size_sd = c(ftCr$sd[1],ftCr$sd[1],
                              ftClE$sd[1],ftClB$sd[1] ))
# 
# Estimate Potential density for eelgrass ------------------------------
Area = Area_eelgr; 
Y_init = 2/Area; reps = 1000; Time = 200; 
set.seed(123); rr = sample(Nsims,reps)
source("pred2prey.R")
Krnd = matrix(0,nrow = reps,ncol = 2)
ii = which(NBpars$Habitat=="Eelgrass")
Krnd[,1] = pmax(0.01,rnorm(reps,NBpars$mu[ii[1]],NBpars$mu_sd[ii[1]])) * (1/omega[1])
Krnd[,2] = pmax(0.001,rnorm(reps,NBpars$mu[ii[2]],NBpars$mu_sd[ii[2]])) * (1/omega[2])
Ottdens_K_E = numeric(length = reps)
for(i in 1:reps){
  pars = list(
    K = Krnd[i,],
    r = r_adjust * as.numeric(mcmc[rr[i],which(startsWith(vn,"r["))]),
    phi = as.numeric(mcmc[rr[i],which(startsWith(vn,"phi["))]),
    z1 = as.numeric(mcmc[rr[i],which(startsWith(vn,"z1"))]),
    z2 = as.numeric(mcmc[rr[i],which(startsWith(vn,"z2"))]),
    a3 = as.numeric(mcmc[rr[i],which(startsWith(vn,"a3"))]),
    a =  as.numeric(mcmc[rr[i],which(startsWith(vn,"a["))]),
    h = stan.data$h ,
    g = stan.data$g,
    Tcf = stan.data$Tcf
  )
  initstate = c(pars$K[1], pars$K[2], Y_init)
  rslt <- ode(initstate, 1:Time, pred2prey, pars,
              method = "ode45")
  if(nrow(rslt)<Time | rslt[nrow(rslt),4] > 57 ){
    Ottdens_K_E[i] = NA 
  }else{
    Ottdens_K_E[i] = rslt[Time,4]*Area
  }
}
Ottdens_K_E = Ottdens_K_E[which(!is.na(Ottdens_K_E))]
# hist(Ottdens_K,25)
invisible(capture.output(ft <- fitdist(Ottdens_K_E,"gamma")))
# plot(ft)
a = as.numeric(ft$estimate[1]); b = as.numeric(ft$estimate[2])
Mean_K_E = a/b
CI_K_E_lo = qgamma(.05,a,b)
CI_K_E_hi = qgamma(.95,a,b)
noquote(paste0("Mean estimated abundance at K, areas with eelgrass = ",
               format(Mean_K_E,digits = 3)," (",
               format(CI_K_E_lo,digits = 3)," - ",
               format(CI_K_E_hi,digits = 3),")"))
#
# Estimate Potential density for bare areas ------------------------------
Area = Area_bare
Krnd = matrix(0,nrow = reps,ncol = 2)
ii = which(NBpars$Habitat=="Bare")
Krnd[,1] = pmax(0.01,rnorm(reps,NBpars$mu[ii[1]],NBpars$mu_sd[ii[1]])) * (1/omega[1])
Krnd[,2] = pmax(0.01,rnorm(reps,NBpars$mu[ii[2]],NBpars$mu_sd[ii[2]])) * (1/omega[2])
Ottdens_K_B = numeric(length = reps)
for(i in 1:reps){
  pars = list(
    K = Krnd[i,],
    r = r_adjust * as.numeric(mcmc[rr[i],which(startsWith(vn,"r["))]),
    phi = as.numeric(mcmc[rr[i],which(startsWith(vn,"phi["))]),
    z1 = as.numeric(mcmc[rr[i],which(startsWith(vn,"z1"))]),
    z2 = as.numeric(mcmc[rr[i],which(startsWith(vn,"z2"))]),
    a3 = as.numeric(mcmc[rr[i],which(startsWith(vn,"a3"))]),
    a =  as.numeric(mcmc[rr[i],which(startsWith(vn,"a["))]),
    h = stan.data$h ,
    g = stan.data$g,
    Tcf = stan.data$Tcf
  )
  initstate = c(pars$K[1], pars$K[2], Y_init)
  rslt <- ode(initstate, 1:Time, pred2prey, pars,
              method = "ode45")
  if(nrow(rslt)<Time | rslt[nrow(rslt),4] > 57 ){
    Ottdens_K_B[i] = NA 
  }else{
    Ottdens_K_B[i] = rslt[Time,4]*Area
  }
}
Ottdens_K_B = Ottdens_K_B[which(!is.na(Ottdens_K_B))]
# hist(Ottdens_K,25)
invisible(capture.output(ft <- fitdist(Ottdens_K_B,"gamma")))
# plot(ft)
a = as.numeric(ft$estimate[1]); b = as.numeric(ft$estimate[2])
Mean_K_B = a/b
CI_K_B_lo = qgamma(.05,a,b)
CI_K_B_hi = qgamma(.95,a,b)
noquote(paste0("Mean estimated abundance at K, outside eelgrass areas = ",
               format(Mean_K_B,digits = 3)," (",
               format(CI_K_B_lo,digits = 3)," - ",
               format(CI_K_B_hi,digits = 3),")"))
#
savename = paste0("Tomales_PredPrey_rslt_radj-",r_adjust,"_",Sys.Date())
#
# save.image(file=paste0(savename,".rdata"))
#
# Both areas -----------------------------------------------------------------
ii = min(length(Ottdens_K_E),length(Ottdens_K_B))
Ottdens_K_TB = Ottdens_K_E[1:ii] + Ottdens_K_B[1:ii]
invisible(capture.output(ft <- fitdist(Ottdens_K_TB,"gamma")))
# plot(ft)
a = as.numeric(ft$estimate[1]); b = as.numeric(ft$estimate[2])
OttK_TB_rnd = rgamma(5000,a,b)
Mean_K_TB = a/b
CI_K_TB_lo = qgamma(.05,a,b)
CI_K_TB_hi = qgamma(.95,a,b)
noquote(paste0("Mean estimated abundance at K, for all of Tomales Bay = ",
               format(Mean_K_TB,digits = 3)," (",
               format(CI_K_TB_lo,digits = 3)," - ",
               format(CI_K_TB_hi,digits = 3),")"))
# Compare to Elkhorn Slough
invisible(capture.output(ft <- fitdist(Ottdens_K,"gamma")))
a = as.numeric(ft$estimate[1]); b = as.numeric(ft$estimate[2])
OttK_ES_rnd = rgamma(5000,a,b)
noquote(paste0("Mean estimated abundance at K, Elkhorn Slough = ",
               format(Mean_K,digits = 3)," (",
               format(CI_lo,digits = 3)," - ",
               format(CI_hi,digits = 3),")"))


# Compare estimates as violin plots
if(r_adjust==1){
  # Load results for Tomales with r reduced by 50%
  load("Tomales_reps_half_r.rdata")
  K_stats_sum = data.frame(Area = factor(c("Elkhorn Slough","Tomales Bay","Tomales Bay r/2"),
                                         levels = c("Elkhorn Slough","Tomales Bay",
                                                    "Tomales Bay r/2")),
                           Mean = c(Mean_K,Mean_K_TB,Mean_K2_TB),
                           CI_lo = c(CI_lo,CI_K_TB_lo,CI_K_TB2_lo),
                           CI_hi = c(CI_hi,CI_K_TB_hi,CI_K_TB2_hi))
  K_est_reps = data.frame(Area = factor(c(rep("Elkhorn Slough", 5000),
                                          rep("Tomales Bay", 5000),
                                          rep("Tomales Bay r/2", 5000)),
                                        levels = c("Elkhorn Slough","Tomales Bay",
                                                   "Tomales Bay r/2")),
                          Value = c(OttK_ES_rnd,OttK_TB_rnd,OttK_TB2_rnd))
  
}else{
  # Load results for Drakes with r reduced by 50%
  load("Drakes_reps_half_r.rdata")
  OttK_TB2_rnd = OttK_TB_rnd
  Mean_K2_TB = Mean_K_TB
  CI_K_TB2_lo = CI_K_TB_lo 
  CI_K_TB2_hi = CI_K_TB_hi
  K_stats_sum = data.frame(Area = factor(c("Elkhorn Slough","Drakes Estero r/2","Tomales Bay r/2"),
                                         levels = c("Elkhorn Slough","Drakes Estero r/2",
                                                    "Tomales Bay r/2")),
                           Mean = c(Mean_K,Mean_K2_DE,Mean_K2_TB),
                           CI_lo = c(CI_lo,CI_K_DE2_lo,CI_K_TB2_lo),
                           CI_hi = c(CI_hi,CI_K_DE2_hi,CI_K_TB2_hi))
  K_est_reps = data.frame(Area = factor(c(rep("Elkhorn Slough", 5000),
                                          rep("Drakes Estero r/2", 5000),
                                          rep("Tomales Bay r/2", 5000)),
                                        levels = c("Elkhorn Slough","Drakes Estero r/2",
                                                   "Tomales Bay r/2")),
                          Value = c(OttK_ES_rnd,OttK_DE2_rnd,OttK_TB2_rnd))
  save(OttK_TB2_rnd, Mean_K2_TB, CI_K_TB2_lo, CI_K_TB2_hi,
       file = "Tomales_reps_half_r.rdata")
}
plt_K_compare = ggplot(K_est_reps,aes(x=Area,y=Value,group=Area,fill=Area)) +
  geom_violin() +
  geom_errorbar(data=K_stats_sum,aes(x=Area,y=Mean,ymin=CI_lo,ymax=CI_hi),
                width=0.1,size=1.05) +
  geom_point(data=K_stats_sum,aes(x=Area,y=Mean),size=3,shape=21,fill="white") +
  labs(x = "Estuary", y = "Estiamted equilibrium density") +
  theme_classic() + theme(legend.position = "none")
print(plt_K_compare)  
# 
# If r_adjust==1, do another graph to compare Elkhorn, Drakes, Tomales at full r
if(r_adjust==1){
  load("Drakes_reps_full_r.rdata")
  OttK_TB2_rnd = OttK_TB_rnd
  Mean_K2_TB = Mean_K_TB
  CI_K_TB2_lo = CI_K_TB_lo 
  CI_K_TB2_hi = CI_K_TB_hi 
  K_stats_sum2 = data.frame(Area = factor(c("Elkhorn Slough","Drakes Estero","Tomales Bay"),
                                         levels = c("Elkhorn Slough","Drakes Estero",
                                                    "Tomales Bay")),
                           Mean = c(Mean_K,Mean_K2_DE,Mean_K2_TB),
                           CI_lo = c(CI_lo,CI_K_DE2_lo,CI_K_TB2_lo),
                           CI_hi = c(CI_hi,CI_K_DE2_hi,CI_K_TB2_hi))
  K_est_reps2 = data.frame(Area = factor(c(rep("Elkhorn Slough", 5000),
                                          rep("Drakes Estero", 5000),
                                          rep("Tomales Bay", 5000)),
                                        levels = c("Elkhorn Slough","Drakes Estero",
                                                   "Tomales Bay")),
                          Value = c(OttK_ES_rnd,OttK_DE2_rnd,OttK_TB2_rnd))
  #
  plt_K_compare2 = ggplot(K_est_reps2,aes(x=Area,y=Value,group=Area,fill=Area)) +
    geom_violin() +
    geom_errorbar(data=K_stats_sum2,aes(x=Area,y=Mean,ymin=CI_lo,ymax=CI_hi),
                  width=0.1,size=1.05) +
    geom_point(data=K_stats_sum2,aes(x=Area,y=Mean),size=3,shape=21,fill="white") +
    labs(x = "Estuary", y = "Estiamted equilibrium density") +
    theme_classic() + theme(legend.position = "none")
  print(plt_K_compare2)  
}
#
# Sample dynamics ------------------------------------------------------
# Tomales
Area = Area_eelgr; reps = 1000; Time = 75; Y_init = 2/Area;
ii = which(NBpars$Habitat=="Eelgrass")
pars = list(
  K = c(NBpars$mu[ii[1]]*(1/omega[1]), NBpars$mu[ii[2]]*(1/omega[2])) ,
  r = r_adjust * sumstats[which(startsWith(vns,"r[")),1],
  phi = sumstats[which(startsWith(vns,"phi[")),1],
  z1 = sumstats[which(startsWith(vns,"z1")),1],
  z2 = sumstats[which(startsWith(vns,"z2")),1],
  a3 = sumstats[which(startsWith(vns,"a3")),1],
  a =  sumstats[which(startsWith(vns,"a[")),1],
  h = stan.data$h ,
  g = stan.data$g,
  Tcf = stan.data$Tcf
)
initstate = c(pars$K[1], pars$K[2], Y_init)
PPmod <- ode(initstate, 1:Time, pred2prey, pars, method = "ode45")
dfPP = data.frame(Year = 1:Time,
                  Crabs_CPUE = PPmod[,2]*omega[1], 
                  Clams_25m2 = PPmod[,3]*25, 
                  Otters_km2 = PPmod[,4])
#
dfPP = dfPP %>% pivot_longer(cols = c(Crabs_CPUE,Clams_25m2,Otters_km2),
                             names_to = "Species",
                             values_to = "Abundance")
plt_dynam_TB = ggplot(data = dfPP, aes(x=Year,y=Abundance,group=Species,color=Species)) +
  geom_line() + labs(x = "Year from otter colonization",y="Projected abundance") + 
  ggtitle("Model-projected dynamics Tomales Bay (eelgrass)") +
  ylim(0,50) +
  theme_classic()
# Drakes
Area = 5.98; Y_init = 2/Area;
ii = which(NBpars2$Habitat=="Eelgrass")
pars = list(
  K = c(NBpars2$mu[ii[1]]*(1/omega[1]), NBpars2$mu[ii[2]]*(1/omega[2])) ,
  r = r_adjust * sumstats[which(startsWith(vns,"r[")),1],
  phi = sumstats[which(startsWith(vns,"phi[")),1],
  z1 = sumstats[which(startsWith(vns,"z1")),1],
  z2 = sumstats[which(startsWith(vns,"z2")),1],
  a3 = sumstats[which(startsWith(vns,"a3")),1],
  a =  sumstats[which(startsWith(vns,"a[")),1],
  h = stan.data$h ,
  g = stan.data$g,
  Tcf = stan.data$Tcf
)
initstate = c(pars$K[1], pars$K[2], Y_init)
PPmod <- ode(initstate, 1:Time, pred2prey, pars, method = "ode45")
dfPP = data.frame(Year = 1:Time,
                  Crabs_CPUE = PPmod[,2]*omega[1], 
                  Clams_25m2 = PPmod[,3]*25, 
                  Otters_km2 = PPmod[,4])
#
dfPP = dfPP %>% pivot_longer(cols = c(Crabs_CPUE,Clams_25m2,Otters_km2),
                             names_to = "Species",
                             values_to = "Abundance")
plt_dynam_DE = ggplot(data = dfPP, aes(x=Year,y=Abundance,group=Species,color=Species)) +
  geom_line() + labs(x = "Year from otter colonization",y="Projected abundance") + 
  ggtitle("Model-projected dynamics Drakes Estero (eelgrass)") +
  ylim(0,50) +
  theme_classic()
# print(plt_dynam_DE)
#
# Elkorn Slough
Area = 3.45; Y_init = 2/Area;
pars = list(
  K = sumstats[which(startsWith(vns,"K[")),1],
  r = sumstats[which(startsWith(vns,"r[")),1],
  phi = sumstats[which(startsWith(vns,"phi[")),1],
  z1 = sumstats[which(startsWith(vns,"z1")),1],
  z2 = sumstats[which(startsWith(vns,"z2")),1],
  a3 = sumstats[which(startsWith(vns,"a3")),1],
  a =  sumstats[which(startsWith(vns,"a[")),1],
  h = stan.data$h ,
  g = stan.data$g,
  Tcf = stan.data$Tcf
)
initstate = c(pars$K[1], pars$K[2], Y_init)
PPmod <- ode(initstate, 1:Time, pred2prey, pars, method = "ode45")
dfPP = data.frame(Year = 1:Time,
                  Crabs_CPUE = PPmod[,2]*omega[1], 
                  Clams_25m2 = PPmod[,3]*25, 
                  Otters_km2 = PPmod[,4])
#
dfPP = dfPP %>% pivot_longer(cols = c(Crabs_CPUE,Clams_25m2,Otters_km2),
                             names_to = "Species",
                             values_to = "Abundance")
plt_dynam_ES = ggplot(data = dfPP, aes(x=Year,y=Abundance,group=Species,color=Species)) +
  geom_line() + labs(x = "Year from otter colonization",y="Projected abundance") + 
  ggtitle("Model-projected dynamics Elkhorn Slough") +
  ylim(0,50) +
  theme_classic()
#print(plt_dynam_ES)
#
# Combined plot (using cowplot)
combined_plot <- plot_grid(plt_dynam_ES + theme(legend.position = "none"),
                           plt_dynam_DE + theme(legend.position = "none"), 
                           plt_dynam_TB + theme(legend.position = "none"),
                           labels = c('A', 'B', 'C'), label_size = 12, nrow = 3)
legend <- get_legend(plt_dynam_ES + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") )
plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, .1))