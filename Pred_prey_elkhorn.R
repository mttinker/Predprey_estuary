# Elkhorn SLough pred-prey model (2 focal prey sp + 1 "other" nuiscence sp)
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
require(htmlTable)
rstan::rstan_options(javascript=FALSE)
# 
# Load data ----------------------------------------------
df_Cr1 = read_excel("Crab_historical_data.xlsx")
df_Cr2 = read_excel("Crab_data_Elkhorn_Tomales__Drakes_2011-2016.xlsx")
df_crab_biomass_fxn = read_excel("Crab_biomass_params.xlsx")
df_Cl = read_excel("Clam_Sum_elkhorn.xlsx")
df_Ot = read_excel("Elkhorn Slough time series.xlsx")
df_Fdat = read_excel("Foraging_Sum_elkhorn.xlsx")
# Process data ------------------------------------------
Npr = 3     # number prey types (default 3: Crabs, clams, other)
TS = seq(1,max(df_Ot$YearN))
minyr = min(df_Ot$Year)
Years = TS + minyr
NT = max(TS)
#
# Otter counts
ii = which(df_Ot$Year>minyr & df_Ot$Usedat==1)
ixy = df_Ot$Year[ii] - minyr; Ny = length(ixy)
Yobs = df_Ot$USGS_count[ii]
# Crab counts (prey X1)
# Historical
tmp = df_Cr1 %>%
  left_join(df_crab_biomass_fxn,by="Species") %>%
  select(a_par,b_par)
df_Cr1$Biomass = 1.2 * (tmp$a_par * df_Cr1$Size_mm^tmp$b_par)
df_Cr1_sum = df_Cr1 %>% group_by(Record_ID) %>% 
  summarise(year = mean(Year),
            num_cr = sum(Number),
            BM_tot = sum(Biomass),
            num_trp = mean(Num_traps))
df_Cr1_sum$BM_units = round(df_Cr1_sum$BM_tot/100)
# Recent
df_Cr2_all = df_Cr2
df_Cr2 = df_Cr2[which(df_Cr2$Site=="Elkhorn Slough"),]
df_Cr2$Size_mm = as.numeric(df_Cr2$`Size (mm)`)
tmp = df_Cr2 %>%
  left_join(df_crab_biomass_fxn,by="Species") %>%
  select(a_par,b_par)
df_Cr2$Biomass = rep(0,nrow(df_Cr2))
ii = which(df_Cr2$Num>0)  
df_Cr2$Biomass[ii] = df_Cr2$Num[ii] * 1.2 * (tmp$a_par[ii] * df_Cr2$Size_mm[ii]^tmp$b_par[ii])
df_Cr2$Num_traps = rep(1,nrow(df_Cr2))
df_Cr2_sum = df_Cr2 %>% group_by(Record_ID) %>% 
  summarise(year = mean(Year),
            num_cr = sum(Num),
            BM_tot = sum(Biomass),
            num_trp = mean(Num_traps))
ii = which(df_Cr2_sum$BM_tot>0 & df_Cr2_sum$BM_tot<100)
df_Cr2_sum$BM_tot[ii] = df_Cr2_sum$BM_tot[ii]*4
df_Cr2_sum$num_trp[ii] = df_Cr2_sum$num_trp[ii]*4
df_Cr2_sum$BM_units = round(df_Cr2_sum$BM_tot/100)
df_Cr2_sum$Record_ID = as.numeric(as.factor(df_Cr2_sum$Record_ID))
df_Cr_sum = rbind(df_Cr1_sum,df_Cr2_sum)
df_Cr_sum$BM_ptr = df_Cr_sum$BM_units/df_Cr_sum$num_trp
# ggplot(df_Cr_sum,aes(x=as.factor(year),y=BM_ptr)) +
#   geom_boxplot()
ix1 = pmax(1, df_Cr_sum$year - minyr); Nx1 = length(ix1); 
ns1 = df_Cr_sum$num_trp
X1obs = df_Cr_sum$BM_units
invisible(capture.output(ft <- fitdist(X1obs[ix1==1],"nbinom")))
tauX1pri = as.numeric(coef(ft)[1])
# Clam counts (prey X2)
# First, determine which records have dredge counts as well as hole counts
ii = which(!is.na(df_Cl$N_clams))
# Adjusting factor, edible-size clams to holes
hladj = sum(df_Cl$N_holes[ii])/(.7 * sum(df_Cl$N_clams[ii]))
ix2 = c(pmax(1,df_Cl$Year - minyr),pmax(1,df_Cl$Year[ii] - minyr))
Nx2 = length(ix2); ix2p = as.numeric(as.factor(ix2)); Nx2p = max(ix2p)
ns2 = c(df_Cl$Nquad,df_Cl$Nquad[ii])
X2obs = c(df_Cl$N_holes,round(df_Cl$N_clams[ii]*hladj))
invisible(capture.output(ft <- fitdist(X2obs[ix2==1],"nbinom")))
tauX2pri = as.numeric(coef(ft)[1])
# Foraging data
df_Fdat$preyF = factor(df_Fdat$Prey,levels=unique(df_Fdat$Prey))
df_Fdat$preyN = as.numeric(df_Fdat$preyF)
# Proportion effort allocation
ixpY = unique(df_Fdat$Year); ixp = ixpY - minyr; Np = length(ixp)
P_obs = matrix(0,nrow = Np, ncol = Npr) 
P_obs_sd = matrix(0,nrow = Np, ncol = Npr) 
tauP = numeric()
for (i in 1:Np){
  P_obs[i,1:Npr] = df_Fdat$P_effort[df_Fdat$Year==ixpY[i]]
  P_obs_sd[i,1:Npr] = df_Fdat$P_sd[df_Fdat$Year==ixpY[i]]
  p1 = P_obs[i,2]
  sd_p1 = P_obs_sd[i,2]*2
  tp = 1
  sdtst = sqrt((tp*p1*(tp-tp*p1))/(tp^2 *(tp+1) ))  
  while (abs(sdtst - sd_p1) >.005){
    tp = tp+1
    sdtst = sqrt((tp*p1*(tp-tp*p1))/(tp^2 *(tp+1) )) 
  } 
  tauP[i] = tp
}
# Energy rate estimates
ii = which(!is.na(df_Fdat$E_rate))
ixfY = df_Fdat$Year[ii]; ixf = ixfY - minyr; Nf = length(ixf)
ixfsp = df_Fdat$preyN[ii]
tauF = numeric()
ER = numeric()
for(i in 1:Nf){
  ER[i] = df_Fdat$E_rate[ii[i]]
  tauF[i] = ER[i] / (df_Fdat$ER_sd[ii[i]] )^2
}
# Various other model parameter constants
h = df_Fdat$HT_min[1:3]
g = df_Fdat$EG_item[1:3]
Tcf = (60*24*(365)*.35)/1000000 # CF: foraging min per year per m2pred 
Acf = 3.45       # conversion factor, 1 km2 to total sub-tidal area of slough
Y_init = df_Ot$USGS_count[df_Ot$Year==minyr] / Acf
omega_cr = pi * 2^2 # assume crabs caught within circle of 2m radius,
omega_cl = 5 # clam conversion = # holes / # edible clams per quad (~5)
omega = c(omega_cr,omega_cl)
# KestCr = mean(df_Cr$CPUE[df_Cr$YEAR==min(df_Cr$YEAR)])*(1/omega_cr) 
# KestCl = mean(df_Cl$NH_p_q[df_Cl$Year==min(df_Cl$Year)])*(1/omega_cl) 
# Kcf = c(KestCr,KestCl)
#
# Prep to fit model --------------------------------------------------------
#
stan.data <- list(N=NT,Nx1=Nx1,ix1=ix1,ns1=ns1,X1obs=X1obs,Nx2=Nx2,ix2=ix2,
                  ns2=ns2,X2obs=X2obs,Nx2p=Nx2p,ix2p=ix2p,Ny=Ny,ixy=ixy,Yobs=Yobs,
                  Np=Np,ixp=ixp,Pobs=P_obs,Nf=Nf,ixf=ixf,ixfsp=ixfsp,ER=ER,
                  g=g,h=h,Tcf=Tcf,Acf=Acf,omega=omega,ts=TS,
                  tauP=tauP,tauF=tauF,Y_init=Y_init)
# 
parms <- c("tauX1","tauX2","tauY","z1","z2","a0","a","a3","r","K","phi",
           "Yexp","X1exp","X2exp","ERexp","Pexp") # 
#
nburnin = 1000                    # number warm-up (burn-in) samples
nsamples = 10000                 # desired total number samples
fitmodel = c("pred_prey_elksl.stan")    # name of file with stan code
cores = detectCores()
ncore = min(5,cores-1)
Niter = round(nsamples/ncore)
mod <- cmdstan_model(fitmodel)   # compiles model (if necessary)
# Fit model ------------------------------------------------------------------
suppressMessages(                # Suppress messages/warnings (if desired)
  suppressWarnings (
    fit <- mod$sample(
      data = stan.data,
      seed = 123,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      iter_warmup = nburnin,
      iter_sampling = Niter,
      max_treedepth = 15,
      adapt_delta = 0.995
    )
  )
)
# tmp = fit$output(); tmp[[1]][40:60]
source("cmdstan_sumstats.r")
#
#
sampler_diagnostics <- fit$sampler_diagnostics()
str(sampler_diagnostics)
diagtab = fit$sampler_diagnostics(format = "df")
#
# Summarize model ----------------------------------------------------------
#
print(mcmc_trace(fit$draws("phi")))
print(mcmc_trace(fit$draws("tauY")))
#for (p in 1:11){
#  print(mcmc_trace(fit$draws(parms[p])))
#}

mcmc_areas(fit$draws(variables = paste0("K[",seq(1,2),"]")),
           area_method="equal height",
           prob = 0.8) + 
  scale_y_discrete(labels = paste0("K estimate, prey",seq(1,2))) +
  ggtitle("Parameter posterior distributions") +
  labs(x="Parameter value",y="Posterior density") +
  theme_classic()
mcmc_areas(fit$draws(variables = paste0("r[",seq(1,2),"]")),
           area_method="equal height",
           prob = 0.8) + 
  scale_y_discrete(labels = paste0("r estimate, prey",seq(1,2))) +
  ggtitle("Parameter posterior distributions") +
  labs(x="Parameter value",y="Posterior density") +
  theme_classic()
mcmc_areas(fit$draws(variables = paste0("a0[",seq(1,2),"]")),
           area_method="equal height",
           prob = 0.8) + 
  scale_y_discrete(labels = paste0("encounter rate (at K), prey ",seq(1,2))) +
  ggtitle("Parameter posterior distributions") +
  labs(x="Parameter value",y="Posterior density") +
  theme_classic()
mcmc_areas(fit$draws(variables = paste0("phi[",seq(1,2),"]")),
           area_method="equal height",
           prob = 0.8) + 
  scale_y_discrete(labels = paste0("prey preference param ",seq(1,2))) +
  ggtitle("Parameter posterior distributions") +
  labs(x="Parameter value",y="Posterior density") +
  theme_classic()
#
Ypred = sumstats[which(startsWith(vns,"Yexp[")),1]
Ypred_lo = sumstats[which(startsWith(vns,"Yexp[")),4]
Ypred_hi = sumstats[which(startsWith(vns,"Yexp[")),8]
df_Yplt = data.frame(Year=Years,Ypred=Ypred,
                     Ypred_lo=Ypred_lo,Ypred_hi=Ypred_hi)
df_Yobs = data.frame(Year = ixy+minyr, Pred_obs = Yobs)
plt_Yexp = ggplot(df_Yplt,aes(x=Year,y=Ypred)) +
  geom_ribbon(aes(ymin=Ypred_lo,ymax=Ypred_hi),alpha=0.3) +
  geom_line(color = "blue") +
  geom_point(data=df_Yobs,aes(x=Year,y=Pred_obs)) +
  labs(x="Year",y="Otter abundance") +
  ggtitle("Otter trends, observed counts vs. model estimates") +
  theme_classic()
print(plt_Yexp)

X1pred = sumstats[which(startsWith(vns,"X1exp[")),1]
X1pred_lo = sumstats[which(startsWith(vns,"X1exp[")),4]
X1pred_hi = sumstats[which(startsWith(vns,"X1exp[")),8]
df_X1plt = data.frame(Year=Years,X1pred=X1pred,
                      X1pred_lo=X1pred_lo,X1pred_hi=X1pred_hi)
df_X1plt$Year = factor(df_X1plt$Year, levels = as.character(Years))
df_X1obs = data.frame(Year = ix1+minyr, Prey_obs = X1obs/ns1 )
df_X1obs$Year = factor(df_X1obs$Year, levels = as.character(Years))
#
plt_X1exp = ggplot(df_X1plt,aes(x=Year,y=X1pred,group=1)) +
  geom_ribbon(aes(ymin=X1pred_lo,ymax=X1pred_hi),alpha=0.3) +
  geom_line(color = "blue") +
  geom_boxplot(data=df_X1obs,aes(x=Year,y=Prey_obs,group=Year)) +
  # geom_point(data=df_X1obs,aes(x=Year,y=Prey_obs)) +
  labs(x="Year",y="Crab abundance (CPUE)") +
  ylim(0,30) +
  ggtitle("Crab trends, observed vs. model estimates") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#print(plt_X1exp)

X2pred = sumstats[which(startsWith(vns,"X2exp[")),1]
X2pred_lo = sumstats[which(startsWith(vns,"X2exp[")),4]
X2pred_hi = sumstats[which(startsWith(vns,"X2exp[")),8]
df_X2plt = data.frame(Year=Years,X2pred=X2pred,
                      X2pred_lo=X2pred_lo,X2pred_hi=X2pred_hi)
df_X2plt$Year = factor(df_X2plt$Year, levels = as.character(Years))
df_X2obs = data.frame(Year = ix2+minyr, Prey_obs = X2obs/ns2 )
df_X2obs$Year = factor(df_X2obs$Year, levels = as.character(Years))
#
plt_X2exp = ggplot(df_X2plt,aes(x=Year,y=X2pred,group=1)) +
  geom_ribbon(aes(ymin=X2pred_lo,ymax=X2pred_hi),alpha=0.3) +
  geom_line(color = "blue") +
  geom_boxplot(data=df_X2obs,aes(x=Year,y=Prey_obs,group=Year)) +
  # geom_point(data=df_X2obs,aes(x=Year,y=Prey_obs)) +
  labs(x="Year",y="Clam abundance (siphons per quadrat)") +
  ggtitle("Clam trends, observed vs. model estimates") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#print(plt_X2exp)
#  
plot_grid(plt_X1exp,plt_X2exp,nrow = 2)

Preytp = character()
Ppn_ef = numeric()
Yr = numeric()
ErateP1_mn = numeric()
ErateP1_lo = numeric()
ErateP1_hi = numeric()
ErateP2_mn = numeric()
ErateP2_lo = numeric()
ErateP2_hi = numeric()
Erate_mn = numeric()
Erate_lo = numeric()
Erate_hi = numeric()
a1 = mcmc[,which(vn==paste0("a[",1,"]"))]  
a2 = mcmc[,which(vn==paste0("a[",2,"]"))]  
a3 = mcmc[,which(vn=="a3")]  
for(i in 1:NT){
  p1 = mcmc[,which(vn==paste0("Pexp[",i,",",1,"]"))]  
  p2 = mcmc[,which(vn==paste0("Pexp[",i,",",2,"]"))]
  p3 = mcmc[,which(vn==paste0("Pexp[",i,",",3,"]"))]
  Preytp = c(Preytp,c("Crab","Clam","Other"))
  Ppn_ef = c(Ppn_ef,c(mean(p1),mean(p2),mean(p3)))
  Yr = c(Yr,rep(Years[i],3))
  y1 = mcmc[,which(vn==paste0("X1exp[",i,"]"))] / omega[1] 
  y2 = mcmc[,which(vn==paste0("X2exp[",i,"]"))] / omega[2]   
  Erate = (g[1]*p1*a1*y1 + g[2]*p2*a2*y2 + g[3]*p3*a3) / 
     (1 + p1*h[1]*a1*y1 + p2*h[2]*a2*y2 + p3*h[3]*a3 )
  Erate_mn[i] = mean(Erate)
  Erate_lo[i] = quantile(Erate,probs=0.05)
  Erate_hi[i] = quantile(Erate,probs=0.95)
  ErateP1_mn[i] = sumstats[which(vn==paste0("ERexp[",i,",",1,"]")),1]
  ErateP1_lo[i] = sumstats[which(vn==paste0("ERexp[",i,",",1,"]")),5]
  ErateP1_hi[i] = sumstats[which(vn==paste0("ERexp[",i,",",1,"]")),7]
  ErateP2_mn[i] = sumstats[which(vn==paste0("ERexp[",i,",",2,"]")),1]
  ErateP2_lo[i] = sumstats[which(vn==paste0("ERexp[",i,",",2,"]")),5]
  ErateP2_hi[i] = sumstats[which(vn==paste0("ERexp[",i,",",2,"]")),7]  
}
df_PE_plt = data.frame(Year=Yr,Prey=Preytp,Ppn_effort =Ppn_ef)
plt_Peff = ggplot(df_PE_plt, aes(x=Year, y=Ppn_effort, fill=Prey)) + 
  geom_area() + 
  labs(y="Proportion of foraging effort") + theme_classic()
print(plt_Peff) 
df_Erate = data.frame(Year=Years,Erate=Erate_mn,
                      Erate_lo = Erate_lo, Erate_hi = Erate_hi)
plt_Erate = ggplot(df_Erate, aes(x=Year, y=Erate)) + 
  geom_ribbon(aes(ymin=Erate_lo, ymax=Erate_hi),alpha=.3) + 
  geom_line() +
  geom_point(data=df_Fdat[which(df_Fdat$Prey != "Other"),],
             aes(x=Year,y=E_rate,group=Prey,color=Prey)) +
  labs(y="Energy rate foraging") + theme_classic()
print(plt_Erate) 

# Expected dynamics ---------------------------------------
Area = 3.45; reps = 1000; Time = 100
source("pred2prey.R")
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
                  Crabs_10m2 = PPmod[,2]*10, 
                  Clams_10m2 = PPmod[,3]*10, 
                  Otters_km2 = PPmod[,4])
#
dfPP = dfPP %>% pivot_longer(cols = c(Crabs_10m2,Clams_10m2,Otters_km2),
                             names_to = "Species",
                             values_to = "Abundance")
ggplot(data = dfPP, aes(x=Year,y=Abundance,group=Species,color=Species)) +
  geom_line() + theme_classic()
#
OttK_mn = PPmod[Time,4]*Area; OttK_mn

Ottdens_K = numeric(length = reps)
set.seed(123)
rr = sample(Nsims,reps)
# errmat = matrix(0,nrow = Time,ncol = 4); errmat[Time,4] = OttK_mn / Area
for(i in 1:reps){
  pars = list(
    K = as.numeric(mcmc[rr[i],which(startsWith(vn,"K["))]),
    r = as.numeric(mcmc[rr[i],which(startsWith(vn,"r["))]),
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
  if(nrow(rslt)<Time | rslt[nrow(rslt),4] > 50 ){
    Ottdens_K[i] = NA 
  }else{
    Ottdens_K[i] = rslt[Time,4]*Area
  }
}
Ottdens_K = Ottdens_K[which(!is.na(Ottdens_K))]
# hist(Ottdens_K,25)
invisible(capture.output(ft <- fitdist(Ottdens_K,"gamma")))
plot(ft)
a = as.numeric(ft$estimate[1]); b = as.numeric(ft$estimate[2])
Mean_K = a/b
CI_lo = qgamma(.05,a,b)
CI_hi = qgamma(.95,a,b)
noquote(paste0("Mean estimated abundance at K = ",format(Mean_K,digits = 3)," (",
             format(CI_lo,digits = 3)," - ",format(CI_hi,digits = 3),")"))
#
# Save results --------------------------------------------------------------
savename = paste0("elkhorn_PredPrey_rslt_",Sys.Date())
#
fit$save_object(file = paste0(savename,"_fit.RDS"))
# fit = readRDS(paste0(filename))
rm(fit,mod)
save.image(file=paste0(savename,".rdata"))
#
fit = readRDS(paste0(savename,"_fit.RDS"))

ii = which(startsWith(vns,"a["))
ii = c(ii, which(startsWith(vns,"a3")))
ii = c(ii, which(startsWith(vns,"r[")))
ii = c(ii, which(startsWith(vns,"K[")))
ii = c(ii, which(startsWith(vns,"z")))
ii = c(ii, which(startsWith(vns,"phi[")))
ii = c(ii, which(startsWith(vns,"tauY")))
ii = c(ii, which(startsWith(vns,"tauX1")))
ii = c(ii, which(startsWith(vns,"tauX2[2]")))
  
stats_sum = sumstats[ii,c(1,3,4,8,9,10)]
htmlTable(stats_sum)
