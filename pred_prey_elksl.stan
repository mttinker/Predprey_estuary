//
//  This Stan program defines a simple predator-prey model
//  with one predator (Y) and 3 prey species (X1, X2, X3) 
//  where X1 and X2 are primary prey, X3 = nuisance "all other prey""
// - predator population dynamics scaled to standard area of 1 km2
// - prey population dynamics scaled to standard area of 1 m2
// - time units of foraging dynamics = 1 minute of feeding bout,
//    (net consumption then re-scaled to time step interval of 1 year)
// - for predator, energy intake rate (kcal/min) converts to instantaneous 
//    growth rate (r) by parameters z[1] and z[2]
//
functions {
  vector dy_dt(real t,       // time
               vector y,     // system state {prey1, prey2, predator}
               vector K,     // parameter K
               vector r,     // parameter K
               vector phi,   // parameter phi
               real z1,      // prey conversion param 1
               real z2,      // prey conversion param 1
               real a3,      // encounter rate base, alternative prey
               vector a,     // encounter rate params
               vector h,     // handling time (constants)
               vector g,     // energy densities (constants)
               real Tcf) {   // time conversion (minutes to year)
    real p1 = (a[1]*y[1]) / (a[1]*y[1] + phi[1]*a[2]*y[2] + phi[2]*a3) ;
    real p2 = (phi[1]*a[2]*y[2]) / (a[1]*y[1] + phi[1]*a[2]*y[2] + phi[2]*a3) ;
    real p3 = (phi[2]*a3) / (a[1]*y[1] + phi[1]*a[2]*y[2] + phi[2]*a3) ;
    vector[3] dydt;
    //
    dydt[1] = ( (r[1] * (1 - y[1] / K[1])) - 
      ( p1*a[1]*y[3]*Tcf / (1 + p1*h[1]*a[1]*y[1] + p2*h[2]*a[2]*y[2] + p3*h[3]*a3 ))) * y[1] ;
    dydt[2] = ( (r[2] * (1 - y[2] / K[2])) - 
      ( p2*a[2]*y[3]*Tcf / (1 + p1*h[1]*a[1]*y[1] + p2*h[2]*a[2]*y[2] + p3*h[3]*a3 ))) * y[2] ;
    dydt[3] = (z1 * ((g[1]*p1*a[1]*y[1] + g[2]*p2*a[2]*y[2] + g[3]*p3*a3) / 
      (1 + p1*h[1]*a[1]*y[1] + p2*h[2]*a[2]*y[2] + p3*h[3]*a3 ) - z2)) * y[3] ;
    //  
    return dydt;
  }
}
//
// data accepted by model
data {
  int<lower=0> N;         // N years of population dynamics
  int<lower=0> Nx1;       // N data observations for X1
  int<lower=1> ix1[Nx1];  // index of years of observation data for X1
  vector[Nx1] ns1;        // number sample units for x1
  int<lower=0> X1obs[Nx1];// observed data for population X1 (counts)
  int<lower=0> Nx2;       // N data observations for X2
  int<lower=1> ix2[Nx2];  // index of years of observation data for X2
  vector[Nx2] ns2;        // number sample units for x2  
  int<lower=0> Nx2p;      // N time periods for X2 data
  int<lower=1> ix2p[Nx2];// index of period for observation data for X2
  int<lower=0> X2obs[Nx2];// observed data for population X2 (counts)
  int<lower=0> Ny;        // N years of survey data for Y
  int<lower=1> ixy[Ny];   // index of years of survey data
  int<lower=0> Yobs[Ny]; // observed data for population Y (log abundance)
  int<lower=0> Np;        // N obs of foraging data estimates for Y
  int<lower=1> ixp[Np];   // index of years for foraging effort by prey
  simplex[3] Pobs[Np] ;   // proportional allocation of feeding effort {p1, p2, p3}
  int<lower=0> Nf;        // N obs of foraging energy rate estimates for Y
  int<lower=1> ixf[Nf];   // index of years for feeding energy rate estimates
  int<lower=1> ixfsp[Nf]; // index of prey species for feeding energy rate estimates
  real<lower=0> ER[Nf];   // energy recovery rate estimates, by prey
  vector[3] g ;           // energy per item of prey type X
  vector[3] h ;           // handling time per item, prey type X
  real<lower=0> Tcf ;     // time correct factor (CF), min to year/otter per m2 
  real<lower=0> Acf ;     // area CF, 1 km2 to total study area
  vector<lower=0>[2]omega ;// measurement transformation factor, prey X1 and X2
  real<lower=0> ts[N] ;   // timesteps ODE to save (1:N if yearly time steps)
  real<lower=0> tauP[Np]; // precision param, proportional effort allocation
  real<lower=0> tauF[Nf]; // inv. scale param, energy intake rates 
  real<lower=0> Y_init;   // initial population size (per km2) for y
  // COMMENT NEXT LINE (use for coding to avoid error message from ODE fxn)
  // vector<lower=0>[3] y[N] ;
  //
}
// The parameters accepted by the model. 
parameters {
  real<lower=.5> tauX1;      // inv. scale param, observed data for X1
  real<lower=.5> tauX2[Nx2p];// inv. scale params, observed data for X2
  real<lower=0.1,upper=10> tauY;// observer error (precision param) surveys of Y
  real<lower=0.01,upper=.025> z1 ; // param 1 for converting energy to growth rate of Y  
  real<lower=3,upper=7> z2 ; // param 1 for converting energy to growth rate of Y
  vector<lower=0>[2] a0 ;  // per-capita encounter rate for prey type X (when X ~ K)
  real a3 ;               // encounter rate for "alternative" prey
  vector<lower=0.1,upper=2>[2] r ;   // recruitment rates (annual, per HA) for X1 and X2
  vector<lower=0>[2] K ;  // carrying capacity for X1, X2 scaled to sample unit
  vector<lower=0,upper=10>[2] phi ;// prey preference params (X2 and X3 vs. X1)
}
// Derived parameters
transformed parameters {
  vector[2] a = a0 ./ K ;            // per-capita max encounter rate for prey type X 
  vector[3] y0 = append_row(K, Y_init) ; // initial abundances for ODE fxn 
  vector[N] Yexp ;        // estimated predator abundance, entire area
  vector[N] X1exp ;       // estimated prey X1 abundance per sample unit
  vector[N] X2exp ;       // estimated prey X2 abundance per sample unit
  matrix[N,3] Pexp ;      // proportional allocation of foraging effort by prey 
  matrix[N,3] ERexp ;     // energy recovery rate for each prey while feeding
  vector[3] y[N] = ode_rk45(dy_dt, y0, 0, ts,
                    K, r, phi, z1, z2, a3, a, h, g, Tcf) ;
  //
  for (t in 1:N){
    X1exp[t] = y[t][1] * omega[1] ;
    X2exp[t] = y[t][2] * omega[2] ;
    Yexp[t] = y[t][3] * Acf ;
    // Calulate expected proportional allocation of effort to each prey type
    Pexp[t,1] = (a[1]*y[t][1]) / 
       (a[1]*y[t][1] + phi[1]*a[2]*y[t][2] + phi[2]*a3 ) ;
    Pexp[t,2] = (phi[1]*a[2]*y[t][2]) / 
       (a[1]*y[t][1] + phi[1]*a[2]*y[t][2] + phi[2]*a3 ) ;
    Pexp[t,3] = (phi[2]*a3) / 
       (a[1]*y[t][1] + phi[1]*a[2]*y[t][2] + phi[2]*a3 ) ;
    // Calculate expected energy recovery rates (kcal/min) by species      
    ERexp[t,1] = (g[1]*a[1]*y[t][1]) / (1 + h[1]*a[1]*y[t][1]) ;
    ERexp[t,2] = (g[2]*a[2]*y[t][2]) / (1 + h[2]*a[2]*y[t][2]) ;
    ERexp[t,3] = (g[3]*a3) / (1 + h[3]*a3) ;   
  }
}
// The model to be estimated. 
model {
  // Observed data, prey X1, X2 and predator Y
  X1obs ~ neg_binomial_2(X1exp[ix1] .* ns1, tauX1) ;
  X2obs ~ neg_binomial_2(X2exp[ix2] .* ns2, tauX2[ix2p]) ;
  Yobs ~ neg_binomial_2(Yexp[ixy], tauY) ;
  //Yobs ~ gamma(Yexp[ixy] * tauY, tauY) ;
  // Foraging effort (p1, p2, p3)
  for(i in 1:Np){
    Pobs[i] ~ dirichlet(to_vector(Pexp[ixf[i]]) * tauP[i]) ;
  }
  // Consumption rates by year, species
  // MAKE ERexp a matrix!!
  for(i in 1:Nf){
    ER[i] ~ gamma(ERexp[ixf[i], ixfsp[i]] * tauF[i], tauF[i]) ;
  }
  // Priors for parameters
  // omega ~ cauchy(0,1) ;
  r ~ normal(0,1) ; 
  K ~ cauchy(0,1) ;  
  z1 ~ normal(0, .05) ;
  z2 ~ normal(5,1) ;
  phi[1] ~ normal(0,2) ;
  phi[2] ~ normal(0,.1) ;
  a0 ~ normal(0,1) ;
  a3 ~ normal(0,2.5) ;
  tauX1 ~ normal(0,10) ;
  tauX2 ~ normal(0,10) ;  
  tauY ~ normal(0,2.5) ;
}
