# Predator with 2 prey ODE function  (assumes nuisance 3rd prey)
pred2prey <-  function (t, y, p){
  X1 <- y[1] ;
  X2 <- y[2] ;
  Y <- y[3] ;
  #
  K = p$K ;
  r = p$r ;
  phi = p$phi ;
  z1 = p$z1 ;
  z2 = p$z2 ;
  a3 = p$a3 ;
  a =  p$a ;
  h = p$h ;
  g = p$g ;
  Tcf = p$Tcf ;
  #
  p1 = (a[1]*y[1]) / (a[1]*y[1] + phi[1]*a[2]*y[2] + phi[2]*a3) ;
  p2 = (phi[1]*a[2]*y[2]) / (a[1]*y[1] + phi[1]*a[2]*y[2] + phi[2]*a3) ;  
  p3 = (phi[2]*a3) / (a[1]*y[1] + phi[1]*a[2]*y[2] + phi[2]*a3) ;
  #
  dX1.dt <- ( (r[1] * (1 - y[1] / K[1])) - 
                ( p1*a[1]*y[3]*Tcf / (1 + p1*h[1]*a[1]*y[1] + p2*h[2]*a[2]*y[2] + p3*h[3]*a3 ))) * y[1] ;
  dX2.dt <- ( (r[2] * (1 - y[2] / K[2])) - 
                ( p2*a[2]*y[3]*Tcf / (1 + p1*h[1]*a[1]*y[1] + p2*h[2]*a[2]*y[2] + p3*h[3]*a3 ))) * y[2] ;
  dY.dt <-  (z1 * ((g[1]*p1*a[1]*y[1] + g[2]*p2*a[2]*y[2] + g[3]*p3*a3) / 
                     (1 + p1*h[1]*a[1]*y[1] + p2*h[2]*a[2]*y[2] + p3*h[3]*a3 ) - z2)) * y[3] ;
  #  
  return(list(c(dX1.dt, dX2.dt, dY.dt )))
}
