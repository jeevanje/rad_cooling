source("~/Dropbox/Rtools/thermo_tools.R")


#===========#
#Parameters #
#===========#

gamma_const = .0065
RH =0.7
Tav   = 250
kappa0 = L*gamma_const/RH/einf/Tav
p0 = 1e4  # Pa
T0 = 260  # K, ref. temp for abs. coeffs
b_kappa = 6500^-1 # (d ln kappa/dk), 1/m^-1

#===========#
# Functions #
#===========#

kappa_100  = 1e1
kappa_1000 = 1e-5
kappa_1500 = 1 
fac = 1
kappa_h20 = function(k){
				if ( (0 <= k) & (k < 1000) ){
					return(fac*kappa_100*exp(log(kappa_1000/kappa_100)*
										(k-100)/(1000-100)))
				} else if ( (1000 <= k) & (k <= 1500) ) {
					return(fac*kappa_1000*exp(log(kappa_1500/kappa_1000)*
							(k-1000)/(1500-1000)))
				} else {stop('k not between 0 and 1500 cm^-1')
					}
			}

Tstar_h20 = function(k){
				if ( (0 <= k) & (k < 1000) ){
					return(500+2000*(k-100)/(1000-100) )
				} else if ( (1000 <= k) & (k <= 1500) ) {
					return(2500-(k-1000)/(1500-1000)*2000)
				} else {stop('k not between 100 and 1500 cm^-1')
					}
			}

temp_scaling = function(tstar,T){
						exp(-tstar*(1/T-1/T0))
				}
