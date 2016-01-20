
#===========#
#Parameters #
#===========#

RH =0.7
Tav   = 250
gamma = .0065
kappa0 = L*gamma/RH/einf/Tav
p0 = 1e4  # Pa


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
				if ( (100 <= k) & (k < 1000) ){
					return(500+(k-100)/(1000-100)*2000)
				} else if ( (1000 <= k) & (k <= 1500) ) {
					return(2500-(k-1000)/(1500-1000)*2000)
				} else {stop('k not between 100 and 1500 cm^-1')
					}
			}

temp_scaling = function(tstar,T){
						T0 = 260  # K, ref. temp for abs. coeffs
						exp(-tstar*(1/T-1/T0))
				}
