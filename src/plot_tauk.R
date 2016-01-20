library(fields)
library(ncdf)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/git/src/h2o_data.R")

ncpath = "~/Dropbox/rad_cooling/data/H2O_only_300K/verticalstats_slice.nc"
nc	   = open.ncdf(ncpath)
z      = get.var.ncdf(nc,"z")
zint   = zinterp(z)
nz	   = length(z)
dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
nt	   = length(time)
tabs   = get.var.ncdf(nc,"tabs")
rho    = get.var.ncdf(nc,"rho")
qv     = get.var.ncdf(nc,"qv")
p	   = get.var.ncdf(nc,"p")
kmin   = which.min(tabs)
zvec   = 1:kmin
rhov   = qv*rho 
ktp 	= which.min(abs(tabs[1:60]-200))
Ts 	   	= 300
Ttp		= tabs[ktp]
ptp		= p[ktp]
gamma 	= g/Rd*log(Ttp/Ts)/log(ptp/ps)
Tav 	= 1/2*(Ts+200)
H 		= (g/Rd/Tav)^-1
xav 	= L/Rv/Tav
p_est	= ps*exp(-1/H*(Ts/gamma - Tav/gamma*(1-1/xav*(L/Rv/tabs - xav))))
indav	= which.min(abs(tabs[1:60]-Tav))
p_av	= p_est[indav]
wvpinf 	= RH*einf/gamma/L*Tav
 
kvals1 = seq(100,999,length.out=100)     # cm^-1
kvals2 = seq(1001,1500,length.out=100)
kvals  = c(kvals1,kvals2)
nk	   = length(kvals)
kint   = 1/2*( c(kvals[1]-diff(kvals)[1],kvals)+c(kvals,kvals[nk]+diff(kvals)[nk-1]) )
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))
ts_array   = outer(Tstar_vals,tabs,temp_scaling)
alpha_k	   = ( 1 + Rv*Tstar_vals/L + Rv*Tav^2/L/gamma/H )^-1


# Ground truth -- parametrize kappa, use temp and p scalings from Humbert (eqn. 3)
dtaudz = (kappa_vals%o%(rhov*p/p0))*ts_array	#sss
tau	= array(dim=c(nk,nz+1))	 				#ssi
tau[,nz+1] <- 0
for (k in nz:1){
		tau[ ,k] = tau[,k+1] + dzvec[k]*dtaudz[ ,k]
} 
weight = array(dim=c(nk,nz))					# sss
for (k in  1:length(tabs)){
	weight[ ,k] = dtaudz[ ,k]*exp(-0.5*(tau[ ,k]+tau[ ,k+1]))
}

# Approx 1 -- approximate tau integral (eqn. 8)
tau1 = ((alpha_k*kappa_vals)%o%(p_av/p0*exp(-Tav^2/gamma/H*(1/tabs-1/Tav))*wvpinf*exp(-L/Rv/tabs)))*ts_array
weight1 = array(dim=c(nk,nz))
for (k in  1:nz){
	weight1[ ,k] =  (L*gamma/Rv/tabs[k]^2+g/Rd/tabs[k])*tau1[ ,k]*exp(-tau1[ ,k])
}

# Approx 2 -- approximate tau integral, p-scaling at z only
tau2 = ((kappa_vals)%o%(ps/p0*(tabs/Ts)^(g/Rd/gamma)*wvpinf*exp(-L/Rv/tabs)))*ts_array
weight2 = array(dim=c(nk,nz))
for (k in  1:nz){
	weight2[ ,k] =  (L*gamma/Rv/tabs[k]^2+g/Rd/tabs[k])*tau2[ ,k]*exp(-tau2[ ,k])
}

# Begin plot #
cex=1.25
logtaulim = c(-20,10)
kindvec = c(which.min(abs(kvals-200)),which.min(abs(kvals-400)),which.min(abs(kvals-600)))
pdf1x3("~/Dropbox/rad_cooling/git/figures/tauk.pdf")
for (k in kindvec){
	plot(log(tau[k,zvec]),tabs[zvec],
	      ylim = rev(range(tabs)),
		  #xlim = logtaulim,
		  xlab = expression(tau[k]),
		  ylab = "Temperature (K)",
		  main = "Optical depth",
		  cex.lab = cex,
		  cex.axis = cex,
		  cex.main = cex,
		  type="l",
		  lwd=2)
	Napprox=2
	colvec=rev(tim.colors(Napprox))
	for (i in 1:Napprox){
		tau_approx = eval(as.name(paste("tau",i,sep=""))) 
		points(log(tau_approx[k,zvec]),tabs[zvec],type="l",lwd=2,col=colvec[i],lty="solid")
		}
	legend("topright",c("Eqn. (3)","Eqn. (8)","z-scaling"),
			col=c("black",colvec),cex=1.25,lwd=2)
	}
dev.off()
