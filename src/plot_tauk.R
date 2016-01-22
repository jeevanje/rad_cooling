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
indav	= which.min(abs(tabs[1:60]-Tav))
wvpinf 	= RH*einf*Tav/gamma/L
 
kvals1 = seq(100,999,length.out=100)     # cm^-1
kvals2 = seq(1001,1500,length.out=100)
kvals  = c(kvals1,kvals2)
nk	   = length(kvals)
kint   = 1/2*( c(kvals[1]-diff(kvals)[1],kvals)+c(kvals,kvals[nk]+diff(kvals)[nk-1]) )
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))
ts_array   = outer(Tstar_vals,tabs,temp_scaling)

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

# Approx 1 -- approximate tau integral, p-scaling at z only
tau1 = ((kappa_vals)%o%(ps/p0*(tabs/Ts)^(g/Rd/gamma)*wvpinf*exp(-L/Rv/tabs)))*ts_array


#============#
# Begin plot #
#============#
cex=1.25
logtaulim = c(-15,8)
kindvec   = c(which.min(abs(kvals-100)),which.min(abs(kvals-500)),which.min(abs(kvals-900)))
Nk		  = length(kindvec)
colvec 	  = tim.colors(Nk)
pdf1x1("~/Dropbox/rad_cooling/git/figures/tauk.pdf")
for (n in 1:Nk){
	kind = kindvec[n]
	if (n ==1){
		plot(log(tau[kind,zvec]),tabs[zvec],
		      ylim = rev(range(tabs)),
			  xlim = logtaulim,
			  xlab = expression(log~tau[k]),
			  ylab = "Temperature (K)",
			  main = "Optical depth",
			  cex.lab = cex,
			  cex.axis = cex,
			  cex.main = cex,
			  type="l",
			  lwd=2,
			  col=colvec[n])
		} else {
		points(log(tau[kind,zvec]),tabs[zvec],type="l",lwd=2,col=colvec[n],lty="solid")
		}
	points(log(tau1[kind,zvec]),tabs[zvec],type="l",lwd=2,col=colvec[n],lty="dashed")
	#legend("topright",c("Eqn. (3)","Eqn. (9)"),
	#		lty==c("black",colvec),cex=1.25,lwd=2)
	}
dev.off()
