library(fields)
library(ncdf)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/git/src/h2o_data.R")

ncpath = "~/Dropbox/rad_cooling/data/H2O_only_300K/verticalstats_slice.nc"
nc	   = open.ncdf(ncpath)
z      = get.var.ncdf(nc,"z")
nz	   = length(z)
zint   = zinterp(z)
dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
nt	   = length(time)
tabs   = get.var.ncdf(nc,"tabs")
rho    = get.var.ncdf(nc,"rho")
qv     = get.var.ncdf(nc,"qv")
p	   = get.var.ncdf(nc,"p")
lwup   = get.var.ncdf(nc,"lwup")
lwdown = get.var.ncdf(nc,"lwdown")
lwupc  = get.var.ncdf(nc,"lwupc")
lwdownc= get.var.ncdf(nc,"lwdownc")
kmin   = which.min(tabs)
zvec   = 2:kmin
rhov   = qv*rho 
ktp 	= which.min(abs(tabs[1:60]-200))
Ts 	   	= 300
Ttp		= tabs[ktp]
ptp		= p[ktp]
gamma 	= g/Rd*log(Ttp/Ts)/log(ptp/ps)
Tav 	= 1/2*(Ts+200)
wvpinf 	= RH*einf*Tav/gamma/L
 
kvals1 = seq(100,999,length.out=100)     # cm^-1
kvals2 = seq(1001,1500,length.out=100)
kvals  = c(kvals1,kvals2)
nk	   = length(kvals)
kint   = 1/2*( c(kvals[1]-diff(kvals)[1],kvals)+c(kvals,kvals[nk]+diff(kvals)[nk-1]) )
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))
ts_array   = outer(Tstar_vals,tabs,temp_scaling)

# Approx 1 -- parametrize kappa, use temp and p scalings from Humbert (eqn. 3)
dtau1dz = (kappa_vals%o%(rhov*p/p0))*ts_array	#sss
tau1	= array(dim=c(nk,nz+1))	 				#ssi
tau1[,nz+1] <- 0
for (k in nz:1){
		tau1[ ,k] = tau1[,k+1] + dzvec[k]*dtau1dz[ ,k]
} 
weight1 = array(dim=c(nk,nz))					# sss
for (k in  1:length(tabs)){
	weight1[ ,k] = dtau1dz[ ,k]*exp(-0.5*(tau1[ ,k]+tau1[ ,k+1]))
}

# Approx 2 -- approximate tau integral (eqn. 8)
Ts = 300 # K
tau2 = (kappa_vals%o%((exp(-L/Rv/tabs)/kappa0)*ps/p0*(tabs/Ts)^(g/Rd/gamma)))*ts_array
weight2 = array(dim=c(nk,nz))
for (k in  1:nz){
	weight2[ ,k] =  (gamma/tabs[k]^2*(L/Rv + Tstar_vals) + g/Rd/tabs[k])*tau2[ ,k]*exp(-tau2[ ,k])
}

# calculate flux_divi
source   = 1e2*t(pi*outer(tabs,1e2*kvals,planck_k))    # W/m^2/(cm^-1)
for (i in 1:2){
	weight = eval(as.name(paste("weight",i,sep="")))
	assign(paste("flux_div",i,sep=""),diff(kint)%*%(weight*source)) # W/m^3
	}

# analytic formulation
b_kappa = 1.5e-4  # -dTstar/dk, 1/m^-1
b_Tstar = 2.2e-2  # -dlnkappa/dk, K/m^-1
l_k		= 1/(b_kappa + b_Tstar*(1/tabs-1/T0))
#k1 = lk*( log(kappa_h20(0)/kappa0) + log(ps/p0) - g/Rd/gamma*log(Ts/tabs) - L/Rv/tabs ) #m^-1
k1 = l_k*( log(kappa_h20(0)/kappa0) + log(ps/p0) - g/Rd/gamma*log(Ts/tabs) - L/Rv/tabs -Tstar_h20(0)*(1/tabs-1/T0) ) #m^-1
k1[k1 < 0 ] <- 1e-4

flux_div3 = pi*planck_k(tabs,k1)*(g/Rd/tabs + gamma/tabs^2*(L/Rv + Tstar_h20(1e-2*k1)))*l_k   # W/m^3

#DAM output
flux_div_dam   = partialder_i2s(3,z,lwup-lwdown)	# W/m^3
flux_div_damc  = partialder_i2s(3,z,lwupc-lwdownc)


# Begin plot #
pdf1x1("~/Dropbox/rad_cooling/git/figures/ppzf.pdf")
cex=1.25
plot(flux_div_damc[zvec],tabs[zvec],
      ylim = rev(range(tabs)),
	  xlab = expression(partialdiff[z]~F~~"("~W~m^-3~")"),
	  ylab = "Temperature (K)",
	  main = "LW flux divergences, H2O only",
	  cex.lab = cex,
	  cex.axis = cex,
	  cex.main = cex,
	  type="l",
	  lwd=2.5)

Napprox=3
colvec=rev(tim.colors(Napprox))
for (i in 1:Napprox){
	flux_div = eval(as.name(paste("flux_div",i,sep=""))) 
	points(flux_div[zvec],tabs[zvec],type="l",lwd=2,col=colvec[i],lty="solid")
	}
legend("topright",c("RRTM","Eqn. (3)","Eqn. (9)","Eqn. (13)"),
		col=c("black",colvec),cex=1.25,lwd=2)

dev.off()
