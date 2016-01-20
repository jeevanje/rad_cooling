library(fields)
library(ncdf)

source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/my_image_plot.R")
source("~/Dropbox/rad_cooling/src/h2o_data.R")

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

	
kvals1 = seq(100,999,length.out=100)
kvals2 = seq(1001,1500,length.out=100)
kvals  = c(kvals1,kvals2)
nk 	   = length(kvals)
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))
ts_array = outer(Tstar_vals,tabs,temp_scaling)

#========================#
# Approximations for tau #
#========================#

Napprox=2

# Approx 1 -- parametrize kappa, use temp and p scalings from Humbert
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

# Approx 2 -- scale p at z only, parametrize path, constant gamma=6.5 K/km
Ts = 300 # K
tau2 = kappa_vals%o%((exp(-L/Rv/tabs)/kappa0)*ps/p0*(tabs/Ts)^(g/Rd/gamma))
weight2 = array(dim=c(nk,nz))
for (k in  1:nz){
	weight2[ ,k] =  (L*gamma/Rv/tabs[k]^2+g/Rd/tabs[k])*tau2[ ,k]*exp(-tau2[ ,k])
}

source = 1e3*t(pi*outer(tabs,1e2*kvals,planck_k))  # W/m^2/(10 cm^-1)
for (n in 1:Napprox){
	weight = eval(as.name(paste("weight",n,sep="")))
	flux_div = weight*source		# W/m^3/(10 cm^-1)
	cooling = array(dim=dim(flux_div))
	for (i in 1:nk){
		cooling[i,]  = flux_div[i,]/rho/Cp*86400  # K/day/(10 cm^-1) 
		}
	assign(paste("flux_div",n,sep=""),flux_div) 
	assign(paste("cooling",n,sep=""),cooling)   

	}

cex=1.4
plot_it = function(field,main="",zlim=range(field)){
			my.image.plot(kvals,tabs[kmin:1],field[ ,kmin:1],
			  ylim = rev(range(tabs[1:kmin])),
			  zlim = zlim,
			  xlab = "wavenumber (cm^-1)",
			  ylab = "Temperature (K)",
			  main = main,
			  cex.lab = cex,
			  cex.axis = cex,
			  cex.legend = cex,
			  cex.main = cex)
			}

zlim_cool = c(0,0.15)

#======#
# Plot #
#======#
pdf1x2("~/Dropbox/rad_cooling/git/figures/cooling.pdf")
plot_it(cooling1,"Cooling rate ( K/day/(10 cm^-1) )",zlim=zlim_cool)

plot_it(cooling2,"",zlim=zlim_cool)
dev.off()