library(fields)
library(ncdf)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/src/h2o_data.R")

ncpath = "~/Dropbox/rad_cooling/data/H2O_only_300K/verticalstats_slice.nc"
nc	   = open.ncdf(ncpath)
z      = get.var.ncdf(nc,"z")
zint   = zinterp(z)
nz	   = length(z)
tabs   = get.var.ncdf(nc,"tabs")
rho    = get.var.ncdf(nc,"rho")
rh     = get.var.ncdf(nc,"rh")
p	   = get.var.ncdf(nc,"p")
qv	   = get.var.ncdf(nc,"qv")
kmin   = which.min(tabs)
zvec   = 1:kmin
dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
rhov   = qv*rho

# Compute WVP
wvp_dam	 = numeric(nz)	 				#ssi
wvp_est	 = numeric(nz)	 				#ssi
wvp_dam[nz] = rhov[nz]*dzvec[nz]
for (k in (nz-1):1){
		wvp_dam[k] = wvp_dam[k+1] + dzvec[k]*rhov[k]
} 
wvp_est = 1/kappa0*exp(-L/Rv/tabs)


cex.lab=1.2
cex.axis=0.9
cex.main=1.2

# Begin plot
pdf1x2("~/Dropbox/rad_cooling/git/figures/wvp_check.pdf")
plot(wvp_dam[zvec],tabs[zvec],
      ylim = rev(range(tabs)),
	  xlab = expression("WVP  ("~kg~m^-2~")"),
	  ylab = "Temperature (K)",
	  main = "Water vapor path",
	  cex.lab = cex.lab,
	  cex.axis = cex.axis,
	  cex.main = cex.main,
	  type="l",
	  lwd=2)
			
points(wvp_est[zvec],tabs[zvec],type="l",lwd=2,lty="dashed")

plot(wvp_dam[zvec],tabs[zvec],
      ylim = rev(range(tabs)),
	  xlab = expression("WVP  ("~kg~m^-2~")"),
	  ylab = "Temperature (K)",
	  main = "Water vapor path",
	  cex.lab = cex.lab,
	  cex.axis = cex.axis,
	  cex.main = cex.main,
	  type="l",
	  log="x",
	  lwd=2)

points(wvp_est[zvec],tabs[zvec],type="l",lwd=2,lty="dashed")

dev.off()
