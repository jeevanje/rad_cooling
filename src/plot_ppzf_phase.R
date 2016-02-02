library(fields)
library(ncdf)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/my_image_plot.R")
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
 
kvals1 = seq(100,999,length.out=100)     # cm^-1
kvals2 = seq(1001,1500,length.out=100)
kvals  = c(kvals1,kvals2)
nk	   = length(kvals)
b_kappa = 1.5e-4  # -dTstar/dk, 1/m^-1
b_Tstar = 2.2e-2  # -dlnkappa/dk, K/m^-1
l_k		= 1/(b_kappa + b_Tstar*(1/tabs-1/T0))
Ts = 300
k1 = l_k*( log(kappa_h20(0)/kappa0) + log(ps/p0) - g/Rd/gamma*log(Ts/tabs) - L/Rv/tabs -Tstar_h20(0)*(1/tabs-1/T0) ) #m^-1
k1[k1 < 0 ] <- 1e-4
ppzf = array(dim=c(nk,nz))
for (i in 1:nk){
		ppzf[i,] = pi*planck_k(tabs,1e2*kvals[i])*(g/Rd/tabs + gamma/tabs^2*(L/Rv + Tstar_h20(1e-2*k1)))*l_k   # W/m^3
		}
	
# Begin plot #
pdf1x1("~/Dropbox/rad_cooling/git/figures/ppzf_phase.pdf")
cex=1.25
my.image.plot(kvals,tabs[kmin:1], ppzf[ ,kmin:1],zlim=range(ppzf[ ,kmin:1]),
      ylim = rev(range(tabs)),
	  main = expression(partialdiff[z]~F~~"("~W~m^-3~")"),
	  ylab = "Temperature (K)",
	  xlab = expression("Wavenumber"~~"("~cm^-1~")"),
	  cex.lab = cex,
	  cex.axis = cex,
	  cex.main = 1.5,
	  cex.legend = cex)
points(1e-2*k1[zvec],tabs[zvec],type="l",lwd=3,col="black",lty="solid")
text(300,260,expression(k[1]*"("*T*")"),cex=1.25)
text(540,243,expression(k[alt]*"("*T*")"),cex=1,col="black")
abline(a=211.5,b=0.075,lty="dashed",lwd=2,col="black") 
dev.off()
