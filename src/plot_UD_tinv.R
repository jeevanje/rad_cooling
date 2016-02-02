library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/src/h2o_data.R")

SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"

#===========#
# Get data  #
#===========#

for (i in 1:N){
	SST  = SSTlist[i]
	ncpath = paste(datadir,"/zeroO3_",SST,"K/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	z    = get.var.ncdf(nc,"z")
	nz   = length(z)
	zmax = 25e3
	kmax = which.min(abs(zmax-z))
	zvec = 3:kmax
	time = get.var.ncdf(nc,"time")
	nt = length(time)
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)	
	p = apply(get.var.ncdf(nc,start=c(1,nt-10),"p"),1,mean)	
	U = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwup"),1,mean)	
	D = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwdown"),1,mean)	
	assign(paste("tabs",SST,sep=""),tabs)
	assign(paste("p",SST,sep=""),p)
	assign(paste("U",SST,sep=""),U)	
	assign(paste("D",SST,sep=""),D)	
	}



pdf1x3("~/Dropbox/rad_cooling/git/figures/UD_tinv.pdf")
cex=1.75
lty="solid"
z=z[zvec]
xlim=c(0,525)
colvec = tim.colors(length(SSTlist))
tabslim = c(310,170)

#===========#
# height    #
#===========#
plot(1,type="n",
	xlab=expression("U, D  ("~W/m^2~")"), 
	ylab = "z (km)",main="LW Fluxes",
	xlim=xlim, ylim=1e-3*range(z),
	cex.axis=cex, cex.main=cex, cex.lab=cex)
for (i in 1:N){
	Ts  = SSTlist[i]
	U =  eval(as.name(paste("U",Ts,sep="")))[zvec]
	D =  eval(as.name(paste("D",Ts,sep="")))[zvec]
	points(U,1e-3*z,type="l",lwd=2, col = colvec[i],lty=lty)
	points(D,1e-3*z,type="l",lwd=2, col = colvec[i],lty=lty)
		}


#===========#
# pressure  #
#===========#
plot(1,type="n",
	xlab=expression("U, D  ("~W/m^2~")"), 
	ylab = "Pressure (hPa)",
	main="LW Fluxes",
	xlim=xlim, 			
	ylim=rev(range(1e-2*p)),
	cex.axis=cex, cex.main=cex, cex.lab=cex)
for (i in 1:N){
	Ts  = SSTlist[i]
	U =  eval(as.name(paste("U",Ts,sep="")))[zvec]
	D =  eval(as.name(paste("D",Ts,sep="")))[zvec]
	p =  eval(as.name(paste("p",Ts,sep="")))[zvec]
	points(U,1e-2*p,type="l",lwd=2, col = colvec[i],lty=lty)
	points(D,1e-2*p,type="l",lwd=2, col = colvec[i],lty=lty)
		}


#==============#
# temperature  #
#==============#

plot(1,type="n",
	xlab=expression("U, D  ("~W/m^2~")"), 
	ylab = "Temperature (K)",
	main="LW Fluxes",
	xlim=xlim, 			
	ylim=tabslim,
	cex.axis=cex, cex.main=cex, cex.lab=cex)
for (i in 1:N){
	Ts  = SSTlist[i]
	U =  eval(as.name(paste("U",Ts,sep="")))[zvec]
	D =  eval(as.name(paste("D",Ts,sep="")))[zvec]
	tabs =  eval(as.name(paste("tabs",Ts,sep="")))[zvec]
	points(U,tabs,type="l",lwd=2, col = colvec[i],lty=lty)
	points(D,tabs,type="l",lwd=2, col = colvec[i],lty=lty)
		}
			
dev.off()
