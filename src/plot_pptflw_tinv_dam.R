library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/src/h2o_data.R")

SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
colvec = tim.colors(length(SSTlist))
tabslim = c(310,170)

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
	lwup = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwup"),1,mean)	
	lwdown = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwdown"),1,mean)	
	lapse  = -partialder_i2s(3,z,s2i(3,z,tabs))
	ppzf =  partialder_i2s(3,z,lwup-lwdown)  # W/m^3
	assign(paste("tabs",SST,sep=""),tabs)
	assign(paste("lapse",SST,sep=""),lapse)
	assign(paste("p",SST,sep=""),p)
	assign(paste("pptf",SST,sep=""),ppzf/lapse)	
	}



pdf(file="~/Dropbox/rad_cooling/git/figures/pptflw_tinv_dam.pdf",width=12,height=5)
par(mfrow=c(1,3),mar=c(5,6,5,3))
cex=2
xlim=c(0,5)
lty="solid"
z=z[zvec]
main = "LW flux divergence"
xlab=expression(-partialdiff[T]*F~~"("~W/m^2/K~")") 

#===========#
# height    #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	pptf =  eval(as.name(paste("pptf",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,1e-3*z,
			xlab=xlab,
			ylab = "z (km)",
			main=main,
			xlim=xlim,
		 	col=colvec[i],lty=lty,
		 	cex.axis=cex,
		 	cex.main=cex,
		 	cex.lab=cex,
		 	type="l",
		 	lwd=2
			 )
		} else {
		points(pptf,1e-3*z,type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}


#===========#
# pressure  #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	pptf =  eval(as.name(paste("pptf",Ts,sep="")))[zvec]
	p =  eval(as.name(paste("p",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,1e-2*p,
			xlab=xlab,
			ylab = "Pressure (hPa)",
			main=main,
			xlim=xlim,
			ylim=rev(range(1e-2*p)),
		 	col=colvec[i],lty=lty,
		 	cex.axis=cex,
		 	cex.main=cex,
		 	cex.lab=cex,
		 	type="l",
		 	lwd=2
			 )
		} else {
		points(pptf,1e-2*p,type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}

#==============#
# temperature  #
#==============#

for (i in 1:N){
	Ts  = SSTlist[i]
	pptf =  eval(as.name(paste("pptf",Ts,sep="")))[zvec]
	tabs =  eval(as.name(paste("tabs",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,tabs,
			xlab=xlab,
			main=main,
			 ylim = tabslim, xlim=xlim,
			 ylab = "Temperature (K)",
			 col=colvec[i],type="l", lty=lty,
			 cex.main=cex,
			 cex.axis=cex,
			 cex.lab=cex
		 )
		} else {
		points(pptf,tabs,type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}
			
dev.off()
