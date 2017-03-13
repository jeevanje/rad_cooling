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
	ncpath = paste(datadir,"/prod_data_4_22_16/",SST,"k/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	z    = get.var.ncdf(nc,"z")
	nz   = length(z)
	time = get.var.ncdf(nc,"time")
	nt = length(time)
	nt_avg = 20
	start = c(1,nt-nt_avg)
	tabs = apply(get.var.ncdf(nc,start=start,"tabs"),1,mean)	
	p = apply(get.var.ncdf(nc,start=start,"p"),1,mean)	
	cloud = apply(get.var.ncdf(nc,start=start,"cloud"),1,mean)	
	klcl = which.max(cloud[1:15])
	zmax = 22.5e3
	kmax = which.min(abs(zmax-z))
	zvec = klcl:kmax
#	zvec = 1:kmax
	lwup = apply(get.var.ncdf(nc,start=start,"lwup"),1,mean)	
	lwdown = apply(get.var.ncdf(nc,start=start,"lwdown"),1,mean)	
	F = lwup - lwdown
	lapse  = -partialder_i2s(3,z,s2i(3,z,tabs))
	ppzf =  partialder_i2s(3,z,F)  # W/m^3
	assign(paste("klcl",SST,sep=""),klcl)
	assign(paste("zvec",SST,sep=""),zvec)
	assign(paste("tabs",SST,sep=""),tabs)
	assign(paste("lapse",SST,sep=""),lapse)
	assign(paste("p",SST,sep=""),p)
	assign(paste("F",SST,sep=""),F)	
	assign(paste("pptf",SST,sep=""),ppzf/lapse)	
	}



colvec = tim.colors(length(SSTlist))
cex=2
xlim=c(0,5)
tabslim = c(max(SSTlist),170)
lty="solid"
main = "LW flux divergence"
xlab=expression(-partialdiff[T]*F[LW]~~"("~W/m^2/K~")") 

pdf(file="~/Dropbox/rad_cooling/git/figures/pptflw_tinv_dam.pdf",width=12,height=5)
par(mfrow=c(1,3),mar=c(5,6,5,3))

#===========#
# height    #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	col = colvec[i]
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	pptf =  eval(as.name(paste("pptf",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,1e-3*z[zvec],
			xlab=xlab,
			ylab = "z (km)",
			main=main,
			xlim=xlim,
			ylim = c(0,1e-3*zmax),
		 	col=colvec[i],lty=lty,
		 	cex.axis=cex,
		 	cex.main=cex,
		 	cex.lab=cex,
		 	type="l",
		 	lwd=2
			 )
		} else {
		points(pptf,1e-3*z[zvec],type="l",lwd=2, col = col,lty=lty)
		}
	legend("topright",legend=SSTlist,lwd=2,cex=1.5,col=colvec,title=expression(T[s] ~~ "(K)"))
	}


#===========#
# pressure  #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	col = colvec[i]
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	pptf =  eval(as.name(paste("pptf",Ts,sep="")))[zvec]
	p =  eval(as.name(paste("p",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,1e-2*p,
			xlab=xlab,
			ylab = "Pressure (hPa)",
			main=main,
			xlim=xlim,
			ylim=c(1000,0),
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
	col = colvec[i]
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	pptf =  eval(as.name(paste("pptf",Ts,sep="")))[zvec]
	tabs =  eval(as.name(paste("tabs",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,tabs,
			xlab=xlab,
			main=main,
			 ylim = tabslim, xlim=xlim,
			 ylab = "Temperature (K)",
			 col=colvec[i],type="l", lty=lty,lwd=2,
			 cex.main=cex,
			 cex.axis=cex,
			 cex.lab=cex
		 )
		} else {
		points(pptf,tabs,type="l",lwd=2, col = colvec[i],lty=lty)
		}
	abline(v=0,col="gray",lty = "dashed",lwd=1.5)
	}
			
dev.off()
