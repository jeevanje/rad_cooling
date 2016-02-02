library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")

# Get data
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
dqdts_vec = numeric(N)
Qvec	= numeric(N)
for (i in 1:N){
	SST  = SSTlist[i]
	ncpath = paste(datadir,"/zeroO3_",SST,"K/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	time = get.var.ncdf(nc,"time")
	z = get.var.ncdf(nc,"z")
	nz = length(z)
	nt = length(time)
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)	
	cloud = apply(get.var.ncdf(nc,start=c(1,nt-10),"cloud"),1,mean)	
	kanvil = which.max(cloud)
#	klcl = which.max(cloud[1:15])
	klcl = which.min(abs(z-1000))
	lapse = - partialder_i2s(3,z,s2i(3,z,tabs))
	p = apply(get.var.ncdf(nc,start=c(1,nt-10),"p"),1,mean)	
    ktp = which.min(abs(tabs-200))
    Ttp=tabs[ktp]
    ptp=p[ktp]
    gamma = g/Rd*log(Ttp/SST)/log(ptp/ps)
	kstar = which.min(abs(lapse[klcl:kanvil]-gamma)) + klcl -1
	lwup = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwup"),1,mean)	
	lwdown = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwdown"),1,mean)	
	F = lwup-lwdown
	ppzf = partialder_i2s(3,z,lwup-lwdown)
	Q = F[nz]-F[klcl]
	Qvec[i] = Q
	dqdts_vec[i] = ppzf[klcl]/lapse[klcl] + Q/gamma*1e-4 
	assign(paste("lapse",SST,sep=""),lapse)
	assign(paste("klcl",SST,sep=""),klcl)
	assign(paste("kstar",SST,sep=""),kstar)
	assign(paste("F",SST,sep=""),F)
	assign(paste("ppzf",SST,sep=""),ppzf)
	assign(paste("tabs",SST,sep=""),tabs)	
	assign(paste("p",SST,sep=""),p)	
	assign(paste("gamma",SST,sep=""),gamma)		
	assign(paste("Q",SST,sep=""),Q)			
	}


