library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")

# Get data
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
dqdtsLW_vec = numeric(N)
dqdtsSW_vec = numeric(N)
dqdtsNet_vec = numeric(N)
QSWvec	= numeric(N)
QLWvec	= numeric(N)
QNetvec	= numeric(N)
precipvec	= numeric(N)
SHFvec	= numeric(N)
zlcl = 1000
for (i in 1:N){
	SST  = SSTlist[i]
	# globalstats
	globalstatspath = paste(datadir,"/zeroO3_",SST,"K/data/globalstats.nc",sep="")
	nc_gs = open.ncdf(globalstatspath)
	time_gs = get.var.ncdf(nc_gs,"time")
	nt_gs = length(time_gs)
	precip = -mean(get.var.ncdf(nc_gs,"dlsurf"))*L
	SHF = mean(get.var.ncdf(nc_gs,"qsurf"))

	# verticalstats
	ncpath = paste(datadir,"/zeroO3_",SST,"K/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	time = get.var.ncdf(nc,"time")
	z = get.var.ncdf(nc,"z")
	nz = length(z)
	nt = length(time)
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)	
	cloud = apply(get.var.ncdf(nc,start=c(1,nt-10),"cloud"),1,mean)	
	kanvil = which.max(cloud)
	klcl = which.max(cloud[1:15])
#	klcl = which.min(abs(z-zlcl))
	lapse = - partialder_i2s(3,z,s2i(3,z,tabs))
	p = apply(get.var.ncdf(nc,start=c(1,nt-10),"p"),1,mean)	
    ktp = which.min(abs(tabs-200))
    Ttp=tabs[ktp]
    ptp=p[ktp]
    gamma = g/Rd*log(Ttp/SST)/log(ptp/ps)
	kstar = which.min(abs(lapse[klcl:kanvil]-gamma)) + klcl -1
	lwup = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwup"),1,mean)	
	lwdown = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwdown"),1,mean)	
	swup = apply(get.var.ncdf(nc,start=c(1,nt-10),"swup"),1,mean)	
	swdown = apply(get.var.ncdf(nc,start=c(1,nt-10),"swdown"),1,mean)	
	Flw = lwup-lwdown
	Fsw = swup-swdown
	Fnet = Flw+Fsw
	ppzflw = partialder_i2s(3,z,Flw)
	ppzfsw = partialder_i2s(3,z,Fsw)
	ppzfnet = ppzflw+ppzfsw
	Qlw = Flw[nz]-Flw[klcl]
	Qsw = Fsw[nz]-Fsw[klcl]
	Qnet = Fnet[nz]-Fnet[klcl]
	QLWvec[i] = Qlw
	QSWvec[i] = Qsw
	QNetvec[i] = Qnet
	dqdtsLW_vec[i] = ppzflw[klcl]/lapse[klcl] 
	dqdtsSW_vec[i] = ppzfsw[klcl]/lapse[klcl] 
	dqdtsNet_vec[i] = ppzfnet[klcl]/lapse[klcl] 
	precipvec[i] = precip
	SHFvec[i] = SHF
	}

pdf(file="~/Dropbox/rad_cooling/git/figures/Qnet_varsst.pdf",width=12,height=5)
par(mfrow=c(1,3),mar=c(5,6,5,3))
cex=2

for (channel in c("SW","LW","Net")) {
	Qvec = eval(as.name(paste("Q",channel,"vec",sep="")))
	dqdts_vec = eval(as.name(paste("dqdts",channel,"_vec",sep="")))
	
	# Construct line segments from predicted slopes
	dT = 7
	x0 = SSTlist - dT/2
	x1 = SSTlist + dT/2
	y0 = Qvec - dqdts_vec*dT/2
	y1 = Qvec + dqdts_vec*dT/2
	if (channel == "Net"){
        ylim=c(min(y0[1],precipvec[1]),max(y1[N],precipvec[N]))
		} else {
		ylim=c(y0[1],y1[N])
		}		
	plot(type="b",SSTlist,Qvec,xlab=expression(T[s]~~"(K)"),
		ylab=bquote(Q[.(channel)]~~"("~W/m^2~")"),
		main=paste(channel," cooling vs. SST",sep=""),
		xlim=c(x0[1],x1[N]),
		ylim=ylim,
		cex.lab=cex,cex.main=cex,cex.axis=cex,pch = 16,cex=2,
		lty="dashed",lwd=2)
	segments(x0,y0,x1,y1,col="red",lty="solid",lwd=2.5)	
	if (channel == "Net"){
		points(SSTlist,precipvec,pch=8, col="blue",cex=1.5)     
		legend("topleft",c(expression("CRM"~~ Q),"CRM P","Theory"),
                lwd=c(NA,NA,2),pch=c(16,8,NA),col=c("black","blue","red"),cex=1.5)
		}
	}
	
dev.off()


