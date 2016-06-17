library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")

# Get data
N		= 5
SSTlist = 270 + 10*(1:N)
datadir = "~/Dropbox/rad_cooling/data/prod_data_4_22_16"
dqdtsLW_vec = numeric(N)
dqdtsSW_vec = numeric(N)
dqdtsNet_vec = numeric(N)
QSWvec	= numeric(N)
QLWvec	= numeric(N)
QNetvec	= numeric(N)
precipvec	= numeric(N)
SHFvec	= numeric(N)
zlcl = 125
nt_avg_gs = 20*4  # 20 days x 4x daily sampling
nt_avg_vs = 20  # 20 days 

for (i in 1:N){
	SST  = SSTlist[i]
    # precip
    globalstatspath = paste(datadir,"/",SST,"k/data/globalstats.nc",sep="")
    nc_gs = open.ncdf(globalstatspath)
    time_gs = get.var.ncdf(nc_gs,"time")
    nt_gs = length(time_gs)
	start = nt_gs-nt_avg_gs
    precip = -mean(get.var.ncdf(nc_gs,start=start,"dlsurf"))*L
	precipvec[i] = precip
    SHF = -mean(get.var.ncdf(nc_gs,start=start,"qsurf"))
	precipvec[i] = precip
	SHFvec[i] = SHF

	# verticalstats
	ncpath = paste(datadir,"/",SST,"k/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	time = get.var.ncdf(nc,"time")
	z = get.var.ncdf(nc,"z")
	nz = length(z)
	nt = length(time)
	start = c(1,nt-nt_avg_vs)
	tabs = apply(get.var.ncdf(nc,start=start,"tabs"),1,mean)	
	cloud = apply(get.var.ncdf(nc,start=start,"cloud"),1,mean)	
	if (SST != 330){
		klcl = which.max(cloud[1:15])
		} else {
		klcl = klcl320	
	}
#	klcl = which.min(abs(z-zlcl))
	lapse = - partialder_i2s(3,z,s2i(3,z,tabs))
	lwup = apply(get.var.ncdf(nc,start=start,"lwup"),1,mean)	
	lwdown = apply(get.var.ncdf(nc,start=start,"lwdown"),1,mean)	
	swup = apply(get.var.ncdf(nc,start=start,"swup"),1,mean)	
	swdown = apply(get.var.ncdf(nc,start=start,"swdown"),1,mean)	
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

	# Set y-limits
	if (channel == "Net"){
        ylim=c(min(y0[1],precipvec[1]),max(y1[N],precipvec[N]))
		} else {
		ylim=c(min(y0,y1),max(y0,y1))
		}		

	# Plot Q
	plot(type="b",SSTlist,Qvec,xlab=expression(T[s]~~"(K)"),
		ylab=bquote(Q[.(channel)]~~"("~W/m^2~")"),
		main=paste(channel," cooling vs. SST",sep=""),
		xlim=c(x0[1],x1[N]),
		ylim=ylim,
		cex.lab=cex,cex.main=cex,cex.axis=cex,pch = 16,cex=cex,
		lty="dashed",lwd=2)

	# Plot slopes
	slopevec=1:5
	segments(x0[slopevec],y0[slopevec],x1[slopevec],y1[slopevec],col="red",lty="solid",lwd=2.5)	

	# Add precip and legend
	if (channel == "Net"){
		points(SSTlist,precipvec,pch=8, col="blue",cex=cex)     
		legend("topleft",c(expression("CRM"~~ Q),"CRM P","Eqn. (9)"),
                lwd=c(NA,NA,2),pch=c(16,8,NA),
                col=c("black","blue","red"),cex=1.5)
		}

	}
	
dev.off()


