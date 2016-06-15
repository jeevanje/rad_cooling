library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")

# Open nc files
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
for (i in 1:N){
	SST  = SSTlist[i]
	ncpath = paste(datadir,"/prod_data_4_22_16/",SST,"k/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	time = get.var.ncdf(nc,"time")
	nt = length(time)
	navg = 20 # days to average
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-navg),"tabs"),1,mean)	
	rho = apply(get.var.ncdf(nc,start=c(1,nt-navg),"rho"),1,mean)	
	qv = apply(get.var.ncdf(nc,start=c(1,nt-navg),"qv"),1,mean)	
	assign(paste("nc",SST,sep=""),nc)
	assign(paste("nt",SST,sep=""),nt)
	assign(paste("tabs",SST,sep=""),tabs)	
	assign(paste("rhov",SST,sep=""),qv*rho)	
	}
	
z    = get.var.ncdf(nc280,"z")
nz   = length(z)
zmax = 25e3
kmax = which.min(abs(zmax-z))
zvec = 1:kmax
colvec = tim.colors(length(SSTlist))
tabslim = c(310,170)
rhovlim = c(1e-8,4e-2)
lwd =2.5
cex = 1.5

# Begin plots
pdf(file="~/Dropbox/rad_cooling/git/figures/rhov.pdf",width=9,height=5)
par(mfrow=c(1,2),mar=c(5,5,5,3))

# linear plot
plot(1,type="n",
     xlab=expression(rho[v]~~"("~kg/m^3~")"), 
     main="Vapor density",
     ylim = tabslim, 
     xlim = rhovlim,
     ylab = "Temperature (K)",
	 log  = "",
     cex.main=cex,
     cex.axis=cex,
     cex.lab=cex
                 )
for (i in 1:N){
	SST  = SSTlist[i]
	tabs = eval(as.name(paste("tabs",SST,sep="")))
	rhov = eval(as.name(paste("rhov",SST,sep="")))
	points(rhov[zvec],tabs[zvec],type="l",lwd=lwd, col = colvec[i])
	}
legend("topright",legend=SSTlist,lty="solid",col=colvec,lwd=2,title="SST (K)",cex=1.25)

	

# Log plot
plot(1,type="n",
     xlab=expression(rho[v]~~"("~kg/m^3~")"), 
     main="Vapor density",
     ylim = tabslim, 
     xlim = rhovlim,
     ylab = "Temperature (K)",
	 log  = "x",
     cex.main=cex,
     cex.axis=cex,
     cex.lab=cex
                 )
for (i in 1:N){
	SST  = SSTlist[i]
	tabs = eval(as.name(paste("tabs",SST,sep="")))
	rhov = eval(as.name(paste("rhov",SST,sep="")))
	points(rhov[zvec],tabs[zvec],type="l",lwd=lwd, col = colvec[i])
	}	
	
dev.off()
