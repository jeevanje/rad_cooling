library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")

# Open nc files
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/fat/data"
for (i in 1:N){
	SST  = SSTlist[i]
	ncpath = paste(datadir,"/zeroO3_",SST,"K/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	time = get.var.ncdf(nc,"time")
	nt = length(time)
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)	
	rho = apply(get.var.ncdf(nc,start=c(1,nt-10),"rho"),1,mean)	
	qv = apply(get.var.ncdf(nc,start=c(1,nt-10),"qv"),1,mean)	
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

# Begin plot
pdf(file="~/Dropbox/rad_cooling/git/figures/rhov.pdf",width=9,height=5)
par(mfrow=c(1,2),mar=c(5,5,5,3))
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
