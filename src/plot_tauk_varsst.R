library(fields)
library(ncdf)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/git/src/h2o_data.R")

# Theory data 
kvals1 = seq(100,999,length.out=100)     # cm^-1
kvals2 = seq(1001,1500,length.out=100)
kvals  = c(kvals1,kvals2)
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))

# DAM data, for computing gamma and then tau
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
for (i in 1:N){
        Ts  = SSTlist[i]
        ncpath = paste(datadir,"/zeroO3_",Ts,"K/data/verticalstats.nc",sep="")
        nc = open.ncdf(ncpath)
        z    = get.var.ncdf(nc,"z")
        time = get.var.ncdf(nc,"time")
        nt = length(time)
        tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)   
        kmin=which.min(tabs)
        p = apply(get.var.ncdf(nc,start=c(1,nt-10),"p"),1,mean)   
        ktp = which.min(abs(tabs-200))
        Ttp=tabs[ktp]
        ptp=p[ktp]
        gamma = g/Rd*log(Ttp/Ts)/log(ptp/ps)
		zvec   = 1:kmin
		Tav 	= 1/2*(Ts+200)
		wvpinf 	= RH*einf*Tav/gamma/L
		ts_array   = outer(Tstar_vals,tabs,temp_scaling)
		tau = ((kappa_vals)%o%(ps/p0*(tabs/Ts)^(g/Rd/gamma)*wvpinf*exp(-L/Rv/tabs)))*ts_array
        assign(paste("tabs",Ts,sep=""),tabs)  
        assign(paste("gamma",Ts,sep=""),gamma)  
        assign(paste("tau",Ts,sep=""),tau)  
        assign(paste("zvec",Ts,sep=""),zvec)  
		}


#============#
# Begin plot #
#============#
cex=1.25
logtaulim = c(-14,9)
tabslim  = c(SSTlist[N],200)
kindvec   = c(which.min(abs(kvals-100)),which.min(abs(kvals-500)),which.min(abs(kvals-900)))
Nk		  = length(kindvec)
colvec 	  = tim.colors(N)
pdf1x1("~/Dropbox/rad_cooling/git/figures/tauk_varsst.pdf")
plot(1,type="n", xlab=expression(ln~ tau[k]),ylab="Temperature (K)",
	xlim=logtaulim,ylim=tabslim,cex.lab=cex,cex.main=cex,cex.axis=cex)
for (n in 1:N){  
	Ts	= SSTlist[n]
	tau = eval(as.name(paste("tau",Ts,sep="")))
	zvec = eval(as.name(paste("zvec",Ts,sep="")))
	tabs = eval(as.name(paste("tabs",Ts,sep="")))
	for (m in 1:Nk){
		kind = kindvec[m]
		points(log(tau[kind,zvec]),tabs[zvec],type="l",lwd=2,col=colvec[n],lty="solid")
		}
	}
	#legend("topright",c("Eqn. (3)","Eqn. (9)"),
	#		lty==c("black",colvec),cex=1.25,lwd=2)
dev.off()
