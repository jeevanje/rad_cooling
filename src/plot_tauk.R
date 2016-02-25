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
nk	   = length(kvals)
kint   = 1/2*( c(kvals[1]-diff(kvals)[1],kvals)+c(kvals,kvals[nk]+diff(kvals)[nk-1]) )
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))

# varsstdata for computing gamma and then tau
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
for (i in 1:N){
        Ts  = SSTlist[i]
        ncpath = paste(datadir,"/zeroO3_",Ts,"K/data/verticalstats.nc",sep="")
        nc = open.ncdf(ncpath)
        z    = get.var.ncdf(nc,"z")
		zint   = zinterp(z)
		nz	   = length(z)
		dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
        time = get.var.ncdf(nc,"time")
        nt = length(time)
        tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)   
        kmin=which.min(tabs)
        p = apply(get.var.ncdf(nc,start=c(1,nt-10),"p"),1,mean)   
        rho = apply(get.var.ncdf(nc,start=c(1,nt-10),"rho"),1,mean)   
        qv = apply(get.var.ncdf(nc,start=c(1,nt-10),"qv"),1,mean)   
		rhov   = qv*rho 
        ktp = which.min(abs(tabs-200))
        Ttp=tabs[ktp]
        ptp=p[ktp]
        gamma = g/Rd*log(Ttp/Ts)/log(ptp/ps)
        zvec   = 1:kmin
        Tav     = 1/2*(Ts+200)
        wvpinf  = RH*einf*Tav/gamma/L
        ts_array   = outer(Tstar_vals,tabs,temp_scaling)
        tau1 = ((kappa_vals)%o%(ps/p0*(tabs/Ts)^(g/Rd/gamma)*wvpinf*exp(-L/Rv/tabs)))
		dtaudz = (kappa_vals%o%(rhov*p/p0))	#sss
		tau	= array(dim=c(nk,nz+1))	 				#ssi
		tau[,nz+1] <- 0
		for (k in nz:1){
				tau[ ,k] = tau[,k+1] + dzvec[k]*dtaudz[ ,k]
		} 
        assign(paste("tabs",Ts,sep=""),tabs)  
        assign(paste("gamma",Ts,sep=""),gamma)  
        assign(paste("tau",Ts,sep=""),tau)  
        assign(paste("tau1",Ts,sep=""),tau1)  
        assign(paste("zvec",Ts,sep=""),zvec)  
                }

# taucheck data
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
Tav 	= 1/2*(Ts+200)
wvpinf 	= RH*einf*Tav/gamma/L
ts_array   = outer(Tstar_vals,tabs,temp_scaling)

# Ground truth -- parametrize kappa, use temp and p scalings from Humbert (eqn. 3)
dtaudz = (kappa_vals%o%(rhov*p/p0))	#sss
tau	= array(dim=c(nk,nz+1))	 				#ssi
tau[,nz+1] <- 0
for (k in nz:1){
		tau[ ,k] = tau[,k+1] + dzvec[k]*dtaudz[ ,k]
} 

# Approx 1 -- approximate tau integral, p-scaling at z only
tau1 = ((kappa_vals)%o%(ps/p0*(tabs/Ts)^(g/Rd/gamma)*wvpinf*exp(-L/Rv/tabs)))

cex=1.25
logtaulim = c(-15,8)
tabslim  = c(SSTlist[N],200)
kindvec   = c(which.min(abs(kvals-100)),which.min(abs(kvals-500)),which.min(abs(kvals-900)))
Nk		  = length(kindvec)
colvec 	  = tim.colors(Nk)

#=============#
# Begin plots #
#=============#
pdf("~/Dropbox/rad_cooling/git/figures/tauk.pdf",width=10,height=5)
par(mfrow=c(1,2),oma=c(0,1,0,1))

# tau check
plot(1,type="n",
      ylim = rev(range(tabs)),
	  xlim = logtaulim,
	  xlab = expression(ln~tau[k]),
	  ylab = "Temperature (K)",
	  main = "Optical depth",
	  cex.lab = cex,
	  cex.axis = cex,
	  cex.main = cex,
		)
for (n in 1:Nk){
	kind = kindvec[n]
	points(log(tau[kind,zvec]),tabs[zvec],type="l",lwd=2,col=colvec[n],lty="solid")
	points(log(tau1[kind,zvec]),tabs[zvec],type="l",lwd=2,col=colvec[n],lty="dashed")
	#legend("topright",c("Eqn. (3)","Eqn. (9)"),
	#		lty==c("black",colvec),cex=1.25,lwd=2)
	}


# varsst
plot(1,type="n", xlab=expression(ln~ tau[k]),ylab="Temperature (K)",main="Optical Depth",
        xlim=logtaulim,ylim=tabslim,cex.lab=cex,cex.main=cex,cex.axis=cex)
for (n in 1:N){  
        Ts      = SSTlist[n]
        tau = eval(as.name(paste("tau1",Ts,sep="")))
        zvec = eval(as.name(paste("zvec",Ts,sep="")))
        tabs = eval(as.name(paste("tabs",Ts,sep="")))
        for (m in 1:Nk){
                kind = kindvec[m]
                points(log(tau[kind,zvec]),tabs[zvec],type="l",lwd=2,col=colvec[n],lty="solid")
                }
        }
text(-9.7,280,expression(k==900~~cm^-1),cex=1)
text(5.7,240,expression(100~~cm^-1),cex=1)
text(-3.0,220,expression(500~~cm^-1),cex=1)
dev.off()
