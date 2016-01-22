library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/my_image_plot.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/git/src/h2o_data.R")


# Theory data
kvals1 = seq(100,999,length.out=100)   # cm^-1
kvals2 = seq(1001,1500,length.out=100) # cm^-1
kvals  = c(kvals1,kvals2)				# cm^-1
nk	   = length(kvals)
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))
lk = 6500 # m^-1

# Get dam data
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
for (i in 1:N){
	Ts  = SSTlist[i]
	ncpath = paste(datadir,"/zeroO3_",Ts,"K/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	z    = get.var.ncdf(nc,"z")
	zint    = zinterp(z)
	nz   = length(z)
	dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
	time = get.var.ncdf(nc,"time")
	nt = length(time)
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)	
	kmin=which.min(tabs)
	rho = apply(get.var.ncdf(nc,start=c(1,nt-10),"rho"),1,mean)	
	qv = apply(get.var.ncdf(nc,start=c(1,nt-10),"qv"),1,mean)	
	p = apply(get.var.ncdf(nc,start=c(1,nt-10),"p"),1,mean)	
	rhov = qv*rho
	lapse = - partialder_i2s(3,z,s2i(3,z,tabs))
	ktp = which.min(abs(tabs-200))
	Ttp=tabs[ktp]
	ptp=p[ktp]
	gamma = g/Rd*log(Ttp/Ts)/log(ptp/ps)
	ts_array = outer(Tstar_vals,tabs,temp_scaling)
	dtaudz = (kappa_vals%o%(rhov*p/p0))*ts_array	#sss
	tau	= array(dim=c(nk,nz+1))	 				#ssi
	tau[,nz+1] <- 0
	for (k in nz:1){
		tau[ ,k] = tau[,k+1] + dzvec[k]*dtaudz[ ,k]
		} 
	weight = array(dim=c(nk,nz))					# 1/m,sss
	for (k in  1:nz){
		weight[ ,k] = dtaudz[ ,k]*exp(-0.5*(tau[ ,k]+tau[ ,k+1]))
		}
	k1=numeric(nz)
	for (k in 1:nz) {
		k1[k]=kvals[which.max(weight[ ,k])]
		}
	assign(paste("tabs",Ts,sep=""),tabs)	
	assign(paste("lapse",Ts,sep=""),lapse)	
	assign(paste("weight",Ts,sep=""),weight)	
	assign(paste("Tweight",Ts,sep=""),Tweight)	
	assign(paste("k1",Ts,sep=""),k1)	
	assign(paste("kmin",Ts,sep=""),kmin)	
	assign(paste("zvec",Ts,sep=""),1:kmin)	
	assign(paste("gamma",Ts,sep=""),gamma)	
	}



cex=1.25
plot_it = function(field,tabs,main,zlim=range(field),zvec){
			my.image.plot(kvals,rev(tabs[zvec]),field[ ,rev(zvec)],
			  ylim = rev(range(tabs[zvec])),
			  zlim = zlim,
			  xlab = "wavenumber (cm^-1)",
			  ylab = "Temperature (K)",
			  main = main,
			  cex.lab = cex,
			  cex.axis = cex,
			  cex.main = cex)
			}

# Begin plot
pdf1x1("~/Dropbox/rad_cooling/git/figures/k1_varsst.pdf")
k1lim = c(0,5e-4)
cex=1.5
colvec=tim.colors(N)
for (i in 1:N){
    Ts  = SSTlist[i]
	tabs =  eval(as.name(paste("tabs",Ts,sep="")))
	k1 =  eval(as.name(paste("k1",Ts,sep="")))
	gamma =  eval(as.name(paste("gamma",Ts,sep="")))
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))	
	b_kappa = 1.5e-4  # -dTstar/dk, 1/m^-1
	b_Tstar = 2.2e-2  # -dlnkappa/dk, K/m^-1
	l_k		= 1/(b_kappa + b_Tstar*(1/tabs-1/T0))
	k1_theory = l_k*( log(kappa_h20(0)/kappa0) + log(ps/p0) - g/Rd/gamma*log(Ts/tabs) - L/Rv/tabs -Tstar_h20(0)*(1/tabs-1/T0) ) #m^-1
	k1_theory[k1_theory < 1e4 ] <- 1e4
	if (i == 1){
		plot(k1[zvec],tabs[zvec],ylim=c(310,200),xlim=c(0,1000),
			ylab = "Temperature (K)",
			xlab = expression(k[1]~~"("~cm^-1~")"),
			main = expression(k[1]*"("*T*")"),
			type="l", col=colvec[i], lwd=2,
			cex.axis = cex,
			cex.main = cex,
			cex.lab  = cex			
			)
		} else {
		points(k1[zvec],tabs[zvec], type="l", col=colvec[i], lwd=2)
		}
	points(1e-2*k1_theory[zvec],tabs[zvec], type="l", col=colvec[i], lwd=2,lty="dashed")
	}
dev.off()

