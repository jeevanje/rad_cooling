library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/git/src/h2o_data.R")

SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
kvals1 = seq(100,999,length.out=100)     # cm^-1
kvals2 = seq(1001,1500,length.out=100)
kvals  = c(kvals1,kvals2)
nk	   = length(kvals)
kint   = 1/2*( c(kvals[1]-diff(kvals)[1],kvals)+c(kvals,kvals[nk]+diff(kvals)[nk-1]) )
kappa_vals = c(kappa_h20(kvals1),kappa_h20(kvals2))
Tstar_vals = c(Tstar_h20(kvals1),Tstar_h20(kvals2))

for (i in 1:N){
	Ts  = SSTlist[i]
	ncpath = paste(datadir,"/zeroO3_",Ts,"K/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	z    = get.var.ncdf(nc,"z")
	nz   = length(z)
	zint   = zinterp(z)
	dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
	time = get.var.ncdf(nc,"time")
	nt = length(time)
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)	
	rho = apply(get.var.ncdf(nc,start=c(1,nt-10),"rho"),1,mean)	
	qv = apply(get.var.ncdf(nc,start=c(1,nt-10),"qv"),1,mean)	
	p = apply(get.var.ncdf(nc,start=c(1,nt-10),"p"),1,mean)	
	rhov = qv*rho
	lapse = - partialder_i2s(3,z,s2i(3,z,tabs))
	kmin=which.min(tabs)
	ktp = which.min(abs(tabs-200))
	Ttp=tabs[ktp]
	ptp=p[ktp]
	gamma = g/Rd*log(Ttp/Ts)/log(ptp/ps)
	ts_array = outer(Tstar_vals,tabs,temp_scaling)

	# Approx 1 -- numerically compute tau
	dtau1dz = (kappa_vals%o%(rhov*p/p0))*ts_array	#sss
	tau1	= array(dim=c(nk,nz+1))	 				#ssi
	tau1[,nz+1] <- 0
	for (k in nz:1){
			tau1[ ,k] = tau1[,k+1] + dzvec[k]*dtau1dz[ ,k]
	} 
	weight1 = array(dim=c(nk,nz))					# sss
	for (k in  1:nz){
		weight1[ ,k] = dtau1dz[ ,k]*exp(-0.5*(tau1[ ,k]+tau1[ ,k+1]))
	}
	source   = 1e2*t(pi*outer(tabs,1e2*kvals,planck_k))    # W/m^2/(cm^-1)
	ppzf1=diff(kint)%*%(weight1*source) # W/m^3

	# analytic formulation
	b_kappa = 1.5e-4  # -dTstar/dk, 1/m^-1
	b_Tstar = 2.2e-2  # -dlnkappa/dk, K/m^-1
	l_k		= 1/(b_kappa + b_Tstar*(1/tabs-1/T0))
	k1 = l_k*( log(kappa_h20(0)/kappa0) + log(ps/p0) - g/Rd/gamma*log(Ts/tabs) - L/Rv/tabs -Tstar_h20(0)*(1/tabs-1/T0) ) #m^-1
	k1[k1 < 0 ] <- 1e-4
	#gamma = .0065
	ppzf2 = pi*planck_k(tabs,k1)*(g/Rd/tabs + gamma/tabs^2*(L/Rv + Tstar_h20(1e-2*k1)))*l_k   # W/m^3
	assign(paste("tabs",Ts,sep=""),tabs)	
	assign(paste("p",Ts,sep=""),p)	
	assign(paste("lapse",Ts,sep=""),lapse)	
	assign(paste("kmin",Ts,sep=""),kmin)	
	assign(paste("zvec",Ts,sep=""),2:kmin)	
	assign(paste("gamma",Ts,sep=""),gamma)	
	assign(paste("ppzf1",Ts,sep=""),ppzf1)	
	assign(paste("ppzf2",Ts,sep=""),ppzf2)	
	}




# Plot parameters
colvec = tim.colors(length(SSTlist))
tabslim = c(310,170)
cex=2
xlim=c(0,0.03)    # W/m^3
lty="solid"

pdf1x3("~/Dropbox/rad_cooling/git/figures/ppzf_tinv_theory.pdf")
#===========#
# height    #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	ppzf =  eval(as.name(paste("ppzf2",Ts,sep="")))
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	if (i == 1){
		plot(ppzf[zvec],1e-3*z[zvec],
			xlab=expression(partialdiff[z]*F~~"("~W/m^3~")"), 
			ylab = "z (km)",
			main="LW Flux Convergence",
			xlim=xlim,
		 	col=colvec[i],lty=lty,
		 	cex.axis=cex,
		 	cex.main=cex,
		 	cex.lab=cex,
		 	type="l",
		 	lwd=2
			 )
		} else {
		points(ppzf[zvec],1e-3*z[zvec],type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}


#===========#
# pressure  #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	ppzf =  eval(as.name(paste("ppzf2",Ts,sep="")))
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	p =  eval(as.name(paste("p",Ts,sep="")))
	if (i == 1){
		plot(ppzf[zvec],1e-2*p[zvec],
			xlab=expression(partialdiff[z]*F~~"("~W/m^3~")"), 
			ylab = "Pressure (hPa)",
			main="LW Flux Convergence",
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
		points(ppzf[zvec],1e-2*p[zvec],type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}


#==============#
# temperature  #
#==============#
for (i in 1:N){
	Ts  = SSTlist[i]
	tabs =  eval(as.name(paste("tabs",Ts,sep="")))
	ppzf2 =  eval(as.name(paste("ppzf2",Ts,sep="")))
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	if (i == 1){
		plot(ppzf2[zvec],tabs[zvec],
			xlab=expression(partialdiff[z]*F~~"("~W/m^3~")"), 
			main="LW Flux Convergence",
			 ylim = tabslim, xlim=xlim,
			 ylab = "Temperature (K)",
			 col=colvec[i],type="l", lty=lty,
			 cex.main=cex,
			 cex.axis=cex,
			 cex.lab=cex
		 )
		} else {
		points(ppzf2[zvec],tabs[zvec],type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}

#==============#
# temperature  #
# benchmark    #	
#==============#
# for (i in 1:N){
	# Ts  = SSTlist[i]
	# tabs =  eval(as.name(paste("tabs",Ts,sep="")))
	# ppzf1 =  eval(as.name(paste("ppzf1",Ts,sep="")))
	# zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	# if (i == 1){
		# plot(ppzf1[zvec],tabs[zvec],
			# xlab=expression(partialdiff[z]*F~~"("~W/m^3~")"), 
			# main="LW Flux Convergence",
			 # ylim = tabslim, xlim=xlim,
			 # ylab = "Temperature (K)",
			 # col=colvec[i],type="l", lty=lty,
			 # cex.main=cex,
			 # cex.axis=cex,
			 # cex.lab=cex
		 # )
		# } else {
		# points(ppzf1[zvec],tabs[zvec],type="l",lwd=2, col = colvec[i],lty=lty)
		# }
	# }

		
dev.off()
