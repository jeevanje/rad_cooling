library(ncdf)
library(fields)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")
source("~/Dropbox/rad_cooling/src/h2o_data.R")

# Get dam data
SSTlist = c(280,290,300,310)
N       = length(SSTlist)
datadir = "~/Dropbox/rad_cooling/data"
for (i in 1:N){
	SST  = SSTlist[i]
	ncpath = paste(datadir,"/zeroO3_",SST,"K/data/verticalstats.nc",sep="")
	nc = open.ncdf(ncpath)
	z    = get.var.ncdf(nc,"z")
	zint    = zinterp(z)
	nz   = length(z)
	dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
	time = get.var.ncdf(nc,"time")
	nt = length(time)
	tabs = apply(get.var.ncdf(nc,start=c(1,nt-10),"tabs"),1,mean)	
	lwup = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwup"),1,mean)	
	lwdown = apply(get.var.ncdf(nc,start=c(1,nt-10),"lwdown"),1,mean)	
	F = lwup-lwdown
	ppzf = partialder_i2s(3,z,F)
	lapse = -partialder_s2i(3,z,tabs)
	lapse[1] = lapse[2]  # fudge
	assign(paste("nc",SST,sep=""),nc)
	assign(paste("nt",SST,sep=""),nt)
	assign(paste("tabs",SST,sep=""),tabs)	
	assign(paste("F",SST,sep=""),F)			
	assign(paste("ppzf",SST,sep=""),ppzf)			
	assign(paste("lapse",SST,sep=""),lapse)			
	}
	
zmax = 25e3
kmax = which.min(abs(zmax-z))
zvec = 1:kmax

# Begin plot
pdf(file="~/Dropbox/rad_cooling/git/figures/f_pptf.pdf",width=10,height=5)

# plotting parameters
par(mfrow=c(1,2),mar=c(6,5,5,3.5))
cex=2
colvec = tim.colors(length(SSTlist))
tabslim = c(310,170)
Flim=c(0,300)
lty="solid"
lwuplim=c(200,525)
pptflim=c(0,6)

#=============#
# DAM fluxes  #
#=============#

for (i in 1:N){
	SST  = SSTlist[i]
	tabs = eval(as.name(paste("tabs",SST,sep="")))
	F = eval(as.name(paste("F",SST,sep="")))	
	if (i == 1){
		plot(F[zvec],tabs[zvec],
			xlab=expression(F~~"("~W/m^2~")"), 
			ylab = "Temperature [K]", 
			main=expression("Longwave flux (DAM)"),
			type="l", lwd=2,
			log="",
			ylim = tabslim, xlim=Flim,
   			col=colvec[i],
   			cex.axis = cex,	cex.lab  = cex, cex.main = cex )
		} else {
		points(F[zvec],tabs[zvec],type="l",lwd=2, col = colvec[i])
		}	
	}


#===========#
# DAM pptf  #
#===========#
 
for (i in 1:N){
	Ts  = SSTlist[i]
	ppzf =  eval(as.name(paste("ppzf",Ts,sep="")))
	tabs =  eval(as.name(paste("tabs",Ts,sep="")))
	lapse =  eval(as.name(paste("lapse",Ts,sep="")))
	pptf = ppzf/lapse
	if (i == 1){
		plot(pptf[zvec],tabs[zvec],
			xlab=expression(partialdiff[T]*F~~"("~W/m^2/K~")"), 
		    main=expression("T-derivative"),
			 ylim = tabslim, 
			 xlim=pptflim,
			 ylab = "Temperature (K)",
			 col=colvec[i],
			 type="l", lty=lty,lwd=2,
			 cex.main=cex,
			 cex.axis=cex,
			 cex.lab=cex
		 )
		} else {
		points(pptf[zvec],tabs[zvec],type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}

dev.off()
