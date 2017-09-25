library(ncdf)
library(fields)
load("../data/crm.Rdata")

SSTlist = c(280,290,300,310)
N       = length(SSTlist)
colvec = tim.colors(length(SSTlist))
cex=2
xlim=c(0,5)
tabslim = c(max(SSTlist),170)
lty="solid"
xlim=c(-2,0)
main = "SW flux divergence"
xlab = expression(-partialdiff[T]*F[SW]~~"("~W/m^2/K~")") 

pdf(file="../figures/pptfsw_tinv_dam.pdf",width=12,height=5)
par(mfrow=c(1,3),mar=c(5,6,5,3))

#===========#
# height    #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	col = colvec[i]
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	pptf =  eval(as.name(paste("pptfsw",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,1e-3*z[zvec],
			xlab=xlab,
			ylab = "z (km)",
			main=main,
			xlim=xlim,
			ylim = c(0,1e-3*zmax),
		 	col=colvec[i],lty=lty,
		 	cex.axis=cex,
		 	cex.main=cex,
		 	cex.lab=cex,
		 	type="l",
		 	lwd=2
			 )
		} else {
		points(pptf,1e-3*z[zvec],type="l",lwd=2, col = col,lty=lty)
		}
	legend("topleft",legend=SSTlist,lwd=2,cex=1.5,col=colvec,title=expression(T[s] ~~ "(K)"))
	}


#===========#
# pressure  #
#===========#

for (i in 1:N){
	Ts  = SSTlist[i]
	col = colvec[i]
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	pptf =  eval(as.name(paste("pptfsw",Ts,sep="")))[zvec]
	p =  eval(as.name(paste("p",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,1e-2*p,
			xlab=xlab,
			ylab = "Pressure (hPa)",
			main=main,
			xlim=xlim,
			ylim=c(1000,0),
		 	col=colvec[i],lty=lty,
		 	cex.axis=cex,
		 	cex.main=cex,
		 	cex.lab=cex,
		 	type="l",
		 	lwd=2
			 )
		} else {
		points(pptf,1e-2*p,type="l",lwd=2, col = colvec[i],lty=lty)
		}
	}

#==============#
# temperature  #
#==============#

for (i in 1:N){
	Ts  = SSTlist[i]
	col = colvec[i]
	zvec =  eval(as.name(paste("zvec",Ts,sep="")))
	pptf =  eval(as.name(paste("pptfsw",Ts,sep="")))[zvec]
	tabs =  eval(as.name(paste("tabs",Ts,sep="")))[zvec]
	if (i == 1){
		plot(pptf,tabs,
			xlab=xlab,
			main=main,
			 ylim = tabslim, xlim=xlim,
			 ylab = "Temperature (K)",
			 col=colvec[i],type="l", lty=lty,lwd=2,
			 cex.main=cex,
			 cex.axis=cex,
			 cex.lab=cex
		 )
		} else {
		points(pptf,tabs,type="l",lwd=2, col = colvec[i],lty=lty)
		}
	abline(v=0,col="gray",lty = "dashed",lwd=1.5)
	}
			
dev.off()
