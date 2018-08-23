library(ncdf)
library(fields)
load("../data/crm.Rdata")

SSTlist = c(280,290,300,310)
caselist= c("lw","sw")
namelist= c("LW","SW")
N       = length(SSTlist)
colvec  = tim.colors(length(SSTlist))
cex     = 2.25
cex_leg = 1.75
tabslim = c(max(SSTlist),170)
pptflw_lim = c(0,5)
pptfsw_lim = c(-2,0)
pptfnet_lim = c(0,3)
lty     = "solid"

pdf(file="../figures/pptfc.pdf",width=12,height=5)
par(mfrow=c(1,3),mar=c(5,6,5,3))

for (i in 1:2){
    case = caselist[i]
    name = namelist[i]
    main = paste(name," clear-sky flux divergence",sep="")
    xlab = bquote(-partialdiff[T]*F[.(name)]^{cs}~~"("~W/m^2/K~")") 
    xlim = eval(as.name(paste("pptf",case,"_lim",sep="")))

    for (n in 1:N){
        Ts   = SSTlist[n]
        col  = colvec[n]
        zvec =  eval(as.name(paste("zvec",Ts,sep="")))
        pptf =  eval(as.name(paste("pptf",case,Ts,sep="")))[zvec]
        pptfc=  eval(as.name(paste("pptf",case,"c",Ts,sep="")))[zvec]
        tabs =  eval(as.name(paste("tabs",Ts,sep="")))[zvec]
		if (n == 1){
	    	plot(pptf,tabs,xlab=xlab,main=main,
				ylim = tabslim, 
    	    	xlim = xlim,
	    		ylab = "Temperature (K)",
	    		col  = col,
    	    	type = "l",
    			lty  = lty,
    			lwd  = 2,
	    		cex.main=cex,
	    		cex.axis=cex,
	    		cex.lab=cex
		     	)
     	} else {
	    	points(pptf,tabs,type="l",lwd=2, col = col,lty=lty)
   	 	}
	 	if (i==1){
			 legend("topright",legend=SSTlist,lwd=2,cex=cex_leg,col=colvec,
	 				title=expression(T[s] ~~ "(K)"))
	 	}
     }  # Ts loop
} # channel loop

# clear vs. cloudy
pptf  = pptflw300  + pptfsw300
pptfc = pptflwc300 + pptfswc300
xlab  = bquote(-partialdiff[T]*F[net]~~"("~W/m^2/K~")") 
col   = colvec[3]
plot(pptf[zvec],tabs,xlab=xlab,main="Net flux divergence, Ts=300 K",
	ylim = tabslim, 
	xlim = pptfnet_lim,
	ylab = "Temperature (K)",
	col  = col,
    type = "l",
    lty  = lty,
    lwd  = 2,
	cex.main=cex,
	cex.axis=cex,
	cex.lab=cex
		     	)
points(pptfc[zvec],tabs,type="l",lwd=2, col = col,lty="dashed")
legend("topright",c("all-sky","clear-sky"),
		lty=c("solid","dashed"),lwd=2,cex=cex_leg)	        
			
dev.off()
