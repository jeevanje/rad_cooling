load('../data/f.RData')

pdf('../figures/fswlw_all.pdf',width=10,height=7)
layout(matrix(1:6,nrow=2))
par(mar=c(5,5,4,3))
cex    = 2
cex_leg= 1.35
lwd    = 2
xlim   = c(-3,7)
ylim   = c(300,200)
Ts	   = 290
Tslist = c(Ts,Ts+4)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red")

for (model_k in 1:length(model_names)) {
   model 	= model_names[model_k]
   plot(1,type="n",
 		xlim	 = xlim,
   		ylim	 = ylim,
        xlab	 = expression(-partialdiff[T]*F),
        ylab	 = 'T (K)',
        main     = model,
        xaxs	 = "i",
        cex.lab  = cex,
		cex.axis = cex,
        cex.main = cex)
   for (case_k in 1:2){
 	   Tsindex <- which(Tvals==Tslist[case_k])
 	   col  = colvec[case_k]
	   xsw <- Fvals[,Tsindex,case_k,1,model_k]
	   xlw <- Fvals[,Tsindex,case_k,2,model_k]
	   xsw[which(Farea[,Tsindex,1,model_k]<0.5*max(Farea[,Tsindex,1,model_k]))] <-NA
	   xlw[which(Farea[,Tsindex,1,model_k]<0.5*max(Farea[,Tsindex,1,model_k]))] <-NA
	   points(-diff(xsw)/2,Tvals[-1],type='l',lty="dashed",col=col,lwd=lwd)
	   points(-diff(xlw)/2,Tvals[-1],type='l',lty="solid",col=col,lwd=lwd)
	   abline(v=0,lty="dotted",lwd=1)	
	}
	if (model == "CanAM4"){
		legend(x=2,y=220,legend=c(caselist,"LW","SW"),col=c(colvec,"black","black"),
		lwd=2,lty=c("solid","solid","solid","dashed"),cex=cex_leg,bty="n")
	}
}

dev.off()
