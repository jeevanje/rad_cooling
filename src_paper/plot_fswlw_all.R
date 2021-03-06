load('../data/f.Rdata')
source('./plot_params.R')

pdf('../figures_paper/fswlw_all.pdf',width=10,height=7)
layout(matrix(1:6,nrow=2))
par(mar=c(5,5,4,3))
cex    = 2
cex_leg= 1.35
lwd    = 2
xlim   = c(-3,7)
ylim   = c(300,200)
Tsbase = 290
Tslist = Tsbase + c(0,4)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red")
cutoff  = 0.9  # for upper trop
deltaT  = 0    # for lower trop

for (model_k in 1:length(model_names)) {
   model 	= model_names[model_k]
   plot(1,type="n",
	xlim	 = xlim,
	ylim	 = ylim,
        xlab	 = "",
        ylab	 = 'T (K)',
        main     = model,
        xaxs	 = "i",
        cex.lab  = cex,
	cex.axis = cex,
        cex.main = cex)
   add_pptflab(pptflab)
   for (case_k in 1:2){
	   Ts	= Tslist[case_k]
 	   Tsindex <- which(Tvals==Ts)
 	   col  = colvec[case_k]
	   xsw <- Fvals[,Tsindex,case_k,1,model_k]
	   xlw <- Fvals[,Tsindex,case_k,2,model_k]
	   xsw[which(Farea[,Tsindex,case_k,model_k] < 
	   			cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
	   xlw[which(Farea[,Tsindex,case_k,model_k] < 
	   		cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
	   xsw[which(Tvals > (Ts - deltaT))] <- NA
	   xlw[which(Tvals > (Ts - deltaT))] <- NA
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
