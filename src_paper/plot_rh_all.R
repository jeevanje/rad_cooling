load('~/Dropbox/16rad_cooling/git/data/rh.Rdata')

model_k <- which(model_names=='IPSL-CM5A-LR')

pdf('../figures_paper/rh_all.pdf',width=9,height=7)
layout(matrix(1:6,nrow=2))
par(mar=c(5,5,4,3))
cex    = 2
cex_leg= 1.25
lwd    = 2
RHlim  = c(0,100)
ylim   = c(300,200)
Ts	   = 290
Tslist = c(Ts,Ts+4)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red")
cutoff  = 0.9  # for upper trop
deltaT  = 5    # for lower trop

for (model_k in 1:length(model_names)) {
   model 	= model_names[model_k]
   plot(1,type="n",
	xlim	 = RHlim,
	ylim	 = ylim,
        xlab	 = "RH (%)",
        ylab	 = 'T (K)',
        main     = model,
        xaxs	 = "i",
        cex.lab  = cex,
		cex.axis = cex,
        cex.main = cex)
   for (case_k in 1:2){
	   Ts      = Tslist[case_k]
 	   Tsindex <- which(Tvals==Ts)
 	   col  = colvec[case_k]
	   x <- RHvals[,Tsindex,case_k,model_k]
	   x[which(RHarea[,Tsindex,case_k,model_k] < 
	   		cutoff*max(RHarea[,Tsindex,case_k,model_k]))] <-NA
       x[which(Tvals > (Ts - deltaT))] <- NA
	   points(x,Tvals,type='l',col=col,lwd=lwd)
	}
	if (model == "IPSL-CM5A-LR"){
		legend("topright",legend=caselist,col=colvec,lwd=2,cex=cex_leg)
	}
}

dev.off()
