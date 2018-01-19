load('../data/rh.RData')

model_k <- which(model_names=='IPSL-CM5A-LR')

pdf('../figures/rh_all.pdf',width=9,height=7)
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

for (model_k in 1:length(model_names)) {
   model 	= model_names[model_k]
   plot(1,type="n",
 		xlim	 = RHlim,
   		ylim	 = ylim,
        xlab	 = "RH",
        ylab	 = 'T (K)',
        main     = model,
        xaxs	 = "i",
        cex.lab  = cex,
		cex.axis = cex,
        cex.main = cex)
   for (case_k in 1:2){
 	   Tsindex <- which(Tvals==Tslist[case_k])
 	   col  = colvec[case_k]
	   x <- RHvals[,Tsindex,case_k,model_k]
	   x[which(RHarea[,Tsindex,1,model_k]<0.5*max(RHarea[,Tsindex,1,model_k]))] <-NA
	   points(x,Tvals,type='l',col=col,lwd=lwd)
	}
	if (model == "IPSL-CM5A-LR"){
		legend("topright",legend=caselist,col=colvec,lwd=2,cex=cex_leg)
	}
}

dev.off()
