load('../data/f.RData')
source("./plot_params.R")

Tslist   = c(270,290)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red")
ltyvec  = c("solid","dashed")
fluxlist= c("all-sky","clear-sky")
legend_posvec =c("bottom","topright")
cutoff  = 0.9  # for upper trop
deltaT  = 0    # for lower trop

for (NTs in 1:length(Tslist)){
	Tsbase  = Tslist[NTs]
	file=paste("../figures_paper/fnetcs_",Tsbase,"_all.pdf",sep="")
	pdf(file,width=9,height=7)
	layout(matrix(1:6,nrow=2))
	par(mar=c(5,5,4,3))
	cex    = 2
	cex_leg= 1.25
	lwd    = 2
	xlim   = c(0,6)
	ylim   = c(300,200)
	Ts_case= Tsbase + c(0,4)

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
	   add_pptflab(pptfnetlab)
	   for (case_k in 1:2){
		   case_k = 1	
		   Ts = Ts_case[case_k]
	 	   Tsindex <- which(Tvals==Ts)
	 	   col  = colvec[case_k]
		   for (flux_k in 1:2){	
			   lty = ltyvec[flux_k]
			   xsw <- Fvals[,Tsindex,case_k,2*flux_k -1,model_k]
			   xlw <- Fvals[,Tsindex,case_k,2*flux_k,model_k]
			   x   = xlw + xsw
		   	   x[which(Farea[,Tsindex,case_k,model_k] < 
		   	   		cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
			   x[which(Tvals > (Ts - deltaT))] <- NA
		   	   points(-diff(x)/2,Tvals[-1],type='l',col=col,lwd=lwd,lty=lty)
			}
		# if (model == model_names[3]){
			# legend("topright",legend=caselist,col=colvec,lwd=2,cex=cex_leg)
		# }
		if (model == model_names[3]){
			legend_pos = legend_posvec[NTs]
			legend(legend_pos,legend=fluxlist,lty=ltyvec,lwd=2,cex=cex_leg)
		}
		}  # case_k
	}  # model_k
	dev.off()
}  # Ts