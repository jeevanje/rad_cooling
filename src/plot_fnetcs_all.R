load('../data/f.RData')

Tslist = c(270,290)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red")
ltyvec  = c("solid","dashed")
fluxlist= c("all-sky","clear-sky")
legend_posvec =c("bottom","topright")

for (NTs in 1:length(Tslist)){
	Ts  = Tslist[NTs]
	file=paste("../figures/fnetcs_",Ts,"_all.pdf",sep="")
	pdf(file,width=9,height=7)
	layout(matrix(1:6,nrow=2))
	par(mar=c(5,5,4,3))
	cex    = 2
	cex_leg= 1.25
	lwd    = 2
	xlim   = c(0,6)
	ylim   = c(300,200)
	Ts_case= c(Ts,Ts+4)

	for (model_k in 1:length(model_names)) {
	   model 	= model_names[model_k]
	   plot(1,type="n",
	 		xlim	 = xlim,
	   		ylim	 = ylim,
	        xlab	 = expression(-partialdiff[T]*F[net]),
	        ylab	 = 'T (K)',
	        main     = model,
	        xaxs	 = "i",
	        cex.lab  = cex,
			cex.axis = cex,
	        cex.main = cex)
	   for (case_k in 1:2){
		   case_k = 1	
	 	   Tsindex <- which(Tvals==Ts_case[case_k])
	 	   col  = colvec[case_k]
		   for (flux_k in 1:2){	
			   lty = ltyvec[flux_k]
			   xsw <- Fvals[,Tsindex,case_k,2*flux_k -1,model_k]
			   xlw <- Fvals[,Tsindex,case_k,2*flux_k,model_k]
			   x   = xlw + xsw
		   	   x[which(Farea[,Tsindex,1,model_k]<0.5*max(Farea[,Tsindex,1,model_k]))] <-NA
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