load('../data/f.Rdata')
load('../data/lapse.Rdata')
source("./plot_params.R")

Tsbases = c(280,290)
for (Tsbase in Tsbases){
	Tslist  = Tsbase + c(0,4)
	cutoff  = 0.9  # area cutoff for upper trop
	NT		= length(Tvals)
	
	cex      = 2.25
	cex_leg  = 1.75
	lwd      = 2
	pptflim  = c(0,5)
	tabslim  = c(300,200)
	colvec	 = c("blue","red","green3")
	caselist = c("AMIP","AMIP4K")
	lowcaselist = c("amip","amip4K")
	file = paste('../figures_paper/fnet_all_',Tsbase,'.pdf',sep="")
	
	pdf(file=file,width=9,height=7)
	layout(matrix(1:6,nrow=2))
	par(mar=c(5,5,4,3))
	for (model_k in 1:length(model_names)) {
	   model 	= model_names[model_k]
	   plot(1,type="n",
 		xlim	 = pptflim,
   		ylim	 = tabslim,
	        xlab	 = "",
	        ylab	 = 'T (K)',
	        main     = model,
	        xaxs	 = "i",
	        cex.lab  = cex,
		cex.axis = cex,
	        cex.main = cex)
	   add_pptflab(pptfnetlab)
	   for (case_k in 1:2){
		   Ts	   = Tslist[case_k]
	 	   Tsindex = which(Tvals==Ts)
	 	   col     = colvec[case_k]
		   xlw     = Fvals[,Tsindex,case_k,2,model_k]
		   xsw     = Fvals[,Tsindex,case_k,1,model_k]
		   x       = xlw + xsw
		   x[which(Farea[,Tsindex,case_k,model_k] < 
					cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
		   points(-diff(x)/2,Tvals_s,type='l',col=col,lwd=lwd)
		}  # case_k
	
	    # construct synthetic profile
	    case_k = 1
	    Ts = Tslist[case_k]
	    Tsindex = which(Tvals==Ts)
	    fnet    = Fvals[,Tsindex,case_k,1,model_k] + Fvals[,Tsindex,case_k,2,model_k]
	    fnet[which(Farea[,Tsindex,case_k,model_k] <
	    							cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
	
		# # kext from pptfvar
	    pptFnet  = pptFvals[ ,Tsindex,case_k,1,model_k] + 
	    		   pptFvals[ ,Tsindex ,case_k,2,model_k]
	    pptFnet2 = pptF2vals[ ,Tsindex,case_k,1,model_k] + 
	    		   pptF2vals[ ,Tsindex ,case_k,2,model_k]
	    pptFnetvar = pptFnet2 - pptFnet^2
	    kext     = min(which( (Tvals_s > 240)&(pptFnetvar >5) ))  # slevs
	
	    # kext from gammavar
	    # gamma      = gammavals[ ,Tsindex,model_k]
	    # gamma2     = gamma2vals[ ,Tsindex,model_k]
	    # gammavar   = gamma2 - gamma^2
	    # kext       = min(which((Tvals_s > 240)&(gammavar>0.5e-6))) 
	    if (!is.infinite(kext))  {
	       Text     = Tvals_s[kext]  
	       pptf_bl  = (fnet[kext+1] -fnet[kext])/2
	       pptf_syn = c(diff(fnet[1:(kext+1)])/2,pptf_bl,pptf_bl,
	                   diff(fnet[(kext+1):(NT-2)])/2 )
	       points(-pptf_syn,Tvals_s,type='l',col=colvec[3],
	              lwd=lwd+2,lty="dashed")
	       points(-pptf_bl,Text,pch=16 ,cex=cex)
	       pptf     = round(pptf_bl,digits=2)
	    }
	
		# Legend + lines
		if (model_k == 1){			
			lines(x=c(-pptf,-pptf),y=c(tabslim[1],240),lty="dashed",type="l")
			legend("topright",legend=c(caselist,expression(A*M*I*P[ext])),
					col=colvec,
				    lty=c("solid","solid","dashed"),	
					lwd=c(lwd,lwd,lwd+2),
				    cex=cex_leg
					)
		} else {
			abline(v=-pptf,lty="dashed")		
		}
	} # model_k
	
	dev.off()
}