load('~/Dropbox/16rad_cooling/git/data/f.Rdata')
load("~/Dropbox/16rad_cooling/git/data/lapse.Rdata")
source("./plot_params.R")

model_k <- which(model_names=='IPSL-CM5A-LR')
model_k = 3

pdf('../figures_paper/fnet_ipsl.pdf',width=9,height=7)
par(mfrow=c(2,3))
par(mar=c(5,5,4,3))
cex    = 2.25
cex_leg= 1.75
lwd    = 2
xlim   = c(0,5)
ylim   = c(305,200)
Tspert = c(0,4)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red","green3")
cutoff  = 0.9  # area cutoff for upper trop
deltaT  = 0    # cutoff for lower trop
dt_bl   = 16   # insertion point
NT		= length(Tvals)

for (Tsbase in 25:30*10){
   Q = numeric(2)
   plot(1,type="n",
		xlim	 = xlim,
		ylim	 = ylim,
        xlab	 = "",
        ylab	 = 'T (K)',
        main     = bquote(T[s] == bold(.(Tsbase))*","~bold(.(Tsbase+4))~"K" ),
        xaxs	 = "r",
        cex.lab  = cex,
	    cex.axis = cex,
        cex.main = cex)
   add_pptflab(pptfnetlab)
   for (case_k in 1:2){
	   Ts = Tsbase + Tspert[case_k]
 	   Tsindex <- which(Tvals==Ts)
 	   col     = colvec[case_k]
	   fnet    = Fvals[,Tsindex,case_k,1,model_k] + Fvals[,Tsindex,case_k,2,model_k]
	   fnet[which(Farea[,Tsindex,case_k,model_k] < 
	   			cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
	   fnet[which(Tvals > (Ts - deltaT))] <- NA
	   Q[case_k] = diff(range(fnet,na.rm=TRUE))
	   points(-diff(fnet)/2,Tvals_s,type='l',col=col,lwd=lwd)
   }  # case_k

   # construct synthetic profile
   case_k = 1
   Ts = Tsbase + Tspert[case_k]
   Tsindex <- which(Tvals==Ts)
   fnet    = Fvals[,Tsindex,case_k,1,model_k] + Fvals[,Tsindex,case_k,2,model_k]
   fnet[which(Farea[,Tsindex,case_k,model_k] <
   								 cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA

   # kext from pptfvar 
   pptFnet  = pptFvals[ ,Tsindex,case_k,1,model_k] + 
   			  pptFvals[ ,Tsindex ,case_k,2,model_k]
   pptFnet2 = pptF2vals[ ,Tsindex,case_k,1,model_k] + 
 			  pptF2vals[ ,Tsindex ,case_k,2,model_k]
   pptFnetvar = pptFnet2 - pptFnet^2
   kext     = min(which( (Tvals_s > 240)&(pptFnetvar > 5 ) ))   # s lev

	# kext from gammavar if necesary
   if (is.infinite(kext))  {   
	   gamma      = gammavals[ ,Tsindex,model_k]
	   gamma2     = gamma2vals[ ,Tsindex,model_k]
	   gammavar   = gamma2 - gamma^2
	   kext       = min(which((Tvals_s > 240)&(gammavar>0.5e-6))) 
   }
   if (!is.infinite(kext))  {
       Text     = Tvals_s[kext]  
       pptf_bl  = (fnet[kext+1] -fnet[kext])/2
       pptf_syn = c(diff(fnet[1:(kext+1)])/2,pptf_bl,pptf_bl,
		   diff(fnet[(kext+1):(NT-2)])/2 )
       points(-pptf_syn,Tvals_s,type='l',col=colvec[3],
	      lwd=lwd+2,lty="dashed")
       points(-pptf_bl,Text,pch=16 ,cex=cex)
       pptf	= round(pptf_bl,digits=2)
   }

   if (Tsbase == 250){
      legend("bottomright",legend=c(caselist,expression(A*M*I*P[ext])),
			    col=colvec,
			    lwd=c(lwd,lwd,lwd+2),
			    cex=cex_leg,
			    lty=c("solid","solid","dashed"))
   }
}  # Tsbase

dev.off()
