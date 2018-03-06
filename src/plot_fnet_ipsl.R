load('~/Dropbox/rad_cooling/git/data/f.Rdata')

model_k <- which(model_names=='IPSL-CM5A-LR')

pdf('~/Dropbox/rad_cooling/git/figures/fnet_ipsl.pdf',width=9,height=7)
par(mfrow=c(2,3))
par(mar=c(5,5,4,3))
cex    = 2
cex_leg= 1.5
lwd    = 2
xlim   = c(0,5)
ylim   = c(305,200)
Tspert = c(0,4)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red","green3")
cutoff  = 0.9  # area cutoff for upper trop
deltaT  = 0    # cutoff for lower trop
dt_bl   = 20   # insertion point
NT		= length(Tvals)

for (Tsbase in 25:30*10){
   Q = numeric(2)
   plot(1,type="n",
 		xlim	 = xlim,
   		ylim	 = ylim,
        xlab	 = expression(-partialdiff[T]*F[net]),
        ylab	 = 'T (K)',
        main     = bquote(T[s] == bold(.(Tsbase))*","~bold(.(Tsbase+4))~"K" ),
        xaxs	 = "r",
        cex.lab  = cex,
		cex.axis = cex,
        cex.main = cex)
   for (case_k in 1:2){
	   Ts = Tsbase + Tspert[case_k]
 	   Tsindex <- which(Tvals==Ts)
 	   col  = colvec[case_k]
	   xsw <- Fvals[,Tsindex,case_k,1,model_k]
	   xlw <- Fvals[,Tsindex,case_k,2,model_k]
	   x   = xlw + xsw
	   x[which(Farea[,Tsindex,case_k,model_k] < 
	   			cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
	   x[which(Tvals > (Ts - deltaT))] <- NA
	   Q[case_k] = diff(range(x,na.rm=TRUE))
	   points(-diff(x)/2,Tvals[-1],type='l',col=col,lwd=lwd)

   	   # construct synthetic profile
	   if (case_k == 1){
			k_bl     = Tsindex - floor(dt_bl/2)
			pptf_bl  = (x[k_bl+1] -x[k_bl])/2
	   		pptf_syn = c(diff(x[1:(k_bl+1)])/2,pptf_bl,pptf_bl,diff(x[(k_bl+1):(NT-2)])/2)
	   		points(-pptf_syn,Tvals[-1],type='l',col=colvec[3],lwd=lwd,lty="solid")
	   }
	}  # case
	pptf	= round(pptf_bl,digits=2)
	dqdts   = round(diff(Q),digits=2)/4
#	abline(v=-pptf,lty="dashed")
#	text(3,208,bquote(frac(d*Q,d*T[s])==.(dqdts)),cex=1.25)
#	text(3,240,bquote(-partialdiff[T]*F[ext]==.(-pptf)),cex=1.25)

	if (Tsbase == 250){
		legend("bottomright",legend=c(caselist,expression(A*M*I*P[ext])),
				col=colvec,lwd=2,cex=cex_leg)
	}
}  # Tsbase

dev.off()
