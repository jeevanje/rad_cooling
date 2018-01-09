load('~/Dropbox/rad_cooling/git/data/f.Rdata')

model_k <- which(model_names=='IPSL-CM5A-LR')

pdf('../figures/fnet_ipsl.pdf',width=9,height=7)
par(mfrow=c(2,3))
par(mar=c(5,5,4,3))
cex    = 2
cex_leg= 1.5
lwd    = 2
xlim   = c(0,6)
ylim   = c(305,200)
Tspert = c(0,4)
caselist = c("AMIP","AMIP4K")
colvec	= c("blue","red")

for (Tsbase in 25:30*10){
   plot(1,type="n",
 		xlim	 = xlim,
   		ylim	 = ylim,
        xlab	 = expression(-partialdiff[T]*F[net]),
        ylab	 = 'T (K)',
        main     = bquote(T[s] == .(Tsbase)*","~.(Tsbase+4)~"K" ),
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
	   x[which(Farea[,Tsindex,1,model_k]<0.5*max(Farea[,Tsindex,1,model_k]))] <-NA
	   points(-diff(x)/2,Tvals[-1],type='l',col=col,lwd=lwd)
	}
	if (Tsbase == 250){
		legend("bottomright",legend=caselist,col=colvec,lwd=2,cex=cex_leg)
	}
}

dev.off()
