load('~/Dropbox/rad_cooling/git/data/f.RData')
load('~/Dropbox/rad_cooling/data/Q.Rdata')
source("~/Dropbox/Rtools/gcm_tools.R")
library(fields)

model_k <- which(model_names=='IPSL-CM5A-LR')
dqdtsvec = numeric()

pdf('~/Dropbox/rad_cooling/git/figures/fnet_all.pdf',width=9,height=7)
layout(matrix(1:6,nrow=2))
par(mar=c(5,5,4,3))
cex    = 2
cex_leg= 1.25
lwd    = 2
xlim   = c(0,5)
ylim   = c(300,200)
Tsbase = 290
Tslist = Tsbase + c(0,4)
caselist = c("AMIP","AMIP4K")
lowcaselist = c("amip","amip4K")
colvec	= c("blue","red","green3")
#colvec  = tim.colors(3)
#colvec  = c(colvec[1],colvec[3],colvec[2])
cutoff  = 0.9  # area cutoff for upper trop
deltaT  = 0    # cutoff for lower trop
dt_bl   = 20   # insertion point
NT		= length(Tvals)

for (model_k in 1:length(model_names)) {
   model 	= model_names[model_k]
   Qdata    = eval(as.name(paste("Qdata_",model,sep="")))	
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
	   Ts	= Tslist[case_k]
 	   Tsindex <- which(Tvals==Ts)
 	   col  = colvec[case_k]
	   xlw <- Fvals[,Tsindex,case_k,2,model_k]
	   xsw <- Fvals[,Tsindex,case_k,1,model_k]
	   x   = xlw + xsw
	   x[which(Farea[,Tsindex,case_k,model_k] < 
				cutoff*max(Farea[,Tsindex,case_k,model_k]))] <-NA
	   x[which(Tvals > (Ts - deltaT))] <- NA
	   points(-diff(x)/2,Tvals[-1],type='l',col=col,lwd=lwd)

   	   # construct synthetic profile
	   if (case_k == 1){
			k_bl     = Tsindex - floor(dt_bl/2)
			pptf_bl  = (x[k_bl+1] -x[k_bl])/2
	   		pptf_syn = c(diff(x[1:(k_bl+1)])/2,pptf_bl,pptf_bl,diff(x[(k_bl+1):(NT-2)])/2)
	   		points(-pptf_syn,Tvals[-1],type='l',col=colvec[3],lwd=lwd,lty="solid")
	   }
	}  # case_k

	# Q
	lon      = Qdata["lon"][[1]]
	lat      = Qdata["lat"][[1]]	
	Q_amip   = global_mean(lon,lat,Qdata["Qlw_amip"][[1]] + Qdata["Qsw_amip"][[1]])
	Q_amip4K = global_mean(lon,lat,Qdata["Qlw_amip4K"][[1]] + Qdata["Qsw_amip4K"][[1]])
	dqdts    = round((Q_amip4K - Q_amip)/4.5,digits=2)
	dqdtsvec[model_k] = dqdts
#	text(1,208,bquote(frac(d*Q,d*T[s])==.(dqdts)),cex=1.25)
#	text(3,240,bquote(-partialdiff[T]*F[ext]==.(-pptf)),cex=1.25)
	pptf	 = round(pptf_bl,digits=2)
	abline(v=-pptf,lty="dashed")
	
	if (model_k == 1){
		legend("topright",legend=c(caselist,expression(A*M*I*P[ext])),col=colvec,lwd=2,cex=cex_leg)
	}
}

dev.off()
