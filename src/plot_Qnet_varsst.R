library(ncdf)
library(fields)
load("../data/crm.Rdata")

pdf(file="../figures/Qnet_varsst.pdf",width=12,height=5)
par(mfrow=c(1,3),mar=c(5,6,5,3))
cex=2.5
cex_leg = 1.75

for (channel in c("SW","LW","Net")) {
	Qvec = eval(as.name(paste("Q",channel,"vec",sep="")))
	dqdts_vec = eval(as.name(paste("dqdts",channel,"_vec",sep="")))
	
	# Construct line segments from predicted slopes
	dT = 7
	x0 = SSTlist - dT/2
	x1 = SSTlist + dT/2
	y0 = Qvec - dqdts_vec*dT/2
	y1 = Qvec + dqdts_vec*dT/2

	# Set y-limits
	if (channel == "Net"){
        ylim=c(min(y0[1],precipvec[1]),max(y1[N],precipvec[N]))
		} else {
		ylim=c(min(y0,y1),max(y0,y1))
		}		

	# Plot Q
	plot(type="b",SSTlist,Qvec,xlab=expression(T[s]~~"(K)"),
		ylab=bquote(Q[.(channel)]~~"("~W/m^2~")"),
		main=paste(channel," cooling vs. SST",sep=""),
		xlim=c(x0[1],x1[N]),
		ylim=ylim,
		cex.lab=cex,cex.main=cex,cex.axis=cex,pch = 16,cex=cex,
		lty="dashed",lwd=2)

	# Plot slopes
	slopevec=1:5
		segments(x0[slopevec],y0[slopevec],x1[slopevec],y1[slopevec],col="red",lty="solid",lwd=2.5)	

	# Add precip and legend
	if (channel == "Net"){
		points(SSTlist,precipvec,pch=8, col="blue",cex=cex,type="p")     
		legend("topleft",c(expression("CRM"~~ Q),"CRM LP","Eqn. (9)"),
                lwd=c(NA,NA,2),pch=c(16,8,NA),
                col=c("black","blue","red"),cex=cex_leg)
		}

	}
	
dev.off()


