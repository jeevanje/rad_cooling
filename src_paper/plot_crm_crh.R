library(fields)
source("~/Dropbox/Rtools/my_image_plot.R")
load("../data/crm_wvp.Rdata")

N = length(Tslist)

crhlim = c(0.60,1)
cex    = 2

pdf(file="../figures_paper/crm_crh.pdf",width=13,height=8)
par(mfrow=c(2,3),mar=c(5,5,5,9))
for (i in 1:N){
	Ts  = Tslist[i]
	crh = wvp_vals[ , ,i,1]/wvp_vals[ , ,i,2]
	my.image.plot(1e-3*x,1e-3*y,crh,
		zlim = crhlim,
		main = bquote("CRH,"~~T[s]==.(Ts)),
		xlab = "x (km)",
		ylab = "y (km)",
		cex.main = cex,
		cex.axis = cex,
		cex.lab  = cex,
		cex.legend  = cex
		)
	mtext("CRH",side=3,adj=1,outer=TRUE)
}
dev.off()