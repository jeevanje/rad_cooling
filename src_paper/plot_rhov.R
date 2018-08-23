library(ncdf)
library(fields)
load("../data/crm.Rdata")

SSTlist = c(280,290,300,310)
N       = length(SSTlist)	
zmax = 25e3
kmax = which.min(abs(zmax-z))
zvec = 1:kmax
colvec = tim.colors(length(SSTlist))
tabslim = c(310,170)
rhovlim = c(1e-8,4e-2)
lwd =2.5
cex = 1.5
leg_title = expression(T[s]~~"(K)")

# Begin plots
pdf(file="../figures_paper/rhov.pdf",width=9,height=5)
par(mfrow=c(1,2),mar=c(5,5,5,3))

# linear plot
plot(1,type="n",
     xlab=expression(rho[v]~~"("~kg/m^3~")"), 
     main="Vapor density",
     ylim = tabslim, 
     xlim = rhovlim,
     ylab = "Temperature (K)",
	 log  = "",
     cex.main=cex,
     cex.axis=cex,
     cex.lab=cex
                 )
for (i in 1:N){
	SST  = SSTlist[i]
	tabs = eval(as.name(paste("tabs",SST,sep="")))
	rhov = eval(as.name(paste("rhov",SST,sep="")))
	points(rhov[zvec],tabs[zvec],type="l",lwd=lwd, col = colvec[i])
	}
legend("topright",legend=SSTlist,lty="solid",col=colvec,
		lwd=2,title=leg_title,cex=1.25)

	

# Log plot
plot(1,type="n",
     xlab=expression(rho[v]~~"("~kg/m^3~")"), 
     main="Vapor density",
     ylim = tabslim, 
     xlim = rhovlim,
     ylab = "Temperature (K)",
	 log  = "x",
     cex.main=cex,
     cex.axis=cex,
     cex.lab=cex
                 )
for (i in 1:N){
	SST  = SSTlist[i]
	tabs = eval(as.name(paste("tabs",SST,sep="")))
	rhov = eval(as.name(paste("rhov",SST,sep="")))
	points(rhov[zvec],tabs[zvec],type="l",lwd=lwd, col = colvec[i])
	}	
	
dev.off()
