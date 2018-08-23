library(ncdf4)
library(fields)

datadir = "../data/rfm_data/"
nk     = 1491
nz     = 76
ssts = c(280,290,300,310)
Nsst = length(ssts)

# Get data
tauvals = array(dim=c(nk,nz,Nsst))
Tvals   = array(dim=c(nz,Nsst))
pvals   = array(dim=c(nz,Nsst))
for (n in 1:Nsst){
    sst = ssts[n]
    ncpath = paste(datadir,sst,"K.nc",sep="")
    nc	   = nc_open(ncpath)
    k	   = ncvar_get(nc,"k")
    z	   = ncvar_get(nc,"z")
    p	   = ncvar_get(nc,"p")
    tau	   = ncvar_get(nc,"opt")
    tabs   = ncvar_get(nc,"tabs")
    tauvals[ , ,n] <- tau
    Tvals[ ,n]     <- tabs
    pvals[ ,n]     <- p
}

# plot
tausvals = c(1e2,1,1e-2)   # in 280 K case
Ntaus    = length(tausvals)
mvals   = numeric(Ntaus)
kmin    = 1e4
mmin    = min(which(k>kmin))
lntaulim= c(-10,10)
tabslim	= c(320,200)
plim	= c(1000,100)
colvec  = tim.colors(Nsst)
lwd	= 2
cex	= 1.5
T1vec = numeric(Nsst)

for (n in 1:Ntaus){
    taus     = tausvals[n]
    mvals[n] = which.min(abs(tauvals[(k>kmin)&(k<10e4),1,1]-taus)) + mmin - 1
}

pdf("../figures_paper/tauk_tinv.pdf",width = 5,height =5)
par(mar=c(5,5,5,3))
plot(1,type="n",
    xlab = expression(ln~tau[lambda]),
    ylab = "Temperature (K)",
    main = "Optical Depth",
    xlim = lntaulim,
    ylim = tabslim,
    cex.axis = cex,
    cex.lab	 = cex,
    cex.main = cex)
for (m in mvals){    
    for (n in 1:Nsst){
        col = colvec[n]
    	points(log(tauvals[m,-nz,n]),Tvals[-nz,n],type="l",lwd=lwd,col=col)
		if (m==mvals[1]){
			k1 = which.min(abs(tauvals[m, ,n]-1))
			T1 = Tvals[k1,n]
			T1vec[n] = T1
		}
    }
	lambda = round(1e6/k[m],digits=1)
	text(log(tauvals[m,1,4])+0.25,313,bquote(.(lambda)~mu*m),cex=0.8)

}
abline(v=0,lty="dashed")
legend("topright",legend=ssts,col=colvec,title=expression(T[s]~~"(K)"),lwd=lwd)
dev.off()

