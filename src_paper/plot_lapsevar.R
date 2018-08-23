load("../data/lapse.Rdata")
load("../data/f.Rdata")  # for Farea_s -- didn't save gammaarea_s!

model_k = 3
case_k  = 1  # AMIP for Farea_s
model	= model_names[model_k]
Tvals_s = seq(151,349,by=2)  # scalar levels 
cutoff  = 0.9  # area cutoff 

varlim  = c(0,5)  # (K/km)^2
tabslim = c(300,220)
lwd     = 2
cex 	= 1.5

file = paste("../figures_paper/lapsevar_",model,".pdf",sep="")
pdf(file=file,width=9,height=7,bg="white")
par(mfrow=c(2,3))
par(mar=c(5,5,4,3))

for (Ts in 25:30*10){
    Tsindex    = which(Tvals==Ts)
#    kmax       = length(Tvals_s)
    kmax       = which(Tvals==(Ts+2))
    kvec       = 1:kmax
    gamma      = gammavals[ ,Tsindex,model_k]
    gamma2     = gamma2vals[ ,Tsindex,model_k]
    gammavar   = gamma2 - gamma^2
    kext       = min(which((Tvals_s > 240)&(gammavar>0.5e-6))) 
    Text       = Tvals_s[kext]  
    gammavar_ext = gammavar[kext]
    gammavar[which(Farea_s[,Tsindex,case_k,model_k] < 
				cutoff*max(Farea_s[,Tsindex,case_k,model_k]))] <-NA

    plot(1e6*gammavar[kvec],Tvals_s[kvec],
	xlim = varlim,
	ylim = tabslim,
	type = "l",
	lwd  = lwd,
    xlab = expression(V*a*r*"("*Gamma*")    ("*K^2/km^2*")"),
    ylab = "Temperature (K)",
	main = bquote(T[s]==.(Ts)~"K"),
        cex.axis = cex,
        cex.main = cex,
        cex.lab  = cex 
     ) 
    abline(v=0,lty="dashed",lwd=1)
    abline(h=Text,lty="longdash",col="black",lwd=1.5)
    text(0.5*max(varlim),Text-5,bquote(T==.(Text)~K),cex=cex,col="black")
    points(1e6*gammavar_ext,Text,pch=16 ,cex=1.75)
}
dev.off()