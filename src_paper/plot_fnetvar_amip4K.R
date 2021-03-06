load("../data/f.Rdata")
source("./plot_params.R")

case_k  = 2
model_k = 3
model   = model_names[model_k]
NT_s	= length(Tvals_s)
cutoff  = 0.9  # area cutoff 

pptFnetvals  = pptFvals[ , ,case_k,1,model_k] + pptFvals[ , ,case_k,2,model_k]
pptFnet2vals = pptF2vals[ , ,case_k,1,model_k] + pptF2vals[ , ,case_k,2,model_k]

varlim  = c(-5,25)
tabslim = c(300,220)
lwd     = 2
cex 	= 1.5
file    = paste("../figures_paper/fnetvar_",model,"_amip4K.pdf",sep="") 
pdf(file=file,width=9,height=7,bg="white")
par(mfrow=c(2,3))
par(mar=c(5,5,4,3))

for (Ts in (25:30*10)+4){
    Tsindex    = which(Tvals==Ts)
    kmax       = which(Tvals==(Ts+2))
#    kvec       = 1:kmax
    kvec       = 1:NT_s     
    pptFnet    = pptFnetvals[ ,Tsindex]
    pptFnet2   = pptFnet2vals[ ,Tsindex]
    pptFnetvar = pptFnet2 - pptFnet^2
    kext       = min(which( (Tvals_s > 240)&(pptFnetvar > 5) ))   # s lev 
    Text       = Tvals_s[kext]  
    pptFnetvar_ext = pptFnetvar[kext]
	pptFnetvar[which(Farea_s[,Tsindex,case_k,model_k] < 
				cutoff*max(Farea_s[,Tsindex,case_k,model_k]))] <-NA

    plot(pptFnetvar[kvec],Tvals_s[kvec],
		xlim=varlim,
		ylim=tabslim,
		type="l",
		lwd=lwd,
		xlab = pptfnetvarlab,
		ylab = "Temperature (K)",
		main = bquote(T[s]==.(Ts)~"K"),
		cex.axis = cex,
		cex.main = cex,
		cex.lab  = cex 
	) 
    abline(v=0,lty="dotted",col="black",lwd=1)
    abline(h=Text,lty="longdash",col="black",lwd=1.5)
    text(varlim-8,Text-5,bquote(T[ext]==.(Text)~K),cex=cex,col="black")
    points(pptFnetvar_ext,Text,pch=16 ,cex=1.75)

}
dev.off()