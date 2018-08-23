pptflab       = expression(-partialdiff[T]*F~~"("~W/m^2/K~")") 
pptfnetlab    = expression(-partialdiff[T]*F[Net]~~"("~W/m^2/K~")") 
pptfnetvarlab = expression(V*a*r*"("*-partialdiff[T]*F[Net]*")"~~"("~W^2/m^4/K^2~")") 

add_pptflab   = function(pptflab){
				mtext(pptflab,side=1,line=3.5,cex=1.5)
} 