all: \
../figures/rhov.pdf \
../figures/pptflw_tinv_dam.pdf \
../figures/pptfsw_tinv_dam.pdf \
../figures/Qnet_varsst.pdf \
../figures/fnet_ipsl.pdf \
../figures/fnet_all.pdf \
../figures/fswlw_all.pdf \

../figures/%.pdf: plot_%.R
	R CMD BATCH --no-restore --no-save $<
