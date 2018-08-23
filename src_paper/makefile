# Variables
figpath = ../figures_paper
datapath= ../data

# Use 'make [class]' to make figures in that class. Only list one target per 
# R script so scripts are not run more than once.   
# $* is auto variable for the stem (what is matched by %)

fig_names    := rhov pptfnet_tinv_dam Qnet_varsst fnet_ipsl\
                crm_crh pptflw_tinv_dam pptfsw_tinv_dam pptfc tauk_tinv\
                rh_all fswlw_all
pdf_names    := $(addprefix $(figpath)/, $(addsuffix .pdf, $(fig_names)))
script_names := $(addprefix plot_, $(addsuffix .R, $(fig_names)))

.PHONY: generic 
generic: $(pdf_names) 

$(figpath)/%.pdf: plot_%.R $(datapath)/crm.Rdata $(datapath)/crm_wvp.Rdata\
                  $(datapath)/f.Rdata $(datapath)/lapse.Rdata $(datapath)/rh.Rdata
	Rscript plot_$*.R


.PHONY: specific
specific: $(figpath)/fnet_all_290.pdf $(figpath)/lapsevar_IPSL-CM5A-LR.pdf\
          $(figpath)/fnetvar_IPSL-CM5A-LR.pdf $(figpath)/fnetcs_270_all.pdf 

$(figpath)/fnet_all_290.pdf: plot_fnet_all.R $(datapath)/f.Rdata $(datapath)/lapse.Rdata
	Rscript $<

$(figpath)/lapsevar_IPSL-CM5A-LR.pdf: plot_lapsevar.R $(datapath)/lapse.Rdata
	Rscript $<

$(figpath)/fnetvar_IPSL-CM5A-LR.pdf: plot_fnetvar.R $(datapath)/f.Rdata
	Rscript $<

$(figpath)/fnetcs_270_all.pdf: plot_fnetcs_all.R  $(datapath)/f.Rdata
	Rscript $<

.PHONY: clean
clean:
	rm ${figpath}/*
