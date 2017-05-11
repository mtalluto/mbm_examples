### GUISANE
g_mbm=mbm/mbmtools.py 1_guisane/code/2_guisane_gp.py
g_dat=1_guisane/dat
g_res=1_guisane/res
g_img=1_guisane/img
g_code=1_guisane/code
g_py=cd 1_guisane; code/2_guisane_gp.py --dir=../mbm

.PHONY: g_beta g_alpha g_all g_figs

g_alpha: $(g_res)/richness/params.csv $(g_res)/simpson/params.csv $(g_res)/funcAlpha/params.csv \
$(g_res)/phyloAlpha/params.csv
g_beta: $(g_res)/funBeta/params.csv $(g_res)/phyBeta/params.csv $(g_res)/taxBeta/params.csv
g_all: g_alpha g_beta
g_figs: $(g_img)/fig1.pdf $(g_img)/fig2.pdf

# guisane figures
$(g_img)/fig1.pdf: g_alpha $(g_code)/fig1.r $(g_code)/mbm_support.r $(g_dat)/alpha_fit.csv
	$(g_code)/fig1.r
$(g_img)/fig2.pdf: g_alpha $(g_code)/fig2.r $(g_code)/mbm_support.r $(g_dat)/alpha_fit.csv
	$(g_code)/fig2.r
$(g_img)/fig3.pdf: g_beta $(g_code)/fig3.r $(g_code)/mbm_support.r $(g_dat)/beta_fit.csv
	$(g_code)/fig3.r

# guisane alpha diversity
$(g_res)/richness/params.csv: $(g_dat)/alpha_fit.csv $(g_dat)/alpha_valid.csv $(g_mbm)
	$(g_py) --richness
$(g_res)/simpson/params.csv: $(g_dat)/alpha_fit.csv $(g_dat)/alpha_valid.csv $(g_mbm)
	$(g_py) --simpson
$(g_res)/funcAlpha/params.csv: $(g_dat)/alpha_fit.csv $(g_dat)/alpha_valid.csv $(g_mbm)
	$(g_py) --fa
$(g_res)/phyloAlpha/params.csv: $(g_dat)/alpha_fit.csv $(g_dat)/alpha_valid.csv $(g_mbm)
	$(g_py) --pa

## guisane beta diversity
$(g_res)/taxBeta/params.csv: $(g_dat)/betaTaxo_fit.csv $(g_dat)/betaTaxo_valid.csv $(g_mbm)
	$(g_py) --tb
$(g_res)/funBeta/params.csv: $(g_dat)/betaFun_fit.csv $(g_dat)/betaFun_valid.csv $(g_mbm)
	$(g_py) --fb
$(g_res)/phyBeta/params.csv: $(g_dat)/betaPhylo_fit.csv $(g_dat)/betaPhylo_valid.csv $(g_mbm)
	$(g_py) --pb


# guisane data processing
$(g_dat)/alpha_fit.csv: 1_guisane/code/1_compute_indices.r
	cd 1_guisane; code/1_compute_indices.r
$(g_dat)/*_valid.csv $(g_dat)/betaTaxo%fit.csv $(g_dat)/betaFun%fit.csv $(g_dat)/betaPhylo%fit.csv: \
$(g_dat)/alpha_fit.csv
	noop