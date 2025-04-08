
meta <- read.csv('clinical_agro_meta.csv', header=T)
amr.summary <- read.csv('results/amrfinder.out.summary.csv', header=T)

head(meta)
head(amr.summary)


combined <-  merge(amr.summary, meta[, c('amr_strain', 'sample_type')], by.x = 'strain', by.y = 'amr_strain')

# Need to update the
meta[which( ! meta$amr_strain %in% combined$strain), 'strain'] 

head(combined)

table(combined$desc, combined$sample_type)
table(meta$sample_type)

############################################### 
# ROARY presence-absence table header compare #
###############################################

strain_list <- c("AP_261","AP_262","Agrobacterium_fabrum_C58","Agrobacterium_genomosp_CFBP_5494","Agrobacterium_pusense_17_1009","Agrobacterium_pusense_76","Agrobacterium_pusense_973_Rhizob3","Agrobacterium_pusense_CC001","Agrobacterium_pusense_CC134","Agrobacterium_pusense_CCGM10","Agrobacterium_pusense_CCGM11","Agrobacterium_pusense_CFBP5496","Agrobacterium_pusense_CFBP5875","Agrobacterium_pusense_CNPSo_3352","Agrobacterium_pusense_FDAARGOS_618","Agrobacterium_pusense_FDAARGOS_619","Agrobacterium_pusense_FDAARGOS_633","Agrobacterium_pusense_GD03648","Agrobacterium_pusense_GD03663","Agrobacterium_pusense_GD03874","Agrobacterium_pusense_GD03880","Agrobacterium_pusense_GD03943","Agrobacterium_pusense_GD03981","Agrobacterium_pusense_GD03982","Agrobacterium_pusense_GD04154","Agrobacterium_pusense_GV2","Agrobacterium_pusense_HU21","Agrobacterium_pusense_HU37","Agrobacterium_pusense_HU39","Agrobacterium_pusense_HU55","Agrobacterium_pusense_HU57","Agrobacterium_pusense_HU59","Agrobacterium_pusense_IRBG74","Agrobacterium_pusense_KCJK7997","Agrobacterium_pusense_LMG_25623","Agrobacterium_pusense_MGBC108980","Agrobacterium_pusense_NRCPB10","Agrobacterium_pusense_RpEC2071","Agrobacterium_pusense_S1","Agrobacterium_pusense_S2","Agrobacterium_pusense_SCN18_30_10_14_R3_B_60_7","Agrobacterium_pusense_ST15_13_056","Agrobacterium_pusense_SX41","Agrobacterium_pusense_ZB01","Agrobacterium_sp_33MFTa1_1","Agrobacterium_sp_33MFTa1_1_ASM432849v1","Agrobacterium_sp_Bin_012_Pasture","Agrobacterium_sp_FDAARGOS_525","Agrobacterium_sp_InxBP2","Agrobacterium_sp_MS2","Agrobacterium_sp_RC10_4_1","Agrobacterium_sp_S2","Agrobacterium_sp_UBA11112","Agrobacterium_sp_UBA11724","Agrobacterium_tumefaciens_MH_0_3_111223_53_1677754698","Agrobacterium_tumefaciens_MH_0_5_111223_17_1677754708","Agrobacterium_tumefaciens_S33","Agrobacterium_tumefaciens_UBA2572","Agrobacterium_tumefaciens_iMGMC_16_1678881632","BWH_1514_R_rad","BWH_80_R_rad","DTU_2020_1000670_1_SI_DNK_HVI_038","Inst1_Iso1_2014","Inst2_Iso1_2014","Inst2_Iso2_2017","Inst3_Iso1_2009","Inst3_Iso2_2014","Inst3_Iso3_2014","Inst3_Iso4_2015","Inst3_Iso5_2019","Inst6_Iso10_2018","Inst6_Iso11_2018","Inst6_Iso9_2017","Inst7_Iso1_2012","Inst7_Iso2_2016","Rhizobium_sp_GHKF11","Rhizobium_sp_H41","Rhizobium_sp_LMB_1","Rhizobium_sp_P007","Rhizobium_sp_RAS22","Rhizobium_sp_S2_005_001_R1_24","Rhizobium_sp_S41","Rhizobium_sp_SORGH_AS260","Rhizobium_sp_SORGH_AS285","Rhizobium_sp_SORGH_AS_755","Rhizobium_sp_UR51a","Rhizobium_sp_Y9","Rhizobium_sp_ZX09")

meta[which(! meta$scoary_strain %in% strain_list), 'strain']