#install.packages("dplyr")
#install.packages("genoPlotR")
#install.packages("magrittr")
# Library
library("magrittr")
library("dplyr")
library("genoPlotR")

# Global variable
setwd("E:\\Windows10\\desktop\\RBFOX1_GeneAcross\\")
vertebrate_tree = "(((((((((((pHuman_T2T:0.001,pHuman_hg38:0.001):6.40000000,(pChimpanzee:2.39810000,pBonobo:2.39810000):4.00190000):2.20000000,pGorilla:8.60000000):34.32000000,(pCommon_marmoset:17.65000000,pMa_s_night_monkey:17.65000000):25.27000000):26.23223000,pPhilippine_tarsier:69.15223000):18.04777000,((mHouse_mouse:10.25000000,mNorway_rat:10.25000000):0.00000000,mNile_rat:10.25000000):76.95000000):6.80000000,((((mBottlenose_dolphin:17.63872400,mVaquita:17.63872400):16.66752600,mBlue_whale:34.30625000):23.48185000,((mCow:15.93000000,mSheep:15.93000000):7.00315000,mRed_deer:22.93315000):34.85495000):23.48106500,(mGreater_horseshoe_bat:62.00000000,mPale_spear_nosed_bat:62.00000000):19.26916500):12.73083500):4.90000000,mNine_banded_armadillo:98.90000000):81.16610000,mPlatypus:180.06610000):138.93390000,(((((((bZebra_finch:56.38808600,bLanced_tailed_manakin:56.38808600):10.22391401,bBudgerigar:66.61200000):0.00000000,bGyrfalcon:66.61200000):13.44800000,bAnna_s_hummingbird:80.06000000):11.58900500,bChicken:91.64900500):169.72099500,rRed_eared_Slider:261.37000000):19.73000000,(rCentral_bearded_dragon:166.89372500,rSand_lizard:166.89372500):114.20627500):37.90000000):33.70449000,aEuropean_common_frog:352.70449000);"
tree <- newick2phylog(vertebrate_tree)
names <- c("pHuman_T2T","pHuman_hg38","pChimpanzee","pBonobo","pGorilla","pCommon_marmoset","pMa_s_night_monkey","pPhilippine_tarsier","mHouse_mouse","mNorway_rat","mNile_rat","mBottlenose_dolphin","mVaquita","mBlue_whale","mCow","mSheep","mRed_deer","mGreater_horseshoe_bat","mPale_spear_nosed_bat","mNine_banded_armadillo","mPlatypus","bZebra_finch","bLanced_tailed_manakin","bBudgerigar","bGyrfalcon","bAnna_s_hummingbird","bChicken","rRed_eared_Slider","rCentral_bearded_dragon","rSand_lizard","aEuropean_common_frog")


# Function for transparency
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}

# Read data seg
read.table('seg/priHomSapT2T.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg1
read.table('seg/priHomSap38.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg2
read.table('seg/priPanTro.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg3
read.table('seg/priPanPan.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg4
read.table('seg/priGorGor.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg5
read.table('seg/priCalJac.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg6
read.table('seg/priAotNan.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg7
read.table('seg/priCarSyr.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg8
read.table('seg/rodMusMus.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg9
read.table('seg/rodRatNor.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg10
read.table('seg/rodArvNil.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg11
read.table('seg/dolTurTru.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg12
read.table('seg/dolPhoSin.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg13
read.table('seg/whaBalMus.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg14
read.table('seg/rumBosTau.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg15
read.table('seg/rumOviAri.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg16
read.table('seg/rumCerEla.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg17
read.table('seg/batRhiFer.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg18
read.table('seg/batPhyDis.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg19
read.table('seg/mamDasNov.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg20
read.table('seg/mamOrnAna.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg21
read.table('seg/aviTaeGut.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg22
read.table('seg/aviChiLan.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg23
read.table('seg/aviMelUnd.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg24
read.table('seg/aviFalRus.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg25
read.table('seg/aviCalAnn.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg26
read.table('seg/aviGalGal.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg27
read.table('seg/turTraScr.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg28
read.table('seg/lizPogVit.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg29
read.table('seg/lizLacAgi.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg30
read.table('seg/ampRanTem.RBFOX1.seg',sep='\t',header=T) %>% dna_seg -> dna_seg31


dna_seg1[order(dna_seg1$start),] -> dna_seg1
dna_seg2[order(dna_seg2$start),] -> dna_seg2
dna_seg3[order(dna_seg3$start),] -> dna_seg3
dna_seg4[order(dna_seg4$start),] -> dna_seg4
dna_seg5[order(dna_seg5$start),] -> dna_seg5
dna_seg6[order(dna_seg6$start),] -> dna_seg6
dna_seg7[order(dna_seg7$start),] -> dna_seg7
dna_seg8[order(dna_seg8$start),] -> dna_seg8
dna_seg9[order(dna_seg9$start),] -> dna_seg9
dna_seg10[order(dna_seg10$start),] -> dna_seg10
dna_seg11[order(dna_seg11$start),] -> dna_seg11
dna_seg12[order(dna_seg12$start),] -> dna_seg12
dna_seg13[order(dna_seg13$start),] -> dna_seg13
dna_seg14[order(dna_seg14$start),] -> dna_seg14
dna_seg15[order(dna_seg15$start),] -> dna_seg15
dna_seg16[order(dna_seg16$start),] -> dna_seg16
dna_seg17[order(dna_seg17$start),] -> dna_seg17
dna_seg18[order(dna_seg18$start),] -> dna_seg18
dna_seg19[order(dna_seg19$start),] -> dna_seg19
dna_seg20[order(dna_seg20$start),] -> dna_seg20
dna_seg21[order(dna_seg21$start),] -> dna_seg21
dna_seg22[order(dna_seg22$start),] -> dna_seg22
dna_seg23[order(dna_seg23$start),] -> dna_seg23
dna_seg24[order(dna_seg24$start),] -> dna_seg24
dna_seg25[order(dna_seg25$start),] -> dna_seg25
dna_seg26[order(dna_seg26$start),] -> dna_seg26
dna_seg27[order(dna_seg27$start),] -> dna_seg27
dna_seg28[order(dna_seg28$start),] -> dna_seg28
dna_seg29[order(dna_seg29$start),] -> dna_seg29
dna_seg30[order(dna_seg30$start),] -> dna_seg30
dna_seg31[order(dna_seg31$start),] -> dna_seg31

dna_segs <- list(dna_seg1,
                 dna_seg2,
                 dna_seg3,
                 dna_seg4,
                 dna_seg5,
                 dna_seg6,
                 dna_seg7,
                 dna_seg8,
                 dna_seg9,
                 dna_seg10,
                 dna_seg11,
                 dna_seg12,
                 dna_seg13,
                 dna_seg14,
                 dna_seg15,
                 dna_seg16,
                 dna_seg17,
                 dna_seg18,
                 dna_seg19,
                 dna_seg20,
                 dna_seg21,
                 dna_seg22,
                 dna_seg23,
                 dna_seg24,
                 dna_seg25,
                 dna_seg26,
                 dna_seg27,
                 dna_seg28,
                 dna_seg29,
                 dna_seg30,
                 dna_seg31)

read.table('comp/comparison1.comp',sep='\t',header=T) %>% comparison -> comparison1
read.table('comp/comparison2.comp',sep='\t',header=T) %>% comparison -> comparison2
read.table('comp/comparison3.comp',sep='\t',header=T) %>% comparison -> comparison3
read.table('comp/comparison4.comp',sep='\t',header=T) %>% comparison -> comparison4
read.table('comp/comparison5.comp',sep='\t',header=T) %>% comparison -> comparison5
read.table('comp/comparison6.comp',sep='\t',header=T) %>% comparison -> comparison6
read.table('comp/comparison7.comp',sep='\t',header=T) %>% comparison -> comparison7
read.table('comp/comparison8.comp',sep='\t',header=T) %>% comparison -> comparison8
read.table('comp/comparison9.comp',sep='\t',header=T) %>% comparison -> comparison9
read.table('comp/comparison10.comp',sep='\t',header=T) %>% comparison -> comparison10
read.table('comp/comparison11.comp',sep='\t',header=T) %>% comparison -> comparison11
read.table('comp/comparison12.comp',sep='\t',header=T) %>% comparison -> comparison12
read.table('comp/comparison13.comp',sep='\t',header=T) %>% comparison -> comparison13
read.table('comp/comparison14.comp',sep='\t',header=T) %>% comparison -> comparison14
read.table('comp/comparison15.comp',sep='\t',header=T) %>% comparison -> comparison15
read.table('comp/comparison16.comp',sep='\t',header=T) %>% comparison -> comparison16
read.table('comp/comparison17.comp',sep='\t',header=T) %>% comparison -> comparison17
read.table('comp/comparison18.comp',sep='\t',header=T) %>% comparison -> comparison18
read.table('comp/comparison19.comp',sep='\t',header=T) %>% comparison -> comparison19
read.table('comp/comparison20.comp',sep='\t',header=T) %>% comparison -> comparison20
read.table('comp/comparison21.comp',sep='\t',header=T) %>% comparison -> comparison21
read.table('comp/comparison22.comp',sep='\t',header=T) %>% comparison -> comparison22
read.table('comp/comparison23.comp',sep='\t',header=T) %>% comparison -> comparison23
read.table('comp/comparison24.comp',sep='\t',header=T) %>% comparison -> comparison24
read.table('comp/comparison25.comp',sep='\t',header=T) %>% comparison -> comparison25
read.table('comp/comparison26.comp',sep='\t',header=T) %>% comparison -> comparison26
read.table('comp/comparison27.comp',sep='\t',header=T) %>% comparison -> comparison27
read.table('comp/comparison28.comp',sep='\t',header=T) %>% comparison -> comparison28
read.table('comp/comparison29.comp',sep='\t',header=T) %>% comparison -> comparison29
read.table('comp/comparison30.comp',sep='\t',header=T) %>% comparison -> comparison30

comparisons <- list(comparison1,
                    comparison2,
                    comparison3,
                    comparison4,
                    comparison5,
                    comparison6,
                    comparison7,
                    comparison8,
                    comparison9,
                    comparison10,
                    comparison11,
                    comparison12,
                    comparison13,
                    comparison14,
                    comparison15,
                    comparison16,
                    comparison17,
                    comparison18,
                    comparison19,
                    comparison20,
                    comparison21,
                    comparison22,
                    comparison23,
                    comparison24,
                    comparison25,
                    comparison26,
                    comparison27,
                    comparison28,
                    comparison29,
                    comparison30)

#B8DEE8

comparisons[[1]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[2]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[3]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[4]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[5]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[6]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[7]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[8]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[9]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[10]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[11]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[12]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[13]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[14]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[15]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[16]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[17]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[18]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[19]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[20]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[21]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[22]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[23]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[24]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[25]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[26]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[27]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[28]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[29]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)
comparisons[[30]]$col <- makeTransparent('#B8DEE8', alpha = 1.00)


names(dna_segs) <- names


left_offset = rep(c(0), times = length(dna_segs))

pdf("0.RBOXP1_alpha1.pdf",17,20)
plot_gene_map(dna_segs, offsets = left_offset,
              tree=tree, tree_width=2,
              comparisons=comparisons)
#dna_seg_scale=TRUE,
#scale=FALSE)
dev.off()

png(file="0.RBOXP1_alpha1.png", bg="transparent", width = 1600, height = 1600)
plot_gene_map(dna_segs, offsets = left_offset,
              tree=tree, tree_width=2,
              comparisons=comparisons)
dev.off()



# Primates
names_primate <- c("pHuman_T2T","pHuman_hg38","pChimpanzee","pBonobo","pGorilla","pCommon_marmoset","pMa_s_night_monkey","pPhilippine_tarsier")#,"mHouse_mouse","mNorway_rat","mNile_rat","mBottlenose_dolphin","mVaquita","mBlue_whale","mCow","mSheep","mRed_deer","mGreater_horseshoe_bat","mPale_spear_nosed_bat","mNine_banded_armadillo","mPlatypus","bZebra_finch","bLanced_tailed_manakin","bBudgerigar","bGyrfalcon","bAnna_s_hummingbird","bChicken","rRed_eared_Slider","rCentral_bearded_dragon","rSand_lizard","aEuropean_common_frog")

dna_segs <- list(dna_seg1,
                 dna_seg2,
                 dna_seg3,
                 dna_seg4,
                 dna_seg5,
                 dna_seg6,
                 dna_seg7,
                 dna_seg8)

names(dna_segs) <- names_primate

comparisons <- list(comparison1,
                    comparison2,
                    comparison3,
                    comparison4,
                    comparison5,
                    comparison6,
                    comparison7)

comparisons[[1]]$col <- makeTransparent('blue', alpha = 0.05)
comparisons[[2]]$col <- makeTransparent('blue', alpha = 0.05)
comparisons[[3]]$col <- makeTransparent('blue', alpha = 0.05)
comparisons[[4]]$col <- makeTransparent('blue', alpha = 0.05)
comparisons[[5]]$col <- makeTransparent('blue', alpha = 0.05)
comparisons[[6]]$col <- makeTransparent('blue', alpha = 0.05)
comparisons[[7]]$col <- makeTransparent('blue', alpha = 0.05)

left_offset = rep(c(0), times = length(dna_segs))

pdf("0.RBOXP1primate_alpha0.05.pdf",10,10)
plot_gene_map(dna_segs, offsets = left_offset,
              #tree=tree, tree_width=2,
              comparisons=comparisons)
#dna_seg_scale=TRUE,
#scale=FALSE)
dev.off()
