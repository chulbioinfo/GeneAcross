#install.packages("dplyr")
#install.packages("genoPlotR")
#install.packages("magrittr")
# Library
library("magrittr")
library("dplyr")
library("genoPlotR")

# Global variable
setwd("C:\\Users\\swear0712\\Desktop\\RBFOX1\\")
#vertebrate_tree = "(((((((((((pHuman_T2T:0.001,pHuman_hg38:0.001):6.40000000,(pChimpanzee:2.39810000,pBonobo:2.39810000):4.00190000):2.20000000,pGorilla:8.60000000):34.32000000,(pCommon_marmoset:17.65000000,pMa_s_night_monkey:17.65000000):25.27000000):26.23223000,pPhilippine_tarsier:69.15223000):18.04777000,((mHouse_mouse:10.25000000,mNorway_rat:10.25000000):0.00000000,mNile_rat:10.25000000):76.95000000):6.80000000,((((mBottlenose_dolphin:17.63872400,mVaquita:17.63872400):16.66752600,mBlue_whale:34.30625000):23.48185000,((mCow:15.93000000,mSheep:15.93000000):7.00315000,mRed_deer:22.93315000):34.85495000):23.48106500,(mGreater_horseshoe_bat:62.00000000,mPale_spear_nosed_bat:62.00000000):19.26916500):12.73083500):4.90000000,mNine_banded_armadillo:98.90000000):81.16610000,mPlatypus:180.06610000):138.93390000,(((((((bZebra_finch:56.38808600,bLanced_tailed_manakin:56.38808600):10.22391401,bBudgerigar:66.61200000):0.00000000,bGyrfalcon:66.61200000):13.44800000,bAnna_s_hummingbird:80.06000000):11.58900500,bChicken:91.64900500):169.72099500,rRed_eared_Slider:261.37000000):19.73000000,(rCentral_bearded_dragon:166.89372500,rSand_lizard:166.89372500):114.20627500):37.90000000):33.70449000,aEuropean_common_frog:352.70449000);"
#tree <- newick2phylog(vertebrate_tree)
#names <- c("pHuman_T2T","pHuman_hg38","pChimpanzee","pBonobo","pGorilla","pCommon_marmoset","pMa_s_night_monkey","pPhilippine_tarsier","mHouse_mouse","mNorway_rat","mNile_rat","mBottlenose_dolphin","mVaquita","mBlue_whale","mCow","mSheep","mRed_deer","mGreater_horseshoe_bat","mPale_spear_nosed_bat","mNine_banded_armadillo","mPlatypus","bZebra_finch","bLanced_tailed_manakin","bBudgerigar","bGyrfalcon","bAnna_s_hummingbird","bChicken","rRed_eared_Slider","rCentral_bearded_dragon","rSand_lizard","aEuropean_common_frog")

primate_tree = "(((((pHuman_T2T:0.001,pHuman_hg38:0.001):6.40000000,(pChimpanzee:2.39810000,pBonobo:2.39810000):4.00190000):2.20000000,pGorilla:8.60000000):34.32000000,(pCommon_marmoset:17.65000000,pMa_s_night_monkey:17.65000000):25.27000000):26.23223000,pPhilippine_tarsier:69.15223000):18.04777000;"
tree <- newick2phylog(primate_tree)

names <- c("pHuman_T2T","pHuman_hg38","pChimpanzee","pBonobo","pGorilla","pCommon_marmoset","pMa_s_night_monkey","pPhilippine_tarsier")


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


dna_seg1[order(dna_seg1$start),] -> dna_seg1
dna_seg2[order(dna_seg2$start),] -> dna_seg2
dna_seg3[order(dna_seg3$start),] -> dna_seg3
dna_seg4[order(dna_seg4$start),] -> dna_seg4
dna_seg5[order(dna_seg5$start),] -> dna_seg5
dna_seg6[order(dna_seg6$start),] -> dna_seg6
dna_seg7[order(dna_seg7$start),] -> dna_seg7
dna_seg8[order(dna_seg8$start),] -> dna_seg8


dna_segs <- list(dna_seg1,
                 dna_seg2,
                 dna_seg3,
                 dna_seg4,
                 dna_seg5,
                 dna_seg6,
                 dna_seg7,
                 dna_seg8)

read.table('comp/comparison1.comp',sep='\t',header=T) %>% comparison -> comparison1
read.table('comp/comparison2.comp',sep='\t',header=T) %>% comparison -> comparison2
read.table('comp/comparison3.comp',sep='\t',header=T) %>% comparison -> comparison3
read.table('comp/comparison4.comp',sep='\t',header=T) %>% comparison -> comparison4
read.table('comp/comparison5.comp',sep='\t',header=T) %>% comparison -> comparison5
read.table('comp/comparison6.comp',sep='\t',header=T) %>% comparison -> comparison6
read.table('comp/comparison7.comp',sep='\t',header=T) %>% comparison -> comparison7


comparisons <- list(comparison1,
                    comparison2,
                    comparison3,
                    comparison4,
                    comparison5,
                    comparison6,
                    comparison7)

#B8DEE8
color_code = makeTransparent('grey', alpha = 0.05)

comparisons[[1]]$col <- color_code
comparisons[[2]]$col <- color_code
comparisons[[3]]$col <- color_code
comparisons[[4]]$col <- color_code
comparisons[[5]]$col <- color_code
comparisons[[6]]$col <- color_code
comparisons[[7]]$col <- color_code

names(dna_segs) <- names


left_offset = rep(c(0), times = length(dna_segs))

left_offset = rep(c(0), times = length(dna_segs))

pdf("0.RBOXP1_primates.pdf",5,8)
plot_gene_map(dna_segs, offsets = left_offset,
              comparisons=comparisons)
#dna_seg_scale=TRUE,
#scale=FALSE)
dev.off()
