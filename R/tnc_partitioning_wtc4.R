# # Estimate the partitioning of the non-structural C among different organs using data from WTC-4 experiment 
# # on similar-sized seedlings of a related species (Eucalyptus parramattensis).
# 
# #  read the harvest WTC-4 TNC data
# tnc.wtc4 <- read.xls("raw_data/Samples_WTC4_heatwave_10.10.2017.xlsx",1)
# keeps = c("Tissue","Chamber","Whole.organ.total..g.")
# tnc.wtc4 = tnc.wtc4[ , keeps, drop = FALSE]
# tnc.wtc4 <- summaryBy(Whole.organ.total..g. ~ Tissue, data=tnc.wtc4, FUN=c(mean,standard.error))
# tnc.wtc4$Tissue = as.character(tnc.wtc4$Tissue)
# # tnc.wtc4[6:9,1] = c("foliage","wood","root","total")
# # tnc.wtc4[6,2:3] = tnc.wtc4[4,2:3] 
# # tnc.wtc4[7,2:3] = tnc.wtc4[1,2:3] + tnc.wtc4[5,2:3]
# # tnc.wtc4[8,2:3] = tnc.wtc4[2,2:3] + tnc.wtc4[3,2:3]
# # tnc.wtc4[9,2:3] = colSums(tnc.wtc4[6:8,2:3])
#   
#   
# tnc.partitioning = data.frame(organs=factor(c("foliage","wood","root","total")),
#                   ambient=numeric(4),
#                   warmed=numeric(4))
# 
# # There was no statistically significant difference in TNC partitioning across the treatments as John tested.
# tnc.partitioning$ambient[1:3] = c(245, 370, 110) # Data directly from John's analysis
# tnc.partitioning$warmed[1:3] = c(380, 650, 150) # Data directly from John's analysis
# 
# tnc.partitioning[4,2:3] = colSums(tnc.partitioning[2:3],na.rm = FALSE)
# tnc.partitioning$amb_ratio = tnc.partitioning[,2] / tnc.partitioning[4,2]
# tnc.partitioning$warm_ratio = tnc.partitioning[,3] / tnc.partitioning[4,3]
# 
# tnc.partitioning$average_ratio = (tnc.partitioning$amb_ratio + tnc.partitioning$warm_ratio) / 2
# 
# tnc.partitioning.ratio = data.frame(matrix(ncol = 4, nrow = 252))
# names(tnc.partitioning.ratio) = c("Date","foliage","wood","root")
# tnc.partitioning.ratio$Date = as.Date(as.Date("2013-09-17"):as.Date("2014-05-26"))
# 
# # Consider a linear change over time for TNC partitioning
# tnc.partitioning.ratio[1,2:4] <- c(0.75,0.16,0.09) # TNC partitioning according to Sink limited Pot experiment (Free Seedling)
# tnc.partitioning.ratio[nrow(tnc.partitioning.ratio),2:4] <- tnc.partitioning$average_ratio[1:3]
# 
# for (i in 2:ncol(tnc.partitioning.ratio)) {
#   tnc.partitioning.ratio[1:nrow(tnc.partitioning.ratio),i] = seq(tnc.partitioning.ratio[1,i], tnc.partitioning.ratio[nrow(tnc.partitioning.ratio),i], length.out = nrow(tnc.partitioning.ratio))
# } 
# 
# tnc.partitioning = tnc.partitioning.ratio


# Estimate the partitioning of the non-structural C among different organs using data from WTC-4 experiment
# on similar-sized seedlings of a related species (Eucalyptus parramattensis).

#  read the harvest WTC-4 TNC data
tnc.wtc4 <- read.xls("raw_data/Samples_WTC4_heatwave_10.10.2017.xlsx",1)
# wtc4.treatment = read.csv("raw_data/WTC_TEMP-PARRA_CM_TREE-HEIGHT-DIAMETER_20151028-20161124_L1.csv")
# names(tnc.wtc4)[3] = "chamber"
# tnc.wtc4 = merge(tnc.wtc4, wtc4.treatment[,c("chamber","T_treatment")], by="chamber")
keeps = c("Tissue","Chamber","Sugar..mg.g.","Starch..mg.g.")
tnc.wtc4 = tnc.wtc4[ , keeps, drop = FALSE]
tnc.wtc4$tnc..mg.g. = tnc.wtc4$Sugar..mg.g. + tnc.wtc4$Starch..mg.g.

# There was no statistically significant difference in TNC partitioning across the treatments as John tested.
tnc.wtc4 <- summaryBy(tnc..mg.g. ~ Tissue, data=tnc.wtc4, FUN=c(mean,standard.error))
tnc.wtc4$Tissue = as.character(tnc.wtc4$Tissue)
tnc.wtc4$tnc.ratio = tnc.wtc4$tnc..mg.g..mean / sum(tnc.wtc4$tnc..mg.g..mean) * 100
tnc.wtc4$tnc.ratio.SE = tnc.wtc4$tnc..mg.g..standard.error / tnc.wtc4$tnc..mg.g..mean * tnc.wtc4$tnc.ratio
tnc.wtc4[c(6:7),1] = as.character(c("wood", "root"))
tnc.wtc4[6,c(2:5)] = tnc.wtc4[2,c(2:5)] + tnc.wtc4[3,c(2:5)]
tnc.wtc4[7,c(2:5)] = tnc.wtc4[1,c(2:5)] + tnc.wtc4[5,c(2:5)]
tnc.wtc4 = tnc.wtc4[c(4,6,7),c(1,4,5)]



# TNC partitioning according to Duan's experiment (used in Pot experiment)
carbohydrates = read.csv("raw_data/Duan_carbohydrates.csv")
carbohydrates = subset(carbohydrates,CO2 == 400 & Temp == "Amb" & Water == "Well watered")
carbohydrates$tnc = carbohydrates$StarchW + carbohydrates$SolSugW # Unit = mg of tnc per g of dry weight biomass
carbohydrates$tnc = carbohydrates$tnc / 10 # Unit = % of dry weight biomass

# unit conversion from g of tnc per g of dry weight biomass to gC in tnc per gC in plant biomass
# 1 g of tnc has 0.4 gC and 1 g of dry weight biomass has 0.48 gC
carbohydrates$tnc = carbohydrates$tnc * (0.4/c1) # Unit = gC in tnc / gC in plant biomass
#-----------------------------------------------------------------------------------------
##### Total TNC calculation considering tree organ biomass partitioning
leaf.tnc = subset(carbohydrates,Organ == "Leaf") # Unit = % of dry weight leafmass
stem.tnc = subset(carbohydrates,Organ == "Stem") # Unit = % of dry weight stemmass
root.tnc = subset(carbohydrates,Organ == "Root") # Unit = % of dry weight rootmass

tnc.pot = data.frame(leaf.tnc$tnc,stem.tnc$tnc,root.tnc$tnc)
names(tnc.pot) <- c("leaf.tnc.C","stem.tnc.C","root.tnc.C") 

tnc.pot[7,] = colMeans(tnc.pot)
tnc.pot[8,] = tnc.pot[7,] / rowSums(tnc.pot[7,]) * 100
# tnc.pot = sapply(tnc.pot,function(x)sd(x)/sqrt(length(x)))


# Consider a linear change over time for TNC partitioning
tnc.partitioning = data.frame(matrix(ncol = 4, nrow = 252))
names(tnc.partitioning) = c("Date","foliage","wood","root")
tnc.partitioning$Date = as.Date(as.Date("2013-09-17"):as.Date("2014-05-26"))

tnc.partitioning[1,2:4] <- tnc.pot[8,]/100 # TNC partitioning according to Duan's experiment (used in Pot experiment)
tnc.partitioning[nrow(tnc.partitioning),2:4] <- tnc.wtc4[,2]/100

for (i in 2:ncol(tnc.partitioning)) {
  tnc.partitioning[1:nrow(tnc.partitioning),i] = seq(tnc.partitioning[1,i], tnc.partitioning[nrow(tnc.partitioning),i], length.out = nrow(tnc.partitioning))
}

