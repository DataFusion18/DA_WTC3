# Estimate the partitioning of the non-structural C among different organs using data from WTC-4 experiment 
# on similar-sized seedlings of a related species (Eucalyptus parramattensis).

#  read the harvest WTC-4 TNC data
tnc.wtc4 <- read.xls("raw_data/Samples_WTC4_heatwave_10.10.2017.xlsx",1)
keeps = c("Tissue","Chamber","Whole.organ.total..g.")
tnc.wtc4 = tnc.wtc4[ , keeps, drop = FALSE]
tnc.wtc4 <- summaryBy(Whole.organ.total..g. ~ Tissue, data=tnc.wtc4, FUN=c(mean,standard.error))
tnc.wtc4$Tissue = as.character(tnc.wtc4$Tissue)
# tnc.wtc4[6:9,1] = c("foliage","wood","root","total")
# tnc.wtc4[6,2:3] = tnc.wtc4[4,2:3] 
# tnc.wtc4[7,2:3] = tnc.wtc4[1,2:3] + tnc.wtc4[5,2:3]
# tnc.wtc4[8,2:3] = tnc.wtc4[2,2:3] + tnc.wtc4[3,2:3]
# tnc.wtc4[9,2:3] = colSums(tnc.wtc4[6:8,2:3])
  
  
tnc.partitioning = data.frame(organs=factor(c("foliage","wood","root","total")),
                  ambient=numeric(4),
                  warmed=numeric(4))

# There was no statistically significant difference in TNC partitioning across the treatments as John tested.
tnc.partitioning$ambient[1:3] = c(245, 370, 110) # Data directly from John's analysis
tnc.partitioning$warmed[1:3] = c(380, 650, 150) # Data directly from John's analysis

tnc.partitioning[4,2:3] = colSums(tnc.partitioning[2:3],na.rm = FALSE)
tnc.partitioning$amb_ratio = tnc.partitioning[,2] / tnc.partitioning[4,2]
tnc.partitioning$warm_ratio = tnc.partitioning[,3] / tnc.partitioning[4,3]

tnc.partitioning$average_ratio = (tnc.partitioning$amb_ratio + tnc.partitioning$warm_ratio) / 2

tnc.partitioning.ratio = data.frame(matrix(ncol = 4, nrow = 252))
names(tnc.partitioning.ratio) = c("Date","foliage","wood","root")
tnc.partitioning.ratio$Date = as.Date(as.Date("2013-09-17"):as.Date("2014-05-26"))

# Consider a linear change over time for TNC partitioning
tnc.partitioning.ratio[1,2:4] <- c(0.75,0.16,0.09) # TNC partitioning according to Sink limited Pot experiment (Free Seedling)
tnc.partitioning.ratio[nrow(tnc.partitioning.ratio),2:4] <- tnc.partitioning$average_ratio[1:3]

for (i in 2:ncol(tnc.partitioning.ratio)) {
  tnc.partitioning.ratio[1:nrow(tnc.partitioning.ratio),i] = seq(tnc.partitioning.ratio[1,i], tnc.partitioning.ratio[nrow(tnc.partitioning.ratio),i], length.out = nrow(tnc.partitioning.ratio))
} 

tnc.partitioning = tnc.partitioning.ratio
