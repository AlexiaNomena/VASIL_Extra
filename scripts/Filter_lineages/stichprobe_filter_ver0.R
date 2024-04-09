rm(list = ls())
setwd("Desktop/percent_study/")
threshold = 3 #percent

# load the daily lineage frequencies in order to filter for variants who pass
# the threshold in the time frame

data_variant_percentage = read.csv("data/Daily_Lineage_Freq.csv")
data_variant_max = as.data.frame(apply(data_variant_percentage,2,max))

colnames(data_variant_max) = c("max", "variants")

dim(data_variant_max)

data_variant_max$variants = rownames(data_variant_max)
data_filtered = data_variant_max[data_variant_max$max > threshold, ]

variants = data_filtered$variants
length(variants)


# load the stichproben file and delete all variants, that are below the threshold

data_stichprobe = read.csv("data/Stichprobe_RKI.tsv", sep = "\t")
head(data_stichprobe)

data_stichprobe_filtered = data_stichprobe[data_stichprobe$lineage %in% variants,]

date_max = max(data_stichprobe_filtered$date)
date_min = min(data_stichprobe_filtered$date)

date_min
date_max

# order of dates might be corrupted. If dates are not in order, VASIL crashes
data_stichprobe_filtered <- data_stichprobe_filtered[order(data_stichprobe_filtered$date),]

write.table(data_stichprobe_filtered, file=paste('data/Stichprobe_filtered_', threshold,'percent.tsv'), quote=FALSE, sep='\t')
#write.table(data_stichprobe_filtered, file=paste('data/Stichprobe_filtered_', threshold,'_percent.csv'), quote=FALSE, sep=',')

# adjust cases to timeline 
data_cases = read.csv("data/caseAscertainmentTable_reportedCasesRatio.csv")
data_cases_filtered = data_cases[(data_cases$date>=date_min)&(data_cases$date<=date_max),]
write.table(data_cases_filtered, file=paste('data/caseAscertainmentTable_reportedCasesRatio_filtered_',threshold, 'percent.csv'), quote=FALSE, sep=',')

