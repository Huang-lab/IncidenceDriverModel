##### calc_drivers_select_on_pancan23.R #####
# Kuan-lin Huang @ MSSM

# Load pancan23 dataset
load('/Users/khuang/Box\ Sync/Huang_lab/Huang_lab_data/pancan23_driver_dataset/pancan23_gam.RData')
# data source: https://www.cell.com/cancer-cell/fulltext/S1535-6108(17)30261-1

# compile into a 
pancan23_d = data.frame(matrix(unlist(pancan23$gam), nrow=6456, byrow=T))
row.names(pancan23_d) = row.names(pancan23$gam)
colnames(pancan23_d) = colnames(pancan23$gam)
pancan23_d[1:5,1:5] 

summary(rowSums(pancan23_d)) # driver per sample distribution, agree with manuscript 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  
# 0.000   1.000   3.000   5.809   6.000 363.000 

driver_count = data.frame(rowSums(pancan23_d))
write.table(driver_count,quote=F, sep = '\t', row.names = T,col.names = F, file="out/pancan23_select_driver_count_by_sample.tsv")

