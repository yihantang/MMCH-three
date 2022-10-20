a = all_merge_filter_filter_filter_extract_all[all_merge_filter_filter_filter_extract_all$WES_alt_ADF>0 &all_merge_filter_filter_filter_extract_all$WES_alt_ADR>0,]

group_pos <- GRanges(seqnames = a$CHROM, ranges = IRanges(start = a$POS, end = a$POS))
targetedoverlap = as.data.frame(findOverlapPairs(group_pos, Targeted_bed))
names(targetedoverlap)[1] = "CHROM"
names(targetedoverlap)[2] = "POS"
names(targetedoverlap)[7] = "Targeted_bed_start"
names(targetedoverlap)[8] = "Targeted_bed_end"
names(targetedoverlap)[11] = "Targeted_gene"
group1_Targeted_targeted <- targetedoverlap[c(1,2,7,8,11)]

b = a %>%
  left_join(.,group1_Targeted_targeted,by = c("CHROM","POS"))
c<- add_column(b,Targeted_normalized, .after = "Targeted_AF")

#############
#add read info
all_reads_rmdupliccate <- read.delim("/Volumes/Lab_Gillis/YiHanTang/MMCH_three/192samples/use_position/reads/all_reads_rmdupliccate.txt")

read1 = all_reads_rmdupliccate 
read1$name = paste(read1[,207],read1[,208],read1[,209],read1[,210],read1[,211],sep = "_")
c$name = paste(c[,14],c[,1],c[,2],c[,3],c[,4],sep = "_")

###calculate ref-alt median
Ref_ALt = read1[,1:101] - read1[,104:204]

apply(Ref_ALt[1,1:10],1,mean)
apply(Ref_ALt[1,1:10],1,median)

x <- as.data.frame(apply(Ref_ALt[,1:10],1,median))

x2 <- as.data.frame(apply(Ref_ALt[,11:20],1,median))
x3 <- as.data.frame(apply(Ref_ALt[,21:30],1,median))
x4 <- as.data.frame(apply(Ref_ALt[,31:40],1,median))
x5 <- as.data.frame(apply(Ref_ALt[,41:50],1,median))
x6 <- as.data.frame(apply(Ref_ALt[,51:60],1,median))
x7 <- as.data.frame(apply(Ref_ALt[,61:70],1,median))
x8 <- as.data.frame(apply(Ref_ALt[,71:80],1,median))
x9 <- as.data.frame(apply(Ref_ALt[,81:90],1,median))
x10 <- as.data.frame(apply(Ref_ALt[,91:101],1,median))

combine = cbind(x,x2,x3,x4,x5,x6,x7,x8,x9,x10)
combine$name =read1$name
names(combine)[1] = "median1"
names(combine)[2] = "median2"
names(combine)[3] = "median3"
names(combine)[4] = "median4"
names(combine)[5] = "median5"
names(combine)[6] = "median6"
names(combine)[7] = "median7"
names(combine)[8] = "median8"
names(combine)[9] = "median9"
names(combine)[10] = "median10"
