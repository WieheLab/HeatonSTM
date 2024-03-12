library(dplyr)
library(tidyr)
library(strinr)

pre <- read.csv('Pre_Var_Pos_counts.csv')

pre$pos1 <- substr(pre$Var1, 1,1)
sum(pre$Freq)

pos1_counts <- aggregate(pre$Freq,  list(pre$pos1),  FUN=sum)
 
pre_sub <- subset(pre, !grepl(x = pre$Var1, pattern = '-|O'))
pre_sub <- select(pre_sub, Var1, Freq)

colnames(pre_sub) <- c('AA_Sequence', 'Count')

pre_sub$Frequency <- pre_sub$Count/sum(pre_sub$Count)

pre_sub <- pre_sub[order(pre_sub$Frequency, decreasing = T),]

write.csv(pre_sub, file = 'Pre_VariablePositions_Clone_Counts.csv', row.names = F, quote = F)


pre_sub2 <- subset(pre, !grepl(x = pre$Var1, pattern = '-'))
pre_sub2 <- select(pre_sub2, Var1, Freq)

colnames(pre_sub2) <- c('AA_Sequence', 'Count')

pre_sub2$Frequency <- pre_sub2$Count/sum(pre_sub2$Count)

pre_sub2 <- pre_sub2[order(pre_sub2$Frequency, decreasing = T),]
pre_sub2$AA_Sequence <- gsub('O', '*', pre_sub2$AA_Sequence)

write.csv(pre_sub2, file = 'Pre_VariablePositions_Clone_Counts_with_Stop.csv', row.names = F, quote = F)


pre_sub3 <- select(pre, Var1, Freq)

colnames(pre_sub3) <- c('AA_Sequence', 'Count')

pre_sub3$Frequency <- pre_sub3$Count/sum(pre_sub3$Count)
pre_sub3 <- pre_sub3[order(pre_sub3$Frequency, decreasing = T),]

write.csv(pre_sub3, file = 'Pre_VariablePositions_Clone_Counts_with_Stop_Gaps.csv', row.names = F, quote = F)


# POST
setwd('/Volumes/mb488/CIVICS/Nick_Heaton/Sb_Sequencing/Post/')

post <- read.csv('Post_Var_Pos_counts.csv')

post$pos1 <- substr(post$Var1, 1,1)
sum(post$Freq)

pos1_counts <- aggregate(post$Freq,  list(post$pos1),  FUN=sum)

post_sub <- subset(post, !grepl(x = post$Var1, pattern = '-|O'))
post_sub <- select(post_sub, Var1, Freq)

colnames(post_sub) <- c('AA_Sequence', 'Count')

post_sub$Frequency <- post_sub$Count/sum(post_sub$Count)

post_sub <- post_sub[order(post_sub$Frequency, decreasing = T),]

write.csv(post_sub, file = 'Post_VariablePositions_Clone_Counts.csv', row.names = F, quote = F)


post_sub2 <- subset(post, !grepl(x = post$Var1, pattern = '-'))
post_sub2 <- select(post_sub2, Var1, Freq)

colnames(post_sub2) <- c('AA_Sequence', 'Count')

post_sub2$Frequency <- post_sub2$Count/sum(post_sub2$Count)

post_sub2 <- post_sub2[order(post_sub2$Frequency, decreasing = T),]
post_sub2$AA_Sequence <- gsub('O', '*', post_sub2$AA_Sequence)

write.csv(post_sub2, file = 'Post_VariablePositions_Clone_Counts_with_Stop.csv', row.names = F, quote = F)


post_sub3 <- select(post, Var1, Freq)

colnames(post_sub3) <- c('AA_Sequence', 'Count')

post_sub3$Frequency <- post_sub3$Count/sum(post_sub3$Count)
post_sub3 <- post_sub3[order(post_sub3$Frequency, decreasing = T),]

write.csv(post_sub3, file = 'Post_VariablePositions_Clone_Counts_with_Stop_Gaps.csv', row.names = F, quote = F)





