library(Biostrings)
library(seqinr)
library(ggseqlogo)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

codons <- read.csv('codon_table.txt', sep = '\t')
codons <- select(codons, X.Codon, Letter)

file_lines <- readLines("Pre/Sb-sequencing-Pre-MACS_S1_L001_R1_001_trimmed_assembled_BC1_Pre_MACS.trimmed.fasta")

df <- matrix(file_lines, ncol = 2, byrow = TRUE) %>%
  as.data.frame() %>%
  rename(sample = V1) %>%
  mutate(sample = str_remove(sample, ">")) %>%
  separate(V2, into = paste0("x", 1:1553), sep = 1:1553)

get_ref_positions <- select(df, x156:x294)
get_ref_positions <- get_ref_positions[1,]
pos <- colnames(get_ref_positions)[which(get_ref_positions!= '-')]

gdf_sub <- select(df, pos)
gdf_first_aa <- paste(gdf_sub$x156, gdf_sub$x158, gdf_sub$x167, sep = '')
gdf_second_aa <- paste(gdf_sub$x171, gdf_sub$x182, gdf_sub$x182, sep = '')
gdf_third_aa <- paste(gdf_sub$x236, gdf_sub$x243, gdf_sub$x250, sep = '')
gdf_fourth_aa <- paste(gdf_sub$x280,gdf_sub$x287,gdf_sub$x294, sep = '')

aas <- data.frame(gdf_first_aa, gdf_second_aa, gdf_third_aa, gdf_fourth_aa)
aas <- aas[-1,]

colnames(codons) <- c('X.Codon', 'First_AA')
aas <- merge(aas, codons, by.x = 'gdf_first_aa', by.y = 'X.Codon', all.x = T)
colnames(codons) <- c('X.Codon', 'Second_AA')
aas <- merge(aas, codons, by.x = 'gdf_second_aa', by.y = 'X.Codon', all.x = T)
colnames(codons) <- c('X.Codon', 'Third_AA')
aas <- merge(aas, codons, by.x = 'gdf_third_aa', by.y = 'X.Codon', all.x = T)
colnames(codons) <- c('X.Codon', 'Fourth_AA')
aas <- merge(aas, codons, by.x = 'gdf_fourth_aa', by.y = 'X.Codon', all.x = T)

aas$First_AA <- ifelse(is.na(aas$First_AA), '-', aas$First_AA)
aas$Second_AA <- ifelse(is.na(aas$Second_AA), '-', aas$Second_AA)
aas$Third_AA <- ifelse(is.na(aas$Third_AA), '-', aas$Third_AA)
aas$Fourth_AA <- ifelse(is.na(aas$Fourth_AA), '-', aas$Fourth_AA)

aas$combo <- paste(aas$First_AA, aas$Second_AA, aas$Third_AA, aas$Fourth_AA, sep = '')


pre_counts <- data.frame(table(aas$combo))

write.csv(pre_counts, file = 'Pre_Var_Pos_counts.csv', row.names = F, quote = F)

# Post 

post <- readDNAStringSet('Post/subset.fasta')
post_ref <- post[1] %>% as.character

gsub('-', '', substr(post_ref, 86,164)) %>% nchar()
file_lines <- readLines("Post/Sb-sequencing-Post-MACS_S2_L001_R1_001_trimmed_assembled_BC2_Post_MACS.trimmed.fasta")

df <- matrix(file_lines, ncol = 2, byrow = TRUE) %>%
  as.data.frame() %>%
  rename(sample = V1) %>%
  mutate(sample = str_remove(sample, ">")) %>%
  separate(V2, into = paste0("x", 1:512), sep = 1:512)

get_ref_positions <- select(df, x86:x164)
get_ref_positions <- get_ref_positions[1,]
pos <- colnames(get_ref_positions)[which(get_ref_positions!= '-')]

gdf_sub <- select(df, pos)
gdf_first_aa <- paste(gdf_sub$x86, gdf_sub$x89, gdf_sub$x93, sep = '')
gdf_second_aa <- paste(gdf_sub$x96, gdf_sub$x99, gdf_sub$x103, sep = '')
gdf_third_aa <- paste(gdf_sub$x131, gdf_sub$x134, gdf_sub$x139, sep = '')
gdf_fourth_aa <- paste(gdf_sub$x155,gdf_sub$x157,gdf_sub$x164, sep = '')

aas <- data.frame(gdf_first_aa, gdf_second_aa, gdf_third_aa, gdf_fourth_aa)
aas <- aas[-1,]

colnames(codons) <- c('X.Codon', 'First_AA')
aas <- merge(aas, codons, by.x = 'gdf_first_aa', by.y = 'X.Codon', all.x = T)
colnames(codons) <- c('X.Codon', 'Second_AA')
aas <- merge(aas, codons, by.x = 'gdf_second_aa', by.y = 'X.Codon', all.x = T)
colnames(codons) <- c('X.Codon', 'Third_AA')
aas <- merge(aas, codons, by.x = 'gdf_third_aa', by.y = 'X.Codon', all.x = T)
colnames(codons) <- c('X.Codon', 'Fourth_AA')
aas <- merge(aas, codons, by.x = 'gdf_fourth_aa', by.y = 'X.Codon', all.x = T)

aas$First_AA <- ifelse(is.na(aas$First_AA), '-', aas$First_AA)
aas$Second_AA <- ifelse(is.na(aas$Second_AA), '-', aas$Second_AA)
aas$Third_AA <- ifelse(is.na(aas$Third_AA), '-', aas$Third_AA)
aas$Fourth_AA <- ifelse(is.na(aas$Fourth_AA), '-', aas$Fourth_AA)

aas$combo <- paste(aas$First_AA, aas$Second_AA, aas$Third_AA, aas$Fourth_AA, sep = '')

post_counts <- data.frame(table(aas$combo))

write.csv(post_counts, file = 'Post_Var_Pos_counts.csv', row.names = F, quote = F)

## seq logos
pre_mat <- read.csv('Pre/get_counts/Sb_Pre.counts.matrix_out_aa.csv', header = T, stringsAsFactors = )
pre_mat <- pre_mat[-c(1,2),]
rownames(pre_mat) <- pre_mat$Reference.Position
pre_mat <- pre_mat[,-1]
pre_mat <- pre_mat[-nrow(pre_mat),]
rows <- rownames(pre_mat)
pre_mat <- apply(pre_mat, 2, as.numeric)

pre_mat2 <- as.matrix(data.frame(pre_mat))
rownames(pre_mat2) <- rows

pre_mat2_sub <- pre_mat2[1:20,3:20]

ggseqlogo( pre_mat2_sub, method = 'prob')

post_mat <- read.csv('Post/get_counts/Sb_Post.counts.matrix_out_aa.csv', header = T, stringsAsFactors = )
post_mat <- post_mat[-c(1,2),]
rownames(post_mat) <- post_mat$Reference.Position
post_mat <- post_mat[,-1]
post_mat <- post_mat[-nrow(post_mat),]
rows <- rownames(post_mat)
post_mat <- apply(post_mat, 2, as.numeric)

post_mat2 <- as.matrix(data.frame(post_mat))
rownames(post_mat2) <- rows

post_mat2_sub <- post_mat2[1:20,3:20]

ggseqlogo( post_mat2_sub, method = 'prob')

mats <- list(Pre = pre_mat2_sub, Post = post_mat2_sub)

pre_logo <- ggseqlogo(pre_mat2_sub, method = 'prob') + scale_x_continuous(breaks=1:18, labels = 186:203) + ggtitle('Pre') + 
   theme_classic() + theme(plot.title = element_text(hjust = .5))

post_logo <- ggseqlogo(post_mat2_sub, method = 'prob') + scale_x_continuous(breaks=1:18, labels = 186:203) + ggtitle('Post') + 
  theme_classic() + theme(plot.title = element_text(hjust = .5))

png('Pre_Post_Logo_Plots.png', height = 6, width = 10, units = 'in', res = 600)
ggarrange(pre_logo, post_logo, nrow=2, common.legend = TRUE, legend="none")
dev.off()
