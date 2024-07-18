source(here::here('R/func.R'))
patient <- 'C161'
sex <- 'male'

fit_dir <- here(paste0('processed_data/copynumber/',patient,'/fits/'))
if(!dir.exists(fit_dir)) dir.create(fit_dir, recursive=T)

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/scna/',patient,'_1000kbp_withXY.rds'))) 
names(obj_list) <- gsub('_aligned','',names(obj_list))
sample_info <- fread(here('processed_data/sample_info.txt'))
valid_samples <- intersect(sample_info$Sample_ID, names(obj_list))
map <- sample_info[Sample_ID %in% valid_samples,c('Sample_ID','Real_Sample_ID'),with=F]
map <- map[order(Sample_ID),]
obj_list <- obj_list[map$Sample_ID]
names(obj_list) <- map$Real_Sample_ID

## remove results from previous run
purity_file=here(paste0('processed_data/copynumber/',patient,'/fits/purity_ploidy.txt'))
if(file.exists(purity_file)) file.remove(purity_file)
samples <- names(obj_list)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fit purity and ploidy for each sample
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Lun1, A1
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=4, purity=0.29, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Lun2, A2
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=4, purity=0.18, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Lun3-A, A5
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=4, purity=0.18, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Plu1-A, A6
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=4, purity=0.34, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Di1-A, B1 
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=4, purity=0.20, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per3-A, B4 
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=4, purity=0.21, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per2-A, B6 
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=4, purity=0.17, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per4-A, B7
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=4, purity=0.16, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Liv1-A, H1
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=4, purity=0.23, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN2, L2
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=4, purity=0.44, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN3, L3
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=4, purity=0.23, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN4, L4
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=4, purity=0.31, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld1-A, Ld1 
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=4, purity=0.32, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld2-A, Ld2 
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=4, purity=0.17, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Normal1, N1
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ov1-A, O1
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=4, purity=0.22, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ov2-A, O2
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=4, purity=0.27, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT2, P2
refit(obj_list[[18]], samplename=samples[18], sex=sex, ploidy=4, purity=0.32, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT3, P3
refit(obj_list[[19]], samplename=samples[19], sex=sex, ploidy=4, purity=0.33, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT4, P4
refit(obj_list[[20]], samplename=samples[20], sex=sex, ploidy=4, purity=0.40, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT5, P5
refit(obj_list[[21]], samplename=samples[21], sex=sex, ploidy=4, purity=0.28, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT6, P6
refit(obj_list[[22]], samplename=samples[22], sex=sex, ploidy=4, purity=0.37, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT7, P7
refit(obj_list[[23]], samplename=samples[23], sex=sex, ploidy=4, purity=0.32, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT8-A, P8
refit(obj_list[[24]], samplename=samples[24], sex=sex, ploidy=4, purity=0.37, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT1, P9
refit(obj_list[[25]], samplename=samples[25], sex=sex, ploidy=4, purity=0.42, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Pa1-A, Pa1
refit(obj_list[[26]], samplename=samples[26], sex=sex, ploidy=4, purity=0.21, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Sp1-A, S1
refit(obj_list[[27]], samplename=samples[27], sex=sex, ploidy=4, purity=0.29, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## SB1-A, SB1
refit(obj_list[[28]], samplename=samples[28], sex=sex, ploidy=4, purity=0.35, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## SB2-A, SB2
refit(obj_list[[29]], samplename=samples[29], sex=sex, ploidy=4, purity=0.28, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## St1-A, St1
refit(obj_list[[30]], samplename=samples[30], sex=sex, ploidy=4, purity=0.20, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## St2-A, St2
refit(obj_list[[31]], samplename=samples[31], sex=sex, ploidy=4, purity=0.33, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process SCNA data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## load angular distance matrix
## subset AD matrices for the samples in common with the SCNA samples
ad <- read_distance_matrix(here(paste0('processed_data/angular_distance_matrices/',patient,'.txt')))
bs_ad <- readRDS(here(paste0('processed_data/angular_distance_matrices_bootstrapped/',patient,'.rds')))
fits <- fread(here(paste0('processed_data/copynumber/',patient,'/fits/purity_ploidy.txt')))
chr <- fread(here('original_data/misc/chr_lengths_b37.txt'))
groups <- sample_info[Patient_ID==patient,c('Real_Sample_ID','group'),with=F]
setnames(groups,'Real_Sample_ID','label')

heatmap_samples <- intersect(colnames(ad), samples)
ad <- ad[heatmap_samples, heatmap_samples]
tree_ad <- nj(ad)
tree_ad <- root(tree_ad,outgroup=grep(paste0('^N'),tree_ad$tip.label),resolve.root=TRUE) 

## process the 1Mb binned SCNA data
info <- process_scna_data_for_patient(heatmap_samples, obj_list, fits, sex)
tree_scna <- nj(info$dm)
tree_scna <- root(tree_scna,outgroup=grep(paste0('^N'),tree_scna$tip.label),resolve.root=TRUE) 
write_tsv(info$bins,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_bins.txt')))
write_tsv(info$segs,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_segments.txt')))
write_tsv(info$mat,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_matrix.txt')))
write_distance_matrix(info$dm,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_distance_matrix.txt')))

## generate SCNA heatmap
p <- scna_segment_heatmap(info$dm, info$segs, chr, sex, groups, group_cols)
ggsave(here(paste0('figures/copynumber/heatmaps/',patient,'_cnv_segment_heatmap.pdf')),width=11,height=4.25)

## test distance matrix similarity
set.seed(42)
message('Testing distance matrix similarity ...')
p2 <- compare_matrices(info$dm, ad, patient, R=1e4)
ggsave(here(paste0('figures/copynumber/distance_matrix_comparisons/',patient,'_cnv_bins_euclidean_matrix_comparison.pdf')),width=11,height=3)

## get bootstrapped SCNA tree
set.seed(42)
message('Resampling chromosomes for bootstrapped SCNA tree ...')
bs_info <- get_bootstrapped_scna_tree(info$bins, info$dm, patient)

## get bootstrapped AD matrices, subset them for the common samples, convert to multiphylo
for(i in 1:length(bs_ad)) bs_ad[[i]] <- bs_ad[[i]][heatmap_samples, heatmap_samples]
bstrees <- lapply(bs_ad, nj)
bstrees <- TreeTools::as.multiPhylo(bstrees)

## test tree similarity
set.seed(42)
message('Testing tree similarity ...')
res <- test_tree_similarity(info$dm, ad, nperm=1e4, title=patient)
write_tsv(res$data, here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_polyg_tree_similarity_test.txt')))

## plot a rotated version of the AD tree for the main figure
tree_ad_rotated <- ape::rotate(tree_ad, node=47)
p <- ggtree(tree_ad_rotated,layout='ape')
p <- p %<+% groups
p <- p + geom_tiplab(fontface=1,size=4,hjust=0,angle=0,aes(color=group))
p <- p + scale_color_manual(values=group_cols,name='Tissue')
p <- p + theme(legend.position='bottom') + guides(color='none')
p <- p + labs(title='C161 AD')
ggsave(here(paste0('figures/copynumber/tree_comparisons/',patient,'_angular_distance_unrooted_rotated.pdf')),plot=p)

## plot unrooted bootstrapped tree comparisons
tree_ad_conf <- addConfidences(tree_ad, bstrees)
tree_ad_conf$node.label <- round(100*tree_ad_conf$node.label)
tree_ad_conf <- scale_branch('Normal1', tree_ad_conf, sf=6)
p_bstree_ad_conf <- plot_bootstrapped_tree_unrooted(tree_ad_conf, groups, paste0(patient,' Poly-G'))
tree_scna <- addConfidences(tree_scna, bs_info$bstrees)
tree_scna$node.label <- round(100*tree_scna$node.label)
tree_scna <- scale_branch('Normal1', tree_scna, sf=5)
p_bstree_scna <- plot_bootstrapped_tree_unrooted(tree_scna, groups, paste0(patient,' SCNA'))
p <- plot_grid(p_bstree_scna, p_bstree_ad_conf, res$plot, nrow=1)
ggsave(here(paste0('figures/copynumber/tree_comparisons/',patient,'_cnv_bins_euclidean_nj_tree_comparison.pdf')),plot=p,width=20,height=8)

## plot bootstrapped AD tree
tree_ad_conf <- addConfidences(tree_ad_conf, bstrees)
tree_ad_conf$node.label <- round(100*tree_ad_conf$node.label)
p_bstree_ad_conf <- plot_bootstrapped_tree(tree_ad_conf, groups, paste0(patient,' Poly-G'))
tree_scna <- bs_info$tree
tree_scna <- addConfidences(tree_scna, bs_info$bstrees)
tree_scna$node.label <- round(100*tree_scna$node.label)
p_bstree_scna <- plot_bootstrapped_tree(tree_scna, groups, paste0(patient,' SCNA'))
p <- plot_grid(p_bstree_scna, p_bstree_ad_conf, ncol=2)
ggsave(here(paste0('figures/copynumber/bootstrapped_trees/',patient,'_scna_polyg_trees_bootstrapped.pdf')),plot=p, width=11,height=8)



