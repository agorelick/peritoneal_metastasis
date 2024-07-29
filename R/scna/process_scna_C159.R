source(here::here('R/func.R'))
patient <- 'C159'
sex <- 'male'

fit_dir <- here(paste0('processed_data/copynumber/',patient,'/fits/'))
if(!dir.exists(fit_dir)) dir.create(fit_dir, recursive=T)

## because hacked to include X,Y, this is now as a list of 'obj', where each 'obj' has 1 sample
obj_list <- readRDS(here(paste0('original_data/scna/',patient,'_1000kbp_withXY.rds'))) 
names(obj_list) <- paste0('C159',gsub('_aligned','',names(obj_list)))
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

## Lun1-A, A1
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=3, purity=0.39, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Lun2c-A, A1w
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=3, purity=0.29, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Lun2a-A, A2a
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=3, purity=0.35, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Lun2b-A, A2b
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=3, purity=0.17, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Lun2d-A, A2w
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=3, purity=0.06, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## AG1b-A, AR1w
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=3, purity=0.54, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per1-A, Ab1
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=3, purity=0.71, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## AG1a-A, Ad1
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=3, purity=0.48, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## AG2-A, Ad2
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=3, purity=0.20, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per3-A, B1
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=3, purity=0.37, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per4-A, B3
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=3, purity=0.52, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Di3b-A, B3w
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=3, purity=0.45, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per5-A, B4
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=3, purity=0.59, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Pa1b-A, D1w
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=3, purity=0.60, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Di1-A, Di1
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=3, purity=0.75, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Di2-A, Di2
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=3, purity=0.56, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Di3a-A, Di3
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=3, purity=0.28, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Liv4b-A, H1w
refit(obj_list[[18]], samplename=samples[18], sex=sex, ploidy=3, purity=0.70, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Liv1-A, H3
refit(obj_list[[19]], samplename=samples[19], sex=sex, ploidy=3, purity=0.53, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Liv2-A, H6b
refit(obj_list[[20]], samplename=samples[20], sex=sex, ploidy=3, purity=0.59, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Liv3a-A, H7a
refit(obj_list[[21]], samplename=samples[21], sex=sex, ploidy=3, purity=0.60, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Liv3b-A, H7b
refit(obj_list[[22]], samplename=samples[22], sex=sex, ploidy=3, purity=0.69, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Liv4a-A, H8
refit(obj_list[[23]], samplename=samples[23], sex=sex, ploidy=3, purity=0.40, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld1-A, LD1
refit(obj_list[[24]], samplename=samples[24], sex=sex, ploidy=3, purity=0.63, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld3b-A, LD1w
refit(obj_list[[25]], samplename=samples[25], sex=sex, ploidy=3, purity=0.32, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

##  Ld2-A, LD2
refit(obj_list[[26]], samplename=samples[26], sex=sex, ploidy=3, purity=0.55, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld4b-A, LD2w
refit(obj_list[[27]], samplename=samples[27], sex=sex, ploidy=3, purity=0.14, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld3a-A, LD3
refit(obj_list[[28]], samplename=samples[28], sex=sex, ploidy=3, purity=0.43, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld5b-A, LD3w
refit(obj_list[[29]], samplename=samples[29], sex=sex, ploidy=3, purity=0.39, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld4a-A, LD4
refit(obj_list[[30]], samplename=samples[30], sex=sex, ploidy=3, purity=0.45, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Ld5a-A, LD5
refit(obj_list[[31]], samplename=samples[31], sex=sex, ploidy=3, purity=0.42, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN1, LN1
refit(obj_list[[32]], samplename=samples[32], sex=sex, ploidy=3, purity=0.51, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN7, LN13
refit(obj_list[[33]], samplename=samples[33], sex=sex, ploidy=3, purity=0.12, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN2, LN2
refit(obj_list[[34]], samplename=samples[34], sex=sex, ploidy=3, purity=0.45, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN3, LN4
refit(obj_list[[35]], samplename=samples[35], sex=sex, ploidy=3, purity=0.21, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN4, LN5
refit(obj_list[[36]], samplename=samples[36], sex=sex, ploidy=3, purity=0.77, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN5, LN6
refit(obj_list[[37]], samplename=samples[37], sex=sex, ploidy=3, purity=0.15, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Me1-A, Me1
refit(obj_list[[38]], samplename=samples[38], sex=sex, ploidy=3, purity=0.69, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Normal1, N1
refit(obj_list[[39]], samplename=samples[39], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT9-A, P10
refit(obj_list[[40]], samplename=samples[40], sex=sex, ploidy=3, purity=0.78, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT10-A, P1w
refit(obj_list[[41]], samplename=samples[41], sex=sex, ploidy=3, purity=0.31, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT1, P2
refit(obj_list[[42]], samplename=samples[42], sex=sex, ploidy=3, purity=0.60, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT11-A, P2w
refit(obj_list[[43]], samplename=samples[43], sex=sex, ploidy=3, purity=0.62, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT2, P3
refit(obj_list[[44]], samplename=samples[44], sex=sex, ploidy=3, purity=0.56, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT12-A, P3w
refit(obj_list[[45]], samplename=samples[45], sex=sex, ploidy=3, purity=0.66, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT3, P4
refit(obj_list[[46]], samplename=samples[46], sex=sex, ploidy=3, purity=0.70, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT13-A, P4w
refit(obj_list[[47]], samplename=samples[47], sex=sex, ploidy=3, purity=0.82, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT5, P5
refit(obj_list[[48]], samplename=samples[48], sex=sex, ploidy=3, purity=0.58, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT7, P6
refit(obj_list[[49]], samplename=samples[49], sex=sex, ploidy=3, purity=0.56, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT4, P7
refit(obj_list[[50]], samplename=samples[50], sex=sex, ploidy=3, purity=0.45, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT6, P8
refit(obj_list[[51]], samplename=samples[51], sex=sex, ploidy=3, purity=0.57, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT8-A, P9
refit(obj_list[[52]], samplename=samples[52], sex=sex, ploidy=3, purity=0.38, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Pa1a-A, Pa1
refit(obj_list[[53]], samplename=samples[53], sex=sex, ploidy=3, purity=0.60, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per2-A, Pg1
refit(obj_list[[54]], samplename=samples[54], sex=sex, ploidy=3, purity=0.39, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Sp1b-A, S1w
refit(obj_list[[55]], samplename=samples[55], sex=sex, ploidy=3, purity=0.45, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## SB1-A, SB1
refit(obj_list[[56]], samplename=samples[56], sex=sex, ploidy=3, purity=0.38, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Sp1a-A, Sp1
refit(obj_list[[57]], samplename=samples[57], sex=sex, ploidy=3, purity=0.42, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## TD1, TD1
refit(obj_list[[58]], samplename=samples[58], sex=sex, ploidy=3, purity=0.23, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## X1a-A, X1a
refit(obj_list[[59]], samplename=samples[59], sex=sex, ploidy=3, purity=0.37, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## X1b-A, X1b
refit(obj_list[[60]], samplename=samples[60], sex=sex, ploidy=3, purity=0.39, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## X1c-A, X1c
refit(obj_list[[61]], samplename=samples[61], sex=sex, ploidy=3, purity=0.55, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## X1d-A, X2w
refit(obj_list[[62]], samplename=samples[62], sex=sex, ploidy=3, purity=0.57, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))


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
ggsave(here(paste0('figures_and_tables/copynumber/heatmaps/',patient,'_cnv_segment_heatmap.pdf')),width=11,height=7)

## test distance matrix similarity
set.seed(42)
message('Testing distance matrix similarity ...')
p2 <- compare_matrices(info$dm, ad, patient, R=1e4)
ggsave(here(paste0('figures_and_tables/copynumber/distance_matrix_comparisons/',patient,'_cnv_bins_euclidean_matrix_comparison.pdf')),width=11,height=3,plot=p2)

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

## plot unrooted bootstrapped tree comparisons
tree_ad_conf <- addConfidences(tree_ad, bstrees)
tree_ad_conf$node.label <- round(100*tree_ad_conf$node.label)
p_bstree_ad_conf <- plot_bootstrapped_tree_unrooted(tree_ad_conf, groups, paste0(patient,' Poly-G'), size=4,bs_size=4)
scna_tree <- bs_info$tree
scna_tree <- addConfidences(scna_tree, bs_info$bstrees)
scna_tree$node.label <- round(100*scna_tree$node.label)
p_bstree_scna <- plot_bootstrapped_tree_unrooted(scna_tree, groups, paste0(patient,' SCNA'), size=4, bs_size=4)
p <- plot_grid(p_bstree_scna, p_bstree_ad_conf, res$plot, nrow=1)
ggsave(here(paste0('figures_and_tables/copynumber/tree_comparisons/',patient,'_cnv_bins_euclidean_nj_tree_comparison.pdf')),plot=p,width=11,height=3)

## plot bootstrapped AD tree
tree_ad_conf <- addConfidences(tree_ad_conf, bstrees)
tree_ad_conf$node.label <- round(100*tree_ad_conf$node.label)
p_bstree_ad_conf <- plot_bootstrapped_tree(tree_ad_conf, groups, paste0(patient,' Poly-G'))
scna_tree <- bs_info$tree
scna_tree <- addConfidences(scna_tree, bs_info$bstrees)
scna_tree$node.label <- round(100*scna_tree$node.label)
p_bstree_scna <- plot_bootstrapped_tree(scna_tree, groups, paste0(patient,' SCNA'))
p <- plot_grid(p_bstree_scna, p_bstree_ad_conf, ncol=2)
ggsave(here(paste0('figures_and_tables/copynumber/bootstrapped_trees/',patient,'_scna_polyg_trees_bootstrapped.pdf')),plot=p, width=11,height=8)



