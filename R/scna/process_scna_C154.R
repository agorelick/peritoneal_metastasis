source(here::here('R/func.R'))
patient <- 'C154'
sex <- 'female'

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

## Per1
refit(obj_list[[1]], samplename=samples[1], sex=sex, ploidy=2, purity=0.58, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## TD1
refit(obj_list[[2]], samplename=samples[2], sex=sex, ploidy=2, purity=0.5, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per2
refit(obj_list[[3]], samplename=samples[3], sex=sex, ploidy=2, purity=0.42, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## TD2
refit(obj_list[[4]], samplename=samples[4], sex=sex, ploidy=2, purity=0.51, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## Per3
refit(obj_list[[5]], samplename=samples[5], sex=sex, ploidy=2, purity=0.4, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## TD3
refit(obj_list[[6]], samplename=samples[6], sex=sex, ploidy=2, purity=0.65, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## TD4
refit(obj_list[[7]], samplename=samples[7], sex=sex, ploidy=2, purity=0.52, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN1
refit(obj_list[[8]], samplename=samples[8], sex=sex, ploidy=2, purity=0.57, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN2a
refit(obj_list[[9]], samplename=samples[9], sex=sex, ploidy=2, purity=0.66, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN2b
refit(obj_list[[10]], samplename=samples[10], sex=sex, ploidy=2, purity=0.30, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## N1
refit(obj_list[[11]], samplename=samples[11], sex=sex, ploidy=2, purity=1, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT2
refit(obj_list[[12]], samplename=samples[12], sex=sex, ploidy=2, purity=0.25, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT3
refit(obj_list[[13]], samplename=samples[13], sex=sex, ploidy=2, purity=0.37, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT4
refit(obj_list[[14]], samplename=samples[14], sex=sex, ploidy=2, purity=0.35, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT7
refit(obj_list[[15]], samplename=samples[15], sex=sex, ploidy=2, purity=0.40, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT8
refit(obj_list[[16]], samplename=samples[16], sex=sex, ploidy=2, purity=0.45, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## PT9
refit(obj_list[[17]], samplename=samples[17], sex=sex, ploidy=2, purity=0.32, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))

## LN5
refit(obj_list[[18]], samplename=samples[18], sex=sex, ploidy=2, purity=0.64, save=T, output_dir=here(paste0('processed_data/copynumber/',patient)))



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
tree_scna <- ape::rotate(tree_scna, 20)
write_tsv(info$bins,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_bins.txt')))
write_tsv(info$segs,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_segments.txt')))
write_tsv(info$mat,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_matrix.txt')))
write_distance_matrix(info$dm,here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_distance_matrix.txt')))

## generate SCNA heatmap
p <- scna_segment_heatmap(info$dm, info$segs, chr, sex, groups, group_cols)
ggsave(here(paste0('figures_and_tables/copynumber/heatmaps/',patient,'_cnv_segment_heatmap.pdf')),width=11,height=4.25)

## test distance matrix similarity
set.seed(42)
message('Testing distance matrix similarity ...')
p2 <- compare_matrices(info$dm, ad, patient, R=1e4)
ggsave(here(paste0('figures_and_tables/copynumber/distance_matrix_comparisons/',patient,'_cnv_bins_euclidean_matrix_comparison.pdf')),width=11,height=3)

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
tree_ad_conf$edge.length[tree_ad_conf$edge.length > 0.2] <- 0.2
tree_ad_conf$node.label <- round(100*tree_ad_conf$node.label)
p_bstree_ad_conf <- plot_bootstrapped_tree_unrooted(tree_ad_conf, groups, paste0(patient,' Poly-G'))
scna_tree <- bs_info$tree
scna_tree$edge.length[scna_tree$edge.length > 3] <- 3
scna_tree <- addConfidences(scna_tree, bs_info$bstrees)
scna_tree$node.label <- round(100*scna_tree$node.label)
scna_tree <- ape::rotate(scna_tree, 20)
p_bstree_scna <- plot_bootstrapped_tree_unrooted(scna_tree, groups, paste0(patient,' SCNA'))
p <- plot_grid(p_bstree_scna, p_bstree_ad_conf, res$plot, nrow=1)
ggsave(here(paste0('figures_and_tables/copynumber/tree_comparisons/',patient,'_cnv_bins_euclidean_nj_tree_comparison.pdf')),plot=p,width=11,height=3)

## plot bootstrapped AD tree
tree_ad_conf <- addConfidences(tree_ad, bstrees)
tree_ad_conf$node.label <- round(100*tree_ad_conf$node.label)
p_bstree_ad_conf <- plot_bootstrapped_tree(tree_ad_conf, groups, paste0(patient,' Poly-G'))
p_bstree_ad_conf <- ggtree::rotate(p_bstree_ad_conf, 23)
p_bstree_ad_conf <- ggtree::rotate(p_bstree_ad_conf, 25)
scna_tree <- bs_info$tree
scna_tree <- addConfidences(scna_tree, bs_info$bstrees)
scna_tree$node.label <- round(100*scna_tree$node.label)
p_bstree_scna <- plot_bootstrapped_tree(scna_tree, groups, paste0(patient,' SCNA'))
p <- plot_grid(p_bstree_scna, p_bstree_ad_conf, ncol=2)
ggsave(here(paste0('figures_and_tables/copynumber/bootstrapped_trees/',patient,'_scna_polyg_trees_bootstrapped.pdf')),plot=p, width=11,height=8)





