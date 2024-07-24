source(here::here('R/func.R'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make annotated poly-G trees for each patient
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sample_info <- fread(here('processed_data/sample_info.txt'))
sample_info <- sample_info[cohort!='lung',]
valid_patients <- unique(sample_info[cohort %in% c('science','natgen','peritoneal') & grepl('E[a-c]3',Patient_ID)==F & grepl('CRC',Patient_ID)==F,(Patient_ID)])
valid_patients <- sort(c(valid_patients,'E3'))
trash <- lapply(valid_patients, make_tree, sample_info=sample_info, collapsed=F, show.depth=T, show.timing=T, outdir=here('figures/polyg_phylogenies'), show.bsvals=T)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig1a
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sample_info <- fread(here('processed_data/sample_info.txt'))
sample_info <- sample_info[cohort=='peritoneal' | Patient_ID %in% c('C38','C89'),] 

## merge the E3 multi-primary data into a single patient
ea3 <- sample_info[Patient_ID=='Ea3',]
eb3 <- sample_info[Patient_ID=='Eb3',]
ec3 <- sample_info[Patient_ID=='Ec3',]
E3 <- rbind(ea3, eb3, ec3)
E3$Patient_ID <- 'E3'
E3 <- E3[!duplicated(Sample_ID),]
sample_info <- sample_info[!Patient_ID %in% c('Ea3','Eb3','Ec3'),]
sample_info <- rbind(sample_info, E3)

tabulate <- function(info) {
    lesions <- sum(info$in_collapsed==T)
    samples <- nrow(info)
    list(lesions=lesions, samples=samples)
}
tbl <- sample_info[,tabulate(.SD),by=c('Patient_ID','group','tissue_type')]

f=function(id) {
    s <- strsplit(id,'')[[1]]
    s <- as.integer(s)
    s <- s[!is.na(s)]
    out <- as.integer(paste(s,collapse=''))
    list(id=id, num=out)
}
l <- lapply(unique(tbl$Patient_ID), f)
info <- rbindlist(l)
info[grep('C',id),group:='C']
info[grep('E',id),group:='E']
info <- info[order(group,num),]

tbl_samples <- data.table::dcast(Patient_ID ~ tissue_type, value.var='samples', data=tbl)
tbl_samples[is.na(tbl_samples)] <- 0
tbl_samples <- as.data.table(reshape2::melt(tbl_samples,id.var='Patient_ID'))
tbl_samples[!is.na(value) & value > 0, label:=value]
tbl_samples <- merge(tbl_samples, sample_info[!duplicated(tissue_type),c('tissue_type','group'),with=F], by.x='variable', by.y='tissue_type',all.x=T)
tbl_samples[value==0, group:=NA]
tbl_samples$variable <- factor(tbl_samples$variable, levels=rev(c('Normal','Primary','Lymph node','Tumor deposit','Peritoneum','Lung','Liver','Ovary (hematogenous)','Lymph node (distant)')))
tbl_samples$Patient_ID <- factor(tbl_samples$Patient_ID, levels=(info$id))

tbl_lesions <- data.table::dcast(Patient_ID ~ tissue_type, value.var='lesions', data=tbl)
tbl_lesions[is.na(tbl_lesions)] <- 0
tbl_lesions <- as.data.table(reshape2::melt(tbl_lesions,id.var='Patient_ID'))
tbl_lesions[!is.na(value) & value > 0, label:=value]
tbl_lesions <- merge(tbl_lesions, sample_info[!duplicated(tissue_type),c('tissue_type','group'),with=F], by.x='variable', by.y='tissue_type',all.x=T)
tbl_lesions[value==0, group:=NA]
tbl_lesions$variable <- factor(tbl_lesions$variable, levels=rev(c('Normal','Primary','Lymph node','Tumor deposit','Peritoneum','Lung','Liver','Ovary (hematogenous)','Lymph node (distant)')))
tbl_lesions$Patient_ID <- factor(tbl_lesions$Patient_ID, levels=(info$id))

p_heatmap_lesions <- ggplot(tbl_lesions[!is.na(group)], aes(x=Patient_ID, y=variable)) + 
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0), position='left') +
    geom_tile(aes(alpha=value/max(tbl_samples$value), fill=group)) +
    scale_alpha_continuous(range=c(0.1,1)) + 
    scale_fill_manual(values=group_cols,'Tissue',na.value='white') + 
    geom_text(aes(label=label)) +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Lesions') +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),legend.position='none')

p_heatmap_samples <- ggplot(tbl_samples[!is.na(group)], aes(x=Patient_ID, y=variable)) + 
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0), position='left') +
    geom_tile(aes(alpha=value/max(tbl_samples$value), fill=group)) +
    scale_alpha_continuous(range=c(0.1,1)) + 
    scale_fill_manual(values=group_cols,'Tissue',na.value='white') + 
    geom_text(aes(label=label)) +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Samples') +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),legend.position='none')

sample_info[met_timing %in% c('metachronous','metachronous after synchronous'), met_timing:='metachronous']
sample_info[group %in% c('Lung','Liver','Distant (other)'), group:='Distant (any)']
cols <- group_cols[c('Locoregional','Peritoneum','Distant (other)')]
names(cols)[3] <- 'Distant (any)'
tbl_timed <- sample_info[met_timing %in% c('synchronous','metachronous'),tabulate(.SD),by=c('Patient_ID','group','tissue_type','met_timing')]
tbl_timed[met_timing=='synchronous', lesions:=-1*lesions]
tbl_timed$Patient_ID <- factor(tbl_timed$Patient_ID, levels=(info$id))

p_timing <- ggplot(tbl_timed, aes(x=Patient_ID, y=lesions)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(limits=c(-35,20),breaks=seq(-30, 20, by=10)) +
    geom_bar(stat='identity', aes(fill=group)) +
    geom_hline(yintercept=0, color='black', linewidth=0.5) +
    scale_fill_manual(values=cols,'Tissue',na.value='white') +
    theme_ang(base_size=12) + guides(fill='none') +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
    labs(y='Syn | Meta', x=NULL)

## get number of lesions and samples per patient, group, and tissue type
vertical_tabulate <- function(info) {
    deep <- sum(info$vertical=='deep')
    luminal <- sum(info$vertical=='mucosal/luminal')
    list(deep=deep, luminal=luminal)
}
depth <- sample_info[group=='Primary',vertical_tabulate(.SD),by=c('Patient_ID')]
depth[,deep := -1*deep]
depth <- data.table::melt(depth, id.var=c('Patient_ID'))
depth[variable=='deep',variable:='Deep']
depth[variable=='luminal',variable:='Mucosal/luminal']
depth_cols <- c('#59b27a','#007639'); names(depth_cols) <- c('Mucosal/luminal','Deep')
depth$variable <- factor(depth$variable, levels=c('Mucosal/luminal','Deep'))
depth$Patient_ID <- factor(depth$Patient_ID, levels=info$id)
p_depth <- ggplot(depth, aes(x=Patient_ID, y=value)) +
    scale_x_discrete(expand=c(0,0)) +
    geom_bar(stat='identity',aes(fill=variable)) +
    geom_hline(yintercept=0, color='black', linewidth=0.5) +
    scale_fill_manual(values=depth_cols,'PT depth',na.value='white') +
    theme_ang(base_size=12) + guides(fill='none') +
    labs(y='Deep     Luminal',x=NULL) +
    theme( axis.text.x=element_blank(),axis.line.x=element_blank(),axis.ticks.x=element_blank()) 

patient_info <- fread(here('original_data/misc/patient_info.txt'))
E3 <- patient_info[PATIENT_ID %in% c('Ea3','Eb3','Ec3')]
E3 <- E3[1,]
E3$PATIENT_ID <- 'E3'
E3$Location <- 'Sig/Desc/Trans'
patient_info <- patient_info[!PATIENT_ID %in% c('Ea3','Eb3','Ec3')]
patient_info <- rbind(patient_info, E3)
patient_info <- patient_info[PATIENT_ID %in% info$id]
patient_info[grepl('^y',Stage)==T,stage_after_tx:='Yes']
patient_info[grepl('^y',Stage)==F,stage_after_tx:='No']
patient_info[is.na(Mstage),Mstage:=0]
patient_info[,Tstage:=paste0('T',Tstage)]
patient_info[,Nstage:=paste0('N',Nstage)]
patient_info[,Mstage:=paste0('M',Mstage)]
patient_info$PATIENT_ID <- factor(patient_info$PATIENT_ID, levels=info$id)
patient_info[MSI_status=='', MSI_status:=NA]
patient_info[PATIENT_ID=='E3', Location:='Sig/Desc/Trans']
patient_info$Location <- factor(patient_info$Location, levels=c('Appendix','Caecum','Ascending','Transverse','Descending','Sigmoid','Rectum','Sig/Desc/Trans'))
location_cols <- c('#971E13','#EA5F59','#F3BBB9','#EDEDED','#BAD0E2','#6694BF','#316BA6','#bfbfbf')
names(location_cols) <- levels(patient_info$Location)

tmp <- patient_info[,c('PATIENT_ID','Age'),with=F]
tmp$Age <- as.integer(tmp$Age)
tmp$notAge <- as.integer(80-tmp$Age)
tmp <- data.table::melt(tmp, id.var='PATIENT_ID')
tmp$variable <- factor(tmp$variable, levels=c('notAge','Age'))

agecols <- c('black','#bfbfbf'); names(agecols) <- c('Age','notAge')
p_age <- ggplot(tmp, aes(x=PATIENT_ID, y=value)) +
    scale_y_continuous(breaks=seq(0,80,by=40),limits=c(0,80),expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    geom_bar(stat='identity',aes(fill=variable)) +
    scale_fill_manual(values=agecols) + 
    theme_ang(base_size=12) +
    guides(fill='none') +
    theme(axis.text.x=element_blank(),axis.line.x=element_blank(),axis.ticks.x=element_blank()) +
    labs(y='Age',x=NULL) 

p_sex <- ggplot(patient_info, aes(x=PATIENT_ID, y=1)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    geom_tile(aes(fill=Sex)) + 
    scale_fill_brewer(palette='Set1',name='Sex') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), 
          axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank()) + 
    labs(y='Sex',x=NULL) 

p_location <- ggplot(patient_info, aes(x=PATIENT_ID, y=1)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    geom_tile(aes(fill=Location), color='black', linewidth=0.25, width=0.8, height=0.8) + 
    scale_fill_manual(values=location_cols,name='Location') +
    theme_ang(base_size=12) +
    theme(
          axis.text.x=element_blank(),axis.line.x=element_blank(),axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank(), legend.position='none')+ 
    labs(y='Location',x=NULL) 

stage <- data.table::melt(patient_info[,c('PATIENT_ID','Tstage','Nstage','Mstage'),with=F], id.vars='PATIENT_ID')
stage$variable <- gsub('stage','',stage$variable)
stage$variable <- factor(stage$variable, levels=c('M','N','T'))
Tcols <- brewer.pal(3,'Greens')[c(2,3)]
names(Tcols) <- c('T3','T4')
Ncols <- brewer.pal(3,'Reds')
names(Ncols) <- c('N0','N1','N2')
Mcols <- brewer.pal(3,'Blues')[c(2,3)]
names(Mcols) <- c('M0','M1')
stage_cols <- c(Tcols, Ncols, Mcols)

p_stage <- ggplot(stage, aes(x=PATIENT_ID, y=variable)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    geom_tile(aes(fill=value),color='white',linewidth=0.25) + 
    scale_fill_manual(values=stage_cols,name='Stage') +
    theme_ang(base_size=12) +
    theme( axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(), legend.position='none') + 
    labs(x=NULL,y='TNM\nStage') 

p_msi <- ggplot(patient_info, aes(x=PATIENT_ID, y=1)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    geom_tile(aes(fill=MSI_status)) + 
    scale_fill_brewer(palette='Pastel1',name='MSI Status') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
          axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank()) + 
    labs(y='MSI Status',x=NULL) 

p_cohort <- ggplot(patient_info, aes(x=PATIENT_ID, y=1)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    geom_tile(aes(fill=Cohort)) + 
    scale_fill_brewer(palette='Pastel2',name='Cohort') +
    theme_ang(base_size=12) +
    theme(
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          axis.text.y=element_blank(),axis.line.y=element_blank(),axis.ticks.y=element_blank()) + 
    labs(y='Cohort',x=NULL) 

p <- plot_grid(p_heatmap_lesions, p_heatmap_samples, p_stage, p_location, p_timing, ncol=1, align='v', axis='lr', rel_heights=c(2,2,1,0.5,2))
ggsave(here('figures/fig_1a.pdf'),width=9, height=9)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 1c-d. C161 SCNA tree
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sample_info <- fread(here('processed_data/sample_info.txt'))
groups <- sample_info[Patient_ID=='C161',c('Real_Sample_ID','group'),with=F]
setnames(groups,'Real_Sample_ID','label')

## get polyG/SCNA distance matrices for trees
dm_scna <- read_distance_matrix(here(paste0('processed_data/copynumber/C161/C161_cnv_distance_matrix.txt')))
dm_polyg <- read_distance_matrix(here('processed_data/angular_distance_matrices/C161.txt'))
common_samples <- intersect(colnames(dm_scna), colnames(dm_polyg))
dm_polyg <- dm_polyg[common_samples, common_samples]
dm_scna <- dm_scna[common_samples, common_samples]

## SCNA tree
tree_scna <- nj(dm_scna)
tree_scna <- root(tree_scna,outgroup=grep(paste0('^N'),tree_scna$tip.label),resolve.root=TRUE) 
tree_scna$edge.length[60] <- 0.35*tree_scna$edge.length[60] # truncate the normal branch
p1 <- ggtree(tree_scna,layout='ape')
p1 <- p1 %<+% groups
p1 <- p1 + geom_tiplab(fontface=1,size=3,hjust=0,angle=0,aes(color=group))
p1 <- p1 + scale_color_manual(values=group_cols,name='Tissue')
p1 <- p1 + theme(legend.position='bottom') + guides(color='none')
p1 <- p1 + labs(title='Fig 1c. C161 SCNA')

## poly-G tree
tree_ad <- nj(dm_polyg)
tree_ad <- root(tree_ad,outgroup=grep(paste0('^N'),tree_ad$tip.label),resolve.root=TRUE) 
tree_ad$edge.length[60] <- 0.5*tree_ad$edge.length[60] # truncate the normal branch
tree_ad_rotated <- ape::rotate(tree_ad, node=47)
p2 <- ggtree(tree_ad_rotated,layout='ape')
p2 <- p2 %<+% groups
p2 <- p2 + geom_tiplab(fontface=1,size=3,hjust=0,angle=0,aes(color=group))
p2 <- p2 + scale_color_manual(values=group_cols,name='Tissue')
p2 <- p2 + theme(legend.position='bottom') + guides(color='none')
p2 <- p2 + labs(title='Fig 1d. C161 poly-G')
p <- plot_grid(p1, p2, nrow=1)
ggsave(here('figures/fig_1cd.pdf'),width=9, height=6)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 1e. Violin plot similarity between SCNA and poly-G trees 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

patients <- c('C146','C154','C157','C159','C161','C186')

extract_polyg_scna_tree_similarity <- function(patient) {
    message(patient)
    similarity_results <- here(paste0('processed_data/copynumber/',patient,'/',patient,'_cnv_polyg_tree_similarity_test.txt'))
    d <- fread(similarity_results, select=c('perm','similarity'))
    setnames(d,'similarity','value')
    d[perm==0,similarity:='Observed']
    d[perm >0,similarity:='Permuted']
    d$patient <- patient
    obs <- d$value[d$similarity=='Observed']
    perms <- d$value[d$similarity=='Permuted']
    x <- sum(perms >= obs); n <- length(perms)
    pval <- (x+1)/(n+1)
    d$pval <- pval
    d
}

l <- lapply(patients, extract_polyg_scna_tree_similarity)
d <- rbindlist(l)
pvals <- d[!duplicated(patient),]
pvals$qval <- p.adjust(pvals$pval, method='BH')
pvals$label <- prettyNum(pvals$qval, digits=1)
pvals <- pvals[order(value,decreasing=F),]
d$patient <- factor(d$patient, levels=pvals$patient)
p <- ggplot(d, aes(x=patient, y=value)) +
    geom_point(data=d[similarity=='Observed'],size=4,pch=21,color='black',aes(fill=similarity)) +
    geom_violin(data=d[similarity=='Permuted'], aes(fill=similarity)) +
    geom_boxplot(data=d[similarity=='Permuted'],outlier.shape=NA,width=0.1) +
    geom_text(data=pvals,aes(x=patient,y=value,label=label),hjust=-0.25,vjust=-0.25) +
    theme_ang(base_size=12) +
    theme(legend.position='bottom') +
    scale_fill_brewer(palette='Accent',name='Comparison') +
    labs(x='Patient', y='Quartet similarity',title='Fig 1e. SCNA/poly-G tree similarity')
ggsave(here('figures/fig_1e.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 2E. RDS compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T] 
res1 <- get_met_specific_distances(si1, ad_table, comparison='rds', distance='node', return_tree=F)
res1 <- res1[type=='Met',]
res1$group <- 'Peritoneum'
res1$tissue_type <- 'Peritoneum'

si2.1 <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary') | tissue_type=='Lymph node') & in_collapsed==T]
res2.1 <- get_met_specific_distances(si2.1, ad_table, comparison='rds', distance='node', return_tree=F)
res2.1 <- res2.1[type=='Met',]
res2.1$group <- 'Locoregional'
res2.1$tissue_type <- 'Lymph node'
si2.2 <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary') | tissue_type=='Tumor deposit') & in_collapsed==T]
res2.2 <- get_met_specific_distances(si2.2, ad_table, comparison='rds', distance='node', return_tree=F)
res2.2 <- res2.2[type=='Met',]
res2.2$group <- 'Locoregional'
res2.2$tissue_type <- 'Tumor deposit'
res2 <- rbind(res2.1, res2.2)

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
res3 <- get_met_specific_distances(si3, ad_table, comparison='rds', distance='node', return_tree=F)
res3 <- res3[type=='Met',]
res3$group <- 'Liver'
res3$tissue_type <- 'Liver'

## supplement the liver data with rds from kim et al wxs-based trees
si4 <- sample_info[grepl('^CRC',Patient_ID) & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
res4 <- get_met_specific_distances(si4, kim_node_distances, comparison='rds', distance='node', return_tree=F)
res4 <- res4[type=='Met',]
res4$group <- 'Liver'
res4$tissue_type <- 'Liver'

res <- rbind(res1, res2, res3, res4)
res$group <- factor(res$group, levels=c('Locoregional','Peritoneum','Liver'))
tst <- dunn_test_ES(res, RDS ~ group)

tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
tmp_cols <- tmp_cols[names(tmp_cols) %in% res$tissue_type]
tmp_cols <- tmp_cols[c('Lymph node','Tumor deposit','Peritoneum','Liver')]

p <- ggplot(res, aes(x=group, y=RDS)) + 
    scale_y_continuous(breaks=seq(0,1,by=0.25)) +
    geom_beeswarm(cex = 3, aes(fill=tissue_type), size=4, pch=21, stroke=0.5) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method = "kruskal.test", label.y = 1.25, geom = "label") + 
    stat_pvalue_manual(tst, label='label', y.position=c(1.1,1.2,1.3)) +
    scale_fill_manual(values=tmp_cols,name='Tissue type') + 
    theme_ang(base_size=12) +
    labs(x=NULL,y='RDS',title='Fig 2e')
ggsave(here('figures/fig_2e.pdf'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 2F. pairwise AD compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

si2 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Locoregional') & in_collapsed==T]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Locoregional'; res2$group2 <- 'Locoregional'; res2$class <- 'LR:LR'

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'

res <- rbind(res1, res2, res3)
res$class <- factor(res$class, levels=c('LR:LR','Per:Per','Liv:Liv'))
tst <- dunn_test_ES(res, distance ~ class)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method = "kruskal.test", label.y = 1.8, geom = "label") + 
    stat_pvalue_manual(tst, label='label', y.position=c(1.5,1.7,1.9)) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Fig 2f')
ggsave(here('figures/fig_2f.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 2G. all-timing intra-lesion heterogeneity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')]
res_per <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_pt <- res_per[group1=='Primary' & group2=='Primary']
res_per <- res_per[group1=='Metastasis' & group2=='Metastasis']
res_per[res_per=='Metastasis'] <- 'Peritoneum'
res_per <- subset_for_intralesion(res_per)
res_per_with_pt <- rbind(res_pt, res_per, fill=T)
patients_with_geq2_per_mets <- names(which(rowSums(xtabs(~ patient + group1, data=res_per_with_pt) > 0) > 1))
res_per_with_pt <- res_per_with_pt[patient %in% patients_with_geq2_per_mets]

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
res_liv <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_pt <- res_liv[group1=='Primary' & group2=='Primary']
res_liv <- res_liv[group1=='Metastasis' & group2=='Metastasis']
res_liv[res_liv=='Metastasis'] <- 'Liver'
res_liv <- subset_for_intralesion(res_liv)
res_liv_with_pt <- rbind(res_pt, res_liv, fill=T)
patients_with_geq2_liv_mets <- names(which(rowSums(xtabs(~ patient + group1, data=res_liv_with_pt) > 0) > 1))

si <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Locoregional')]
res_lr <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_lr <- res_lr[group1=='Metastasis' & group2=='Metastasis',]
res_lr <- res_lr[group1=='Metastasis' & group2=='Metastasis']
res_lr[res_lr=='Metastasis'] <- 'Locoregional'
res_lr <- subset_for_intralesion(res_lr)

res <- rbind(res_per_with_pt, res_liv)
stat.test <- mywilcox2(res[group1 %in% c('Peritoneum','Liver')], distance ~ group1, paired=F)
res$group1 <- factor(res$group1, levels=unique(res$group1))

p <- ggplot(res[group1 %in% c('Peritoneum','Liver')], aes(x=group1, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Fig 2g')
ggsave(here('figures/fig_2g.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3a. RDS compared between tissue types with only untreated/synchronous mets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

si2.1 <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary') | tissue_type=='Lymph node') & in_collapsed==T]
si2.1 <- si2.1[group %in% c('Normal','Primary') | (met_treated=='untreated' & met_timing=='synchronous')]
res2.1 <- get_met_specific_distances(si2.1, ad_table, comparison='rds', distance='node', return_tree=F)
res2.1 <- res2.1[type=='Met',]
res2.1$group <- 'Locoregional'
res2.1$tissue_type <- 'Lymph node'
si2.2 <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary') | tissue_type=='Tumor deposit') & in_collapsed==T]
si2.2 <- si2.2[group %in% c('Normal','Primary') | (met_treated=='untreated' & met_timing=='synchronous')]
res2.2 <- get_met_specific_distances(si2.2, ad_table, comparison='rds', distance='node', return_tree=F)
res2.2 <- res2.2[type=='Met',]
res2.2$group <- 'Locoregional'
res2.2$tissue_type <- 'Tumor deposit'
res2 <- rbind(res2.1, res2.2)

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T] 
si1 <- si1[group %in% c('Normal','Primary') | (met_treated=='untreated' & met_timing=='synchronous')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='rds', distance='node', return_tree=F)
res1 <- res1[type=='Met',]
res1$group <- 'Peritoneum'
res1$tissue_type <- 'Peritoneum'

# verified that this does not exclude post-HIPEC-only liv mets
si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated=='untreated' & met_timing=='synchronous')] 
res3 <- get_met_specific_distances(si3, ad_table, comparison='rds', distance='node', return_tree=F)
res3 <- res3[type=='Met',]
res3$group <- 'Liver'
res3$tissue_type <- 'Liver'

## supplement the liver data with rds from kim et al wxs-based trees
si4 <- sample_info[grepl('^CRC',Patient_ID) & group %in% c('Normal','Primary','Liver') & in_collapsed==T & Patient_ID %in% liv_patients]
si4 <- si4[group %in% c('Normal','Primary') | (met_treated=='untreated' & met_timing=='synchronous')]
res4 <- get_met_specific_distances(si4, kim_node_distances, comparison='rds', distance='node', return_tree=F)
res4 <- res4[type=='Met',]
res4$group <- 'Liver'
res4$tissue_type <- 'Liver'

res <- rbind(res1, res2, res3, res4)
res$group <- factor(res$group, levels=c('Locoregional','Peritoneum','Liver'))
res$category <- 'untreated and synchronous'
out <- res[,c('patient','RDS','group','category'),with=F]

tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
tmp_cols <- tmp_cols[names(tmp_cols) %in% res$tissue_type]
tmp_cols <- tmp_cols[c('Lymph node','Tumor deposit','Peritoneum','Liver')]

tst <- dunn_test_ES(res, RDS ~ group)
tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
tmp_cols <- tmp_cols[names(tmp_cols) %in% res$tissue_type]
tmp_cols <- tmp_cols[c('Lymph node','Tumor deposit','Peritoneum','Liver')]

p <- ggplot(res, aes(x=group, y=RDS)) + 
    scale_y_continuous(breaks=seq(0,1,by=0.25)) +
    geom_beeswarm(cex = 3, aes(fill=tissue_type), size=4, pch=21, stroke=0.5) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method = "kruskal.test", label.y = 1.25, geom = "label") + 
    stat_pvalue_manual(tst, label='label', y.position=c(1.1,1.2,1.3)) +
    scale_fill_manual(values=tmp_cols,name='Tissue type') + 
    theme_ang(base_size=12) +
    labs(x=NULL,y='RDS (untreated+synchronous)',title='Fig 3a')
ggsave(here('figures/fig_3a.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3b. PM RDS: synchronous/untreated vs metachronous/treated
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

si1 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% c('treated','treated after untreated') & met_treated_type=='systemic chemo' & met_timing %in% c('metachronous','metachronous after synchronous')))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='rds', distance='node', return_tree=F)
res1 <- res1[type=='Met',]
res1$group <- 'Peritoneum'; 
res1$class <- 'Metachronous/\ntreated'

si2 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated' & met_timing=='synchronous'))]
res2 <- get_met_specific_distances(si2, ad_table, comparison='rds', distance='node', return_tree=F)
res2 <- res2[type=='Met',]
res2$group <- 'Peritoneum'; 
res2$class <- 'Synchronous/\nuntreated'

res_per <- rbind(res2, res1)
res_per$class <- factor(res_per$class, levels=unique(res_per$class))
stat.test <- mywilcox2(res_per, RDS ~ class, paired=F, include_n=F)
stat.test$y.position <- 1.1

p <- ggplot(res_per, aes(x=class, y=RDS)) + 
    scale_y_continuous(breaks=seq(0,1,by=0.25)) +
    geom_beeswarm(cex = 3, aes(fill=group), size=4, pch=21, color='black', stroke=0.25) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='RDS',title='Fig 3b')
ggsave(here('figures/fig_3b.pdf'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3c. PM inter-lesion AD: synchronous/untreated vs metachronous/treated
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

si1 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% c('treated','treated after untreated') & met_treated_type=='systemic chemo' & met_timing %in% c('metachronous','metachronous after synchronous')))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group <- 'Peritoneum'; 
res1$class <- 'Metachronous/\ntreated'

si2 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated' & met_timing=='synchronous'))]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group <- 'Peritoneum'; 
res2$class <- 'Synchronous/\nuntreated'

res_per <- rbind(res2, res1)
res_per$class <- factor(res_per$class, levels=unique(res_per$class))
stat.test <- mywilcox2(res_per, distance ~ class, paired=F, include_n=F)

p <- ggplot(res_per, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=42),pch=21,color='black',size=4,aes(fill=group)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Between-lesion AD',title='Fig 3c')
ggsave(here('figures/fig_3c.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3d. Liver RDS synchronous+untreated vs metachronous chemo
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1.1 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T] 
si1.1 <- si1.1[group %in% c('Normal','Primary') | (met_timing=='synchronous' & met_treated=='untreated')]
res1.1 <- get_met_specific_distances(si1.1, ad_table, comparison='rds', distance='node', return_tree=F)
res1.1 <- res1.1[type=='Met',]
res1.1$group <- 'Liver'
res1.1$class <- 'Synchronous, untreated'

## supplement the liver data with rds from kim et al wxs-based trees
si1.2 <- sample_info[grepl('^CRC',Patient_ID) & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si1.2 <- si1.2[group %in% c('Normal','Primary') | (met_timing=='synchronous' & met_treated=='untreated')]
res1.2 <- get_met_specific_distances(si1.2, kim_node_distances, comparison='rds', distance='node', return_tree=F)
res1.2 <- res1.2[type=='Met',]
res1.2$group <- 'Liver'
res1.2$class <- 'Synchronous, untreated'

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si2.1 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T] 
si2.1 <- si2.1[group %in% c('Normal','Primary') | (met_timing=='metachronous' & met_treated_type=='systemic chemo')]
res2.1 <- get_met_specific_distances(si2.1, ad_table, comparison='rds', distance='node', return_tree=F)
res2.1 <- res2.1[type=='Met',]
res2.1$group <- 'Liver'
res2.1$class <- 'Metachronous, systemic chemo'

## supplement the liver data with rds from kim et al wxs-based trees
si2.2 <- sample_info[grepl('^CRC',Patient_ID) & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si2.2 <- si2.2[group %in% c('Normal','Primary') | (met_timing=='metachronous' & met_treated_type=='systemic chemo')]
res2.2 <- get_met_specific_distances(si2.2, kim_node_distances, comparison='rds', distance='node', return_tree=F)
res2.2 <- res2.2[type=='Met',]
res2.2$group <- 'Liver'
res2.2$class <- 'Metachronous, systemic chemo'

res <- rbind(res1.1, res1.2, res2.1, res2.2)
res$class <- factor(res$class, levels=c('Synchronous, untreated','Metachronous, systemic chemo'))
stat.test <- mywilcox2(res, RDS ~ class, paired=F)

p <- ggplot(res, aes(x=class, y=RDS)) + 
    scale_y_continuous(breaks=seq(0,1,by=0.25)) +
    geom_beeswarm(cex = 3, aes(fill=group), size=4, pch=21, stroke=0.5) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.05) +
    scale_fill_manual(values=group_cols,name='Tissue type') + 
    theme_ang(base_size=12) +
    labs(x=NULL,y='RDS',title='Fig 3d. Met-specific RDS vs treatment')
ggsave(here('figures/fig_3d.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3e. Liver inter-lesion AD synchronous+untreated vs metachronous chemo
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Liv; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T] 
si1 <- si1[group %in% c('Normal','Primary') | (met_treated_type!='systemic chemo' & met_timing=='synchronous')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Liver'; res1$group2 <- 'Liver'; res1$class <- 'Synchronous, no systemic chemo'

## subset the sample_info for N, PT, and Liv; get met-spec node distance from peritoneum to normal
si2 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T] 
si2 <- si2[group %in% c('Normal','Primary') | (met_treated_type=='systemic chemo' & met_timing=='metachronous')]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Liver'; res2$group2 <- 'Liver'; res2$class <- 'Metachronous, systemic chemo'

res <- rbind(res1, res2)
res$class <- factor(res$class, levels=c('Synchronous, no systemic chemo','Metachronous, systemic chemo'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.25) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=10) +
    labs(x=NULL,y='Between-lesion AD',title='Fig 3e')
ggsave(here('figures/fig_3e.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3h. PM vs Liver (synchronous, no chemo) intra-lesion AD
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')] 
si1 <- si1[group %in% c('Normal','Primary') | (met_treated_type!='systemic chemo' & met_timing=='synchronous')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si2 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')] 
si2 <- si2[group %in% c('Normal','Primary') | (met_treated_type!='systemic chemo' & met_timing=='synchronous')]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Liver'; res2$group2 <- 'Liver'; res2$class <- 'Liv:Liv'

res <- rbind(res1, res2)
res <- subset_for_intralesion(res)
res$group1 <- factor(res$group1, levels=c('Peritoneum','Liver'))
stat.test <- mywilcox2(res, distance ~ group1, paired=F, include_n=F)

p <- ggplot(res, aes(x=group1, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=42),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.5) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=10) +
    labs(x=NULL,y='Within-lesion AD\n(synchronous, untreated)',title='Fig 3h')
ggsave(here('figures/fig_3h.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 3i. PM vs Liver (untreated, any timing) intra-lesion AD
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')] 
si1 <- si1[group %in% c('Normal','Primary') | (met_treated=='untreated')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si2 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
si2 <- si2[group %in% c('Normal','Primary') | (met_treated=='untreated')]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Liver'; res2$group2 <- 'Liver'; res2$class <- 'Liv:Liv'

res <- rbind(res1, res2)
res <- subset_for_intralesion(res)
res$group1 <- factor(res$group1, levels=c('Peritoneum','Liver'))
stat.test <- mywilcox2(res, distance ~ group1, paired=F, include_n=F)

p <- ggplot(res, aes(x=group1, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=42),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.53) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=10) +
    labs(x=NULL,y='Within-lesion AD\n(untreated, any timing)',title='Fig 3i')
ggsave(here('figures/fig_3i.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 4b. sync metastases types vs t-stage
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m <- read_distance_matrix(here('original_data/misc/lemmens_ijc_2011.txt'))
m <- m[rownames(m)!='Tx',]
m <- t(m)
prop <- copy(m)
for(i in 1:4) prop[,i] <- prop[,i] / sum(prop[,i])
dat <- as.data.table(reshape2::melt(m))
names(dat) <- c('mets','tstage','n')
prop <- as.data.table(reshape2::melt(prop))
names(prop) <- c('mets','tstage','prop')
dat$prop <- prop$prop
dat$mets <- factor(dat$mets, levels=rev(c('PC','PC+','Liv')))
cols <- c('#4C86C6','#FAB31D','#af7d14')
names(cols) <- c('Liv','PC','PC+')

get_label_pos <- function(dat) {
    dat <- dat[order(mets,decreasing=T),]
    dat$pos <- (cumsum(dat$prop) - 0.5*dat$prop)
    dat
}
dat2 <- dat[,get_label_pos(.SD), by=c('tstage')]

p <- ggplot(dat2, aes(x=tstage, y=prop)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,by=0.25)) +
    geom_bar(stat='identity', aes(fill=mets),color='black',linewidth=0.25) +
    geom_text(data=dat2[prop > 0],aes(label=n,y=pos)) +
    labs(x='T-stage',y='Fraction of patients',title='Fig 4b') + 
    scale_fill_manual(values=cols, name='Synchronous\nmetastases') + 
    theme_ang(base_size=12) +
    theme(legend.position='right')
ggsave(here('figures/fig_4b.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 4e. permutation-test for association between PMs and deep/luminal PTs
# ED Fig 6. same thing for LN, TD, Liver mets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test <- function(i, res) {
    if(i > 0) {
        shuffled <- sample(res$vertical, replace=F)
        res$vertical <- shuffled
    }
    res_deep <- res[vertical=='deep',]
    res_deep <- res_deep[,c('sample1','sample2','vertical','distance'),with=F]
    res_deep <- res_deep[order(sample1, distance, decreasing=F),]
    res_deep <- res_deep[!duplicated(sample1),]
    res_luminal <- res[vertical=='mucosal/luminal',]
    res_luminal <- res_luminal[,c('sample1','sample2','vertical','distance'),with=F]
    res_luminal <- res_luminal[order(sample1, distance, decreasing=F),]
    res_luminal <- res_luminal[!duplicated(sample1),]
    res_merged <- merge(res_deep[,c('sample1','distance'),with=F], 
                        res_luminal[,c('sample1','distance'),with=F], 
                        by='sample1')
    res_merged[,ratio:=distance.x / distance.y]
    ratio <- mean(res_merged$ratio)
    list(i=i, ratio=ratio)
}

run_test <- function(res, R, ncpus) {
    patient <- unique(res$patient)
    group <- unique(res$group)
    cohort <- unique(res$cohort)
    message(paste0(patient,', ',group,', ',cohort)) 

    ## run test for the patient
    l <- mclapply(0:R, test, res, mc.cores=ncpus)
    l <- rbindlist(l)
    obs <- l$ratio[l$i==0]
    exp <- l$ratio[l$i > 0]

    ## get two-sided p-value
    n <- R + 1
    x1 <- sum(exp <= obs) + 1
    p1 <- x1 / n
    x2 <- sum(exp >= obs) + 1
    p2 <- x2 / n
    p_twosided <- 2*min(c(p1,p2))

    ## get numbers for this patient
    tmp <- res[!duplicated(sample2),]
    n_deep <- sum(tmp$vertical=='deep')
    n_lum <- sum(tmp$vertical=='mucosal/luminal')
    n_met <- length(unique(res$sample1))

    list(patient=patient, group=group, cohort=cohort, pval=p_twosided, obs=obs, exp=exp, 
         exp_center=median(exp,na.rm=T), exp_mean=mean(exp,na.rm=T), n_deep=n_deep, n_lum=n_lum, n_met=n_met)
}

get_results <- function(l) l[names(l)!='exp']
get_exp <- function(l) {
    tmp <- l[names(l)!='exp']
    out <- data.table(expected=l$exp)
    for(f in names(tmp)) {
        out[[f]] <- tmp[[f]]
    }
    out
}

si <- sample_info[in_collapsed==T]
res <- get_met_specific_distances(si, ad_table, comparison='primary', distance='angular', return_tree=F)
res <- merge(res, si[,c('Patient_ID','Real_Sample_ID','vertical','cohort'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res <- merge(res, si[,c('Patient_ID','Real_Sample_ID','group'),with=F], by.x=c('patient','sample1'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res <- res[group %in% c('Locoregional','Peritoneum','Liver')]
res$plotgroup <- as.character(NA)
res[(patient %in% c('C38','C89') | cohort=='peritoneal') & group=='Locoregional' & grepl('LN',sample1), plotgroup:='Lymph node']
res[(patient %in% c('C38','C89') | cohort=='peritoneal') & group=='Locoregional' & grepl('TD',sample1), plotgroup:='Tumor deposit']
res[(patient %in% c('C38','C89') | cohort=='peritoneal') & group=='Peritoneum', plotgroup:='Peritoneum']
res[cohort %in% c('science','natgen','peritoneal') & group=='Liver', plotgroup:='Liver']
res <- res[!is.na(plotgroup),]
res[,id:=paste0(patient,', ',plotgroup)]
res[,group:=plotgroup]
res_split <- split(res, by='id')
RNGkind("L'Ecuyer-CMRG") 
set.seed(42)

## this takes a long time to run
required_file <- here('processed_data/misc/tissue_type_vs_depth_permutation_tests.txt')
if(!file.exists(required_file)) {
    l <- lapply(res_split, run_test, R=10000, ncpus=8)
    results <- rbindlist(lapply(l, get_results))
    results <- results[n_deep > 0 & n_lum > 0,]
    results[,lfc:=log2(obs / exp_center)]
    getq <- function(results) {
        results$qval <- p.adjust(results$pval, method='BH')
        results
    }
    results <- results[,getq(.SD),by=c('group')]
    results[,nlog10p:=-log10(pval)]
    results[,nlog10q:=-log10(qval)]
    results[,frac_deep:=n_deep / (n_deep+n_lum)]
    write_tsv(results,required_file)
} else {
    message('Using pre-processed file: ',required_file)
    results <- fread(required_file)
}

results$group <- factor(results$group, levels=c('Lymph node','Tumor deposit','Liver','Peritoneum'))
results[pval < 0.01 & qval >= 0.01, significance:='n.s.']
results[qval < 0.01, significance:='p-adj < 0.01']
results[is.na(significance),significance:='n.s.']
cols <- c('#bfbfbf','red')
names(cols) <- c('n.s.','p-adj < 0.01')

p <- ggplot(results[group=='Peritoneum'], aes(x=lfc, y=nlog10q)) +
    scale_y_continuous(limits=c(0,3.0)) +
    scale_x_continuous(limits=c(-1.5,1.5)) +
    geom_vline(xintercept=0, color='black', linewidth=0.5) +
    geom_point(pch=21, size=2.5, aes(fill=significance), color='black', stroke=0.25) + 
    geom_text_repel(data=results[pval < 0.05 & group=='Peritoneum'], 
                    aes(label=patient, color=significance)) + 
    scale_color_manual(values=cols, name='Significance') + 
    scale_fill_manual(values=cols, name='Significance') + 
    theme_bw(base_size=12) +
    labs(x='Effect size', y='-log10(q-value)', title='Fig 4e')
ggsave(here('figures/fig_4e.pdf'),width=5, height=3.5)

p <- ggplot(results[group!='Peritoneum'], aes(x=lfc, y=nlog10q)) +
    scale_y_continuous(limits=c(0,3.0)) +
    scale_x_continuous(limits=c(-1.5,1.5)) +
    geom_vline(xintercept=0, color='black', linewidth=0.5) +
    geom_point(pch=21, size=2.5, aes(fill=significance), color='black', stroke=0.25) + 
    geom_text_repel(data=results[pval < 0.05 & group!='Peritoneum'],
                    aes(label=patient, color=significance)) + 
    scale_color_manual(values=cols, name='Significance') + 
    scale_fill_manual(values=cols, name='Significance') + 
    facet_wrap(facets=~group, scale='free', ncol=4) +
    theme_bw(base_size=12) +
    labs(x='Effect size', y='-log10(q-value)', title='ED Fig 6') +
    theme(legend.position='none')
ggsave(here('figures/ed_fig_6.pdf'), width=9, height=3.5)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 4f. AD between met and deep/luminal PTs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get patients with any peritoneal mets and with any liver mets
per_patients <- unique(sample_info[group=='Peritoneum','Patient_ID',with=F][[1]])
liv_patients <- unique(sample_info[group=='Liver','Patient_ID',with=F][[1]])
dm_patients <- unique(sample_info[group %in% c('Liver','Lung','Distant (other)'),'Patient_ID',with=F][[1]])
non_per_liv_patients <- liv_patients[!liv_patients %in% per_patients]
non_per_dm_patients <- dm_patients[!dm_patients %in% per_patients]

si_per <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Peritoneum')) & in_collapsed==T]
res_per <- get_met_specific_distances(si_per, ad_table, comparison='primary', distance='angular', return_tree=F)
res_per$group1 <- 'Peritoneum'
si_lr <- sample_info[Patient_ID %in% per_patients & cohort=='peritoneal' & (group %in% c('Normal','Primary','Locoregional')) & in_collapsed==T ]
res_lr <- get_met_specific_distances(si_lr, ad_table, comparison='primary', distance='angular', return_tree=F)
res_lr$group1 <- 'Locoregional'
si_liv1 <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
res_liv1 <- get_met_specific_distances(si_liv1, ad_table, comparison='primary', distance='angular', return_tree=F)
res_liv1$group1 <- 'Liver (PM patients)'
si_liv2 <- sample_info[Patient_ID %in% non_per_liv_patients & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
res_liv2 <- get_met_specific_distances(si_liv2, ad_table, comparison='primary', distance='angular', return_tree=F)
res_liv2$group1 <- 'Liver (non-PM patients)'
si_dm1 <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Liver','Lung','Distant (other)')) & in_collapsed==T ]
res_dm1 <- get_met_specific_distances(si_dm1, ad_table, comparison='primary', distance='angular', return_tree=F)
res_dm1$group1 <- 'Distant (PM patients)'
si_dm2 <- sample_info[Patient_ID %in% non_per_dm_patients & (group %in% c('Normal','Primary','Liver','Lung','Distant (other)')) & in_collapsed==T ]
res_dm2 <- get_met_specific_distances(si_dm2, ad_table, comparison='primary', distance='angular', return_tree=F)
res_dm2$group1 <- 'Distant (non-PM patients)'
si_dm3 <- sample_info[Patient_ID %in% dm_patients & (group %in% c('Normal','Primary','Liver','Lung','Distant (other)')) & in_collapsed==T ]
res_dm3 <- get_met_specific_distances(si_dm3, ad_table, comparison='primary', distance='angular', return_tree=F)
res_dm3$group1 <- 'Distant (All patients)'

info <- fread(here('processed_data/pm_cohort_per_and_liver_case_samples.txt'))
info <- info[in_collapsed==T,]
info <- info[order(Patient_ID,timing_order,group),]
info[,id:=paste(Patient_ID, Real_Sample_ID)]

get_earlier_and_later_pms_for_each_liver_met <- function(liv_sample, info) {
    message(liv_sample)
    current_patient <- info[id==liv_sample,(Patient_ID)]
    tmp <- info[Patient_ID==current_patient]
    sample_pos <- tmp[id==liv_sample,(timing_order)]    
    earlier_pms <- tmp[timing_order < sample_pos & group=='Peritoneum',(id)]
    sync_pms <- tmp[timing_order == sample_pos & group=='Peritoneum',(id)]
    later_pms <- tmp[timing_order > sample_pos & group=='Peritoneum',(id)]
    d1 <- data.table(sample=earlier_pms, pos='liv after pm')
    d2 <- data.table(sample=sync_pms, pos='liv sync with pm')
    d3 <- data.table(sample=later_pms, pos='liv before pm')
    out <- rbind(d1, d2, d3)
    out$liv_sample <- liv_sample
    out$patient <- current_patient
    out
}

liv_samples <- unique(info[group=='Liver',(id)])
l <- lapply(liv_samples, get_earlier_and_later_pms_for_each_liver_met, info)
ll <- rbindlist(l)

summarize_liv_sample <- function(ll) {
    n_pms_prior_to_liv <- sum(ll$pos=='liv after pm')
    n_pms_same_time_as_liv <- sum(ll$pos=='liv sync with pm')
    n_pms_after_liv <- sum(ll$pos=='liv before pm')
    list(n_pms_prior_to_liv=n_pms_prior_to_liv, 
         n_pms_same_time_as_liv=n_pms_same_time_as_liv, 
         n_pms_after_liv=n_pms_after_liv)
}
n <- ll[,summarize_liv_sample(.SD), by=c('patient','liv_sample')]

bucket1 <- n[n_pms_prior_to_liv == 0 & n_pms_same_time_as_liv==0,(liv_sample)] # Liver before PMs
bucket2 <- n[n_pms_prior_to_liv == 0 & n_pms_same_time_as_liv > 0,(liv_sample)] # Liver same time as PMs
bucket3 <- n[n_pms_prior_to_liv > 0,(liv_sample)] # Liver after PMs

pd <- copy(res_liv1)
pd[,id:=paste(patient, sample1)]
pd[id %in% bucket1, group1:='Liver before PMs']
pd[id %in% bucket2, group1:='Liver same time as PMs']
pd[id %in% bucket3, group1:='Liver after PMs']
pd[group1 %in% c('Liver before PMs'),group1:='Liver before PMs']
pd[group1 %in% c('Liver same time as PMs','Liver after PMs'),group1:='Liver same time/after PMs']
pd[,id:=NULL]
res <- rbind(res_per, res_lr, pd, res_liv2)

res[,patient_type:=group1]
res[grepl('Liver',patient_type), group1:='Liver']
res[grepl('Distant',patient_type), group1:='Distant']
res <- merge(res, sample_info[,c('Patient_ID','Real_Sample_ID','vertical','cohort'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res$vertical <- factor(res$vertical, levels=c('mucosal/luminal','deep'))
res <- res[!is.na(vertical)]
res[grepl('LN',sample1),group1:=gsub('Locoregional','Lymph node',group1)]
res[grepl('TD',sample1),group1:=gsub('Locoregional','Tumor deposit',group1)]
res[group1 %in% c('Lymph node','Tumor deposit'), patient_type:=group1]

pd1 <- res[patient_type %in% c('Peritoneum','Lymph node','Tumor deposit')]
pd1$patient_type <- factor(pd1$patient_type, levels=c('Peritoneum','Lymph node','Tumor deposit'))
stat.test1 <- mywilcox2(pd1, distance ~ vertical, facet_field='patient_type', paired=F, include_n=F)
pd2 <- res[patient_type %in% c('Liver same time/after PMs','Liver before PMs','Liver (non-PM patients)')]
pd2$patient_type <- factor(pd2$patient_type, levels=c('Liver same time/after PMs','Liver before PMs','Liver (non-PM patients)'))
stat.test2 <- mywilcox2(pd2, distance ~ vertical, facet_field='patient_type', paired=F, include_n=F)

p1 <- ggplot(pd1, aes(x=vertical, y=distance)) +
    scale_y_continuous(limits=c(0,2.25),breaks=seq(0,2.35,by=0.5)) +
    geom_point(aes(fill=vertical),position=position_jitter(width=0.1,height=0,seed=2), pch=21, color='white', stroke=0.25, size=3.0) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA,width=0.5) +
    scale_fill_manual(values=cols_depth2) +
    guides(fill='none', color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    facet_wrap(facets=~patient_type, nrow=1, ncol=3) +
    stat_pvalue_manual(stat.test1, label = "label", tip.length = 0.02, y.position=1.75) +
    labs(x='Invasion depth of primary tumor', y='Angular distance',title='Fig 4f')
ggsave(here('figures/fig_4f.pdf'), width=6, height=4.5)


p2 <- ggplot(pd2, aes(x=vertical, y=distance)) +
    scale_y_continuous(limits=c(0,2.25),breaks=seq(0,2.35,by=0.5)) +
    geom_point(aes(fill=vertical),
               position=position_jitter(width=0.1,height=0,seed=2), 
               pch=21, color='white', stroke=0.25, size=3.0) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA,width=0.5) +
    scale_fill_manual(values=cols_depth2) +
    guides(fill='none', color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    facet_wrap(facets=~patient_type, nrow=1, ncol=3) +
    stat_pvalue_manual(stat.test2, label = "label", tip.length = 0.02, y.position=2.15) +
    labs(x='Invasion depth of primary tumor', y='Angular distance',title='Fig 4g')
ggsave(here('figures/fig_4g.pdf'), width=6, height=4.5)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2024-07-23
# Fig 5d
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_patient_met_level <- function(patient, si, ad_table, query_group, comparitor_group, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary',query_group, comparitor_group)]

    ## get the original set of valid-samples for the true number of distant/pt/ln samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))
    if(query_group!='Liver') {
        n_comparitor <- length(unique(si$Real_Sample_ID[si$group %in% comparitor_group]))
        n_query <- length(unique(si$Real_Sample_ID[si$group %in% query_group]))
    } else {
        n_comparitor <- length(unique(si$Real_Sample_ID[si$group %in% c('Distant','Liver')])) - 1
        n_query <- length(unique(si$Real_Sample_ID[si$group=='Liver']))
    }

    ## load the bootstrapped data for this patient
    bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/',patient,'.rds'))
    bs <- readRDS(bs_file)

    test_patient_bs <- function(i, patient, si, ad_table, bs) {
        if(i > 0) {
            ## for the bootstrap replicates, obtain new ad and valid_samples
            ad <- bs[[i]]
            valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
        }
        ad <- ad[valid_samples, valid_samples]
        ad <- as.data.table(reshape2::melt(ad))    
        groups <- si[,c('Real_Sample_ID','group'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        ad <- ad[,c('Var1','Var2','group1','group2','value'),with=F] 

        ## get the min distances between PMs and DMs
        comparitor <- ad[group1 %in% query_group & group2 %in% comparitor_group]
        comparitor <- comparitor[Var1!=Var2,]
        comparitor <- comparitor[order(Var1, value, decreasing=F),]
        comparitor <- comparitor[!duplicated(Var1),]
        pt <- ad[group1 %in% query_group & group2=='Primary']
        pt <- pt[order(Var1, value, decreasing=F),]
        pt <- pt[!duplicated(Var1),]

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        out <- merge(comparitor[,c('Var1','Var2','group2','value'),with=F], pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T)
        names(out) <- c('sample','comparitor_sample','comparitor_type','comparitor_distance','pt_sample','pt_distance')
        out$i <- i
        out$patient <- patient
        out
    }

    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(comparitor_distance / pt_distance)]
    closest_samples <- res[!duplicated(sample),c('sample','comparitor_sample'),with=F]
    collapse_ln <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_ln(.SD),by=c('patient','sample')] 
    res <- merge(res, closest_samples, by='sample', all.x=T)
    res$n_comparitor <- n_comparitor
    res$n_query <- n_query
    res$n_pt <- n_pt
    res
}

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Shared (95%)','Shared (80%)','n.s.','Distinct (80%)','Distinct (95%)')

# Fig 5d-i, common/distinct origin of PMs and distant mets
si <- copy(sample_info)
si[group %in% c('Liver','Lung','Distant (other)'), group:='Distant']
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
liv_patients <- unique(si[group=='Distant','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, liv_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]
RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, query_group='Peritoneum', comparitor_group='Distant', ncpus=ncpus)
res <- rbindlist(l)
res[,query:='Peritoneum']
res[,comparitor:='Distant']
write_tsv(res,here('processed_data/min_distance_ratios_pm_to_dm_vs_pt.txt'))

# Fig 5d-ii, common/distinct origin of LNs and distant mets
si <- copy(sample_info)
si[group %in% c('Liver','Lung','Distant (other)'), group:='Distant']
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
si[Patient_ID %in% per_patients & tissue_type %in% 'Lymph node', group:='Lymph node']
ln_patients <- unique(si[group=='Lymph node','Patient_ID',with=F][[1]])
dm_patients <- unique(si[group=='Distant','Patient_ID',with=F][[1]])
valid_patients <- intersect(ln_patients, dm_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]
RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, query_group='Lymph node', comparitor_group='Distant', ncpus=ncpus)
res <- rbindlist(l)
res[,query:='Lymph node']
res[,comparitor:='Distant']
write_tsv(res,here('processed_data/min_distance_ratios_ln_to_dm_vs_pt.txt'))

# Fig 5d-iii, common/distinct origin of TDs and distant mets
si <- copy(sample_info)
si[group %in% c('Liver','Lung','Distant (other)'), group:='Distant']
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
si[Patient_ID %in% per_patients & tissue_type %in% 'Tumor deposit', group:='Tumor deposit']
td_patients <- unique(si[group=='Tumor deposit','Patient_ID',with=F][[1]])
dm_patients <- unique(si[group=='Distant','Patient_ID',with=F][[1]])
valid_patients <- intersect(td_patients, dm_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]
RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, query_group='Tumor deposit', comparitor_group='Distant', ncpus=ncpus)
res <- rbindlist(l)
res[,query:='Tumor deposit']
res[,comparitor:='Distant']
write_tsv(res,here('processed_data/min_distance_ratios_td_to_dm_vs_pt.txt'))

# Fig 5d-iv, common/distinct origin of Liver mets and other DMs (which may include Liver mets) 
si <- copy(sample_info)
count_mets <- function(qc) {
    n <- length(unique(qc$Real_Sample_ID))
    list(n=n)
}
counts <- si[in_collapsed==T,count_mets(.SD),by=c('Patient_ID','group')]
counts <- data.table::dcast(Patient_ID ~ group, value.var='n', data=counts)
counts[is.na(counts)] <- 0
counts[,Distant_non_liver:=Lung + `Distant (other)`]
valid_patients <- unique(counts[Liver >= 2 | (Liver==1 & Distant_non_liver >= 1), (Patient_ID)])
valid_patients <- valid_patients[grepl('CRC',valid_patients)==F]
si[group %in% c('Lung','Distant (other)'), group:='Distant']
RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, query_group='Liver', comparitor_group=c('Liver','Distant'), ncpus=ncpus)
res <- rbindlist(l)
res[,query:='Liver']
res[,comparitor:='Distant']
write_tsv(res,here('processed_data/min_distance_ratios_liv_to_dm_vs_pt.txt'))


## compare distributions of min-dist-ratios for any difference between PM and LNs
res_pm <- fread(here('processed_data/min_distance_ratios_pm_to_dm_vs_pt.txt'))
res_ln <- fread(here('processed_data/min_distance_ratios_ln_to_dm_vs_pt.txt'))
res_td <- fread(here('processed_data/min_distance_ratios_td_to_dm_vs_pt.txt'))
res_liv <- fread(here('processed_data/min_distance_ratios_liv_to_dm_vs_pt.txt'))
res <- rbind(res_pm, res_ln, res_liv, res_td, fill=T)
res[,comparitor_pt_ratio:=n_comparitor/n_pt]
res$i <- 1:nrow(res)
res[,group:=query]
res <- res[,c('patient','sample','obs','2.5%','10%','90%','97.5%','query','group','n_query','n_comparitor','comparitor_pt_ratio'),with=F]
res[group %in% c('Lymph node','Tumor deposit'), group:='Locoregional']

## plot on the lesion-level
pd <- res[order(obs,decreasing=F),]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`90%` < 0, origin:='Shared (80%)']
pd[`97.5%` < 0, origin:='Shared (95%)']
pd[`10%` > 0, origin:='Distinct (80%)']
pd[`2.5%` > 0, origin:='Distinct (95%)']
pd$origin <- factor(pd$origin, levels=c('Shared (95%)','Shared (80%)','n.s.','Distinct (80%)','Distinct (95%)'))
to_percentile <- function(pd) {
    pd <- pd[order(obs,decreasing=F),]
    pd$percentile <- (1:nrow(pd)) / nrow(pd)
    pd
}
pd <- pd[,to_percentile(.SD),by=group]
pd$group <- factor(pd$group, levels=c('Peritoneum','Locoregional','Liver'))
p1 <- ggplot(pd, aes(x=percentile, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', linewidth=0.5) +
    facet_wrap(facets=~group, ncol=1, scale='free_y') +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio',title='Fig 5d')
ggsave(here('figures/fig_5d.pdf'), height=6.5, width=6)



stat.test1 <- mywilcox2(res[group %in% c('Peritoneum','Locoregional')], obs ~ group, paired=F)
stat.test1$y.position <- 3.25
stat.test2 <- mywilcox2(res[group %in% c('Peritoneum','Liver')], obs ~ group, paired=F)
stat.test2$y.position <- 3.5
stat.test3 <- mywilcox2(res[group %in% c('Locoregional','Liver')], obs ~ group, paired=F)
stat.test3$y.position <- 3.75
res$group <- factor(res$group, levels=c('Peritoneum','Locoregional','Liver'))
stat.test <- rbind(stat.test1, stat.test2, stat.test3)
tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')

p <- ggplot(res, aes(x=group, y=obs)) +
    geom_hline(yintercept=0, color='#bfbfbf', linetype='dashed', linewidth=0.25) + 
    geom_point(position=position_jitter(width=0.1, height=0, seed=42), pch=21, size=2.25, aes(fill=query), stroke=0.25) +
    geom_boxplot(fill=NA, outlier.shape=NA, width=0.4) +
    scale_fill_manual(values=tmp_cols) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = NA) +
    theme_ang(base_size=10) +
    guides(fill='none') +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))  +
    labs(x='Met type', y='log2 min-distance-ratio',title='Fig 5f')
ggsave(here('figures/fig_5f.pdf'))

res$group_factor <- factor(res$group, levels=c('Liver','Peritoneum','Locoregional'))
m1 <- lm(obs ~ group_factor + comparitor_pt_ratio, data=res); summary(m1)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 5e. Liver/PM origins vs relative timing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## test bucket 1,2,3 for shared origins with PMs
info <- fread(here('processed_data/pm_cohort_per_and_liver_case_samples.txt'))
info <- info[in_collapsed==T,]
info <- info[order(Patient_ID,timing_order,group),]
info[,id:=paste(Patient_ID, Real_Sample_ID)]
valid_patients <- unique(info$Patient_ID)

# Fig 5e, common/distinct origin of Liver/PMs vs met timing
si <- copy(sample_info)
RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, query_group='Peritoneum', comparitor_group='Liver', ncpus=ncpus)
res <- rbindlist(l)
res[,id:=paste(patient,comparitor_sample)]
res[id %in% bucket1, bucket:='Liver before PMs']
res[id %in% bucket2, bucket:='Liver same time as PMs']
res[id %in% bucket3, bucket:='Liver after PMs']
res <- res[order(obs),]
cols <- c('steelblue','#d8d8d8','salmon')
names(cols) <- c('Liver after PMs','Liver same time as PMs','Liver before PMs')

stat.test <- mywilcox2(res, obs ~ bucket, paired=F, include_n=F)
cols <- c('steelblue','#d8d8d8','salmon')
names(cols) <- c('Liver after PMs','Liver same time as PMs','Liver before PMs')
res$bucket <- factor(res$bucket, levels=names(cols))

p3 <- ggplot(res, aes(x=bucket, y=obs)) +
    geom_point(aes(color=bucket),position=position_jitter(width=0.1,height=0,seed=2), pch=16, size=4) + 
    geom_hline(yintercept=0, linetype='dashed', linewidth=0.5) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA,width=0.5) +
    scale_color_manual(values=cols, name='Timing') +
    theme_ang(base_size=12) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = NA, y.position=c(2.1,2,2.2)) +
    labs(x='Time of liver met. diagnosis', y='Origin ratio', title='Fig 5e')
ggsave(here('figures/fig_5e.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ED Fig 2. multi-primary tumor origins
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to get mean lengths 
get_markerlengths  <- function(dir) { 

    subject <- str_split(dir, "/") %>%
        purrr::map(tail, n = 1) %>%
        unlist() %>%
        str_split("_") %>%
        purrr::map(1) %>%
        unlist()

    ## path to poly-G raw data directory (marker length files)
    marker_dir <- list.files(dir, pattern = "repre_repli", full.names = TRUE)

    ## load marker lengths and get the mean length of each marker in each sample
    get_markers <- function(marker_dir) {
        marker <- suppressMessages(read_tsv(marker_dir))

        # getting the marker name
        int_file <- tail(str_split(marker_dir, "/")[[1]], 1) 
        marker_name <- head(str_split(int_file, "_")[[1]], 1)
        cols <- NCOL(marker) +1 

        marker  %>% 
            rowid_to_column("length")  %>% 
            pivot_longer(cols=c(2:cols), names_to="sample", values_to="Frequency")  %>% 
            group_by(sample)  %>% 
            mutate(Frequency=Frequency/sum(Frequency), sample=str_remove(sample, "_[1-3]$")) %>%
            summarize(length = sum(length * Frequency)) %>%
            mutate(marker=marker_name) 
    }

    marker_files <- list.files(marker_dir, full.names = TRUE)
    markers <- lapply(marker_files, get_markers)
    markers <- bind_rows(markers)

    # subtract marker length of normal from each marker   
    root_sample <- str_remove(str_subset(unique(markers$sample), "N"), "_[1-3]$")[1] 

    markers <- markers %>% 
        mutate(sample=str_remove(sample, "_[1-3]$"))  %>% 
        group_by(marker) %>% 
        filter(any(sample == root_sample)) %>% 
        mutate(length = length - length[sample == root_sample]) %>% 
        ungroup 

    markers$subject <- subject
    return(markers)
}

# directories of all patients with suspected multiple primary tumors
multi_tumor_dirs <- c(
                      here("processed_data/polyG/peritoneal/results/sample_exclusion_0.01_rep_cut_0.11/E3_R"),
                      here("processed_data/polyG/peritoneal/results/sample_exclusion_0.05_rep_cut_0.11/E10_R"),
                      here("processed_data/polyG/peritoneal/results/sample_exclusion_0.1_rep_cut_0.11/E11_R"),
                      here("processed_data/polyG/peritoneal/results/sample_exclusion_0.1_rep_cut_0.11/E15_R"))

# getting mean lengths per patient
l <- lapply(multi_tumor_dirs, get_markerlengths)
markerlengths <- bind_rows(l)

# functions to get markers that are present in all samples
get_sampler_marker  <- function(sample_i, marker_tbl) {

    marker_tbl  %>% 
        filter(sample==sample_i)  %>% 
        pull(marker)  %>% 
        unique()
}

get_minimum_marker_table  <- function(subject_i, markerlengths){

    marker_tbl <-  markerlengths  %>% 
        filter(subject==subject_i)  %>% 
        group_by(sample)  %>% 
        add_count(sample)  %>% 
        ungroup()  %>% 
        mutate(minimum_markers = max(n)*0.7) %>% 
        filter(n>minimum_markers) 

    samples  <- unique(marker_tbl$sample)
    marker_list  <- lapply(samples, get_sampler_marker, marker_tbl)
    common_markers  <- Reduce(intersect, marker_list)
    n_common_markers  <- length(common_markers)

    markerlengths  %>% 
        filter(subject==subject_i, marker %in% common_markers)   %>% 
        group_by(sample)  %>% 
        add_count(sample)  %>% 
        filter(n==n_common_markers)  %>% 
        select(-n)
}

# finding all shared markers per patient
subjects  <- unique(markerlengths$subject)
minimum_marker_tbls  <- lapply(subjects, get_minimum_marker_table, markerlengths)
minimum_markerlengths <- bind_rows(minimum_marker_tbls) 

# find samples
samples <- minimum_markerlengths$sample %>% unique
combos_wide <- combn(samples, m= 2) %>% as.data.frame()

## make combo table longer
combos <- data.frame(a = as.character(combos_wide[1, ]), b = as.character(combos_wide[2, ])) %>%
    mutate(subject_a =  str_extract(a, "E[:digit:]+"), subject_b =  str_extract(b, "E[:digit:]+")) %>% 
    filter(subject_a==subject_b) %>% 
    select(1:2)

# function to get L1 and cor for a specific combination of two samples
get_l1_r_for_combination <- function(i, combos, markerlengths) {
    sample_a <- combos$a[i]
    sample_b <- combos$b[i]

    markerlengths_a  <- markerlengths %>% 
        dplyr::filter(sample==sample_a) %>% 
        arrange(marker)  %>% 
        pull(length)

    markerlengths_b  <- markerlengths %>% 
        dplyr::filter(sample==sample_b) %>% 
        arrange(marker)  %>% 
        pull(length)

    n_markers_a  <- length(markerlengths_a)
    n_markers_b  <- length(markerlengths_b)
    if (n_markers_a != n_markers_b) (error)
    l1  <- sum(abs(markerlengths_a - markerlengths_b))/n_markers_a
    r <- suppressWarnings(cor(markerlengths_a, markerlengths_b))

    list(a=sample_a, b=sample_b, l1=l1, r=r, marker=n_markers_a)
}

# create a table for the cor  and L1 of marker lengths 
non_boot_l1 <- lapply(1:nrow(combos), get_l1_r_for_combination, combos, minimum_markerlengths)
non_boot_l1 <- bind_rows(non_boot_l1)

# get table with only polyG correlation
d <- as.data.table(non_boot_l1) %>%
    select(a, b, r)

d[grepl('^E10',a) & grepl('^E10',b), patient:='E10']
d[grepl('^E11',a) & grepl('^E11',b), patient:='E11']
d[grepl('^E15',a) & grepl('^E15',b), patient:='E15']
d[grepl('^E3',a) & grepl('^E3',b), patient:='E3']

info <- fread(here('processed_data/sample_info.txt'))
info[Patient_ID %in% c('Ea3','Eb3','Ec3'), Patient_ID:='E3']
info <- info[!duplicated(Sample_ID),]
info <- info[Patient_ID %in% c('E3','E10','E11','E15'),]
info <- info[,c('Sample_ID','Real_Sample_ID','group'),with=F]
d <- merge(d, info, by.x='a', by.y='Sample_ID', all.x=T)
d <- merge(d, info, by.x='b', by.y='Sample_ID', all.x=T)
d <- d[!is.na(group.x) & !is.na(group.y) & group.x=='Primary' & group.y=='Primary']
d[,c('group.x','group.y'):=NULL]
d[,pt_a:=strtrim(Real_Sample_ID.x,3)]
d[,pt_b:=strtrim(Real_Sample_ID.y,3)]
d[pt_a==pt_b,same_or_different:='Same']
d[pt_a!=pt_b,same_or_different:='Different']

order_combo <- function(i, d) {
    out <- d[i,]
    t1 <- out$Real_Sample_ID.x
    t2 <- out$Real_Sample_ID.y
    sorted <- sort(c(t1,t2),decreasing=F)
    out$t1 <- sorted[1]
    out$t2 <- sorted[2]
    out
}
n <- nrow(d)
d <- rbindlist(lapply(1:n, order_combo, d))
d[,comparison:=paste(t1,'vs',t2)]

res <- d[patient!='E3',]
e3 <- d[patient=='E3']
e3_a_vs_b <- e3[patient=='E3' & pt_a %in% c('PTa','PTb') & pt_b %in% c('PTa','PTb'),]
e3_a_vs_b$patient <- 'E3 (A vs B)'
e3_a_vs_c <- e3[patient=='E3' & pt_a %in% c('PTa','PTc') & pt_b %in% c('PTa','PTc'),]
e3_a_vs_c$patient <- 'E3 (A vs C)'
e3_b_vs_c <- e3[patient=='E3' & pt_a %in% c('PTb','PTc') & pt_b %in% c('PTb','PTc'),]
e3_b_vs_c$patient <- 'E3 (B vs C)'
res <- rbind(res, e3_a_vs_b, e3_a_vs_c, e3_b_vs_c)
res$patient <- factor(res$patient)
res$same_or_different <- factor(res$same_or_different, levels=c('Same','Different'))
stat.test <- mywilcox2(res, r ~ same_or_different, paired=F, facet_field='patient', include_n=F)
cols <- c('#bfbfbf','steelblue')
names(cols) <- c('Same','Different')
p <- ggplot(res, aes(x=same_or_different, y=r)) +
    scale_y_continuous(breaks=seq(0,1,by=0.25), limits=c(0,1.15)) +
    geom_point(position=position_jitter(width=0.1, height=0, seed=42), pch=21, size=3.5, 
               aes(fill=same_or_different),color='white', stroke=0.25) +
    scale_fill_manual(values=cols) + 
    geom_boxplot(fill=NA,outlier.shape=NA,color='black',width=0.4) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.08) +
    facet_wrap(facets=~patient, scale='free_x') +
    guides(fill='none') +
    theme_bw(base_size=12) +
    labs(x='Samples from same/different PT', y='Coalescence ratio', title='ED Fig 2. Multi-PT CRs')
ggsave(here('figures/ed_fig_2.pdf'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ED Fig 3. Mouse data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m <- fread(here('original_data/misc/SDI_sliding_window_nBins_10_nIterations_10_RGB_sum_min_50_sameRGBexcluded_FALSE.csv'))
m[organ=='CAECUM', group:='Primary']
m[organ=='LIVER', group:='Liver']
m[organ=='PERI', group:='Peritoneum']
m$group <- factor(m$group, levels=c('Primary','Peritoneum','Liver'))
pixel_size <- 1.25e-6 # 1 pixel = 1.25um
pixel_area <- pixel_size^2 # um^2
m[,region_size_mm2:=(n_pixels * pixel_area) * (1e3)^2] # 1m^2 = (1000mm)^2
m$region_size_mm2 <- round(m$region_size_mm2, 3)
m[region_size_mm2 < 1, region_size_mm2:=1]

tst <- dunn_test_ES(m, SDI ~ group)
p <- ggplot(m, aes(x=group, y=SDI))  +
    geom_point(position=position_jitter(width=0.15, height=0, seed=42), aes(size=region_size_mm2, fill=group), pch=21, color='black', stroke=0.25) + 
    geom_boxplot(fill=NA, outlier.shape=NA, width=0.5) + 
    theme_ang(base_size=12) + 
    scale_size_area(breaks=c(1,5,10)) +   
    scale_fill_manual(values=group_cols, name='Tissue type') +
    stat_compare_means(method = "kruskal.test", label.y = 1.1, size=3, geom = "label") +
    stat_pvalue_manual(tst, label='label', y.position=c(0.95,1.0,1.05), size=3, tip.length=0) +
    labs(x='Sample type', y='SDI', title='ED Fig 3')
ggsave(here('figures/ed_fig_3.pdf'),width=6, height=4.5)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ED Fig XX (to add) correlation between primary tumor size and number of regions sampled
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d <- fread(here('original_data/misc/primary_sizes_and_num_regions.txt'))
tst <- cor.test(d$pt_size_cm, d$pt_regions_sampled, method='pearson')
fit <- lm(pt_regions_sampled ~ pt_size_cm, data=d)
coefs <- as.data.frame(summary(fit)$coef)
label1 <- paste0('Pearson R=',round(tst$estimate,2),', p=',prettyNum(tst$p.value,digits=2))
label2 <- paste0('fit line intercept=',round(coefs[1,1], 2),', slope=',round(coefs[2,1], 2))
label <- paste(label1,label2,sep='\n')
p <- ggplot(d, aes(x=pt_size_cm, y=pt_regions_sampled)) + 
    geom_point(pch=16,size=4,color='#008C45') + 
    geom_smooth(method='lm') +
    geom_text(data=d[1,], x=2.5, y=14,label=label,hjust=0) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x='Primary tumor longest dimension [cm]',y='N regions sampled') 
ggsave(here('figures/ed_fig_XX.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ED Fig XX (to add) intra-lesion in E15 between PTa and PTb
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
ad <- ad_table[['E15']]
da <- ad[grepl('PTa',rownames(ad)), grepl('PTa',colnames(ad))]
da[lower.tri(da,diag=T)] <- NA
da <- as.data.table(reshape2::melt(da))
db <- ad[grepl('PTb',rownames(ad)), grepl('PTb',colnames(ad))]
db[lower.tri(db,diag=T)] <- NA
db <- as.data.table(reshape2::melt(db))
da$primary <- 'PTa'
db$primary <- 'PTb'
d <- rbind(da, db)
d <- d[!is.na(value)]
d$group <- 'Primary'
stat.test <- mywilcox2(d, value ~ primary, paired=F)

p1 <- ggplot(d, aes(x=primary, y=value)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=16,size=4,aes(color=group)) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.3) +
    scale_color_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x='Primary tumor',y='Intra-lesion angular distance') 
ggsave(here('figures/ed_fig_XX.pdf'))


if(FALSE) { 
    ## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
    ad <- ad_table[['E10']]
    da <- ad[grepl('PTa',rownames(ad)), grepl('PTa',colnames(ad))]
    da[lower.tri(da,diag=T)] <- NA
    da <- as.data.table(reshape2::melt(da))
    db <- ad[grepl('PTb',rownames(ad)), grepl('PTb',colnames(ad))]
    db[lower.tri(db,diag=T)] <- NA
    db <- as.data.table(reshape2::melt(db))
    da$primary <- 'PTa'
    db$primary <- 'PTb'
    d <- rbind(da, db)
    d <- d[!is.na(value)]
    d$group <- 'Primary'
    stat.test <- mywilcox2(d, value ~ primary, paired=F)

    p2 <- ggplot(d, aes(x=primary, y=value)) + 
        geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=16,size=4,aes(color=group)) + 
        geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
        stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.3) +
        scale_color_manual(values=group_cols) + 
        guides(fill='none') +
        theme_ang(base_size=12) +
        labs(x='Primary tumor',y='Intra-lesion angular distance') 
    #ggsave(here('figures/ed_fig_XX.pdf'))


    ## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
    ad <- ad_table[['E11']]
    da <- ad[grepl('PTa',rownames(ad)), grepl('PTa',colnames(ad))]
    da[lower.tri(da,diag=T)] <- NA
    da <- as.data.table(reshape2::melt(da))
    db <- ad[grepl('PTb',rownames(ad)), grepl('PTb',colnames(ad))]
    db[lower.tri(db,diag=T)] <- NA
    db <- as.data.table(reshape2::melt(db))
    da$primary <- 'PTa'
    db$primary <- 'PTb'
    d <- rbind(da, db)
    d <- d[!is.na(value)]
    d$group <- 'Primary'
    stat.test <- mywilcox2(d, value ~ primary, paired=F)

    p3 <- ggplot(d, aes(x=primary, y=value)) + 
        geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=16,size=4,aes(color=group)) + 
        geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
        stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.3) +
        scale_color_manual(values=group_cols) + 
        guides(fill='none') +
        theme_ang(base_size=12) +
        labs(x='Primary tumor',y='Intra-lesion angular distance') 
    #ggsave(here('figures/ed_fig_XX.pdf'))


    ## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
    ad_file <- here(paste0('processed_data/angular_distance_matrices/E3.txt'))
    ad <- read_distance_matrix(ad_file)
    da <- ad[grepl('PTa',rownames(ad)), grepl('PTa',colnames(ad))]
    da[lower.tri(da,diag=T)] <- NA
    da <- as.data.table(reshape2::melt(da))
    db <- ad[grepl('PTb',rownames(ad)), grepl('PTb',colnames(ad))]
    db[lower.tri(db,diag=T)] <- NA
    db <- as.data.table(reshape2::melt(db))
    dc <- ad[grepl('PTc',rownames(ad)), grepl('PTc',colnames(ad))]
    dc[lower.tri(dc,diag=T)] <- NA
    dc <- as.data.table(reshape2::melt(dc))
    da$primary <- 'PTa'
    db$primary <- 'PTb'
    dc$primary <- 'PTc'
    d <- rbind(da, db, dc)
    d <- d[!is.na(value)]
    d$group <- 'Primary'
    stat.test <- mywilcox2(d, value ~ primary, paired=F)

    p3 <- ggplot(d, aes(x=primary, y=value)) + 
        geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=16,size=4,aes(color=group)) + 
        geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
        stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.3) +
        scale_color_manual(values=group_cols) + 
        guides(fill='none') +
        theme_ang(base_size=12) +
        labs(x='Primary tumor',y='Intra-lesion angular distance') 
    #ggsave(here('figures/ed_fig_XX.pdf'))

}




