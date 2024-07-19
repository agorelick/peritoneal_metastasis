source(here::here('R/func.R'))
options(repr.plot.width=7, repr.plot.height=7)
ncpus = 4


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make annotated poly-G trees for each patient
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_tree <- function(this.patient, sample_info, collapsed, outdir, show.depth=F, show.timing=F, show.bsvals=F, tree.layout='ape') { 
    set.seed(42)
    ad_file <- here(paste0('processed_data/angular_distance_matrices/',this.patient,'.txt'))
    ad <- read_distance_matrix(ad_file)

    if(!dir.exists(outdir)) dir.create(outdir, recursive=T)
    message(this.patient)
    if(this.patient!='E3') {
        si <- sample_info[Patient_ID %in% this.patient]
    } else {
        si <- sample_info[grepl('E[a-c]3',Patient_ID),]
        si <- si[!duplicated(Sample_ID),]   
    }
    if(collapsed) si <- si[in_collapsed==T,]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    ad <- ad[valid_samples, valid_samples]

    ## get the NJ tree
    get_tree <- function(mat) { 
        tree <- nj(mat)
        tree <- ape::root(tree, outgroup=grep('^N',tree$tip.label))
        tree$tip.label <- gsub('Normal','N',tree$tip.label)
        tree
    }
    tree <- get_tree(ad)
  
    ## truncating normal branch if it's too long
    truncated_normal <- F
    normal_branch <- which(tree$edge[,2]==grep('^N',tree$tip.label))
    other_lengths <- tree$edge.length[-normal_branch]

    if(tree$edge.length[normal_branch] > 1.5*max(other_lengths)) {
        truncated_normal <- T
        tree$edge.length[normal_branch] <- 1.5*max(other_lengths)
    }  
    bsvals <- F 
    bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/',this.patient,'.rds'))

    if(file.exists(bs_file) & show.bsvals==T) {
        message('loading bootstrapped angular distance matrices ...')
        bs_patient <- readRDS(bs_file)
        ## subset each bs_matrix for the valid samples
        f=function(bs_patient, valid_samples) { 
            bs_patient[valid_samples, valid_samples]
        }
        bs_list <- lapply(bs_patient, f, valid_samples)
        tree_list <- lapply(bs_list, get_tree)
        bstrees <- TreeTools::as.multiPhylo(tree_list)
        tree <- addConfidences(tree, bstrees) 
        tree$node.label <- round(100*tree$node.label)
        bsvals <- T
    }

    if(tree.layout=='ape' & show.bsvals==T) tree$node.label[tree$node.label < 50] <- NA
    groups <- si[,c('Real_Sample_ID','group','vertical','met_timing'),with=F]
    groups[met_timing=='metachronous after synchronous', met_timing:='metachronous']
    setnames(groups,'Real_Sample_ID','label')
    groups$label <- gsub('Normal','N',groups$label)
    cols <- c(group_cols,'blue')
    names(cols)[length(cols)] <- 'bsval'

    p <- ggtree(tree, layout=tree.layout) 
    p <- p %<+% groups
    pd <- as.data.table(p$data)

    if(bsvals==T) { 
        pd[isTip==F,group:='bsval']
    }

    if(tree.layout=='ape') {
        p <- p + geom_text_repel(data=pd, aes(x=x, y=y, label=label, color=group), min.segment.length=0.2,max.overlaps=100) 
    } else {
        p <- p + geom_text(data=pd, aes(x=x, y=y, label=label, color=group), hjust=-0.25, vjust=0.5)
    }
    p <- p + scale_color_manual(values=cols) + guides(color='none') 
    p <- p + guides(color='none')

    if(show.depth) {
        pd <- as.data.table(p$data)
        pd <- pd[!is.na(vertical) & vertical %in% c('deep','mucosal/luminal'),]
        depth_cols <- c('#2e388e','#be1e2d'); names(depth_cols) <- c('mucosal/luminal','deep')
        p <- p + geom_tippoint(data=pd, aes(fill=vertical), pch=21, stroke=0.25, color='black', size=2.5) + 
            scale_fill_manual(values=depth_cols, name='PT region depth')
    }

    if(show.timing) {
        pd <- as.data.table(p$data)
        pd <- pd[!is.na(met_timing) & met_timing %in% c('synchronous','metachronous'),]
        timing_shapes <- c(1,4); names(timing_shapes) <- c('synchronous','metachronous')
        p <- p + geom_tippoint(data=pd, aes(pch=met_timing), size=2.5) + scale_shape_manual(values=timing_shapes,name='Metastasis timing')
    }

    tree_title <- paste0(this.patient,'. BS values 50%+ shown')
    if(truncated_normal) tree_title <- paste0(tree_title,'. Normal branch truncated.')
    p <- p + theme(legend.position='bottom') + ggtitle(tree_title)
    outfile <- file.path(outdir,paste0(this.patient,'.pdf'))

    if(tree.layout=='ape') 
        ggsave(outfile, plot=p, width=10, height=8)
    else
        ggsave(outfile, plot=p, width=8, height=10)
}

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
    #theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank())+
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
ggsave(here('figures/fig_1a_overview.pdf'),width=9, height=9)



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
ggsave(here('figures/fig_1cd_C161_tree_comparison.pdf'),width=9, height=6)



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
ggsave(here('figures/fig_1e_scna_polyg_similarity.pdf'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reload polyG and kim et al data for following
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sample_info <- fread(here('processed_data/sample_info.txt'))
sample_info <- sample_info[cohort!='lung',]
sample_info$met_treated_type <- ''
sample_info[met_treated %in% c('hipec','hipec after untreated'), met_treated_type:='hipec']
sample_info[met_treated %in% c('systemic chemo','systemic chemo after untreated'), met_treated_type:='systemic chemo']
sample_info[met_treated %in% c('systemic chemo','hipec'), met_treated:='treated']
sample_info[met_treated %in% c('systemic chemo after untreated','hipec after untreated'), met_treated:='treated after untreated']

## load uncollapsed angular distance matrices into a named list
load_ad_matrix <- function(patient) {
    ad_file <- here(paste0('processed_data/angular_distance_matrices/',patient,'.txt'))
    if(!file.exists(ad_file)) {
        ad <- NULL
    } else {
        ad <- read_distance_matrix(ad_file)
    }
    ad
}
patients <- unique(sample_info[!grepl('CRC',Patient_ID), (Patient_ID)])
ad_table <- lapply(patients, load_ad_matrix)
names(ad_table) <- patients

## get Kim et al WES tree topologies
kim_node_distances <- load_kim_node_distances()



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
ggsave(here('figures/fig_2e_tissue_type_rds.pdf'))


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
ggsave(here('figures/fig_2f_tissue_type_angular_distance.pdf'))



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
ggsave(here('figures/fig_2g_tissue_type_intralesion_angular_distance.pdf'))



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
ggsave(here('figures/fig_3a_tissue_type_rds_untreated_synchronous.pdf'))



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
# script to generate CR table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)

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

# saving table with only polyG correlation
non_boot_l1 %>%
  select(a, b, r) %>% 
  write_tsv(here("processed_data/per_multitumor_cr.tsv"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load coalescence-ratio data for multi-PT patients
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
source(here('R/func.R'))

d <- fread(here('processed_data/per_multitumor_cr.tsv'))
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
#res[,patient:=paste(patient,'(A vs B)')]
e3 <- d[patient=='E3']
e3_a_vs_b <- e3[patient=='E3' & pt_a %in% c('PTa','PTb') & pt_b %in% c('PTa','PTb'),]
e3_a_vs_b$patient <- 'E3 (A vs B)'
e3_a_vs_c <- e3[patient=='E3' & pt_a %in% c('PTa','PTc') & pt_b %in% c('PTa','PTc'),]
e3_a_vs_c$patient <- 'E3 (A vs C)'
e3_b_vs_c <- e3[patient=='E3' & pt_a %in% c('PTb','PTc') & pt_b %in% c('PTb','PTc'),]
e3_b_vs_c$patient <- 'E3 (B vs C)'
res <- rbind(res, e3_a_vs_b, e3_a_vs_c, e3_b_vs_c)

res$patient <- factor(res$patient)#, levels=c('E3','E10','E11','E15'))
stat.test <- mywilcox2(res, r ~ same_or_different, paired=F, facet_field='patient', include_n=F)

p <- ggplot(res, aes(x=same_or_different, y=r)) +
    scale_y_continuous(breaks=seq(0,1,by=0.25), limits=c(0,1.15)) +
    geom_point(position=position_jitter(width=0.1, height=0, seed=42), pch=16, size=3.5, aes(color=same_or_different)) +
    geom_boxplot(fill=NA,outlier.shape=NA,color='black',width=0.4) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.08) +
    facet_wrap(facets=~patient) +
    guides(fill='none') +
    theme_bw(base_size=12) +
    labs(x='Samples from same/different PT', y='Correlation coefficient')
ggsave(here('figures/multi_primary_CRs.pdf'), width=8, height=6)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Untreated version: Fig 2F. pairwise AD compared between Liv/PM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T]
si1 <- si1[group %in% c('Normal','Primary') | (met_treated=='untreated')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated=='untreated')]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
res$category <- 'untreated'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
write_tsv(out,here('processed_data/misc/interlesion_ad_untreated.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Pairwise angular distance vs tissue type (untreated)', 
         subtitle='NatGen, Science, Peritoneum (wilcox test)')
ggsave(here('figures/unk_fig_2f_tissue_type_angular_distance_untreated.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Synchronous+untreated version: Fig 2F. pairwise AD compared between LR/Liv/PM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T]
si1 <- si1[group %in% c('Normal','Primary') | (met_timing=='synchronous' & met_treated=='untreated')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

si2 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Locoregional') & in_collapsed==T]
si2 <- si2[group %in% c('Normal','Primary') | (met_timing=='synchronous' & met_treated=='untreated')]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Locoregional'; res2$group2 <- 'Locoregional'; res2$class <- 'LR:LR'

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si3 <- si3[group %in% c('Normal','Primary') | (met_timing=='synchronous' & met_treated=='untreated')]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'

res <- rbind(res1, res2, res3)
res$class <- factor(res$class, levels=c('LR:LR','Per:Per','Liv:Liv'))
tst <- dunn_test_ES(res, distance ~ class)

table(res$class)
x <- as.data.frame.matrix(xtabs(~ patient + class, data=res))
sum(x$`LR:LR` > 0)
sum(x$`Per:Per` > 0)
sum(x$`Liv:Liv` > 0)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method = "kruskal.test", label.y = 1.8, geom = "label") + 
    stat_pvalue_manual(tst, label='label', y.position=c(1.5,1.7,1.9)) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Pairwise angular distance vs tissue type (untreated+synchronous)', 
         subtitle='NatGen, Science, Peritoneum (Dunn\'s test + Holm\'s)')
ggsave(here('figures/unk_fig_2f_tissue_type_angular_distance_untreated_synchronous_3group.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Synchronous+untreated version: Fig 2F. pairwise AD compared between Liv/PM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T]
si1 <- si1[group %in% c('Normal','Primary') | (met_timing=='synchronous' & met_treated=='untreated')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si3 <- si3[group %in% c('Normal','Primary') | (met_timing=='synchronous' & met_treated=='untreated')]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
res$category <- 'untreated+synchronous'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
write_tsv(out,here('processed_data/misc/interlesion_ad_untreated_synchronous.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Pairwise angular distance vs tissue type (untreated+synchronous)', 
         subtitle='NatGen, Science, Peritoneum (wilcox test)')
ggsave(here('figures/unk_fig_2f_tissue_type_angular_distance_untreated_synchronous.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Systemic chemo version: Fig 2F. pairwise AD compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T]
si1 <- si1[group %in% c('Normal','Primary') | (met_treated_type=='systemic chemo')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated_type=='systemic chemo')]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)
res$category <- 'prior systemic chemo'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
write_tsv(out,here('processed_data/misc/interlesion_ad_systemic_chemo.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
tmp_cols <- tmp_cols[names(tmp_cols) %in% res$tissue_type]
tmp_cols <- tmp_cols[c('Lymph node','Tumor deposit','Peritoneum','Liver')]

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Pairwise angular distance vs tissue type (Prior systemic chemo)', 
         subtitle='NatGen, Science, Peritoneum (Wilcox test)')
ggsave(here('figures/unk_fig_2f_tissue_type_angular_distance_systemic_chemo.pdf'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# any treatment version: Fig 2F. pairwise AD compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T]
si1 <- si1[group %in% c('Normal','Primary') | (met_treated %in% c('treated','treated after untreated'))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated %in% c('treated','treated after untreated'))]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)
res$category <- 'any prior treatment'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
write_tsv(out,here('processed_data/misc/interlesion_ad_any_treatment.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
tmp_cols <- tmp_cols[names(tmp_cols) %in% res$tissue_type]
tmp_cols <- tmp_cols[c('Lymph node','Tumor deposit','Peritoneum','Liver')]

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Pairwise angular distance vs tissue type (Any prior treatment)', 
         subtitle='NatGen, Science, Peritoneum (Wilcox test)')
ggsave(here('figures/unk_fig_2f_tissue_type_angular_distance_treated.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# correlation between primary tumor size and number of regions sampled
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
ggsave(here('figures/unk_primary_size_vs_regions.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# intra-lesion in E15 between PTa and PTb
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

p <- ggplot(d, aes(x=primary, y=value)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=16,size=4,aes(color=group)) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.3) +
    scale_color_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x='Primary tumor',y='Intra-lesion angular distance',subtitle='E15 pairwise intra-lesion angular distance in primary tumors A and B') 
ggsave(here('figures/unk_E15_primaries_intralesion_ad.pdf'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Untreated intra-lesion: pairwise AD compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subset_for_intralesion <- function(res) {
    ## annotate each sample with the type+lesion, so that we can find comparisons of
    ## different samples within the same type+lesion
    bc1 <- parse_barcode(res$sample1)
    bc1[,typelesion:=paste0(type,lesion)]
    bc2 <- parse_barcode(res$sample2)
    bc2[,typelesion:=paste0(type,lesion)]
    res$typelesion1 <- bc1$typelesion
    res$typelesion2 <- bc2$typelesion
    res <- res[typelesion1==typelesion2]
    res
}

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')]
si1 <- si1[group %in% c('Normal','Primary') | (met_treated=='untreated')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'
res1 <- subset_for_intralesion(res1)

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated=='untreated')]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'
res3 <- subset_for_intralesion(res3)

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)
res$category <- 'untreated'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
write_tsv(out,here('processed_data/misc/intralesion_ad_untreated.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Intra-lesion angular distance',title='Pairwise intra-lesion angular distance vs tissue type (untreated)', 
         subtitle='NatGen, Science, Peritoneum (Wilcox test)')
ggsave(here('figures/unk_intralesion_tissue_type_angular_distance_untreated.pdf'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Untreated+synchronous intra-lesion: pairwise AD compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subset_for_intralesion <- function(res) {
    ## annotate each sample with the type+lesion, so that we can find comparisons of
    ## different samples within the same type+lesion
    bc1 <- parse_barcode(res$sample1)
    bc1[,typelesion:=paste0(type,lesion)]
    bc2 <- parse_barcode(res$sample2)
    bc2[,typelesion:=paste0(type,lesion)]
    res$typelesion1 <- bc1$typelesion
    res$typelesion2 <- bc2$typelesion
    res <- res[typelesion1==typelesion2]
    res
}

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')]
si1 <- si1[group %in% c('Normal','Primary') | (met_treated=='untreated' & met_timing=='synchronous')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'
res1 <- subset_for_intralesion(res1)

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated=='untreated' & met_timing=='synchronous')]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'
res3 <- subset_for_intralesion(res3)

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)
res$category <- 'untreated+synchronous'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
#write_tsv(out,here('processed_data/misc/intralesion_ad_untreated_synchronous.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Intra-lesion angular distance',title='Pairwise intra-lesion angular distance vs tissue type (synchronous+untreated)', 
         subtitle='NatGen, Science, Peritoneum (Wilcox test)')
ggsave(here('figures/unk_intralesion_tissue_type_angular_distance_untreated_synchronous.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Systemic chemo version: pairwise AD compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')]
si1 <- si1[group %in% c('Normal','Primary') | (met_treated_type=='systemic chemo')]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'
res1 <- subset_for_intralesion(res1)

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated_type=='systemic chemo')]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'
res3 <- subset_for_intralesion(res3)

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)
res$category <- 'systemic chemo'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
write_tsv(out,here('processed_data/misc/intralesion_ad_systemic_chemo.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
tmp_cols <- tmp_cols[names(tmp_cols) %in% res$tissue_type]
tmp_cols <- tmp_cols[c('Lymph node','Tumor deposit','Peritoneum','Liver')]

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Intra-lesion angular distance',title='Pairwise intra-lesion angular distance vs tissue type (Prior systemic chemo)', 
         subtitle='NatGen, Science, Peritoneum (wilcox test)')
ggsave(here('figures/unk_intralesion_tissue_type_angular_distance_systemic.pdf'))





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Any prior treatment version: pairwise AD compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')]
si1 <- si1[group %in% c('Normal','Primary') | (met_treated %in% c('treated','treated after untreated'))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'
res1 <- subset_for_intralesion(res1)

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
si3 <- si3[group %in% c('Normal','Primary') | (met_treated %in% c('treated','treated after untreated'))]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'
res3 <- subset_for_intralesion(res3)

res <- rbind(res1, res3)
res$class <- factor(res$class, levels=c('Per:Per','Liv:Liv'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)
res$category <- 'any prior treatment'
res[,pair:=paste0(sample1,',',sample2)]
out <- res[,c('patient','pair','distance','group1','category'),with=F]
setnames(out,c('distance','group1'),c('ad','group'))
write_tsv(out,here('processed_data/misc/intralesion_ad_any_treatment.txt'))
stat.test <- mywilcox2(res, distance ~ class, paired=F)

tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
tmp_cols <- tmp_cols[names(tmp_cols) %in% res$tissue_type]
tmp_cols <- tmp_cols[c('Lymph node','Tumor deposit','Peritoneum','Liver')]

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, y.position=1.65) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Intra-lesion angular distance',title='Pairwise intra-lesion angular distance vs tissue type (Any prior treatment)', 
         subtitle='NatGen, Science, Peritoneum (wilcox test)')
ggsave(here('figures/unk_intralesion_tissue_type_angular_distance_any_treated.pdf'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SI Fig 2A. pairwise nodes compared between tissue types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si1 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum') & in_collapsed==T]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='node', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per'

si2 <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Locoregional') & in_collapsed==T]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='node', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Locoregional'; res2$group2 <- 'Locoregional'; res2$class <- 'LR:LR'

si3 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='node', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Liver'; res3$group2 <- 'Liver'; res3$class <- 'Liv:Liv'

si4 <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver') & in_collapsed==T]
res4 <- get_met_specific_distances(si4, kim_node_distances, comparison='pairwise', distance='node', return_tree=F)
res4 <- res4[group1=='Metastasis' & group2=='Metastasis']
res4$group1 <- 'Liver'; res4$group2 <- 'Liver'; res4$class <- 'Liv:Liv'

res <- rbind(res1, res2, res3, res4)
res$class <- factor(res$class, levels=c('LR:LR','Per:Per','Liv:Liv'))
tst <- dunn_test_ES(res, distance ~ class)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method = "kruskal.test", label.y = 2.5, geom = "label") + 
    stat_pvalue_manual(tst, label='label', y.position=c(2.3,2.45,2.6)) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Normalized node distance',title='Pairwise node distance vs tissue type', 
         subtitle='NatGen, Science, Peritoneum + Kim et al. (Dunn\'s test + Holm\'s)')
ggsave(here('figures/ed_fig_2a_tissue_type_node_distance.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 4C, ed_fig_5f. pairwise AD for PMs by timing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; all liver-met patients; all LR-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

## run for Per synchronous vs metachronous vs metachronous-after-synchronous
si1 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='synchronous'))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per (synchronous)'
res1$type <- 'Per_synchronous'
si2 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing %in% 'metachronous'))]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Peritoneum'; res2$group2 <- 'Peritoneum'; res2$class <- 'Per:Per (metachronous)'
res2$type <- 'Per_metachronous'
si3 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing %in% 'metachronous after synchronous'))]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Peritoneum'; res3$group2 <- 'Peritoneum'; res3$class <- 'Per:Per (metachronous after synchronous)'
res3$type <- 'Per_metachronous'
res_per <- rbind(res1, res2, res3)

## AD plot but only for the peritoneal mets
res <- copy(res_per)
res$class <- factor(res$class, levels=unique(res$class))
res[grepl('(synchronous)',class),class3:='synchronous']
res[grepl('(metachronous)',class),class3:='metachronous']
res[grepl('(metachronous after synchronous)',class),class3:='metachronous after synchronous']
res[,timing:=class3]
res[class3=='metachronous after synchronous', timing:='metachronous']
res$class3 <- factor(res$class3, levels=c('synchronous','metachronous','metachronous after synchronous'))

tst <- dunn_test_ES(res, distance ~ class3)
kw.p <- prettyNum(kruskal.test(distance ~ class3, data=res)$p.value, digits=2)

p <- ggplot(res, aes(x=class3, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,size=4,aes(fill=group1),stroke=0.5,color='black') +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    geom_label(data=res_per[1,], label=paste0('Kruskal-Wallis, p=',kw.p), aes(x=1.2,y=1.8)) + 
    stat_pvalue_manual(tst, label='label', y.position=seq(1.5,by=0.1,length.out=nrow(tst)),tip.length=NA) +
    scale_fill_manual(values=group_cols) + 
    scale_color_manual(values=cols_timing) + 
    guides(fill='none',color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x=NULL,y='Angular distance',title='Peritoneal met pairwise angular distance vs met timing',subtitle='All patients with PMs')
ggsave(here('figures/fig_4c_timing_per_angular_distance.pdf'))

collapse <- function(res) {
    mu <- mean(res$distance)
    list(mu=mu)
}
tmp <- res[class3 %in% c('synchronous','metachronous'),collapse(.SD),by=c('patient','timing','class3')]
tmp$group <- 'Peritoneum'
stat.test <- mywilcox2(tmp, mu ~ class3, paired=F)

p <- ggplot(tmp, aes(x=class3, y=mu)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Peritoneal met pairwise AD vs met timing, patient-average',subtitle='All patients with PMs')
ggsave(here('figures/ed_fig_5d_timing_per_angular_distance_average.pdf'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 4d. Pairwwise AD for PMs subset by treatment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; all liver-met patients; all LR-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

## run for Per synchronous vs metachronous vs metachronous-after-synchronous
si1 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated'))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per (untreated)'
res1$type <- 'Per_untreated'
si2 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% c('treated','treated after untreated') & met_treated_type=='systemic chemo'))]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Peritoneum'; res2$group2 <- 'Peritoneum'; res2$class <- 'Per:Per (systemic chemo)'
res2$type <- 'Per_systemics'
si3 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% c('treated','treated after untreated') & met_treated_type=='hipec'))]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Peritoneum'; res3$group2 <- 'Peritoneum'; res3$class <- 'Per:Per (hipec)'
res3$type <- 'Per_hipec'
res_per <- rbind(res1, res2, res3)

## AD plot but only for the peritoneal mets
res <- copy(res_per)
res$class <- factor(res$class, levels=unique(res$class))

tst <- dunn_test_ES(res, distance ~ class)
kw.p <- prettyNum(kruskal.test(distance ~ class, data=res)$p.value, digits=2)

p <- ggplot(res, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,size=4,aes(fill=group1),stroke=0.5,color='black') +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    geom_label(data=res_per[1,], label=paste0('Kruskal-Wallis, p=',kw.p), aes(x=1.2,y=1.8)) + 
    stat_pvalue_manual(tst, label='label', y.position=seq(1.5,by=0.1,length.out=nrow(tst)),tip.length=NA) +
    scale_fill_manual(values=group_cols) + 
    scale_color_manual(values=cols_timing) + 
    guides(fill='none',color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x=NULL,y='Angular distance',title='Peritoneal met pairwise angular distance vs met timing',subtitle='All patients with PMs')
ggsave(here('figures/fig_4d_treatment_pm_angular_distance.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4 way comparison of RDS between synchronous/metachronous, 
# run for Per untreated vs treated vs treated-after-untreated
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

si1 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated' & met_timing=='synchronous'))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='rds', distance='node', return_tree=F)
res1 <- res1[type=='Met',]
res1$group <- 'Peritoneum'; 
res1$class <- 'Untreated, sync'
res1$type <- 'Per_untreated_sync'

si2 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated' & met_timing=='metachronous'))]
res2 <- get_met_specific_distances(si2, ad_table, comparison='rds', distance='node', return_tree=F)
res2 <- res2[type=='Met',]
res2$group <- 'Peritoneum'; 
res2$class <- 'Untreated, meta'
res2$type <- 'Per_untreated_meta'

si3 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated_type=='systemic chemo' & met_timing=='synchronous'))]
res3 <- get_met_specific_distances(si3, ad_table, comparison='rds', distance='node', return_tree=F)
res3 <- res3[type=='Met',]
res3$group <- 'Peritoneum'; 
res3$class <- 'Systemic chemo, sync'
res3$type <- 'Per_treated_sync'

si4 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated_type=='systemic chemo' & met_timing=='metachronous'))]
res4 <- get_met_specific_distances(si4, ad_table, comparison='rds', distance='node', return_tree=F)
res4 <- res4[type=='Met',]
res4$group <- 'Peritoneum'; 
res4$class <- 'Systemic chemo, meta'
res4$type <- 'Per_treated_meta'

si5 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='treated after untreated' & met_treated_type=='systemic chemo' & met_timing=='metachronous after synchronous'))]
res5 <- get_met_specific_distances(si5, ad_table, comparison='rds', distance='node', return_tree=F)
res5 <- res5[type=='Met',]
res5$group <- 'Peritoneum'; 
res5$class <- 'Systemic chemo, meta-after-sync'
res5$type <- 'Per_treated_meta'

si6 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% 'treated' & met_treated_type=='hipec' & met_timing=='metachronous'))]
res6 <- get_met_specific_distances(si6, ad_table, comparison='rds', distance='node', return_tree=F)
res6 <- res6[type=='Met',]
res6$group <- 'Peritoneum'; 
res6$class <- 'HIPEC, meta'
res6$type <- 'Per_treated_meta'

si7 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% c('treated after untreated') & met_treated_type=='hipec' & met_timing=='metachronous after synchronous'))]
res7 <- get_met_specific_distances(si7, ad_table, comparison='rds', distance='node', return_tree=F)
res7 <- res7[type=='Met',]
res7$group <- 'Peritoneum'; 
res7$class <- 'HIPEC after untreated, meta-after-sync'
res7$type <- 'Per_treated_meta_after'

res_per <- rbind(res1, res2, res3, res4, res5, res6, res7) 
res_per$class <- factor(res_per$class, levels=unique(res_per$class))

p <- ggplot(res_per, aes(x=class, y=RDS)) + 
    scale_y_continuous(breaks=seq(0,1,by=0.25)) +
    geom_beeswarm(cex = 3, aes(color=group), size=4, pch=16) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    geom_text_repel(position=position_jitter(width=0.15,height=0,seed=2),color='black',size=4,aes(label=patient)) + 
    scale_color_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x=NULL,y='RDS',title='Met-specific RDS vs treatment and timing', 
         subtitle='All peritoneal met patients')
ggsave(here('figures/unused_fig_treatment_and_timing_per_rds.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pairwise AD in PM vs treatment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89; and all liver-met patients)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

## run for Per untreated vs treated vs treated-after-untreated
si1 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated'))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res1 <- res1[group1=='Metastasis' & group2=='Metastasis']
res1$group1 <- 'Peritoneum'; res1$group2 <- 'Peritoneum'; res1$class <- 'Per:Per (untreated)'
res1$type <- 'Per_untreated'

si2 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% 'treated'))]
res2 <- get_met_specific_distances(si2, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res2 <- res2[group1=='Metastasis' & group2=='Metastasis']
res2$group1 <- 'Peritoneum'; res2$group2 <- 'Peritoneum'; res2$class <- 'Per:Per (treated)'
res2$type <- 'Per_treated'

si3 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T &
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated %in% 'treated after untreated'))]
res3 <- get_met_specific_distances(si3, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res3 <- res3[group1=='Metastasis' & group2=='Metastasis']
res3$group1 <- 'Peritoneum'; res3$group2 <- 'Peritoneum'; res3$class <- 'Per:Per (treated after untreated)'
res3$type <- 'Per_treated'

res_per <- rbind(res1, res2, res3)
res_per$class <- factor(res_per$class, levels=unique(res_per$class))
tst <- dunn_test_ES(res_per, distance ~ class)
kw.p <- prettyNum(kruskal.test(distance ~ class, data=res_per)$p.value, digits=2)

## ED Figure 5e, legacy treated/untreated/treated-after-untreated AD
p <- ggplot(res_per, aes(x=class, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA)+ 
    geom_label(data=res_per[1,], label=paste0('Kruskal-Wallis, p=',kw.p), aes(x=1.2,y=2.2)) + 
    stat_pvalue_manual(tst, label='label', y.position=c(1.6,1.7,1.8), tip.length=NA) + 
    scale_fill_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x=NULL,y='Angular distance',title='Pairwise angular distance vs treatment', 
         subtitle='All peritoneal met patients (Dunn\'s test + Holm\'s)')
ggsave(here('figures/ed_fig_5e_treatment_pm_angular_distance.pdf'))



## run for Per untreated vs treated vs treated-after-untreated
si1 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated'))]
res1 <- get_met_specific_distances(si1, ad_table, comparison='rds', distance='node', return_tree=F)
res1 <- res1[type=='Met',]
res1$group <- 'Peritoneum'; res1$class <- 'Per (untreated)'
res1$type <- 'Per_untreated'

si2 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='treated'))]
res2 <- get_met_specific_distances(si2, ad_table, comparison='rds', distance='node', return_tree=F)
res2 <- res2[type=='Met',]
res2$group <- 'Peritoneum'; res2$class <- 'Per (treated)'
res2$type <- 'Per_treated'

si3 <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & 
                   (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='treated after untreated'))]
res3 <- get_met_specific_distances(si3, ad_table, comparison='rds', distance='node', return_tree=F)
res3 <- res3[type=='Met',]
res3$group <- 'Peritoneum'; res3$class <- 'Per (treated after untreated)'
res3$type <- 'Per_treated'

res_per <- rbind(res1, res2, res3)
res_per$class <- factor(res_per$class, levels=unique(res_per$class))
tst <- dunn_test_ES(res_per, RDS ~ class)
kw.p <- prettyNum(kruskal.test(RDS ~ class, data=res_per)$p.value, digits=2)

p <- ggplot(res_per, aes(x=class, y=RDS)) + 
    scale_y_continuous(breaks=seq(0,1,by=0.25), limits=c(0,1.40)) + 
    geom_beeswarm(cex = 3, aes(color=group), size=4, pch=16) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA)+ 
    geom_text(position=position_jitter(width=0.15,height=0,seed=2),size=3, hjust=1,vjust=1, aes(label=patient)) + 
    geom_label(data=res_per[1,], label=paste0('Kruskal-Wallis, p=',kw.p), aes(x=1.2,y=1.3)) + 
    stat_pvalue_manual(tst, label='label', y.position=c(1.1,1.2,1.3), tip.length=NA) +
    scale_color_manual(values=group_cols) + 
    guides(fill='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x=NULL,y='RDS',title='Met-specific RDS vs treatment', 
         subtitle='All peritoneal met patients (Dunn\'s test + Holm\'s)')
ggsave(here('figures/unused_fig_treatment_per_rds.pdf'))





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# original version: synchronous/untreated intra-lesion heterogeneity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subset_for_intralesion <- function(res) {
    ## annotate each sample with the type+lesion, so that we can find comparisons of
    ## different samples within the same type+lesion
    bc1 <- parse_barcode(res$sample1)
    bc1[,typelesion:=paste0(type,lesion)]
    bc2 <- parse_barcode(res$sample2)
    bc2[,typelesion:=paste0(type,lesion)]
    res$typelesion1 <- bc1$typelesion
    res$typelesion2 <- bc2$typelesion
    res <- res[typelesion1==typelesion2]
    res
}

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='synchronous')]
res_syn <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_syn <- res_syn[group1=='Metastasis' & group2=='Metastasis']
res_syn[res_syn=='Metastasis'] <- 'Peritoneum'
res_syn <- subset_for_intralesion(res_syn)
res_syn$timing <- 'Synchronous'

si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='metachronous')]
res_meta <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_meta <- res_meta[group1=='Metastasis' & group2=='Metastasis']
res_meta[res_meta=='Metastasis'] <- 'Peritoneum'
res_meta <- subset_for_intralesion(res_meta)
res_meta$timing <- 'Metachronous'

## unused fig
res <- rbind(res_syn, res_meta)
res$timing <- factor(res$timing, levels=c('Synchronous','Metachronous'))
p_ad <- ggplot(res, aes(x=timing, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method="wilcox.test", comparisons=list(c("Synchronous","Metachronous"))) +
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Angular distance',title='Intralesion angular distance',subtitle='All patients with PMs')
ggsave(here('figures/unused_fig_timing_per_intralesion_angular_distance.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compare heterogeneity within PT between sync and meta PM patients
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89)
all_sync_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='synchronous',(Patient_ID)])
all_meta_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='metachronous',(Patient_ID)])

## get pairwise inter-lesion AD
si_sync <- sample_info[Patient_ID %in% all_sync_patients & group %in% c('Normal','Primary') & in_collapsed==T]
res_sync <- get_met_specific_distances(si_sync, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_sync$timing <- 'Synchronous PMs only'
si_meta <- sample_info[Patient_ID %in% all_meta_patients & group %in% c('Normal','Primary') & in_collapsed==T]
res_meta <- get_met_specific_distances(si_meta, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_meta$timing <- 'Metachronous PMs only'
res <- rbind(res_sync, res_meta)
res$timing <- factor(res$timing, levels=c('Synchronous PMs only','Metachronous PMs only'))
stat.test <- mywilcox2(res, distance ~ timing, paired=F)

## unused fig
p_inter_sync <- ggplot(res, aes(x=timing, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) + 
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x='PT samples among patients with ...',y='Angular distance',subtitle='Inter-sample AD, PT samples from patients with synch/meta PMs')
ggsave(here('figures/unused_fig_timing_pt_intersample_angular_distance.pdf'))


## compare pairwise AD in synchronous LNs between patients who go on to have (metachronous) PMs
si_sync <- sample_info[Patient_ID %in% all_sync_patients & (group %in% c('Normal','Primary') | (tissue_type=='Lymph node' & met_timing=='synchronous')) & in_collapsed==T]
res_sync <- get_met_specific_distances(si_sync, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_sync <- res_sync[group1=='Metastasis' & group2=='Metastasis']
res_sync$timing <- 'Synchronous PMs only'

si_meta <- sample_info[Patient_ID %in% all_meta_patients & (group %in% c('Normal','Primary') | (tissue_type=='Lymph node' & met_timing=='synchronous')) & in_collapsed==T]
si_meta <- si_meta[in_collapsed==T]
res_meta <- get_met_specific_distances(si_meta, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_meta <- res_meta[group1=='Metastasis' & group2=='Metastasis']
res_meta$timing <- 'Metachronous PMs only'

res <- rbind(res_sync, res_meta)
res$timing <- factor(res$timing, levels=c('Synchronous PMs only','Metachronous PMs only'))
res$group1 <- 'Locoregional'
stat.test <- mywilcox2(res, distance ~ timing, paired=F)


## unused fig
p_LN <- ggplot(res, aes(x=timing, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) + 
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x='Synchronous LN samples among patients with ...',y='Angular distance',subtitle='Inter-lesion AD, LN samples from patients with synchronous/metachronous PMs')
ggsave(here('figures/unused_fig_timing_ln_intersample_angular_distance.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test for a shared origin between sync/meta PMs and synchronous tumor deposit
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89)
all_sync_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='synchronous',(Patient_ID)])
all_meta_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='metachronous',(Patient_ID)])

## check for shared origin between PMs and Synchronous LNs between Sync and Meta PM patients
test_patient_met_level <- function(patient, si, ad_table, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary','Peritoneum', 'Locoregional')]

    if(patient %in% c('Ea3','Eb3','Ec3')) {    
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/E3.rds'))
        bs <- readRDS(bs_file)
    } else {
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/',patient,'.rds'))
        bs <- readRDS(bs_file)
    }

    ## get the original set of valid-samples for the true number of distant/pt/pm samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    #n_td <- length(unique(si$Real_Sample_ID[si$group=='Locoregional']))
    n_td <- length(unique(si$Real_Sample_ID[si$tissue_type=='Tumor deposit']))
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))
    n_pm <- length(unique(si$Real_Sample_ID[si$group=='Peritoneum']))

    test_patient_bs <- function(i, patient, si, ad_table, bs) {
        if(i > 0) {
            ## for the bootstrap replicates, obtain new ad and valid_samples
            ad <- bs[[i]]
            valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
        }
        ad <- ad[valid_samples, valid_samples]
        ad <- as.data.table(reshape2::melt(ad))    
        groups <- si[,c('Real_Sample_ID','group','tissue_type'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        setnames(ad,'tissue_type','tissue_type1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        setnames(ad,'tissue_type','tissue_type2')
        ad <- ad[,c('Var1','Var2','group1','group2','tissue_type1','tissue_type2','value'),with=F] 

        ## get the min distances between PMs and DMs
        pm_td <- ad[group1=='Peritoneum' & group2=='Locoregional']
        pm_td <- pm_td[order(Var1, value, decreasing=F),]
        pm_td <- pm_td[!duplicated(Var1),]
        pm_pt <- ad[group1=='Peritoneum' & group2=='Primary']
        pm_pt <- pm_pt[order(Var1, value, decreasing=F),]
        pm_pt <- pm_pt[!duplicated(Var1),]

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        out <- merge(pm_td[,c('Var1','Var2','tissue_type2','value'),with=F], pm_pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T)
        names(out) <- c('sample','td_sample','td_type','td_distance','pt_sample','pt_distance')
        to_add <- data.table(sample='average', td_sample='average', td_type='average', td_distance=mean(out$td_distance), pt_sample='average', pt_distance=mean(out$pt_distance))
        out <- rbind(out, to_add)
        out$i <- i
        out$patient <- patient
        out
    }

    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(td_distance / pt_distance)]
    
    collapse_pm <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_pm(.SD),by=c('patient','sample')] 
    res$n_td <- n_td
    res$n_pt <- n_pt
    res$n_pm <- n_pm
    res
}


## test common/distinct origin of synchronous PMs and synchronous LNs
si <- copy(sample_info)
si <- si[
         (group %in% c('Normal','Primary') | 
         (tissue_type %in% c('Tumor deposit') & met_timing=='synchronous') | 
         (group %in% 'Peritoneum')) & Patient_ID %in% all_sync_patients]

## get patients with any peritoneal mets and with any synchronous TDs
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
td_patients <- unique(si[tissue_type=='Tumor deposit','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, td_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]

per_patients_sync <- unique(si[met_timing=='synchronous' & group=='Peritoneum','Patient_ID',with=F][[1]])
length(intersect(per_patients_sync, td_patients))
per_patients_meta <- unique(si[met_timing!='synchronous' & group=='Peritoneum','Patient_ID',with=F][[1]])
length(intersect(per_patients_meta, td_patients))


cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Common (95%)','Common (80%)','n.s.','Distinct (80%)','Distinct (95%)')
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res1 <- rbindlist(l)

pd1 <- res1[sample!='average']
pd1 <- pd1[order(obs,decreasing=F),]
pd1[,n_td_div_n_pm:=n_td / n_pm]
pd1[,id:=paste(patient, sample)]
pd1$id <- factor(pd1$id, levels=pd1$id)
pd1$origin <- 'n.s.'
pd1[`97.5%` < 0, origin:='Common (95%)']
pd1[`90%` < 0, origin:='Common (80%)']
pd1[`10%` > 0, origin:='Distinct (80%)']
pd1[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(pd1, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Sync-PM, Sync-TD\n(origin confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    labs(x=NULL, y='log2 min-distance-ratio\n(PM:TD vs PM:PT)')

## test common/distinct origin of metachronous PMs and synchronous LNs
si <- copy(sample_info)
si <- si[
         (group %in% c('Normal','Primary') | 
         (tissue_type %in% c('Tumor deposit') & met_timing=='synchronous') | 
         (group %in% 'Peritoneum')) & Patient_ID %in% all_meta_patients]

## get patients with any peritoneal mets and with any synchronous LN mets
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
td_patients <- unique(si[tissue_type=='Tumor deposit','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, td_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Common (95%)','Common (80%)','n.s.','Distinct (80%)','Distinct (95%)')
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res2 <- rbindlist(l)

pd2 <- res2[sample!='average']
pd2 <- pd2[order(obs,decreasing=F),]
pd2[,n_td_div_n_pm:=n_td / n_pm]
pd2[,id:=paste(patient, sample)]
pd2$id <- factor(pd2$id, levels=pd2$id)
pd2$origin <- 'n.s.'
pd2[`97.5%` < 0, origin:='Common (95%)']
pd2[`90%` < 0, origin:='Common (80%)']
pd2[`10%` > 0, origin:='Distinct (80%)']
pd2[`2.5%` > 0, origin:='Distinct (95%)']

p2 <- ggplot(pd2, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Meta-PM, Sync-TD\n(origin confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    labs(x=NULL, y='log2 min-distance-ratio\n(PM:TD vs PM:PT)')

p <- plot_grid(p1, p2, ncol=1)
ggsave(here('figures/unused_fig_origins_pm_timing_vs_sync_td.pdf'), width=8, height=8)

tbl1 <- table_freq(pd1$origin)
tbl2 <- table_freq(pd2$origin)
tbl <- merge(tbl1, tbl2, by='value', all=T)
tbl[is.na(tbl)] <- 0
names(tbl) <- c('origin','Sync PM','Meta PM')
toadd <- data.table("Common (95%)",0,0)
names(toadd) <- names(tbl)
tbl <- rbind(tbl, toadd)
tbl$origin <- factor(tbl$origin, levels=c('Common (95%)','Common (80%)','n.s.','Distinct (80%)','Distinct (95%)'))
tbl <- tbl[order(origin),]
tmp <- d2m(tbl)
tmp <- tmp[-1,]
p_ca <- CochranArmitageTest(tmp)$p.value
p_chi <- chisq.test(tmp)$p.value
st <- paste0('Chi-Squared P=',round(p_chi,4),'; Cochran-Armitage P=',round(p_ca,4))
prop <- copy(tbl)
prop[,2] <- prop[,2] / sum(prop[,2])
prop[,3] <- prop[,3] / sum(prop[,3])

dat <- data.table::melt(tbl, id.var='origin')
prop <- data.table::melt(prop, id.var='origin')
dat$prop <- prop$value
dat$origin <- factor(dat$origin, levels=rev(levels(dat$origin)))
get_label_pos <- function(dat) {
    dat <- dat[order(origin,decreasing=T),]
    dat$pos <- (cumsum(dat$prop) - 0.5*dat$prop)
    dat
}
dat2 <- dat[,get_label_pos(.SD), by=c('variable')]
p <- ggplot(dat2, aes(x=variable, y=prop)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    geom_bar(stat='identity', aes(fill=origin),color='black',linewidth=0.25) +
    geom_text(data=dat2[value > 0],aes(label=value,y=pos)) +
    labs(x='PM timing',y=paste('Proportion of',nrow(pd1) + nrow(pd2),'PM lesions'),subtitle=st) +
    scale_fill_manual(values=cols, name='Origin with sync. tumor deposits') + 
    theme_ang(base_size=12) +
    theme(legend.position='bottom') +
    coord_flip()
ggsave(here('figures/unused_fig_origins_pm_timing_vs_td_timing_summary.pdf'), width=7, height=2.25)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test for a shared origin between sync/meta PMs and synchronous lymph node mets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89)
all_sync_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='synchronous',(Patient_ID)])
all_meta_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='metachronous',(Patient_ID)])

## check for shared origin between PMs and Synchronous LNs between Sync and Meta PM patients
test_patient_met_level <- function(patient, si, ad_table, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary','Peritoneum', 'Locoregional')]

    if(patient %in% c('Ea3','Eb3','Ec3')) {    
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/E3.rds'))
        bs <- readRDS(bs_file)
    } else {
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/',patient,'.rds'))
        bs <- readRDS(bs_file)
    }

    ## get the original set of valid-samples for the true number of distant/pt/pm samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    #n_ln <- length(unique(si$Real_Sample_ID[si$group=='Locoregional']))
    n_ln <- length(unique(si$Real_Sample_ID[si$tissue_type=='Lymph node']))
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))
    n_pm <- length(unique(si$Real_Sample_ID[si$group=='Peritoneum']))

    test_patient_bs <- function(i, patient, si, ad_table, bs) {
        if(i > 0) {
            ## for the bootstrap replicates, obtain new ad and valid_samples
            ad <- bs[[i]]
            valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
        }
        ad <- ad[valid_samples, valid_samples]
        ad <- as.data.table(reshape2::melt(ad))    
        groups <- si[,c('Real_Sample_ID','group','tissue_type'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        setnames(ad,'tissue_type','tissue_type1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        setnames(ad,'tissue_type','tissue_type2')
        ad <- ad[,c('Var1','Var2','group1','group2','tissue_type1','tissue_type2','value'),with=F] 

        ## get the min distances between PMs and DMs
        pm_ln <- ad[group1=='Peritoneum' & group2=='Locoregional']
        pm_ln <- pm_ln[order(Var1, value, decreasing=F),]
        pm_ln <- pm_ln[!duplicated(Var1),]
        pm_pt <- ad[group1=='Peritoneum' & group2=='Primary']
        pm_pt <- pm_pt[order(Var1, value, decreasing=F),]
        pm_pt <- pm_pt[!duplicated(Var1),]

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        out <- merge(pm_ln[,c('Var1','Var2','tissue_type2','value'),with=F], pm_pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T)
        names(out) <- c('sample','ln_sample','ln_type','ln_distance','pt_sample','pt_distance')
        to_add <- data.table(sample='average', ln_sample='average', ln_type='average', ln_distance=mean(out$ln_distance), pt_sample='average', pt_distance=mean(out$pt_distance))
        out <- rbind(out, to_add)
        out$i <- i
        out$patient <- patient
        out
    }

    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(ln_distance / pt_distance)]
    
    collapse_pm <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_pm(.SD),by=c('patient','sample')] 
    res$n_ln <- n_ln
    res$n_pt <- n_pt
    res$n_pm <- n_pm
    res
}


## test common/distinct origin of synchronous PMs and synchronous LNs
si <- copy(sample_info)
si <- si[
         (group %in% c('Normal','Primary') | 
         (tissue_type %in% c('Lymph node') & met_timing=='synchronous') | 
         (group %in% 'Peritoneum')) & Patient_ID %in% all_sync_patients]

## get patients with any peritoneal mets and with any synchronous LNs
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
ln_patients <- unique(si[tissue_type=='Lymph node','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, ln_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Common (95%)','Common (80%)','n.s.','Distinct (80%)','Distinct (95%)')
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res1 <- rbindlist(l)

pd1 <- res1[sample!='average']
pd1 <- pd1[order(obs,decreasing=F),]
pd1[,n_ln_div_n_pm:=n_ln / n_pm]
pd1[,id:=paste(patient, sample)]
pd1$id <- factor(pd1$id, levels=pd1$id)
pd1$origin <- 'n.s.'
pd1[`97.5%` < 0, origin:='Common (95%)']
pd1[`90%` < 0, origin:='Common (80%)']
pd1[`10%` > 0, origin:='Distinct (80%)']
pd1[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(pd1, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Sync-PM, Sync-LN\n(origin confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    labs(x=NULL, y='log2 min-distance-ratio\n(PM:LN vs PM:PT)')

## test common/distinct origin of metachronous PMs and synchronous LNs
si <- copy(sample_info)
si <- si[
         (group %in% c('Normal','Primary') | 
         (tissue_type %in% c('Lymph node') & met_timing=='synchronous') | 
         (group %in% 'Peritoneum')) & Patient_ID %in% all_meta_patients]

## get patients with any peritoneal mets and with any synchronous LN mets
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
ln_patients <- unique(si[tissue_type=='Lymph node','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, ln_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Common (95%)','Common (80%)','n.s.','Distinct (80%)','Distinct (95%)')
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res2 <- rbindlist(l)

pd2 <- res2[sample!='average']
pd2 <- pd2[order(obs,decreasing=F),]
pd2[,n_ln_div_n_pm:=n_ln / n_pm]
pd2[,id:=paste(patient, sample)]
pd2$id <- factor(pd2$id, levels=pd2$id)
pd2$origin <- 'n.s.'
pd2[`97.5%` < 0, origin:='Common (95%)']
pd2[`90%` < 0, origin:='Common (80%)']
pd2[`10%` > 0, origin:='Distinct (80%)']
pd2[`2.5%` > 0, origin:='Distinct (95%)']

p2 <- ggplot(pd2, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Meta-PM, Sync-LN\n(origin confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    labs(x=NULL, y='log2 min-distance-ratio\n(PM:LN vs PM:PT)')

p <- plot_grid(p1, p2, ncol=1)
ggsave(here('figures/unused_fig_origins_pm_timing_vs_sync_ln.pdf'), width=8, height=8)

tbl1 <- table_freq(pd1$origin)
tbl2 <- table_freq(pd2$origin)
tbl <- merge(tbl1, tbl2, by='value', all=T)
tbl[is.na(tbl)] <- 0
names(tbl) <- c('origin','Sync PM','Meta PM')
toadd <- data.table("Common (95%)",0,0)
names(toadd) <- names(tbl)
tbl <- rbind(tbl, toadd)
tbl$origin <- factor(tbl$origin, levels=c('Common (95%)','Common (80%)','n.s.','Distinct (80%)','Distinct (95%)'))
tbl <- tbl[order(origin),]
tmp <- d2m(tbl)
tmp <- tmp[-1,]
p_ca <- CochranArmitageTest(tmp)$p.value
p_chi <- chisq.test(tmp)$p.value
st <- paste0('Chi-Squared P=',round(p_chi,4),'; Cochran-Armitage P=',round(p_ca,4))
prop <- copy(tbl)
prop[,2] <- prop[,2] / sum(prop[,2])
prop[,3] <- prop[,3] / sum(prop[,3])

dat <- data.table::melt(tbl, id.var='origin')
prop <- data.table::melt(prop, id.var='origin')
dat$prop <- prop$value
dat$origin <- factor(dat$origin, levels=rev(levels(dat$origin)))
get_label_pos <- function(dat) {
    dat <- dat[order(origin,decreasing=T),]
    dat$pos <- (cumsum(dat$prop) - 0.5*dat$prop)
    dat
}
dat2 <- dat[,get_label_pos(.SD), by=c('variable')]
p <- ggplot(dat2, aes(x=variable, y=prop)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    geom_bar(stat='identity', aes(fill=origin),color='black',linewidth=0.25) +
    geom_text(data=dat2[value > 0],aes(label=value,y=pos)) +
    labs(x='PM timing',y=paste('Proportion of',nrow(pd1) + nrow(pd2),'PM lesions'),subtitle=st) +
    scale_fill_manual(values=cols, name='Origin with sync. lymph node mets') + 
    theme_ang(base_size=12) +
    theme(legend.position='bottom') +
    coord_flip()
ggsave(here('figures/unused_fig_origins_pm_timing_vs_ln_timing_summary.pdf'), width=7, height=2.25)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sync/meta-only PM vs any liver inter-lesion pairwise AD
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

collapse <- function(res) {
    dist=mean(res$distance,na.rm=T)
    list(distance=dist)
}

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## get pairwise inter-lesion AD
si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='synchronous')]
si <- si[in_collapsed==T,]
res_per_sync <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_per_sync <- res_per_sync[group1=='Metastasis' & group2=='Metastasis']
res_per_sync[res_per_sync=='Metastasis'] <- 'Peritoneum'
res_per_sync$class <- 'PM (sync only)'

## get pairwise inter-lesion AD
si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='metachronous')]
si <- si[in_collapsed==T,]
res_per_meta <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_per_meta <- res_per_meta[group1=='Metastasis' & group2=='Metastasis']
res_per_meta[res_per_meta=='Metastasis'] <- 'Peritoneum'
res_per_meta$class <- 'PM (meta only)'

## get for any liver mets
si <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
si <- si[in_collapsed==T,]
res_liv <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_liv <- res_liv[group1=='Metastasis' & group2=='Metastasis']
res_liv[res_liv=='Metastasis'] <- 'Liver'
res_liv$class <- 'Liver (any timing)'

res <- rbind(res_per_sync, res_per_meta, res_liv)
res$class <- factor(res$class, levels=unique(res$class))
res[grep('PM',class),group:='Peritoneum']
res[grep('Liver',class),group:='Liver']

p_inter <- ggplot(res, aes(x=class, y=distance)) + 
    #scale_y_continuous(breaks=seq(0,1.5,by=0.5),limits=c(0,1.5),expand=c(0,0)) +
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method="wilcox.test", tip.length=NA, 
                       comparisons=list(
                                        c("PM (sync only)","PM (meta only)"), 
                                        c("PM (sync only)","Liver (any timing)"),
                                        c("PM (meta only)","Liver (any timing)"))) + 
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Inter-lesion angular distance',title='PM (sync x meta) vs Liver inter-lesion angular distance',subtitle='(un-corrected wilcox tests)')
ggsave(here('figures/ed_fig_5a_timing_per_vs_liver_inter_angular_distance.pdf'), width=7, height=6)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sync/meta PM intra-lesion AD vs Liver (any timing)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subset_for_intralesion <- function(res) {
    ## annotate each sample with the type+lesion, so that we can find comparisons of different samples within the same type+lesion
    bc1 <- parse_barcode(res$sample1)
    bc1[,typelesion:=paste0(type,lesion)]
    bc2 <- parse_barcode(res$sample2)
    bc2[,typelesion:=paste0(type,lesion)]
    res$typelesion1 <- bc1$typelesion
    res$typelesion2 <- bc2$typelesion
    res <- res[typelesion1==typelesion2]
    res
}

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])
liv_patients <- unique(sample_info$Patient_ID[sample_info$group=='Liver'])

## subset the sample_info for N, PT, and Per; get met-spec angular distance from peritoneum to normal
si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='synchronous')]
res_per_sync <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_per_sync <- res_per_sync[group1=='Metastasis' & group2=='Metastasis']
res_per_sync[res_per_sync=='Metastasis'] <- 'Peritoneum'
res_per_sync <- subset_for_intralesion(res_per_sync)
res_per_sync$class <- 'PM (sync only)'

si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='metachronous')]
res_per_meta <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_per_meta <- res_per_meta[group1=='Metastasis' & group2=='Metastasis']
res_per_meta[res_per_meta=='Metastasis'] <- 'Peritoneum'
res_per_meta <- subset_for_intralesion(res_per_meta)
res_per_meta$class <- 'PM (meta only)'

si <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
res_liv <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_liv <- res_liv[group1=='Metastasis' & group2=='Metastasis']
res_liv[res_liv=='Metastasis'] <- 'Liver'
res_liv <- subset_for_intralesion(res_liv)
res_liv$class <- 'Liver (any timing)'

res <- rbind(res_per_sync, res_per_meta, res_liv)
res$class <- factor(res$class, levels=unique(res$class))

res[grep('PM',class),group:='Peritoneum']
res[grep('Liver',class),group:='Liver']

p_intra <- ggplot(res, aes(x=class, y=distance)) + 
    #scale_y_continuous(breaks=seq(0,1.5,by=0.5),limits=c(0,1.5),expand=c(0,0)) +
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method="wilcox.test", tip.length=NA, 
                       comparisons=list(
                                        c("PM (sync only)","PM (meta only)"), 
                                        c("PM (sync only)","Liver (any timing)"),
                                        c("PM (meta only)","Liver (any timing)"))) + 
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Intra-lesion angular distance',title='PM (sync x meta) vs Liver intra-lesion angular distance',subtitle='(un-corrected wilcox tests)')
ggsave(here('figures/ed_fig_5b_timing_per_vs_liver_intra_angular_distance.pdf'), width=7, height=6)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ED Fig 2B. i, iii: intralesion pairwise node distance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si <- sample_info[Patient_ID %in% per_patients & group %in% c('Normal','Primary','Peritoneum')]
res_per <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='node', return_tree=F)
res_pt <- res_per[group1=='Primary' & group2=='Primary']
res_per <- res_per[group1=='Metastasis' & group2=='Metastasis']
res_per[res_per=='Metastasis'] <- 'Peritoneum'
res_per <- subset_for_intralesion(res_per)
res_per_with_pt <- rbind(res_pt, res_per, fill=T)
patients_with_geq2_per_mets <- names(which(rowSums(xtabs(~ patient + group1, data=res_per_with_pt) > 0) > 1))
res_per_with_pt <- res_per_with_pt[patient %in% patients_with_geq2_per_mets]

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si <- sample_info[Patient_ID %in% liv_patients & group %in% c('Normal','Primary','Liver')]
res_liv <- get_met_specific_distances(si, ad_table, comparison='pairwise', distance='node', return_tree=F)
res_pt <- res_liv[group1=='Primary' & group2=='Primary']
res_liv <- res_liv[group1=='Metastasis' & group2=='Metastasis']
res_liv[res_liv=='Metastasis'] <- 'Liver'
res_liv <- subset_for_intralesion(res_liv)
res_liv_with_pt <- rbind(res_pt, res_liv, fill=T)
patients_with_geq2_liv_mets <- names(which(rowSums(xtabs(~ patient + group1, data=res_liv_with_pt) > 0) > 1))
res_liv_with_pt <- res_liv_with_pt[patient %in% patients_with_geq2_liv_mets]

res <- rbind(res_per_with_pt, res_liv)
res$group1 <- factor(res$group1, levels=unique(res$group1))

p <- ggplot(res[group1 %in% c('Peritoneum','Liver')], aes(x=group1, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_compare_means(method="wilcox.test", comparisons=list(c("Peritoneum","Liver"))) +
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    labs(x=NULL,y='Node distance / overall average node distance',title='Intralesion met-specific node distance (all met timing)',subtitle='All patients with PMs or liver mets')
ggsave(here('figures/ed_fig_2b_tissue_type_intralesion_node_distance.pdf'))


## pairwise plot
collapse <- function(res) {
    dist=mean(res$distance,na.rm=T)
    list(dist=dist)
}
per_paired <- res_per_with_pt[,collapse(.SD),by=c('patient','group1')]
per_paired <- data.table::dcast(patient ~ group1, value.var='dist', data=per_paired)
per_paired <- per_paired[!is.na(Peritoneum)]
liv_paired <- res_liv_with_pt[,collapse(.SD),by=c('patient','group1')]
liv_paired <- data.table::dcast(patient ~ group1, value.var='dist', data=liv_paired)
liv_paired <- liv_paired[!is.na(Liver)]

per_pairedM <- data.table::melt(per_paired, id.var='patient')
per_pairedM$comparison <- 'Peritoneum vs Primary'
per_pairedM$variable <- factor(per_pairedM$variable, levels=c('Primary','Peritoneum'))
per_pairedM <- per_pairedM[order(variable),]
liv_pairedM <- data.table::melt(liv_paired, id.var='patient')
liv_pairedM$comparison <- 'Liver vs Primary'
liv_pairedM$variable <- factor(liv_pairedM$variable, levels=c('Primary','Liver'))
liv_pairedM <- liv_pairedM[order(variable),]
pairedM <- rbind(per_pairedM, liv_pairedM)
pairedM[,class:=paste0(variable,' (',comparison,')')]
pairedM$class <- factor(pairedM$class, levels=unique(pairedM$class))
pairedM[,patient_comparison:=paste(patient,comparison)]

p1 <- ggplot(pairedM, aes(x=class, y=value)) +
    scale_y_continuous(limits=c(0,2),breaks=seq(0,1.5,by=0.5)) +
    geom_line(aes(group=patient_comparison)) +
    geom_point(pch=21,aes(fill=variable,group=patient), size=4) + 
    stat_compare_means(method="wilcox.test", 
                       comparisons=list(c('Peritoneum (Peritoneum vs Primary)','Primary (Peritoneum vs Primary)'),
                                        c('Liver (Liver vs Primary)','Primary (Liver vs Primary)')), label.y=1.6, paired=T) +
    stat_compare_means(method="wilcox.test", 
                       comparisons=list(c('Peritoneum (Peritoneum vs Primary)','Liver (Liver vs Primary)')),label.y=1.8, paired=F) +
    scale_fill_manual(values=group_cols) +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    guides(fill='none') +
    labs(x=NULL,y='Node distance (patient average)',title='Intralesion node distance vs tissue type (all met timing)',subtitle='All patients with PMs or liver mets')
ggsave(here('figures/ed_fig_3b_tissue_type_intralesion_node_distance_average_wrt_primary.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Maybe add to EDF5: intra-lesion AD subset by timing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si_per_syn <- sample_info[Patient_ID %in% per_patients & 
                  (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='synchronous'))]
res_per_syn <- get_met_specific_distances(si_per_syn, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_per_syn <- res_per_syn[group1=='Metastasis' & group2=='Metastasis']
res_per_syn[res_per_syn=='Metastasis'] <- 'Peritoneum'
res_per_syn$timing <- 'synchronous'
res_per_syn <- subset_for_intralesion(res_per_syn)

si_per_met <- sample_info[Patient_ID %in% per_patients & 
                  (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_timing=='metachronous'))]
res_per_met <- get_met_specific_distances(si_per_met, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_per_met <- res_per_met[group1=='Metastasis' & group2=='Metastasis']
res_per_met[res_per_met=='Metastasis'] <- 'Peritoneum'
res_per_met$timing <- 'metachronous'
res_per_met <- subset_for_intralesion(res_per_met)

res <- rbind(res_per_syn, res_per_met)
res$timing <- factor(res$timing, levels=unique(res$timing))

collapse <- function(res) {
    d <- mean(res$distance)
    list(d=d)
}
res2 <- res[group1=='Peritoneum',collapse(.SD),by=c('patient','timing','group1')]
stat.test <- mywilcox2(res2, d ~ timing, paired=F)

p <- ggplot(res2, aes(x=timing, y=d)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=group1)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    scale_fill_manual(values=group_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x=NULL,y='Angular distance',title='Intralesion angular distance by timing (patient-average)')
ggsave(here('figures/unused_fig_timing_per_intralesion_angular_distance_average.pdf'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No longer used: min. angular distance to PT (luminal vs deep)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

per_patients <- sample_info[group=='Peritoneum',(Patient_ID)] %>% unique

## angular distance from peritoneum to PT subset by PT depth
si <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Peritoneum')) & in_collapsed==T]
res <- get_met_specific_distances(si, ad_table, comparison='primary', distance='angular', return_tree=F)
res <- merge(res, si[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res <- merge(res, si[,c('Patient_ID','Real_Sample_ID','met_timing'),with=F], by.x=c('patient','sample1'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res[group1=='Metastasis',group1:='Peritoneum']
res$vertical <- factor(res$vertical, levels=c('mucosal/luminal','deep'))
res <- res[!is.na(vertical),]
summarize <- function(res) {
    min <- min(res$distance)
    list(d=min)
}
info <- res[,summarize(.SD),by=c('patient','vertical','met_timing')]
info <- data.table::dcast(patient + met_timing ~ vertical, value.var='d', data=info)
info <- info[!is.na(`mucosal/luminal`) & !is.na(deep),]
info <- data.table::melt(info, id.vars=c('patient','met_timing'))
setnames(info,c('variable','value'),c('vertical','d'))
info$group <- 'Peritoneum'
stat.test <- mywilcox2(info[met_timing %in% c('synchronous','metachronous')], d ~ vertical, paired=T, facet_field='met_timing')

p1 <- ggplot(info[met_timing %in% c('synchronous','metachronous')], aes(x=vertical, y=d, group=patient)) +
    scale_y_continuous(limits=c(0,0.8), breaks=seq(0,0.8,by=0.2)) + 
    geom_line() +
    geom_point(aes(fill=group), pch=21, size=4, color='black') +
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    scale_fill_manual(values=group_cols,name='Tissue') +
    facet_wrap(facets=~met_timing, nrow=1) +
    guides(fill='none') +  
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x='Primary tumor invasion depth',y='Min. angular distance',subtitle='Min. angular distance from PMs to deep and luminal/mucosal PT regions') 
ggsave(here('figures/unused_fig_depth_vs_timing_per_min_angular_distance_paired.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# patient-level permutation-based test of proximity between Mets and deep/luminal PT
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

#res[cohort=='peritoneal' & group=='Locoregional' & grepl('TD',sample1), group:='Tumor deposit']
#res[cohort=='science' & group=='Locoregional' & patient!='C12', group:='Lymph node']
#res[cohort=='science' & group=='Locoregional' & patient=='C12', group:='Tumor deposit']
#res[group=='Locoregional', group:='Locoregional, NOS']
res <- res[!is.na(plotgroup),]
res[,id:=paste0(patient,', ',plotgroup)]
res[,group:=plotgroup]

qc <- res[patient=='E15' & group=='Peritoneum',]
l <- run_test(qc, R=10000, ncpus=8)
tmp <- data.table(expected=l$exp)
pval <- 2*(sum(tmp$expected >= l$obs) + 1) / (nrow(tmp) + 1) 
plab <- paste0('P=',round(pval, 3))
p <- ggplot(tmp, aes(x=expected)) +
    geom_histogram(color='white',fill='#bfbfbf',bins=50) +
    geom_vline(xintercept=l$obs,color='red') +
    geom_text(data=tmp[1,], x=1, y=580, label=plab,color='red') +
    theme_ang(base_size=12) +
    labs(y='N', x='Average distance ratio from each met\nto nearest deep and luminal PTs',subtitle='Patient E15')
ggsave(here('figures/tmp_fig_E15_peritoneum_pt_depth_volcano.pdf'))


res_split <- split(res, by='id')
require(parallel)
RNGkind("L'Ecuyer-CMRG") 
set.seed(42)

## this takes a long time to run
required_file <- here('processed_data/misc/fig_3b_depth_angular_distance_by_tissue_type_all_liver_mindistance_volcanoes_averageratio_twosided_validpatients.txt')
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
cols <- c('#bfbfbf','red')
names(cols) <- c('n.s.','p-adj < 0.01')

library(ggrepel)
p <- ggplot(results, aes(x=lfc, y=nlog10q)) +
    scale_y_continuous(limits=c(0,3.0)) +
    scale_x_continuous(limits=c(-1.5,1.5)) +
    geom_vline(xintercept=0, color='black', linewidth=0.5) +
    geom_point(data=results[pval >= 0.01], pch=19, size=2, color='#bfbfbf') + 
    geom_point(data=results[qval >= 0.01 & pval < 0.01], pch=19, size=2, aes(color=significance)) + 
    geom_point(data=results[qval < 0.01], pch=19, size=3, aes(color=significance)) + 
    geom_text_repel(data=results[pval < 0.05], aes(label=patient, color=significance)) + 
    scale_color_manual(values=cols, name='Significance') + 
    facet_wrap(facets=~group, scale='free', ncol=4) +
    theme_bw(base_size=12) +
    labs(x='log2(observed / median expected)', y='-log10(Q-value)',subtitle='Distance ratio from each met to nearest deep and luminal PTs; average across mets per patient')
ggsave(here('figures/fig_3b_depth_angular_distance_by_tissue_type_all_liver_mindistance_volcanoes_averageratio_twosided.pdf'),width=11,height=3)


## are the cases with deep-invading evidence enriched for synchronoud and/or untreated PMs?
x <- results[group=='Peritoneum',]
x[lfc < 0, class:='deep-invading']
x[lfc > 0, class:='luminal']

sync_pm_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='synchronous',(Patient_ID)])
meta_pm_patients <- unique(sample_info[group=='Peritoneum' & met_timing=='metachronous',(Patient_ID)])
x[patient %in% sync_pm_patients & !patient %in% meta_pm_patients, met_timing:='sync']
x[!patient %in% sync_pm_patients & patient %in% meta_pm_patients, met_timing:='meta']
x[patient %in% sync_pm_patients & patient %in% meta_pm_patients, met_timing:='sync+meta']

untreated_pm_patients <- unique(sample_info[group=='Peritoneum' & met_treated=='untreated',(Patient_ID)])
treated_pm_patients <- unique(sample_info[group=='Peritoneum' & met_treated=='treated',(Patient_ID)])
x[patient %in% untreated_pm_patients & !patient %in% treated_pm_patients, met_treated:='untreated']
x[!patient %in% untreated_pm_patients & patient %in% treated_pm_patients, met_treated:='treated']
x[patient %in% untreated_pm_patients & patient %in% treated_pm_patients, met_treated:='untreated+treated']
x$untreated <- grepl('untreated', x$met_treated)
#fisher.test(xtabs(~ class + untreated, data=x))
xtabs(~ class + untreated, data=x)

untreated_pm_patients <- unique(sample_info[group=='Peritoneum' & met_treated=='untreated',(Patient_ID)])

#x[patient %in% untreated_pm_patients & !patient %in% sync_pm_patients, class:='untreated']
#x[patient %in% untreated_pm_patients & patient %in% sync_pm_patients, class:='sync AND untreated']
#x[patient %in% untreated_pm_patients | patient %in% sync_pm_patients, class:='sync OR untreated']



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# updated: 2024-05-17
# angular distance to PT from LN, TD, Liv, and Per mets, subset by PT depth
# NOW INCLUDE ALL PATIENTS WITH LIVER METS AND DEPTH INFO
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
#si_liv <- sample_info[Patient_ID %in% liv_patients & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
#res_liv <- get_met_specific_distances(si_liv, ad_table, comparison='primary', distance='angular', return_tree=F)
#res_liv$group1 <- 'Liver (all liver patients)'
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


## 2024-06-17
## repeat this with two specified buckets of liver mets:
## 1.liver Mets came out before peri Mets or synchronously
## 2.liver Mets came out out at least 3 months after peri Mets

#pm_cohort_per <- unique(sample_info[Patient_ID %in% per_patients & group=='Peritoneum', (Patient_ID)])
#pm_cohort_liv <- unique(sample_info[Patient_ID %in% per_patients & group=='Liver', (Patient_ID)])
#pm_cohort_both <- sort(intersect(pm_cohort_per, pm_cohort_liv))
#info <- sample_info[Patient_ID %in% pm_cohort_both & group %in% c('Peritoneum','Liver'),c('Patient_ID','group','Real_Sample_ID','in_collapsed','treatment','met_timing','met_treated','met_treated_type','resection_procedure'),with=F]
#write_tsv(info, here('processed_data/pm_cohort_per_and_liver_case_samples.txt')) # manually added timing_order field, relative timing (sep by 3 months)
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
    list(n_pms_prior_to_liv=n_pms_prior_to_liv, n_pms_same_time_as_liv=n_pms_same_time_as_liv, n_pms_after_liv=n_pms_after_liv)
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

## potentially new version of Fig 3b
p1 <- ggplot(pd1, aes(x=vertical, y=distance)) +
    scale_y_continuous(limits=c(0,2.25),breaks=seq(0,2.25,by=0.5)) +
    geom_point(aes(color=vertical),position=position_jitter(width=0.1,height=0,seed=2), pch=16, size=4) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA,width=0.5) +
    scale_color_manual(values=cols_depth2) +
    guides(fill='none', color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    facet_wrap(facets=~patient_type, nrow=1, ncol=3) +
    stat_pvalue_manual(stat.test1, label = "label", tip.length = 0.02) +
    labs(x='Invasion depth of primary tumor', y='Angular distance',title='Angular distance from primary tumor to mets')

p2 <- ggplot(pd2, aes(x=vertical, y=distance)) +
    scale_y_continuous(limits=c(0,2.25),breaks=seq(0,2.25,by=0.5)) +
    geom_point(aes(color=vertical),position=position_jitter(width=0.1,height=0,seed=2), pch=16, size=4) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA,width=0.5) +
    scale_color_manual(values=cols_depth2) +
    guides(fill='none', color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    facet_wrap(facets=~patient_type, nrow=1, ncol=3) +
    stat_pvalue_manual(stat.test2, label = "label", tip.length = 0.02) +
    labs(x='Invasion depth of primary tumor', y='Angular distance',title='Angular distance from primary tumor to mets')

p <- plot_grid(p1, p2, nrow=1)
ggsave(here('figures/unused_fig_depth_angular_distance_by_tissue_type_with_liver_buckets.pdf'), width=12, height=4)


## LR among NG/science cohorts
lr_patients <- unique(sample_info[cohort %in% c('science','natgen') & group=='Locoregional',(Patient_ID)])
lr_patients <- lr_patients[!lr_patients %in% c('C38','C89')]
si_lr <- sample_info[Patient_ID %in% lr_patients & (group %in% c('Normal','Primary','Locoregional')) & in_collapsed==T ]
res_lr <- get_met_specific_distances(si_lr, ad_table, comparison='primary', distance='angular', return_tree=F)
res_lr$group1 <- 'Locoregional'
res_lr <- merge(res_lr, sample_info[,c('Patient_ID','Real_Sample_ID','vertical','cohort'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res_lr$vertical <- factor(res_lr$vertical, levels=c('mucosal/luminal','deep'))
res_lr <- res_lr[!is.na(vertical)]
res_lr$class <- 'Locoregional (non-PM patients)'

## PM synchronous+untreated
per_patients <- unique(sample_info[group=='Peritoneum','Patient_ID',with=F][[1]])
si_pm1 <- sample_info[in_collapsed==T & Patient_ID %in% per_patients & (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated=='untreated' & met_timing=='synchronous'))]
res_pm1 <- get_met_specific_distances(si_pm1, ad_table, comparison='primary', distance='angular', return_tree=F)
res_pm1$group1 <- 'Peritoneum'
res_pm1 <- merge(res_pm1, sample_info[,c('Patient_ID','Real_Sample_ID','vertical','cohort'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res_pm1$vertical <- factor(res_pm1$vertical, levels=c('mucosal/luminal','deep'))
res_pm1 <- res_pm1[!is.na(vertical)]
res_pm1$class <- 'Peritoneum (synchronous, untreated)'

per_patients <- unique(sample_info[group=='Peritoneum','Patient_ID',with=F][[1]])
si_pm2 <- sample_info[in_collapsed==T & Patient_ID %in% per_patients & (group %in% c('Normal','Primary') | (group=='Peritoneum' & met_treated_type=='systemic chemo' & met_timing %in% c('metachronous','metachronous after synchronous')))]
res_pm2 <- get_met_specific_distances(si_pm2, ad_table, comparison='primary', distance='angular', return_tree=F)
res_pm2$group1 <- 'Locoregional'
res_pm2 <- merge(res_pm2, sample_info[,c('Patient_ID','Real_Sample_ID','vertical','cohort'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res_pm2$vertical <- factor(res_pm2$vertical, levels=c('mucosal/luminal','deep'))
res_pm2 <- res_pm2[!is.na(vertical)]
res_pm2$class <- 'Peritoneum (adjuvant systemic chemo)'

res <- rbind(res_lr, res_pm1, res_pm2)
res$class <- factor(res$class, levels=c('Locoregional (non-PM patients)','Peritoneum (synchronous, untreated)','Peritoneum (adjuvant systemic chemo)'))
stat.test <- mywilcox2(res, distance ~ vertical, facet_field='class', paired=F, include_n=F)

p3 <- ggplot(res, aes(x=vertical, y=distance)) +
    scale_y_continuous(limits=c(0,2.25),breaks=seq(0,2.25,by=0.5)) +
    geom_point(aes(color=vertical),position=position_jitter(width=0.1,height=0,seed=2), pch=16, size=4) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA,width=0.5) +
    scale_color_manual(values=cols_depth2) +
    guides(fill='none', color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    facet_wrap(facets=~class, nrow=1, ncol=3) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    labs(x='Invasion depth of primary tumor', y='Angular distance', title='SI panels')

ggsave(here('figures/unused_fig_depth_angular_distance_by_tissue_type_with_SI_panels.pdf'), width=9, height=5)




## test bucket 1,2,3 for shared origins with PMs
info <- fread(here('processed_data/pm_cohort_per_and_liver_case_samples.txt'))
info <- info[in_collapsed==T,]
info <- info[order(Patient_ID,timing_order,group),]
info[,id:=paste(Patient_ID, Real_Sample_ID)]
valid_patients <- unique(info$Patient_ID)

test_patient_per_met_level <- function(patient, si, ad_table, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary','Peritoneum','Liver')]

    ## get the original set of valid-samples for the true number of distant/pt/liv samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    n_liv <- length(unique(si$Real_Sample_ID[si$tissue_type=='Liver']))
    n_pm <- length(unique(si$Real_Sample_ID[si$group=='Peritoneum']))
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))

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
        groups <- si[,c('Real_Sample_ID','group','tissue_type'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        setnames(ad,'tissue_type','tissue_type1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        setnames(ad,'tissue_type','tissue_type2')
        ad <- ad[Var1!=Var2,]

        ## get the min distances between PMs and DMs
        pm_liv <- ad[tissue_type1=='Peritoneum' & group2=='Liver']
        pm_liv <- pm_liv[order(Var1, value, decreasing=F),]
        pm_liv <- pm_liv[!duplicated(Var1),]
        pm_pt <- ad[tissue_type1=='Peritoneum' & group2=='Primary']
        pm_pt <- pm_pt[order(Var1, value, decreasing=F),]
        pm_pt <- pm_pt[!duplicated(Var1),]

        ## merge the liv-met level nearest PM and PT, and return the values for this bootstrap replicate
        out <- merge(pm_liv[,c('Var1','Var2','tissue_type2','value'),with=F], pm_pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T, allow.cartesian=T)
        names(out) <- c('sample','liv_sample','liv_type','liv_distance','pt_sample','pt_distance')
        out$i <- i
        out$patient <- patient
        out
    }
    #debugonce(test_patient_bs)
    #l <- test_patient_bs(1, patient, si, ad_table, bs)
    
    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(liv_distance / pt_distance)]
    closest <- res[i==0,c('patient','sample','liv_sample','liv_type'),with=F]
    names(closest) <- c('patient','sample','closest_sample','closest_sample_type')
    collapse_sample <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_sample(.SD),by=c('patient','sample')] 
    res <- merge(res, closest, by=c('patient','sample'), all.x=T)
    res$n_liv <- n_liv
    res$n_pm <- n_pm
    res$n_pt <- n_pt
    res
}

RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_per_met_level, sample_info, ad_table, ncpus)
res <- rbindlist(l)
res[,id:=paste(patient,closest_sample)]
res[id %in% bucket1, bucket:='Liver before PMs']
res[id %in% bucket2, bucket:='Liver same time as PMs']
res[id %in% bucket3, bucket:='Liver after PMs']
res <- res[order(obs),]
write_tsv(res,here('processed_data/min_distance_ratios_pm_to_liv_vs_pm_to_pt_buckets.txt'))


#cols <- c('#cc99a2','#d8d8d8','black')
cols <- c('steelblue','#d8d8d8','salmon')
names(cols) <- c('Liver after PMs','Liver same time as PMs','Liver before PMs')

res <- fread(here('processed_data/min_distance_ratios_pm_to_liv_vs_pm_to_pt_buckets.txt'))
res[,pmid:=paste(patient,sample)]
res <- res[order(obs),]
res$pmid <- factor(res$pmid, levels=res$pmid)
res$origin <- 'n.s.'
res[`90%` < 0, origin:='Shared (80%)']
res[`97.5%` < 0, origin:='Shared (95%)']
res[`10%` > 0, origin:='Distinct (80%)']
res[`2.5%` > 0, origin:='Distinct (95%)']

res$shared <- F
res[grepl('Shared',origin), shared:=T]
res$beforePM <- 'Yes'
res[grepl('Liver before PMs',bucket), beforePM:='No']
res$beforePM <- factor(res$beforePM, levels=c('Yes','No'))
res[bucket=='Liver before PMs', newbucket:='PMs after liver']
res[bucket=='Liver same time as PMs', newbucket:='PMs same time as liver']
res[bucket=='Liver after PMs', newbucket:='PMs before liver']
res$bucket <- factor(res$bucket, levels=c('Liver before PMs','Liver same time as PMs','Liver after PMs'))
stat.test <- mywilcox2(res, obs ~ bucket, paired=F, include_n=F)
cols <- c('steelblue','#d8d8d8','salmon')
names(cols) <- c('Liver after PMs','Liver same time as PMs','Liver before PMs')
#res$newbucket <- factor(res$newbucket, levels=rev(c('PMs after liver','PMs same time as liver','PMs before liver')))
#stat.test <- mywilcox2(res, obs ~ newbucket, paired=F, include_n=F)
#cols <- c('steelblue','#d8d8d8','salmon')
#names(cols) <- c('PMs before liver','PMs same time as liver','PMs after liver')

p3 <- ggplot(res, aes(x=bucket, y=obs)) +
    geom_point(aes(color=bucket),position=position_jitter(width=0.1,height=0,seed=2), pch=16, size=4) + 
    geom_hline(yintercept=0, linetype='dashed', linewidth=0.5) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA,width=0.5) +
    scale_color_manual(values=cols, name='Timing') +
    theme_ang(base_size=12) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = NA) +
    labs(x='Time of liver met. diagnosis', y='Origin ratio')
ggsave(here('figures/unused_fig_distance_ratios_boxplot_with_liv_timing_v2.pdf'))












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# angular distance to PT from LN, TD, Liv, and Per mets, subset by PT depth
# NOW INCLUDE ALL PATIENTS WITH LIVER METS AND DEPTH INFO
# (more verbose option showing samples from specific cohorts)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get patients with any peritoneal mets and with any liver mets
per_patients <- unique(sample_info[group=='Peritoneum','Patient_ID',with=F][[1]])
lr_patients <- unique(sample_info[group=='Locoregional','Patient_ID',with=F][[1]])
liv_patients <- unique(sample_info[group=='Liver','Patient_ID',with=F][[1]])

getN <- function(res) {
    n_per_sample <- function(res) {
        mets <- length(unique(res$sample1))
        pts <- length(unique(res$sample2))
        list(mets=mets, pts=pts)
    }
    cnt <- res[,n_per_sample(.SD),by=patient]
    paste0('(',nrow(cnt),' patients, ',sum(cnt$mets),' Mets, ',sum(cnt$pts),' PTs)')
}

si_per <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Peritoneum')) & in_collapsed==T]
res_per <- get_met_specific_distances(si_per, ad_table, comparison='primary', distance='angular', return_tree=F)
res_per$group1 <- 'Peritoneum'
res_per$cohort <- 'Patients with PMs'

si_lr <- sample_info[Patient_ID %in% lr_patients & cohort=='science' & (group %in% c('Normal','Primary','Locoregional')) & in_collapsed==T ]
res_lr1 <- get_met_specific_distances(si_lr, ad_table, comparison='primary', distance='angular', return_tree=F)
res_lr1$group1 <- 'Locoregional'
res_lr1$cohort <- 'Science cohort'

si_lr <- sample_info[Patient_ID %in% lr_patients & cohort=='natgen' & (group %in% c('Normal','Primary','Locoregional')) & in_collapsed==T ]
res_lr2 <- get_met_specific_distances(si_lr, ad_table, comparison='primary', distance='angular', return_tree=F)
res_lr2$group1 <- 'Locoregional'
res_lr2$cohort <- 'NatGen cohort'

si_lr <- sample_info[Patient_ID %in% lr_patients & cohort=='peritoneal' & (group %in% c('Normal','Primary','Locoregional')) & in_collapsed==T ]
res_lr3 <- get_met_specific_distances(si_lr, ad_table, comparison='primary', distance='angular', return_tree=F)
res_lr3$group1 <- 'Locoregional'
res_lr3$cohort <- 'Peritoneal cohort'

si_ln <- sample_info[Patient_ID %in% lr_patients & cohort=='peritoneal' & ((group %in% c('Normal','Primary') | tissue_type=='Lymph node')) & in_collapsed==T ]
res_ln <- get_met_specific_distances(si_ln, ad_table, comparison='primary', distance='angular', return_tree=F)
res_ln$group1 <- 'Lymph node'
res_ln$cohort <- 'Peritoneal cohort'

si_td <- sample_info[Patient_ID %in% lr_patients & cohort=='peritoneal' & ((group %in% c('Normal','Primary') | tissue_type=='Tumor deposit')) & in_collapsed==T ]
res_td <- get_met_specific_distances(si_td, ad_table, comparison='primary', distance='angular', return_tree=F)
res_td$group1 <- 'Tumor deposit'
res_td$cohort <- 'Peritoneal cohort'

si_liv <- sample_info[Patient_ID %in% liv_patients & cohort=='science' & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
res_liv1 <- get_met_specific_distances(si_liv, ad_table, comparison='primary', distance='angular', return_tree=F)
res_liv1$group1 <- 'Liver'
res_liv1$cohort <- 'Science cohort'

si_liv <- sample_info[Patient_ID %in% liv_patients & cohort=='natgen' & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
res_liv2 <- get_met_specific_distances(si_liv, ad_table, comparison='primary', distance='angular', return_tree=F)
res_liv2$group1 <- 'Liver'
res_liv2$cohort <- 'NatGen cohort'

si_liv <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
res_liv3 <- get_met_specific_distances(si_liv, ad_table, comparison='primary', distance='angular', return_tree=F)
res_liv3$group1 <- 'Liver'
res_liv3$cohort <- 'Patients with PMs'

si_liv <- sample_info[Patient_ID %in% liv_patients & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
res_liv4 <- get_met_specific_distances(si_liv, ad_table, comparison='primary', distance='angular', return_tree=F)
res_liv4$group1 <- 'Liver'
res_liv4$cohort <- 'Patients with liver mets'

res <- rbind(res_lr1, res_lr2, res_lr3, res_liv1, res_liv2, res_liv3,  
             res_ln, res_td, res_per, res_liv4)

res$pos <- 1:nrow(res)
res <- merge(res, sample_info[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res$vertical <- factor(res$vertical, levels=c('mucosal/luminal','deep'))
res <- res[!is.na(vertical)]
res <- res[order(pos),]
res <- res[!(cohort=='Science cohort' & group1=='Locoregional' & patient=='C12')]
check_valid_patient <- function(res) {
    res <- res[!duplicated(sample2),]
    deep <- sum(res$vertical=='deep')
    lum <- sum(res$vertical=='mucosal/luminal')
    list(deep=deep, lum=lum)
}
cnt <- res[,check_valid_patient(.SD),by=patient]
valid <- cnt$patient[cnt$deep > 0 & cnt$lum > 0]
res <- res[patient %in% valid,]

addlabel <- function(res) {
    res$label <- getN(res)
    res
}
res2 <- res[,addlabel(.SD),by=c('cohort','group1')]
res2$cohort2 <- paste0(res$cohort,'\n',res2$label)
res2[,category:=paste0(group1, ', ',cohort2)]
res2[,pos:=NULL]
res2[,cohort2:=NULL]
res2$category <- factor(res2$category, levels=unique(res2$category))
res2[group1 %in% c('Lymph node','Tumor deposit'), group1:='Locoregional']
stat.test <- mywilcox2(res2, distance ~ vertical, facet_field='category', paired=F)

## potentially new version of Fig 3b
p <- ggplot(res2, aes(x=vertical, y=distance)) +
    scale_y_continuous(limits=c(0,2.2),breaks=seq(0,2.2,by=0.5)) +
    geom_point(aes(fill=group1),position=position_jitter(width=0.15,height=0,seed=2), pch=21, size=4, color='black') + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) +
    scale_fill_manual(values=group_cols) + 
    guides(fill='none', color='none') +
    theme_ang(base_size=10) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    facet_wrap(facets=~category, ncol=6) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, size=2.5) +
    labs(x='Invasion depth of primary tumor', y='Angular distance',title='Angular distance from primary tumor to mets')
ggsave(here('figures/unused_fig_depth_angular_distance_by_tissue_type_all_liver_separate_cohorts_patients_with_both_depths.pdf'),width=11,height=8)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dev: permutation test for enriched comingling between deep/luminal and given met type
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get patients with any peritoneal mets and with any liver mets
per_patients <- unique(sample_info[group=='Peritoneum','Patient_ID',with=F][[1]])

si_per <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Peritoneum')) & in_collapsed==T]
res_per <- get_met_specific_distances(si_per, ad_table, comparison='primary', distance='angular', return_tree=F)
res_per$group1 <- 'Peritoneum'
res_per <- merge(res_per, sample_info[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res_per <- merge(res_per, sample_info[,c('Patient_ID','Real_Sample_ID','met_timing'),with=F], by.x=c('patient','sample1'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res_per$vertical <- factor(res_per$vertical, levels=c('mucosal/luminal','deep'))
res_sync <- res_per[group1=='Peritoneum' & met_timing=='synchronous' & !is.na(vertical) & !is.na(met_timing),]
res_meta <- res_per[group1=='Peritoneum' & met_timing=='metachronous' & !is.na(vertical) & !is.na(met_timing),]


test <- function(i, res) {
    if(i > 0) {
        shuffled <- sample(res$vertical, replace=F)
        res$vertical <- shuffled
    }
    mu_deep <- mean(res$distance[res$vertical=='deep'])
    mu_luminal <- mean(res$distance[res$vertical=='mucosal/luminal'])
    ratio <- mu_deep / mu_luminal
    list(i=i, ratio=ratio)
}

run_test <- function(res, R) {
    set.seed(42)
    l <- lapply(0:R, test, res)
    l <- rbindlist(l)
    obs <- l$ratio[l$i==0]
    exp <- l$ratio[l$i > 0]
    x <- sum(exp <= obs) + 1
    n <- R + 1
    pval <- x / n
    list(pva=pval, x=x, n=n, obs=obs, exp=exp)
}

tst_sync <- run_test(res_sync, R=1e4)
tst_meta <- run_test(res_meta, R=1e4)

pvals <- data.table(pval=c(tst_sync$pva, tst_meta$pva), group=c('Synchronous','Metachronous'))
obs <- data.table(obs=c(tst_sync$obs, tst_meta$obs), group=c('Synchronous','Metachronous'))
pvals <- merge(pvals, obs, by='group')
pvals$pval <- paste0('p=',prettyNum(pvals$pval,digits=2))
dat <- rbind(data.table(tst_sync$exp, group='Synchronous'),
             data.table(tst_meta$exp, group='Metachronous'))
names(dat)[1] <- 'ratio'
dat$group <- factor(dat$group, levels=(c('Synchronous','Metachronous')))
pvals$group <- factor(pvals$group, levels=(c('Synchronous','Metachronous')))
obs$group <- factor(obs$group, levels=(c('Synchronous','Metachronous')))

p <- ggplot(dat, aes(x=ratio)) +
    geom_histogram(fill='#bfbfbf',color='white',linewidth=0.75,bins=50) + 
    facet_wrap(facets=~group, ncol=1) +
    theme_ang(base_size=12) +
    geom_vline(data=obs, aes(xintercept=obs), color='red') + 
    geom_text(data=pvals, y=900, aes(x=obs, label=pval), hjust=-0.1, color='red') +
    labs(x='Deep-PT to Met / Luminal-PT to Met',title='Permutation test for enriched deep/PM-timing association', subtitle='All pairwise combos of PTs and PMs.') 
ggsave(here('figures/unused_fig_depth_angular_distance_by_pm_timing_permutation.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# unused fig.
# Min node distance to deep/luminal PT per group
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get patients with any peritoneal mets and with any liver mets
per_patients <- unique(sample_info[group=='Peritoneum','Patient_ID',with=F][[1]])
liv_patients <- unique(sample_info[group=='Liver','Patient_ID',with=F][[1]])

si_per <- sample_info[Patient_ID %in% per_patients & (group %in% c('Normal','Primary','Peritoneum')) & in_collapsed==T]
res_per <- get_met_specific_distances(si_per, ad_table, comparison='primary', distance='node', return_tree=F)
res_per$group1 <- 'Peritoneum'
si_lr <- sample_info[Patient_ID %in% per_patients & cohort=='peritoneal' & (group %in% c('Normal','Primary','Locoregional')) & in_collapsed==T ]
res_lr <- get_met_specific_distances(si_lr, ad_table, comparison='primary', distance='node', return_tree=F)
res_lr$group1 <- 'Locoregional'
si_liv <- sample_info[Patient_ID %in% liv_patients & (group %in% c('Normal','Primary','Liver')) & in_collapsed==T ]
res_liv <- get_met_specific_distances(si_liv, ad_table, comparison='primary', distance='node', return_tree=F)
res_liv$group1 <- 'Liver'
res <- rbind(res_per, res_lr, res_liv)
res <- merge(res, sample_info[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res <- merge(res, sample_info[,c('Patient_ID','Real_Sample_ID','group'),with=F], by.x=c('patient','sample1'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res$vertical <- factor(res$vertical, levels=c('mucosal/luminal','deep'))
res <- res[!is.na(vertical)]
res[,class:=paste0(group,'-',vertical)]
res[class=='Liver-mucosal/luminal',class:='Liv_luminal']
res[class=='Liver-deep',class:='Liv_deep']
res[class=='Peritoneum-mucosal/luminal',class:='Per_luminal']
res[class=='Peritoneum-deep',class:='Per_deep']
res[class=='Locoregional-mucosal/luminal',class:='LR_luminal']
res[class=='Locoregional-deep',class:='LR_deep']
res[grepl('LN',sample1),class:=gsub('LR','LN',class)]
res[grepl('TD',sample1),class:=gsub('LR','TD',class)]
res[grepl('LN',sample1),group:=gsub('Locoregional','Lymph node',group)]
res[grepl('TD',sample1),group:=gsub('Locoregional','Tumor deposit',group)]
res$group <- factor(res$group, levels=c('Lymph node','Tumor deposit','Peritoneum','Liver'))

collapse <- function(res) {
    min <- min(res$distance,na.rm=T)
    list(min=min)
}
collapsed <- res[,collapse(.SD),by=c('patient','group','vertical')]
collapsed <- data.table::melt(collapsed, id.vars=c('patient','group','vertical'))
collapsed <- data.table::dcast(variable + patient + group ~ vertical, value.var='value', data=collapsed)
collapsed <- collapsed[!is.na(deep) & !is.na(`mucosal/luminal`)]
setnames(collapsed,'variable','distance_type')
collapsed <- data.table::melt(collapsed, id.vars=c('patient','group','distance_type'))
names(collapsed) <- c('patient','group','distance_type','vertical','distance')
stat.test <- mywilcox2(collapsed, distance ~ vertical, facet_field='group', paired=T)

tmp_cols <- c(group_cols,'#bc4822','#ef7b55')
names(tmp_cols)[8:9] <- c('Lymph node','Tumor deposit')
p <- ggplot(collapsed, aes(x=vertical, y=distance)) +
    scale_y_continuous(limits=c(0,2.0),breaks=seq(0,2.0,by=0.5)) +
    geom_line(aes(group=patient)) +
    geom_point(pch=21,color='black',size=4,aes(fill=group)) +
    scale_fill_manual(values=tmp_cols) +
    facet_wrap(facets=~group, nrow=1) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02, size=2.75) +
    theme_ang(base_size=12) +
    guides(fill='none') +
    labs(x='Invasion depth of primary tumor', y='Node distance / average overall node distance',
         title='Min met-specific node distance to primary tumor',subtitle='Subset by vertical invasion of the primary tumor')
ggsave(here('figures/unused_fig_depth_node_distance_paired.pdf'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# unused fig, formerly SI Fig 4B (i) node distance from PM to deep/luminal PT (PM)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get patients with any peritoneal mets
per_patients <- unique(sample_info[group=='Peritoneum','Patient_ID',with=F][[1]])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & ((group %in% c('Normal','Primary')) | (group=='Peritoneum' & met_timing=='synchronous'))]
res_per_syn <- get_met_specific_distances(si, ad_table, comparison='primary', distance='node', return_tree=F)
res_per_syn <- merge(res_per_syn, sample_info[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res_per_syn[group1=='Metastasis',group1:='Peritoneum']
res_per_syn$vertical <- factor(res_per_syn$vertical, levels=c('mucosal/luminal','deep'))
res_per_syn[,group:=paste0(group1,'-',vertical)]
res_per_syn$met_timing <- 'synchronous'

si <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & ((group %in% c('Normal','Primary')) | (group=='Peritoneum' & met_timing=='metachronous'))]
res_per_meta <- get_met_specific_distances(si, ad_table, comparison='primary', distance='node', return_tree=F)
res_per_meta <- merge(res_per_meta, sample_info[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
res_per_meta[group1=='Metastasis',group1:='Peritoneum']
res_per_meta$vertical <- factor(res_per_meta$vertical, levels=c('mucosal/luminal','deep'))
res_per_meta[,group:=paste0(group1,'-',vertical)]
res_per_meta$met_timing <- 'metachronous'

res_per <- rbind(res_per_syn, res_per_meta)
res_per <- res_per[!is.na(vertical)]
res_per[vertical=='mucosal/luminal',class:='Per_luminal']
res_per[vertical=='deep',class:='Per_deep']
res_per$vertical <- factor(res_per$vertical, levels=c('mucosal/luminal','deep'))
res_per$met_timing <- factor(res_per$met_timing, levels=c('synchronous','metachronous'))

p <- ggplot(res_per, aes(x=vertical, y=distance)) +
    geom_point(aes(fill=class),color='black',position=position_jitter(width=0.15,height=0,seed=2), pch=21, size=4) + 
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) +
    scale_fill_manual(values=cols_depth) +
    scale_color_manual(values=cols_timing) +
    facet_wrap(facets=~met_timing) + 
    guides(fill='none', color='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    stat_compare_means(method="wilcox.test", label.y=2.25, 
                       comparisons=list(c("mucosal/luminal",'deep'))) +
    labs(x='Invasion depth of primary tumor', y='Node distance / average overall node distance',title='Node distance from primary tumor to peritoneal mets', subtitle='Subset by vertical invasion of the primary tumor\n(excludes meta after syn)')
ggsave(here('figures/unused_fig_depth_node_distance_by_pm_timing.pdf'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# unused fig. formerly Fig 4C. Clinical T-stage vs PM iming
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_for_data <- function(dat,title=NULL) {
    ## dat should be a data.frame 2x2 contingency table, rows = 3,4 (stage), columns = 'synchronous','metachronous' 

    ## test for significant difference between stage vs timing
    tst <- fisher.test(dat)
    n <- rowSums(dat)
    x <- dat[,'metachronous']
    ci <- as.data.table(binom.confint(x, n, conf.level = 0.95, methods = "exact"))

    ## prep for plotting as proportions
    prop <- dat / n
    prop <- cbind(stage=rownames(prop), as.data.table(prop))
    prop$lwr <- ci$lower; prop$upr <- ci$upper
    prop$n_meta <- ci$x; prop$n_syn <- ci$n - ci$x
    prop <- as.data.table(melt(prop, id.vars=c('stage','lwr','upr','n_meta','n_syn')))
    prop[variable=='synchronous',lwr:=NA]
    prop[variable=='synchronous',upr:=NA]
    prop[variable=='synchronous',cnt:=n_syn]
    prop[variable=='metachronous',cnt:=n_meta]
    prop[,n_meta:=NULL]; prop[,n_syn:=NULL]
    prop$variable <- factor(prop$variable, levels=c('synchronous','metachronous'))

    ## get position for bar labels with Ns
    prop <- prop[order(stage,variable,decreasing=T),]
    prop[c(1,3), y_pos:=value*0.5]
    prop$y_pos[2] <- prop$value[1] + prop$value[2]*0.5
    prop$y_pos[4] <- prop$value[3] + prop$value[4]*0.5
    cols <- brewer.pal(5,'Accent')[1:2]
    names(cols) <- c('synchronous','metachronous')
    siglab <- paste0('p=',prettyNum(tst$p.value,digits=2))

    ## make plot
    p <- ggplot(prop, aes(x=stage, y=value)) +
        scale_y_continuous(expand=c(0,0), limits=c(0,1.1), breaks=seq(0,1,by=0.25)) +
        geom_bar(stat='identity', aes(fill=variable),color='black',size=0.25) +
        geom_text(aes(y=y_pos,label=cnt)) +
        scale_fill_manual(values=cols,name='Metastasis timing') +
        geom_errorbar(aes(min=lwr, max=upr), width=0.1) +
        theme_ang(base_size=12) +
        geom_signif(annotation=siglab, y_position=1.025, comparison=list(c('3','4'))) +
        labs(x='T stage', y='Proportion',title=title)
    p
}


## extract timing and T,N stage
clin <- fread(here('original_data/misc/HIPECDatabase_timing_mets_and_T&Nstage.csv'))
names(clin)[1:2] <- c('ID','timing')
clin <- clin[pN!='Nx' & pT!='Tx',]
clin <- clin %>%
    mutate(N=as.integer(str_extract(pN, "[0-9]")), T=as.integer(str_extract(pT, "[0-9]")),
          timing=fct_relevel(timing, list("synchronous", "metachronous")))
dat <- as.data.frame.matrix(xtabs(~ T + timing, data=clin[T %in% c(3,4)]))
p1 <- run_for_data(dat,title='HIPEC database')

dat2 <- read_distance_matrix(here('original_data/misc/Tstage_Jayne.txt'))
dat2 <- as.data.frame(dat2)
rownames(dat2) <- gsub('T','',rownames(dat2))
dat2 <- dat2[c('3','4'),]
colnames(dat2) <- tolower(colnames(dat2))
dat2 <- dat2[,1:2]
p2 <- run_for_data(dat2, title='Jayne et al.')

p <- plot_grid(p1, p2, ncol=2, rel_widths=c(2,2,1))
ggsave(here('figures/ed_fig_6_tstage_timing.pdf'))  ## formerly Fig 4c



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig 5B-C. Common/distinct origin of PM and distant mets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

si <- copy(sample_info)
si[group %in% c('Liver','Lung','Distant (other)'), group:='Distant']

## get patients with any peritoneal mets and with any distant mets
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
liv_patients <- unique(si[group=='Distant','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, liv_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]

test_patient_met_level <- function(patient, si, ad_table, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary','Peritoneum', 'Distant')]

    ## get the original set of valid-samples for the true number of distant/pt/pm samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    n_dm <- length(unique(si$Real_Sample_ID[si$group=='Distant']))
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))
    n_pm <- length(unique(si$Real_Sample_ID[si$group=='Peritoneum']))

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
        groups <- si[,c('Real_Sample_ID','group','tissue_type'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        setnames(ad,'tissue_type','tissue_type1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        setnames(ad,'tissue_type','tissue_type2')
        ad <- ad[,c('Var1','Var2','group1','group2','tissue_type1','tissue_type2','value'),with=F] 

        ## get the min distances between PMs and DMs
        pm_dm <- ad[group1=='Peritoneum' & group2=='Distant']
        pm_dm <- pm_dm[order(Var1, value, decreasing=F),]
        pm_dm <- pm_dm[!duplicated(Var1),]
        pm_pt <- ad[group1=='Peritoneum' & group2=='Primary']
        pm_pt <- pm_pt[order(Var1, value, decreasing=F),]
        pm_pt <- pm_pt[!duplicated(Var1),]

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        out <- merge(pm_dm[,c('Var1','Var2','tissue_type2','value'),with=F], pm_pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T)
        names(out) <- c('sample','dm_sample','dm_type','dm_distance','pt_sample','pt_distance')
        to_add <- data.table(sample='average', dm_sample='average', dm_type='average', dm_distance=mean(out$dm_distance), pt_sample='average', pt_distance=mean(out$pt_distance))
        out <- rbind(out, to_add)
        out$i <- i
        out$patient <- patient
        out
    }

    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(dm_distance / pt_distance)]
    
    collapse_pm <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_pm(.SD),by=c('patient','sample')] 
    res$n_dm <- n_dm
    res$n_pt <- n_pt
    res$n_pm <- n_pm
    res
}

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Shared (95%)','Shared (80%)','n.s.','Distinct (80%)','Distinct (95%)')

RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res <- rbindlist(l)
write_tsv(res,here('processed_data/min_distance_ratios_pm_to_dm_vs_pt.txt'))


res <- fread(here('processed_data/min_distance_ratios_pm_to_dm_vs_pt.txt'))
avg <- res[sample=='average']
avg <- avg[order(obs,decreasing=F),]
avg[,n_dm_div_n_pm:=n_dm / n_pm]
avg$patient <- factor(avg$patient, levels=avg$patient)
avg$origin <- 'n.s.'
avg[`97.5%` < 0, origin:='Shared (95%)']
avg[`90%` < 0, origin:='Shared (80%)']
avg[`10%` > 0, origin:='Distinct (80%)']
avg[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(avg, aes(x=patient, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='PM-DM origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', linewidth=0.5) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio (PM:DM vs PM:PT)') +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) 

# save this for possible use later
tiles <- avg[,c('patient','n_pt','n_dm','n_pm'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm <- tiles$n_dm / max(tiles$n_dm)
tiles$frac_pm <- tiles$n_pm / max(tiles$n_pm)
tiles_n <- data.table::melt(tiles[,c('patient','n_pt','n_dm','n_pm'),with=F], id.var='patient')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('patient','frac_pt','frac_dm','frac_pm'),with=F], id.var='patient')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('patient','variable'), all=T)
names(tiles) <- c('patient','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='pm',type:='Peritoneum']
tiles[type=='dm',type:='Distant']

fit <- lm(obs ~ n_pm + n_dm + n_pt, data=avg)
coefs <- coef(summary(fit))
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_pm',type:='Peritoneum']
info[type=='n_dm',type:='Distant']

p2 <- ggplot(tiles, aes(x=patient, y=type)) +
    geom_tile(aes(fill=norm)) +
    geom_text(aes(label=n)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    labs(x='Patient', y='Number of samples')

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', size=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank()) +
    scale_y_continuous(position='right',limits=c(-0.5,0.25),expand=c(0,0)) + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(2,1), rel_widths=c(3,1.5))
ggsave(here('figures/fig_5c_pm_dm_distinct_origin.pdf'),width=10,height=8)


## repeat on the lesion-level
pd <- res[sample!='average']
pd <- pd[order(obs,decreasing=F),]
pd[,n_dm_div_n_pm:=n_dm / n_pm]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`90%` < 0, origin:='Shared (80%)']
pd[`97.5%` < 0, origin:='Shared (95%)']
pd[`10%` > 0, origin:='Distinct (80%)']
pd[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(pd, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='PM-DM origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', linewidth=0.5) +
    #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio (PM:DM vs PM:PT)')

tiles <- pd[,c('id','n_pt','n_dm','n_pm'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm <- tiles$n_dm / max(tiles$n_dm)
tiles$frac_pm <- tiles$n_pm / max(tiles$n_pm)
tiles_n <- data.table::melt(tiles[,c('id','n_pt','n_dm','n_pm'),with=F], id.var='id')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('id','frac_pt','frac_dm','frac_pm'),with=F], id.var='id')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('id','variable'), all=T)
names(tiles) <- c('id','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='pm',type:='Peritoneum']
tiles[type=='dm',type:='Distant']

p2 <- ggplot(tiles, aes(x=id, y=type)) +
    geom_tile(aes(fill=norm)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    #theme(axis.text.x=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position='none') +
    labs(x='Patient', y='Number of samples')

fit <- lm(obs ~ n_pm + n_dm + n_pt, data=pd)
coefs <- coef(summary(fit))
ps <- summary(fit)$coef[2:4,4]
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_pm',type:='Peritoneum']
info[type=='n_dm',type:='Distant']

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', linewidth=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank()) +
    scale_y_continuous(position='right',limits=c(-0.2,0.05),expand=c(0,0)) + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(3,1.5), rel_widths=c(3,1.5))
ggsave(here('figures/fig_5b_pm_dm_distinct_origin_metlevel_withIDs.pdf'),width=20,height=11)

#x <- pd[!is.na(`10%`)]
#nrow(x)
#length(unique(x$patient))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version of Fig 5B-C with LNs (common/distinct origin of LN and distant mets)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

si <- copy(sample_info)
si[group %in% c('Liver','Lung','Distant (other)'), group:='Distant']

## get patients with any peritoneal mets and with any distant mets
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
si[Patient_ID %in% per_patients & tissue_type %in% 'Lymph node', group:='Lymph node']
ln_patients <- unique(si[group=='Lymph node','Patient_ID',with=F][[1]])
liv_patients <- unique(si[group=='Distant','Patient_ID',with=F][[1]])
valid_patients <- intersect(ln_patients, liv_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]


test_patient_met_level <- function(patient, si, ad_table, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary','Lymph node', 'Distant')]

    ## get the original set of valid-samples for the true number of distant/pt/ln samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    n_dm <- length(unique(si$Real_Sample_ID[si$group=='Distant']))
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))
    n_ln <- length(unique(si$Real_Sample_ID[si$group=='Lymph node']))

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
        groups <- si[,c('Real_Sample_ID','group','tissue_type'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        setnames(ad,'tissue_type','tissue_type1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        setnames(ad,'tissue_type','tissue_type2')
        ad <- ad[,c('Var1','Var2','group1','group2','tissue_type1','tissue_type2','value'),with=F] 

        ## get the min distances between PMs and DMs
        ln_dm <- ad[group1=='Lymph node' & group2=='Distant']
        ln_dm <- ln_dm[order(Var1, value, decreasing=F),]
        ln_dm <- ln_dm[!duplicated(Var1),]
        ln_pt <- ad[group1=='Lymph node' & group2=='Primary']
        ln_pt <- ln_pt[order(Var1, value, decreasing=F),]
        ln_pt <- ln_pt[!duplicated(Var1),]

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        out <- merge(ln_dm[,c('Var1','Var2','tissue_type2','value'),with=F], ln_pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T)
        names(out) <- c('sample','dm_sample','dm_type','dm_distance','pt_sample','pt_distance')
        to_add <- data.table(sample='average', dm_sample='average', dm_type='average', dm_distance=mean(out$dm_distance), pt_sample='average', pt_distance=mean(out$pt_distance))
        out <- rbind(out, to_add)
        out$i <- i
        out$patient <- patient
        out
    }

    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(dm_distance / pt_distance)]
    
    collapse_ln <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_ln(.SD),by=c('patient','sample')] 
    res$n_dm <- n_dm
    res$n_pt <- n_pt
    res$n_ln <- n_ln
    res
}

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Shared (95%)','Shared (80%)','n.s.','Distinct (80%)','Distinct (95%)')

RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res <- rbindlist(l)
write_tsv(res,here('processed_data/min_distance_ratios_ln_to_dm_vs_pt.txt'))

avg <- res[sample=='average']
avg <- avg[order(obs,decreasing=F),]
avg[,n_dm_div_n_ln:=n_dm / n_ln]
avg$patient <- factor(avg$patient, levels=avg$patient)
avg$origin <- 'n.s.'
avg[`97.5%` < 0, origin:='Shared (95%)']
avg[`90%` < 0, origin:='Shared (80%)']
avg[`10%` > 0, origin:='Distinct (80%)']
avg[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(avg, aes(x=patient, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='LN-DM origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio (LN:DM vs LN:PT)') +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) 

# save this for possible use later
tiles <- avg[,c('patient','n_pt','n_dm','n_ln'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm <- tiles$n_dm / max(tiles$n_dm)
tiles$frac_ln <- tiles$n_ln / max(tiles$n_ln)
tiles_n <- data.table::melt(tiles[,c('patient','n_pt','n_dm','n_ln'),with=F], id.var='patient')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('patient','frac_pt','frac_dm','frac_ln'),with=F], id.var='patient')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('patient','variable'), all=T)
names(tiles) <- c('patient','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='ln',type:='Peritoneum']
tiles[type=='dm',type:='Distant']

fit <- lm(obs ~ n_ln + n_dm + n_pt, data=avg)
coefs <- coef(summary(fit))
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_ln',type:='Peritoneum']
info[type=='n_dm',type:='Distant']

p2 <- ggplot(tiles, aes(x=patient, y=type)) +
    geom_tile(aes(fill=norm)) +
    geom_text(aes(label=n)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    labs(x='Patient', y='Number of samples')

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', size=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank())
    #scale_y_continuous(position='right',limits=c(-0.5,0.25),expand=c(0,0)) + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(2,1), rel_widths=c(3,1.5))
ggsave(here('figures/fig_5c_ln_dm_distinct_origin.pdf'),width=10,height=8)


## repeat on the lesion-level
pd <- res[sample!='average']
pd <- pd[order(obs,decreasing=F),]
pd[,n_dm_div_n_ln:=n_dm / n_ln]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`97.5%` < 0, origin:='Shared (95%)']
pd[`90%` < 0, origin:='Shared (80%)']
pd[`10%` > 0, origin:='Distinct (80%)']
pd[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(pd, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='LN-DM origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio (LN:DM vs LN:PT)')

tiles <- pd[,c('id','n_pt','n_dm','n_ln'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm <- tiles$n_dm / max(tiles$n_dm)
tiles$frac_ln <- tiles$n_ln / max(tiles$n_ln)
tiles_n <- data.table::melt(tiles[,c('id','n_pt','n_dm','n_ln'),with=F], id.var='id')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('id','frac_pt','frac_dm','frac_ln'),with=F], id.var='id')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('id','variable'), all=T)
names(tiles) <- c('id','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='ln',type:='Lymph node']
tiles[type=='dm',type:='Distant']

p2 <- ggplot(tiles, aes(x=id, y=type)) +
    geom_tile(aes(fill=norm)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    #theme(axis.text.x=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position='none') +
    labs(x='Patient', y='Number of samples')

fit <- lm(obs ~ n_ln + n_dm + n_pt, data=pd)
coefs <- coef(summary(fit))
ps <- summary(fit)$coef[2:4,4]
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_ln',type:='Lymph node']
info[type=='n_dm',type:='Distant']

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', size=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank()) + 
    scale_y_continuous(position='right') + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(3,1.5), rel_widths=c(3,1.5))
ggsave(here('figures/fig_5b_ln_dm_distinct_origin_metlevel_withIDs.pdf'),width=20,height=11)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version of Fig 5B-C with TDs (common/distinct origin of TD and distant mets)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

si <- copy(sample_info)
si[group %in% c('Liver','Lung','Distant (other)'), group:='Distant']

## get patients with any peritoneal mets and with any distant mets
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
si[Patient_ID %in% per_patients & tissue_type %in% 'Tumor deposit', group:='Tumor deposit']
td_patients <- unique(si[group=='Tumor deposit','Patient_ID',with=F][[1]])
liv_patients <- unique(si[group=='Distant','Patient_ID',with=F][[1]])
valid_patients <- intersect(td_patients, liv_patients)
valid_patients <- valid_patients[grepl('Mm',valid_patients)==F]


test_patient_met_level <- function(patient, si, ad_table, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary','Tumor deposit', 'Distant')]

    ## get the original set of valid-samples for the true number of distant/pt/td samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    n_dm <- length(unique(si$Real_Sample_ID[si$group=='Distant']))
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))
    n_td <- length(unique(si$Real_Sample_ID[si$group=='Tumor deposit']))

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
        groups <- si[,c('Real_Sample_ID','group','tissue_type'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        setnames(ad,'tissue_type','tissue_type1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        setnames(ad,'tissue_type','tissue_type2')
        ad <- ad[,c('Var1','Var2','group1','group2','tissue_type1','tissue_type2','value'),with=F] 

        ## get the min distances between PMs and DMs
        td_dm <- ad[group1=='Tumor deposit' & group2=='Distant']
        td_dm <- td_dm[order(Var1, value, decreasing=F),]
        td_dm <- td_dm[!duplicated(Var1),]
        td_pt <- ad[group1=='Tumor deposit' & group2=='Primary']
        td_pt <- td_pt[order(Var1, value, decreasing=F),]
        td_pt <- td_pt[!duplicated(Var1),]

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        out <- merge(td_dm[,c('Var1','Var2','tissue_type2','value'),with=F], td_pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T)
        names(out) <- c('sample','dm_sample','dm_type','dm_distance','pt_sample','pt_distance')
        to_add <- data.table(sample='average', dm_sample='average', dm_type='average', dm_distance=mean(out$dm_distance), pt_sample='average', pt_distance=mean(out$pt_distance))
        out <- rbind(out, to_add)
        out$i <- i
        out$patient <- patient
        out
    }

    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(dm_distance / pt_distance)]
    
    collapse_td <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_td(.SD),by=c('patient','sample')] 
    res$n_dm <- n_dm
    res$n_pt <- n_pt
    res$n_td <- n_td
    res
}

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Shared (95%)','Shared (80%)','n.s.','Distinct (80%)','Distinct (95%)')

RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res <- rbindlist(l)
write_tsv(res,here('processed_data/min_distance_ratios_td_to_dm_vs_pt.txt'))

avg <- res[sample=='average']
avg <- avg[order(obs,decreasing=F),]
avg[,n_dm_div_n_td:=n_dm / n_td]
avg$patient <- factor(avg$patient, levels=avg$patient)
avg$origin <- 'n.s.'
avg[`97.5%` < 0, origin:='Shared (95%)']
avg[`90%` < 0, origin:='Shared (80%)']
avg[`10%` > 0, origin:='Distinct (80%)']
avg[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(avg, aes(x=patient, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='TD-DM origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio (TD:DM vs TD:PT)') +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) 

# save this for possible use later
tiles <- avg[,c('patient','n_pt','n_dm','n_td'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm <- tiles$n_dm / max(tiles$n_dm)
tiles$frac_td <- tiles$n_td / max(tiles$n_td)
tiles_n <- data.table::melt(tiles[,c('patient','n_pt','n_dm','n_td'),with=F], id.var='patient')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('patient','frac_pt','frac_dm','frac_td'),with=F], id.var='patient')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('patient','variable'), all=T)
names(tiles) <- c('patient','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='td',type:='Peritoneum']
tiles[type=='dm',type:='Distant']

fit <- lm(obs ~ n_td + n_dm + n_pt, data=avg)
coefs <- coef(summary(fit))
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_td',type:='Peritoneum']
info[type=='n_dm',type:='Distant']

p2 <- ggplot(tiles, aes(x=patient, y=type)) +
    geom_tile(aes(fill=norm)) +
    geom_text(aes(label=n)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    labs(x='Patient', y='Number of samples')

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', size=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank())
    #scale_y_continuous(position='right',limits=c(-0.5,0.25),expand=c(0,0)) + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(2,1), rel_widths=c(3,1.5))
ggsave(here('figures/fig_5c_td_dm_distinct_origin.pdf'),width=10,height=8)


## repeat on the lesion-level
pd <- res[sample!='average']
pd <- pd[order(obs,decreasing=F),]
pd[,n_dm_div_n_td:=n_dm / n_td]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`97.5%` < 0, origin:='Shared (95%)']
pd[`90%` < 0, origin:='Shared (80%)']
pd[`10%` > 0, origin:='Distinct (80%)']
pd[`2.5%` > 0, origin:='Distinct (95%)']

p1 <- ggplot(pd, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='TD-DM origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio (TD:DM vs TD:PT)')

tiles <- pd[,c('id','n_pt','n_dm','n_td'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm <- tiles$n_dm / max(tiles$n_dm)
tiles$frac_td <- tiles$n_td / max(tiles$n_td)
tiles_n <- data.table::melt(tiles[,c('id','n_pt','n_dm','n_td'),with=F], id.var='id')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('id','frac_pt','frac_dm','frac_td'),with=F], id.var='id')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('id','variable'), all=T)
names(tiles) <- c('id','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='td',type:='Tumor deposit']
tiles[type=='dm',type:='Distant']

p2 <- ggplot(tiles, aes(x=id, y=type)) +
    geom_tile(aes(fill=norm)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    #theme(axis.text.x=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position='none') +
    labs(x='Patient', y='Number of samples')

fit <- lm(obs ~ n_td + n_dm + n_pt, data=pd)
coefs <- coef(summary(fit))
ps <- summary(fit)$coef[2:4,4]
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_td',type:='Tumor deposit']
info[type=='n_dm',type:='Distant']

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', size=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank()) +
    scale_y_continuous(position='right') + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(3,1.5), rel_widths=c(3,1.5))
ggsave(here('figures/fig_5b_td_dm_distinct_origin_metlevel_withIDs.pdf'),width=20,height=11)






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Version of Fig 5B-C with Liv mets (common/distinct origin of Liver mets with other Liver mets)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

si <- copy(sample_info)

## get patients with any peritoneal mets and with any distant mets
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
si[group %in% c('Liver','Lung','Distant (other)'), group:='Distant']

test_patient_met_level <- function(patient, si, ad_table, ncpus) { 
    message(patient)
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary','Distant')]

    ## get the original set of valid-samples for the true number of distant/pt/liv samples
    ad <- ad_table[[patient]]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    si <- si[Real_Sample_ID %in% valid_samples,]
    n_liv <- length(unique(si$Real_Sample_ID[si$tissue_type=='Liver']))
    n_dm <- length(unique(si$Real_Sample_ID[si$group=='Distant']))
    n_pt <- length(unique(si$Real_Sample_ID[si$group=='Primary']))

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
        groups <- si[,c('Real_Sample_ID','group','tissue_type'),with=F]
        ad <- merge(ad, groups, by.x='Var1', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group1')
        setnames(ad,'tissue_type','tissue_type1')
        ad <- merge(ad, groups, by.x='Var2', by.y='Real_Sample_ID', all.x=T)
        setnames(ad,'group','group2')
        setnames(ad,'tissue_type','tissue_type2')
        ad <- ad[Var1!=Var2,]

        ## get the min distances between PMs and DMs
        liv_dm <- ad[tissue_type1=='Liver' & group2=='Distant']
        liv_dm <- liv_dm[order(Var1, value, decreasing=F),]
        liv_dm <- liv_dm[!duplicated(Var1),]
        liv_pt <- ad[tissue_type1=='Liver' & group2=='Primary']
        liv_pt <- liv_pt[order(Var1, value, decreasing=F),]
        liv_pt <- liv_pt[!duplicated(Var1),]

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        out <- merge(liv_dm[,c('Var1','Var2','tissue_type2','value'),with=F], liv_pt[,c('Var1','Var2','value'),with=F], by='Var1', all=T)
        names(out) <- c('sample','dm_sample','dm_type','dm_distance','pt_sample','pt_distance')
        to_add <- data.table(sample='average', dm_sample='average', dm_type='average', dm_distance=mean(out$dm_distance), pt_sample='average', pt_distance=mean(out$pt_distance))
        out <- rbind(out, to_add)
        out$i <- i
        out$patient <- patient
        out
    }
    #debugonce(test_patient_bs)
    #test_patient_bs(1, patient, si, ad_table, bs)
 
    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(dm_distance / pt_distance)]
    closest <- res[i==0,c('patient','sample','dm_sample','dm_type'),with=F]
    names(closest) <- c('patient','sample','closest_sample','closest_sample_type')
    collapse_liv <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_liv(.SD),by=c('patient','sample')] 
    res <- merge(res, closest, by=c('patient','sample'), all.x=T)
    res$n_liv <- n_liv
    res$n_dm <- n_dm
    res$n_pt <- n_pt
    res
}

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Shared (95%)','Shared (80%)','n.s.','Distinct (80%)','Distinct (95%)')

RNGkind("L'Ecuyer-CMRG") 
set.seed(42)
l <- lapply(valid_patients, test_patient_met_level, si, ad_table, ncpus)
res <- rbindlist(l)
write_tsv(res,here('processed_data/min_distance_ratios_liv_to_dm_vs_liv_to_pt.txt'))






## repeat on the lesion-level
res <- fread(here('processed_data/min_distance_ratios_liv_to_dm_vs_liv_to_pt.txt'))
pd <- res[sample!='average']
pd <- pd[order(obs,decreasing=F),]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`97.5%` < 0, origin:='Shared (95%)']
pd[`90%` < 0, origin:='Shared (80%)']
pd[`10%` > 0, origin:='Distinct (80%)']
pd[`2.5%` > 0, origin:='Distinct (95%)']
pd[,n_dm_nonliv:=n_dm - n_liv]
p1 <- ggplot(pd, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Liv-DM origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', linewidth=0.5) +
    #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio (Liv:DM vs Liv:PT)')

tiles <- pd[,c('id','n_pt','n_liv','n_dm_nonliv'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm_nonliv <- tiles$n_dm_nonliv / max(tiles$n_dm_nonliv)
tiles$frac_liv <- tiles$n_liv / max(tiles$n_liv)
tiles_n <- data.table::melt(tiles[,c('id','n_pt','n_liv','n_dm_nonliv'),with=F], id.var='id')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('id','frac_pt','frac_liv','frac_dm_nonliv'),with=F], id.var='id')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('id','variable'), all=T)
names(tiles) <- c('id','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='liv',type:='Liver']
tiles[type=='dm_nonliv',type:='Distant (non-liver)']

p2 <- ggplot(tiles, aes(x=id, y=type)) +
    geom_tile(aes(fill=norm)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    #theme(axis.text.x=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position='none') +
    labs(x='Patient', y='Number of samples')

fit <- lm(obs ~ n_liv + n_pt + n_dm_nonliv, data=pd)
coefs <- coef(summary(fit))
ps <- summary(fit)$coef[2:4,4]
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_liv',type:='Liver']
info[type=='n_dm_nonliv',type:='Distant (non-liver)']

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', linewidth=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank()) +
    scale_y_continuous(position='right') + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(3,1.5), rel_widths=c(3,1.5))
ggsave(here('figures/min_distance_ratios_liv_to_dm_vs_pt_metlevel_withIDs.pdf'),width=20,height=11)


## compare distributions of min-dist-ratios for any difference between PM and LNs
res_pm <- fread(here('processed_data/min_distance_ratios_pm_to_dm_vs_pt.txt'))
res_pm$group <- 'Peritoneum'
res_pm[, n_test:=n_pm]
res_pm[, n_comparitor:=n_dm]

res_ln <- fread(here('processed_data/min_distance_ratios_ln_to_dm_vs_pt.txt'))
res_ln$group <- 'Locoregional'
res_ln$tissue_type <- 'Lymph node'
res_ln[,n_test:=n_ln]
res_ln[,n_comparitor:=n_dm]

res_liv <- fread(here('processed_data/min_distance_ratios_liv_to_dm_vs_liv_to_pt.txt'))
res_liv$group <- 'Liver'
res_liv[, n_test:=n_liv]
res_liv[, n_comparitor:=n_dm-1]

res_td <- fread(here('processed_data/min_distance_ratios_td_to_dm_vs_pt.txt'))
res_td$group <- 'Locoregional'
res_td$tissue_type <- 'Tumor deposit'
res_td[,n_test:=n_td]
res_td[,n_comparitor:=n_dm]

res <- rbind(res_pm, res_ln, res_liv, res_td, fill=T)
res <- res[sample!='average',]
res[is.na(tissue_type), tissue_type:=group]
res[,comparitor_pt_ratio:=n_comparitor/n_pt]
res$i <- 1:nrow(res)
out <- res[,c('patient','sample','obs','2.5%','10%','90%','97.5%','group','tissue_type','n_ln','n_liv','n_pm','n_dm','n_pt','n_test','n_comparitor','comparitor_pt_ratio'),with=F]
write_tsv(out,here('processed_data/min_distance_ratio_vs_comparitor_ratio_raw_data.txt'))


res <- fread(here('processed_data/min_distance_ratio_vs_comparitor_ratio_raw_data.txt'))
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
    geom_point(position=position_jitter(width=0.1, height=0, seed=42), pch=21, size=2.25, aes(fill=tissue_type), stroke=0.25) +
    geom_boxplot(fill=NA, outlier.shape=NA, width=0.4) +
    scale_fill_manual(values=tmp_cols) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = NA) +
    theme_ang(base_size=10) +
    guides(fill='none') +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))  +
    labs(x='Met type', y='log2 min-distance-ratio\n(sample to DM vs sample to PT)')
ggsave(here('figures/unk_5b_distinct_origin_distributions_3cat.pdf'))


m1 <- lm(obs ~ group + comparitor_pt_ratio, data=res[group %in% c('Peritoneum','Locoregional')]); summary(m1)
m2 <- lm(obs ~ group + comparitor_pt_ratio, data=res[group %in% c('Peritoneum','Liver')]); summary(m2)
m3 <- lm(obs ~ group + comparitor_pt_ratio, data=res[group %in% c('Locoregional','Liver')]); summary(m3)
m1 <- lm(obs ~ I(group!='Liver') + comparitor_pt_ratio, data=res); summary(m1)

res$group_factor <- factor(res$group, levels=c('Liver','Peritoneum','Locoregional'))

m1 <- lm(obs ~ group_factor + comparitor_pt_ratio, data=res); summary(m1)






stat.test1$y.position <- 3.25
stat.test2 <- mywilcox2(res[group %in% c('Peritoneum','Liver')], obs ~ group, paired=F)
stat.test2$y.position <- 3.5
stat.test3 <- mywilcox2(res[group %in% c('Locoregional','Liver')], obs ~ group, paired=F)
stat.test3$y.position <- 3.75





## repeat on the lesion-level
res <- fread(here('processed_data/min_distance_ratio_vs_comparitor_ratio_raw_data.txt'))
pd <- res[sample!='average']
pd <- pd[order(obs,decreasing=F),]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`90%` < 0, origin:='Shared (80%)']
pd[`97.5%` < 0, origin:='Shared (95%)']
pd[`10%` > 0, origin:='Distinct (80%)']
pd[`2.5%` > 0, origin:='Distinct (95%)']
pd[,n_dm_nonliv:=n_dm - n_liv]
pd$origin <- factor(pd$origin, levels=c('Shared (95%)','Shared (80%)','n.s.','Distinct (80%)','Distinct (95%)'))
to_percentile <- function(pd) {
    pd <- pd[order(obs,decreasing=F),]
    pd$percentile <- (1:nrow(pd)) / nrow(pd)
    pd
}
pd <- pd[,to_percentile(.SD),by=group]

p1 <- ggplot(pd, aes(x=percentile, y=obs)) +
    #scale_y_continuous(limits=c(-3.5,4.0),breaks=seq(-3,4,by=1)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    facet_wrap(facets=~group, ncol=1, scale='free_y') +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio')
ggsave(here('figures/unk_5b_distinct_origin_min_dist_ratios_3_cat_stacked.pdf'), height=10, width=6.5)


p1 <- ggplot(pd, aes(x=percentile, y=obs)) +
    scale_y_continuous(limits=c(-3.5,6.5),breaks=seq(-3,4,by=1)) +
    geom_text(aes(label=id, color=origin), y=4, angle=90, hjust=0, vjust=0.5, size=2) + 
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Origin\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    facet_wrap(facets=~group, ncol=1, scale='free_y') +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank()) + 
    guides(color='none') +
    labs(x=NULL, y='log2 min-distance-ratio')
ggsave(here('figures/unk_5b_distinct_origin_min_dist_ratios_3_cat_stacked_labeled.pdf'), height=10, width=10)

## quantification via boxplot
p2 <- ggplot(pd, aes(x=group, y=obs)) +
    geom_hline(yintercept=0, color='#bfbfbf', linetype='dashed', linewidth=0.25) + 
    geom_point(position=position_jitter(width=0.1, height=0, seed=42), pch=21, size=2.25, aes(fill=tissue_type), stroke=0.25) +
    geom_boxplot(fill=NA, outlier.shape=NA, width=0.4) +
    scale_fill_manual(values=tmp_cols) +
    stat_pvalue_manual(stat.test, label = "label", tip.length = NA) +
    theme_ang(base_size=10) +
    guides(fill='none') +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))  +
    labs(x='Met type', y='log2 min-distance-ratio\n(sample to DM vs sample to PT)')
ggsave(here('figures/unk_5b_distinct_origin_min_dist_ratios_3_cat_quantification.pdf')) # height=10, width=6.5)





m <- as.data.frame.matrix(xtabs(~group + origin, data=pd))

test <- function(g1,g2,m) {
    tmp <- m[c(g1,g2),]
    p_ca <- CochranArmitageTest(tmp)$p.value
    p_chi <- chisq.test(tmp)$p.value
    data.table(group1=g1, group2=g2, p_ca=p_ca, p_chi=p_chi)    
}
t1 <- test('Peritoneum','Liver',m)
t2 <- test('Peritoneum','Locoregional',m)
t3 <- test('Locoregional','Liver',m)
tests <- rbind(t1, t2, t3)
tests[,label:=paste0('p=',prettyNum(p_chi,digits=2))]


#tbl$origin <- factor(tbl$origin, levels=c('Common (95%)','Common (80%)','n.s.','Distinct (80%)','Distinct (95%)'))

prop <- copy(m)
for(i in 1:3) prop[i,] <- prop[i,] / sum(prop[i,])
dat <- cbind(group=rownames(m),as.data.table(m))
dat <- data.table::melt(dat, id.var='group')
prop <- cbind(group=rownames(prop),as.data.table(prop))
prop <- data.table::melt(prop, id.var='group')
dat$prop <- prop$value
dat$variable <- factor(dat$variable, levels=rev(colnames(m)))
get_label_pos <- function(dat) {
    dat <- dat[order(variable,decreasing=T),]
    dat$pos <- (cumsum(dat$prop) - 0.5*dat$prop)
    dat
}
dat2 <- dat[,get_label_pos(.SD), by=c('group')]
dat2$group <- factor(dat2$group, levels=c('Peritoneum','Locoregional','Liver'))

p <- ggplot(dat2, aes(x=group, y=prop)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,by=0.25), limits=c(0,1.3)) +
    geom_bar(stat='identity', aes(fill=variable),color='black',linewidth=0.25) +
    geom_text(data=dat2[value > 0],aes(label=value,y=pos)) +
    labs(x='Group',y='Proportion of lesions') +
    scale_fill_manual(values=(cols), name='Origin with sync. tumor deposits') + 
    stat_pvalue_manual(data=tests, label='label', y.position=c(1.1,1.05,1.15), hjust=0.5, size=3, angle=0, tip.length=NA) +
    theme_ang(base_size=12) +
    theme(legend.position='right')
ggsave(here('figures/unk_5b_distinct_origin_min_dist_ratios_3_cat_summary.pdf'), height=7, width=5)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# unused figs of distance from PT (either deep or luminal) to each PT depth
# e.g. pairwise AD from luminal PTs to either deep PT, or other luminal PTs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define applicable patients
all_patients <- unique(sample_info$Patient_ID[!is.na(sample_info$vertical)])

## subset the sample_info for N, PT, and Per; get met-spec node distance from peritoneum to normal
si_deep <- sample_info[Patient_ID %in% all_patients & (group!='Primary' | (group=='Primary' & vertical=='deep')) & in_collapsed==T]
res_deep <- get_met_specific_distances(si_deep, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_deep <- res_deep[group1=='Primary' & group2=='Primary']
res_deep$vertical1 <- 'deep'
res_deep$vertical2 <- 'deep'

si_luminal <- sample_info[Patient_ID %in% all_patients & (group!='Primary' | (group=='Primary' & vertical=='mucosal/luminal')) & in_collapsed==T]
res_luminal <- get_met_specific_distances(si_luminal, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_luminal <- res_luminal[group1=='Primary' & group2=='Primary']
res_luminal$vertical1 <- 'mucosal/luminal'
res_luminal$vertical2 <- 'mucosal/luminal'

si_deep_vs_luminal <- sample_info[Patient_ID %in% all_patients & (group!='Primary' | (group=='Primary' & vertical %in% c('deep','mucosal/luminal'))) & in_collapsed==T]
res_deep_vs_luminal <- get_met_specific_distances(si_deep_vs_luminal, ad_table, comparison='pairwise', distance='angular', return_tree=F)
res_deep_vs_luminal <- res_deep_vs_luminal[group1=='Primary' & group2=='Primary']
res_deep_vs_luminal <- merge(res_deep_vs_luminal, 
                             sample_info[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample1'),
                             by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
setnames(res_deep_vs_luminal,'vertical','vertical1')
res_deep_vs_luminal <- merge(res_deep_vs_luminal, 
                             sample_info[,c('Patient_ID','Real_Sample_ID','vertical'),with=F], by.x=c('patient','sample2'),
                             by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)
setnames(res_deep_vs_luminal,'vertical','vertical2')
res_deep_vs_luminal <- res_deep_vs_luminal[vertical1!=vertical2]

depth_cols <- c('#006230','#7fc5a2')
names(depth_cols) <- c('deep','mucosal/luminal')

## repeat for distance from luminal to either deep or luminal
res <- rbind(res_luminal, res_deep_vs_luminal)
res$from <- 'mucosal/luminal'
res[vertical1=='mucosal/luminal' & vertical2=='mucosal/luminal', to:='mucosal/luminal']
res[vertical1=='deep' | vertical2=='deep', to:='deep']
stat.test <- mywilcox2(res, distance ~ to, paired=F)

p <- ggplot(res, aes(x=to, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=to)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    scale_fill_manual(values=depth_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x='Primary tumor depth',y='Angular distance',title='Pairwise AD from luminal/mucosal PT to PT of given depth')
ggsave(here('figures/unused_fig_depth_angular_distance_from_luminal_to_PT.pdf'))


## repeat for distance from deep to either deep or luminal
res <- rbind(res_deep, res_deep_vs_luminal)
res$from <- 'deep'
res[vertical1=='deep' & vertical2=='deep', to:='deep']
res[vertical1=='mucosal/luminal' | vertical2=='mucosal/luminal', to:='mucosal/luminal']
stat.test <- mywilcox2(res, distance ~ to, paired=F)

p <- ggplot(res, aes(x=to, y=distance)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=2),pch=21,color='black',size=4,aes(fill=to)) +
    geom_boxplot(fill=NA,color='black',outlier.shape=NA) + 
    stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) +
    scale_fill_manual(values=depth_cols) +
    guides(fill='none') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
    labs(x='Primary tumor depth',y='Angular distance',title='Pairwise AD from deep PT to PT of given depth')
ggsave(here('figures/unused_fig_depth_angular_distance_from_deep_to_PT.pdf'))





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# association between synchronous and deep-invading PT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# patient/timing-level version:

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

## get the distance from each PM to all PTs in that patient, then annotate the PT depth
si <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & group %in% c('Normal','Primary','Peritoneum')]
res <- get_met_specific_distances(si, ad_table, comparison='primary', distance='angular', return_tree=F)
depth <- si[group=='Primary',c('Patient_ID','Real_Sample_ID','vertical')]
res <- merge(res, depth, by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)

## annotate each patient as whether they had any synchronous PMs
timing <- si[group=='Peritoneum',c('Patient_ID','Real_Sample_ID','met_timing'),with=F]
res <- merge(res, timing, by.x=c('patient','sample1'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)

## for each patient, get the min distance between any PM and a deep vs luminal PT
get_min_distance <- function(res) {
    min_d <- min(res$distance[res$vertical=='deep'])
    min_l <- min(res$distance[res$vertical=='mucosal/luminal'])
    list(min_d=min_d, min_l=min_l)
}
mins <- res[,get_min_distance(.SD), by=c('patient','met_timing')]
mins <- mins[!is.infinite(min_d) & !is.infinite(min_l)]
mins[,ratio:=log2(min_d/min_l)]
mins <- mins[order(ratio,decreasing=F),]
mins$ID <- 1:nrow(mins)
mins[met_timing=='metachronous after synchronous', met_timing:='metachronous']
wilcox.test(mins$ratio[mins$met_timing=='synchronous'], mu=0)
wilcox.test(mins$ratio[mins$met_timing!='synchronous'], mu=0)





# met-level version:

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

## get the distance from each PM to all PTs in that patient, then annotate the PT depth
si <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & group %in% c('Normal','Primary','Peritoneum')]
res <- get_met_specific_distances(si, ad_table, comparison='primary', distance='angular', return_tree=F)
depth <- si[group=='Primary',c('Patient_ID','Real_Sample_ID','vertical')]
res <- merge(res, depth, by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)

## annotate each patient as whether they had any synchronous PMs
timing <- si[group=='Peritoneum',c('Patient_ID','Real_Sample_ID','met_timing'),with=F]
res <- merge(res, timing, by.x=c('patient','sample1'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)

## for each patient, get the min distance between any PM and a deep vs luminal PT
get_min_distance <- function(res) {
    min_d <- min(res$distance[res$vertical=='deep'])
    min_l <- min(res$distance[res$vertical=='mucosal/luminal'])
    list(min_d=min_d, min_l=min_l)
}
mins <- res[,get_min_distance(.SD), by=c('patient','sample1','met_timing')]
mins <- mins[!is.infinite(min_d) & !is.infinite(min_l)]
mins[,ratio:=log2(min_d/min_l)]
mins <- mins[order(ratio,decreasing=F),]
mins$ID <- 1:nrow(mins)
mins[met_timing=='metachronous after synchronous', met_timing:='synchronous']
wilcox.test(mins$min_d[mins$met_timing=='synchronous'], mins$min_l[mins$met_timing=='synchronous'], paired=T)
wilcox.test(mins$min_d[mins$met_timing=='metachronous'], mins$min_l[mins$met_timing=='metachronous'], paired=T)

wilcox.test(mins$min_d[mins$met_timing=='synchronous'], mins$min_d[mins$met_timing=='metachronous'])
wilcox.test(mins$ratio[mins$met_timing=='synchronous'], mu=0)
wilcox.test(mins$ratio[mins$met_timing!='synchronous'], mu=0)


p <- ggplot(mins, aes(x=ID, y=ratio)) +
    geom_bar(stat='identity', aes(fill=met_timing)) 



# original version

## define applicable patients (peritoneum, with per cohort + C38, C89)
per_patients <- unique(sample_info$Patient_ID[sample_info$group=='Peritoneum'])

## get the distance from each PM to all PTs in that patient, then annotate the PT depth
si <- sample_info[Patient_ID %in% per_patients & in_collapsed==T & group %in% c('Normal','Primary','Peritoneum')]
res <- get_met_specific_distances(si, ad_table, comparison='primary', distance='angular', return_tree=F)
depth <- si[group=='Primary',c('Patient_ID','Real_Sample_ID','vertical')]
res <- merge(res, depth, by.x=c('patient','sample2'), by.y=c('Patient_ID','Real_Sample_ID'), all.x=T)

## for each patient, get the min distance between any PM and a deep vs luminal PT
get_min_distance <- function(res) {
    min_d <- min(res$distance[res$vertical=='deep'])
    min_l <- min(res$distance[res$vertical=='mucosal/luminal'])
    list(min_d=min_d, min_l=min_l)
}
mins <- res[,get_min_distance(.SD), by=c('patient')]

## annotate each patient as whether they had any synchronous PMs
timing <- si[group=='Peritoneum',c('Patient_ID','Real_Sample_ID','met_timing'),with=F]
any_timing <- function(timing) {
    any_sync <- any(timing$met_timing=='synchronous')
    any_meta <- any(timing$met_timing %in% c('metachronous','metachronous after synchronous'))
    list(any_sync=any_sync, any_meta=any_meta)
}
timing <- timing[,any_timing(.SD),by=Patient_ID]
mins <- merge(mins, timing, by.x='patient', by.y='Patient_ID', all.x=T)
mins <- mins[!is.infinite(min_d) & !is.infinite(min_l)]
mins[,ratio:=log2(min_d/min_l)]
mins <- mins[order(ratio,decreasing=F),]
mins$patient <- factor(mins$patient, levels=mins$patient)

p <- ggplot(mins, aes(x=patient, y=ratio)) +
    geom_bar(stat='identity', aes(fill=any_meta)) +
    geom_text(aes(label=patient))
p <- ggplot(mins, aes(x=patient, y=ratio)) +
    geom_bar(stat='identity', aes(fill=any_sync)) +
    geom_text(aes(label=patient))

wilcox.test(mins$ratio[mins$any_sync==T], mu=0)
wilcox.test(mins$ratio[mins$any_sync==F], mu=0)
wilcox.test(mins$ratio[mins$any_meta==T], mu=0)
wilcox.test(mins$ratio[mins$any_meta==F], mu=0)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2024-04-29
# Test for enriched co-mingling between mets of a given type and deep vs luminal PTs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_met_depth <- function(patient, si, ad_table, met_group, ncpus) { 
    message(patient)

    ## load AD matrix, subset metadata for samples with AD data
    ad <- as.data.frame(ad_table[[patient]])
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary',met_group)]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    valid_met_samples <- intersect(si[group %in% met_group,c(Real_Sample_ID)], valid_samples)
    valid_pt_samples <- intersect(si[group %in% 'Primary',c(Real_Sample_ID)], valid_samples)
    si <- si[Real_Sample_ID %in% valid_samples,]

    ## get number of samples for QC
    n_deep <- length(unique(si$Real_Sample_ID[si$group=='Primary' & si$vertical=='deep']))
    n_luminal <- length(unique(si$Real_Sample_ID[si$group=='Primary' & si$vertical=='mucosal/luminal']))
    n_met <- length(unique(si$Real_Sample_ID[si$group %in% met_group]))
    if(n_met < 1 | n_deep < 1 | n_luminal < 1) {
        message('Insufficient samples.')
        return(NULL)
    }

    ## load the bootstrapped data for this patient
    if(patient %in% c('Ea3','Eb3','Ec3')) {
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/E3.rds'))
    } else {
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/',patient,'.rds'))
    }
    bs <- readRDS(bs_file)

    test_patient_bs <- function(i, patient, si, ad_table, bs) {
        if(i > 0) {
            ## for the bootstrap replicates, obtain new ad and valid_samples
            ad <- as.data.frame(bs[[i]])
        }
        ad <- ad[valid_met_samples, valid_pt_samples]
        ad$sample <- rownames(ad)
        ad <- as.data.table(reshape2::melt(ad,id.vars='sample'))
        setnames(ad,'variable','pt')
        groups <- si[,c('Real_Sample_ID','group'),with=F]
        depths <- si[,c('Real_Sample_ID','vertical'),with=F]
        ad <- merge(ad, groups, by.x='sample', by.y='Real_Sample_ID', all.x=T)
        ad <- merge(ad, depths, by.x='pt', by.y='Real_Sample_ID', all.x=T)

        ## get the min distances between PMs and DMs
        dist_deep <- ad[vertical=='deep',c('sample','pt','value'),with=F]
        dist_deep <- dist_deep[order(sample,value,decreasing=F),]
        dist_deep <- dist_deep[!duplicated(sample),]
        names(dist_deep) <- c('sample','nearest_deep','dist_deep')
        dist_luminal <- ad[vertical=='mucosal/luminal',c('sample','pt','value'),with=F]
        dist_luminal <- dist_luminal[order(sample,value,decreasing=F),]
        dist_luminal <- dist_luminal[!duplicated(sample),]
        names(dist_luminal) <- c('sample','nearest_luminal','dist_luminal')

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        dist <- merge(dist_deep, dist_luminal, by='sample', all=T)
        dist$i <- i
        dist$patient <- patient
        dist
    }
    #debugonce(test_patient_bs)
    #test_patient_bs(1, patient, si, ad_table, bs)
    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(dist_deep / dist_luminal)]
    
    collapse_td <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_td(.SD),by=c('patient','sample')] 
    res$n_deep <- n_deep
    res$n_luminal <- n_luminal
    res$n_met <- n_met
    res
}

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#e5e5e5'
names(cols) <- c('Deep (95%)','Deep (80%)','n.s.','Luminal (80%)','Luminal (95%)')

si <- copy(sample_info)
deep_patients <- si[vertical=='deep',(Patient_ID)]
luminal_patients <- si[vertical=='mucosal/luminal',(Patient_ID)]
patients_with_both_deep_and_luminal <- intersect(deep_patients,luminal_patients)

per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])

## PM (any)
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si, ad_table, met_group='Peritoneum', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum'
write_tsv(res,here('processed_data/min_distance_ratios_pm_to_deep_vs_luminal.txt'))

## Locoregional (either)
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si, ad_table, met_group='Locoregional', ncpus=4)
res <- rbindlist(l)
res$group <- 'Locoregional'
write_tsv(res,here('processed_data/min_distance_ratios_lr_to_deep_vs_luminal.txt'))

## Liver
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si, ad_table, met_group='Liver', ncpus=4)
res <- rbindlist(l)
res$group <- 'Liver'
write_tsv(res,here('processed_data/min_distance_ratios_liv_to_deep_vs_luminal.txt'))

## Any distant
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si, ad_table, met_group=c('Liver','Lung','Distant (other)'), ncpus=4)
res <- rbindlist(l)
res$group <- 'Distant (any)'
write_tsv(res,here('processed_data/min_distance_ratios_dm_to_deep_vs_luminal.txt'))


res1 <- fread(here('processed_data/min_distance_ratios_pm_to_deep_vs_luminal.txt'))
res2 <- fread(here('processed_data/min_distance_ratios_lr_to_deep_vs_luminal.txt'))
res3 <- fread(here('processed_data/min_distance_ratios_liv_to_deep_vs_luminal.txt'))
res4 <- fread(here('processed_data/min_distance_ratios_dm_to_deep_vs_luminal.txt'))
res <- rbind(res1, res2, res3, res4)


## repeat on the lesion-level
res$group <- factor(res$group, levels=unique(res$group))
pd <- copy(res)
pd <- pd[order(obs,decreasing=F),]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=unique(pd$id))
pd$origin <- 'n.s.'
pd[`90%` < 0, origin:='Deep (80%)']
pd[`97.5%` < 0, origin:='Deep (95%)']
pd[`10%` > 0, origin:='Luminal (80%)']
pd[`2.5%` > 0, origin:='Luminal (95%)']
pd$origin <- factor(pd$origin, levels=rev(names(cols)))

p1 <- ggplot(pd, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Nearest PT Depth:\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    facet_wrap(facets=~group, ncol=1, scale='free_x') +
    #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), legend.position='right') + 
    labs(x=NULL, y='log2 min-distance-ratio (Met:Deep vs Met:Luminal)', subtitle='Met vs PT-depth min-distance-ratios')

x <- as.data.frame.matrix(xtabs(~ group + origin, data=pd))
for(i in 1:nrow(x)) x[i,] <- x[i,] / sum(x[i,])
x$group <- rownames(x)
xm <- reshape2::melt(x, id.var='group')
xm$variable <- factor(xm$variable, levels=rev(names(cols)))
xm$group <- factor(xm$group, levels=unique(res$group))

x <- as.data.frame.matrix(xtabs(~ group + origin, data=pd))
chisq.test(x[c('Peritoneum','Distant (any)'),])
chisq.test(x[c('Locoregional','Distant (any)'),])
chisq.test(x[c('Locoregional','Peritoneum'),])

p2 <- ggplot(xm, aes(x=group, y=value)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,by=0.25)) +
    geom_bar(stat='identity', aes(fill=variable),color='black',linewidth=0.25) +
    scale_fill_manual(values=cols,name='Nearest PT Depth:\n(confidence)') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position='none') + 
    labs(x='Metastasis type', y='Fraction of lesions')

p <- plot_grid(p1, p2, nrow=1, rel_widths=c(2,1))
ggsave(here('figures/met_deep_luminal_mindistance_pmpatients_only.pdf'),width=14, height=9)






## PM (untreated+synchronous)
si2 <- copy(si)
si2[met_timing=='synchronous' & met_treated=='untreated' & group=='Peritoneum', group:='Peritoneum (untreated+synchronous)']
per_patients <- unique(si2[group=='Peritoneum (untreated+synchronous)','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Peritoneum (untreated+synchronous)', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum (untreated+synchronous)'
write_tsv(res,here('processed_data/min_distance_ratios_pm_untreatedsync_to_deep_vs_luminal.txt'))

## PM (treated+metachronous)
si2 <- copy(si)
#si2[met_timing %in% c('metachronous','metachronous after synchronous') & met_treated %in% c('treated','treated after untreated') & group=='Peritoneum', group:='Peritoneum (treated+metachronous)']
si2[met_timing %in% c('metachronous') & met_treated %in% c('treated') & group=='Peritoneum', group:='Peritoneum (treated+metachronous)']
per_patients <- unique(si2[group=='Peritoneum (treated+metachronous)','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Peritoneum (treated+metachronous)', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum (treated+metachronous)'
write_tsv(res,here('processed_data/min_distance_ratios_pm_treatedmeta_to_deep_vs_luminal.txt'))

## PM (treated(any)+metachronous(any))
si2 <- copy(si)
si2[met_timing %in% c('metachronous','metachronous after synchronous') & met_treated %in% c('treated','treated after untreated') & group=='Peritoneum', group:='Peritoneum (treated(any)+metachronous(any))']
per_patients <- unique(si2[group=='Peritoneum (treated(any)+metachronous(any))','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Peritoneum (treated(any)+metachronous(any))', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum (treated(any)+metachronous(any))'
write_tsv(res,here('processed_data/min_distance_ratios_pm_treatedmetaany_to_deep_vs_luminal.txt'))

## Lymph node
si2 <- copy(si)
si2[tissue_type %in% 'Lymph node', group:=tissue_type]
ln_patients <- unique(si2[(cohort=='peritoneal' | Patient_ID %in% c('C38','C89')) & group=='Lymph node','Patient_ID',with=F][[1]])
valid_patients <- intersect(ln_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Lymph node', ncpus=4)
res <- rbindlist(l)
res$group <- 'Lymph node'
write_tsv(res,here('processed_data/min_distance_ratios_ln_to_deep_vs_luminal.txt'))

## Tumor deposit
si2 <- copy(si)
si2[tissue_type %in% 'Tumor deposit', group:=tissue_type]
td_patients <- unique(si2[(cohort=='peritoneal' | Patient_ID %in% c('C38','C89')) & group=='Tumor deposit','Patient_ID',with=F][[1]])
valid_patients <- intersect(td_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Tumor deposit', ncpus=4)
res <- rbindlist(l)
res$group <- 'Tumor deposit'
write_tsv(res,here('processed_data/min_distance_ratios_td_to_deep_vs_luminal.txt'))


res_pm <- fread(here('processed_data/min_distance_ratios_pm_to_deep_vs_luminal.txt'))
res_pm1 <- fread(here('processed_data/min_distance_ratios_pm_untreatedsync_to_deep_vs_luminal.txt'))
res_pm2 <- fread(here('processed_data/min_distance_ratios_pm_treatedmeta_to_deep_vs_luminal.txt'))
res_pm3 <- fread(here('processed_data/min_distance_ratios_pm_treatedmetaany_to_deep_vs_luminal.txt'))
res_ln <- fread(here('processed_data/min_distance_ratios_ln_to_deep_vs_luminal.txt'))
res_td <- fread(here('processed_data/min_distance_ratios_td_to_deep_vs_luminal.txt'))
#res_lr <- fread(here('processed_data/min_distance_ratios_lr_to_deep_vs_luminal.txt'))
res_liv <- fread(here('processed_data/min_distance_ratios_liv_to_deep_vs_luminal.txt'))
res <- rbind(res_pm, res_pm1, res_pm2, res_ln, res_td, res_liv)
res$group <- factor(res$group, levels=unique(res$group))


## plot results for each lesion showing samples from the same patient
get_plot_for_group <- function(this.group, res) {
    message(this.group)
    pd <- res[group==this.group]
    get_patient_center <- function(pd) {
        center <- median(pd$obs,na.rm=T)
        list(center=center)
    }
    pt <- pd[,get_patient_center(.SD),by=patient]
    pt <- pt[order(center),]
    pd$patient <- factor(pd$patient, levels=pt$patient)
    pd <- pd[order(patient,obs,decreasing=F),]
    pd$id <- paste(pd$patient, pd$sample)
    pd$id <- factor(pd$id, levels=pd$id)
    pd$pos <- as.integer(pd$id)

    collapse_patient <- function(pd) {
        start <- min(pd$pos)-0.5
        end <- max(pd$pos)+0.5
        list(start=start,end=end)
    }
    rects <- pd[,collapse_patient(.SD),by=c('group','patient')]
    rects <- rects[order(patient),]
    rects$fill <- rep(c('white','#f7f7f7'),length.out=nrow(rects))
    rects[,center:=0.5*(start+end)]

    #pd <- pd[order(group,patient,obs,decreasing=F),]
    pd$origin <- 'n.s.'
    pd[`90%` < 0, origin:='Deep (80%)']
    pd[`97.5%` < 0, origin:='Deep (95%)']
    pd[`10%` > 0, origin:='Luminal (80%)']
    pd[`2.5%` > 0, origin:='Luminal (95%)']
    pd$origin <- factor(pd$origin, levels=rev(names(cols)))

    p1 <- ggplot(pd) +
        scale_x_continuous(expand=c(0.02,0.02)) + 
        scale_y_continuous(limits=c(-2.5,2.5),breaks=seq(-2,2,by=1.0)) + 
        geom_rect(data=rects,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),color=NA,fill=rects$fill) +
        geom_text(data=rects,aes(x=center,label=patient),y=2,angle=90) + #hjust=0,vjust=0.5) +
        geom_point(aes(x=pos, y=obs, color=origin)) +
        geom_errorbar(aes(x=pos, y=obs, min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
        scale_color_manual(values=cols,name='Nearest PT Depth:\n(confidence)') +
        theme_ang(base_size=12) +
        geom_hline(yintercept=0, color='black', size=0.5) +
        theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(),legend.position='none') + 
        labs(x=NULL, y='log2 MDR', subtitle=this.group)
    p1
}
plot_list <- lapply(c('Peritoneum (untreated+synchronous)','Peritoneum (treated+metachronous)','Lymph node','Tumor deposit','Liver'), get_plot_for_group, res)
p <- plot_grid(plotlist=plot_list, ncol=1, align='v')
ggsave(here('figures/met_deep_luminal_mindistance_patientnames.pdf'),width=12, height=12)


## repeat on the lesion-level
res <- rbind(res_pm1, res_pm2, res_ln, res_td, res_liv)
res$group <- factor(res$group, levels=unique(res$group))
pd <- copy(res)
pd <- pd[order(obs,decreasing=F),]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`90%` < 0, origin:='Deep (80%)']
pd[`97.5%` < 0, origin:='Deep (95%)']
pd[`10%` > 0, origin:='Luminal (80%)']
pd[`2.5%` > 0, origin:='Luminal (95%)']
pd$origin <- factor(pd$origin, levels=rev(names(cols)))

p1 <- ggplot(pd, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Nearest PT Depth:\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    facet_wrap(facets=~group, ncol=1, scale='free_x') +
    #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), legend.position='right') + 
    labs(x=NULL, y='log2 min-distance-ratio (Met:Deep vs Met:Luminal)', subtitle='Met vs PT-depth min-distance-ratios')

x <- as.data.frame.matrix(xtabs(~ group + origin, data=pd))
for(i in 1:nrow(x)) x[i,] <- x[i,] / sum(x[i,])
x$group <- rownames(x)
xm <- reshape2::melt(x, id.var='group')
xm$variable <- factor(xm$variable, levels=rev(names(cols)))
xm$group <- factor(xm$group, levels=unique(res$group))

p2 <- ggplot(xm, aes(x=group, y=value)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,by=0.25)) +
    geom_bar(stat='identity', aes(fill=variable),color='black',linewidth=0.25) +
    scale_fill_manual(values=cols,name='Nearest PT Depth:\n(confidence)') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position='none') + 
    labs(x='Metastasis type', y='Fraction of lesions')

p <- plot_grid(p1, p2, nrow=1, rel_widths=c(2,1))
ggsave(here('figures/met_deep_luminal_mindistance.pdf'),width=14, height=9)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2024-04-30
# alternative test for enriched co-mingling between mets of a given type and deep vs luminal PTs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

patient <- 'E20'
met_group <- 'Peritoneum'

test_met_depth <- function(patient, si, ad_table, met_group, ncpus) { 
    message(patient)

    ## load AD matrix, subset metadata for samples with AD data
    ad <- as.data.frame(ad_table[[patient]])
    si <- si[Patient_ID %in% patient & in_collapsed==T & group %in% c('Normal','Primary',met_group)]
    valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
    valid_met_samples <- intersect(si[group %in% met_group,c(Real_Sample_ID)], valid_samples)
    valid_pt_samples <- intersect(si[group %in% 'Primary',c(Real_Sample_ID)], valid_samples)
    si <- si[Real_Sample_ID %in% valid_samples,]

    ## get number of samples for QC
    n_deep <- length(unique(si$Real_Sample_ID[si$group=='Primary' & si$vertical=='deep']))
    n_luminal <- length(unique(si$Real_Sample_ID[si$group=='Primary' & si$vertical=='mucosal/luminal']))
    n_met <- length(unique(si$Real_Sample_ID[si$group %in% met_group]))
    if(n_met < 1 | n_deep < 1 | n_luminal < 1) {
        message('Insufficient samples.')
        return(NULL)
    }

    ## load the bootstrapped data for this patient
    if(patient %in% c('Ea3','Eb3','Ec3')) {
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/E3.rds'))
    } else {
        bs_file <- here(paste0('processed_data/angular_distance_matrices_bootstrapped/',patient,'.rds'))
    }
    bs <- readRDS(bs_file)

    test_patient_bs <- function(i, patient, si, ad_table, bs) {
        if(i > 0) {
            ## for the bootstrap replicates, obtain new ad and valid_samples
            ad <- as.data.frame(bs[[i]])
        }
        ad <- ad[valid_met_samples, valid_pt_samples]
        ad$sample <- rownames(ad)
        ad <- as.data.table(reshape2::melt(ad,id.vars='sample'))
        setnames(ad,'variable','pt')
        groups <- si[,c('Real_Sample_ID','group'),with=F]
        depths <- si[,c('Real_Sample_ID','vertical'),with=F]
        ad <- merge(ad, groups, by.x='sample', by.y='Real_Sample_ID', all.x=T)
        ad <- merge(ad, depths, by.x='pt', by.y='Real_Sample_ID', all.x=T)

        ## get the min distances between PMs and DMs
        dist_deep <- ad[vertical=='deep',c('sample','pt','value'),with=F]
        dist_deep <- dist_deep[order(sample,value,decreasing=F),]
        collapse <- function(d) {
            mid <- median(d$value,na.rm=T)
            list(mid=mid)
        }
        dist_deep <- dist_deep[,collapse(.SD),by=sample]
        names(dist_deep) <- c('sample','dist_deep')
        dist_luminal <- ad[vertical=='mucosal/luminal',c('sample','pt','value'),with=F]
        dist_luminal <- dist_luminal[order(sample,value,decreasing=F),]
        dist_luminal <- dist_luminal[,collapse(.SD),by=sample]
        names(dist_luminal) <- c('sample','dist_luminal')

        ## merge the PM-level nearest DM and PT, and return the values for this bootstrap replicate
        dist <- merge(dist_deep, dist_luminal, by='sample', all=T)
        dist$i <- i
        dist$patient <- patient
        dist
    }
    #debugonce(test_patient_bs)
    #test_patient_bs(1, patient, si, ad_table, bs)
    l <- mclapply(0:length(bs), test_patient_bs, patient, si, ad_table, bs, mc.cores=ncpus)
    res <- rbindlist(l,fill=TRUE)
    res[, logratio:=log2(dist_deep / dist_luminal)]
    
    collapse_td <- function(res) {
        obs <- res$logratio[res$i==0]
        boots <- res$logratio[res$i > 0]
        qs <- as.list(quantile(boots,c(0.025,0.1,0.90,0.975)))
        append(list(obs=obs), qs)
    }
    res <- res[,collapse_td(.SD),by=c('patient','sample')] 
    res$n_deep <- n_deep
    res$n_luminal <- n_luminal
    res$n_met <- n_met
    res
}

cols <- brewer.pal(5,'PRGn')
cols[3] <- '#d8d8d8'
names(cols) <- c('Deep (95%)','Deep (80%)','n.s.','Luminal (80%)','Luminal (95%)')

si <- copy(sample_info)
deep_patients <- si[vertical=='deep',(Patient_ID)]
luminal_patients <- si[vertical=='mucosal/luminal',(Patient_ID)]
patients_with_both_deep_and_luminal <- intersect(deep_patients,luminal_patients)


## PM (any)
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si, ad_table, met_group='Peritoneum', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum'
write_tsv(res,here('processed_data/median_distance_ratios_pm_to_deep_vs_luminal.txt'))

## PM (untreated+synchronous)
si2 <- copy(si)
si2[met_timing=='synchronous' & met_treated=='untreated' & group=='Peritoneum', group:='Peritoneum (untreated+synchronous)']
per_patients <- unique(si2[group=='Peritoneum (untreated+synchronous)','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Peritoneum (untreated+synchronous)', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum (untreated+synchronous)'
write_tsv(res,here('processed_data/median_distance_ratios_pm_untreatedsync_to_deep_vs_luminal.txt'))

## PM (treated+metachronous)
si2 <- copy(si)
#si2[met_timing %in% c('metachronous','metachronous after synchronous') & met_treated %in% c('treated','treated after untreated') & group=='Peritoneum', group:='Peritoneum (treated+metachronous)']
si2[met_timing %in% c('metachronous') & met_treated %in% c('treated') & group=='Peritoneum', group:='Peritoneum (treated+metachronous)']
per_patients <- unique(si2[group=='Peritoneum (treated+metachronous)','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Peritoneum (treated+metachronous)', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum (treated+metachronous)'
write_tsv(res,here('processed_data/median_distance_ratios_pm_treatedmeta_to_deep_vs_luminal.txt'))

## PM (treated(any)+metachronous(any))
si2 <- copy(si)
si2[met_timing %in% c('metachronous','metachronous after synchronous') & met_treated %in% c('treated','treated after untreated') & group=='Peritoneum', group:='Peritoneum (treated(any)+metachronous(any))']
per_patients <- unique(si2[group=='Peritoneum (treated(any)+metachronous(any))','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Peritoneum (treated(any)+metachronous(any))', ncpus=4)
res <- rbindlist(l)
res$group <- 'Peritoneum (treated(any)+metachronous(any))'
write_tsv(res,here('processed_data/median_distance_ratios_pm_treatedmetaany_to_deep_vs_luminal.txt'))

## Lymph node
si2 <- copy(si)
si2[tissue_type %in% 'Lymph node', group:=tissue_type]
ln_patients <- unique(si2[(cohort=='peritoneal' | Patient_ID %in% c('C38','C89')) & group=='Lymph node','Patient_ID',with=F][[1]])
valid_patients <- intersect(ln_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Lymph node', ncpus=4)
res <- rbindlist(l)
res$group <- 'Lymph node'
write_tsv(res,here('processed_data/median_distance_ratios_ln_to_deep_vs_luminal.txt'))

## Tumor deposit
si2 <- copy(si)
si2[tissue_type %in% 'Tumor deposit', group:=tissue_type]
td_patients <- unique(si2[(cohort=='peritoneal' | Patient_ID %in% c('C38','C89')) & group=='Tumor deposit','Patient_ID',with=F][[1]])
valid_patients <- intersect(td_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si2, ad_table, met_group='Tumor deposit', ncpus=4)
res <- rbindlist(l)
res$group <- 'Tumor deposit'
write_tsv(res,here('processed_data/median_distance_ratios_td_to_deep_vs_luminal.txt'))

## Locoregional (either)
per_patients <- unique(si[group=='Peritoneum','Patient_ID',with=F][[1]])
valid_patients <- intersect(per_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si, ad_table, met_group='Locoregional', ncpus=4)
res <- rbindlist(l)
res$group <- 'Locoregional'
write_tsv(res,here('processed_data/median_distance_ratios_lr_to_deep_vs_luminal.txt'))

## Liver
liv_patients <- unique(si[group=='Liver','Patient_ID',with=F][[1]])
valid_patients <- intersect(liv_patients, patients_with_both_deep_and_luminal)
l <- lapply(valid_patients, test_met_depth, si, ad_table, met_group='Liver', ncpus=4)
res <- rbindlist(l)
res$group <- 'Liver'
write_tsv(res,here('processed_data/median_distance_ratios_liv_to_deep_vs_luminal.txt'))


## repeat on the lesion-level
res_pm <- fread(here('processed_data/median_distance_ratios_pm_to_deep_vs_luminal.txt'))
res_pm1 <- fread(here('processed_data/median_distance_ratios_pm_untreatedsync_to_deep_vs_luminal.txt'))
res_pm2 <- fread(here('processed_data/median_distance_ratios_pm_treatedmeta_to_deep_vs_luminal.txt'))
res_pm3 <- fread(here('processed_data/median_distance_ratios_pm_treatedmetaany_to_deep_vs_luminal.txt'))
res_ln <- fread(here('processed_data/median_distance_ratios_ln_to_deep_vs_luminal.txt'))
res_td <- fread(here('processed_data/median_distance_ratios_td_to_deep_vs_luminal.txt'))
#res_lr <- fread(here('processed_data/median_distance_ratios_lr_to_deep_vs_luminal.txt'))
res_liv <- fread(here('processed_data/median_distance_ratios_liv_to_deep_vs_luminal.txt'))
res <- rbind(res_pm, res_pm1, res_pm2, res_ln, res_td, res_liv)
res$group <- factor(res$group, levels=unique(res$group))

plot_list <- lapply(c('Peritoneum (untreated+synchronous)','Peritoneum (treated+metachronous)','Lymph node','Tumor deposit','Liver'), get_plot_for_group, res)
p <- plot_grid(plotlist=plot_list, ncol=1, align='v')
ggsave(here('figures/met_deep_luminal_mediandistance_patientnames.pdf'),width=12, height=12)



res <- rbind(res_pm1, res_pm2, res_ln, res_td, res_liv)
res$group <- factor(res$group, levels=unique(res$group))
pd <- copy(res)
pd <- pd[order(obs,decreasing=F),]
pd$id <- paste(pd$patient, pd$sample)
pd$id <- factor(pd$id, levels=pd$id)
pd$origin <- 'n.s.'
pd[`90%` < 0, origin:='Deep (80%)']
pd[`97.5%` < 0, origin:='Deep (95%)']
pd[`10%` > 0, origin:='Luminal (80%)']
pd[`2.5%` > 0, origin:='Luminal (95%)']
pd$origin <- factor(pd$origin, levels=rev(names(cols)))

p1 <- ggplot(pd, aes(x=id, y=obs)) +
    geom_point(aes(color=origin)) +
    geom_errorbar(aes(min=`2.5%`, max=`97.5%`,color=origin),width=NA) +
    scale_color_manual(values=cols,name='Nearest PT Depth:\n(confidence)') +
    theme_ang(base_size=12) +
    geom_hline(yintercept=0, color='black', size=0.5) +
    facet_wrap(facets=~group, ncol=1, scale='free_x') +
    #theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), legend.position='right') + 
    labs(x=NULL, y='log2 median-distance-ratio (Met:Deep vs Met:Luminal)', subtitle='Met vs PT-depth median-distance-ratios')

x <- as.data.frame.matrix(xtabs(~ group + origin, data=pd))
for(i in 1:nrow(x)) x[i,] <- x[i,] / sum(x[i,])
x$group <- rownames(x)
xm <- reshape2::melt(x, id.var='group')
xm$variable <- factor(xm$variable, levels=rev(names(cols)))
xm$group <- factor(xm$group, levels=unique(res$group))

p2 <- ggplot(xm, aes(x=group, y=value)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,by=0.25)) +
    geom_bar(stat='identity', aes(fill=variable),color='black',linewidth=0.25) +
    scale_fill_manual(values=cols,name='Nearest PT Depth:\n(confidence)') +
    theme_ang(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position='none') + 
    labs(x='Metastasis type', y='Fraction of lesions')

p <- plot_grid(p1, p2, nrow=1, rel_widths=c(2,1))
ggsave(here('figures/met_deep_luminal_mediandistance.pdf'),width=14, height=9)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2024-06-20
# Mouse data
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

#stat.test <- mywilcox2(m, SDI ~ group, paired=F)
tst <- dunn_test_ES(m, SDI ~ group)
p <- ggplot(m, aes(x=group, y=SDI))  +
    #geom_beeswarm(aes(size=region_size_mm2, color=group), pch=16) + 
    geom_point(position=position_jitter(width=0.15, height=0, seed=42), aes(size=region_size_mm2, color=group), pch=16) + 
    geom_boxplot(fill=NA, outlier.shape=NA, width=0.5) + 
    theme_ang(base_size=12) + 
    scale_size(range=c(3,12)) +
    #scale_size_area(max_size=10, breaks=seq(1,13,by=3)) +
    scale_color_manual(values=group_cols, name='Tissue type') +
    stat_compare_means(method = "kruskal.test", label.y = 1.05, geom = "label") +
    stat_pvalue_manual(tst, label='label', y.position=c(0.95,1.0,1.05), tip.length=0) 
ggsave(here('figures/LEGO_boxplot_with_ROI_size.pdf'),width=8, height=6)

m <- fread(here('original_data/misc/SDI_sliding_window_nBins_10_nIterations_10_RGB_sum_min_50_sameRGBexcluded_FALSE.csv'))
m[organ=='CAECUM', group:='Primary']
m[organ=='LIVER', group:='Liver']
m[organ=='PERI', group:='Peritoneum']
m$group <- factor(m$group, levels=c('Primary','Peritoneum','Liver'))
pixel_size <- 1.25e-6 # 1 pixel = 1.25um
pixel_area <- pixel_size^2 # um^2
m[,region_size_mm2:=(n_pixels * pixel_area) * (1e3)^2] # 1m^2 = (1000mm)^2
m$region_size_mm2 <- round(m$region_size_mm2, 3)

collapse <- function(m) {
    size <- mean(m$region_size_mm2)
    n <- nrow(m)
    SDI <- sum(m$SDI * m$region_size_mm2) / sum(m$region_size_mm2)
    #SDI <- median(m$SDI)
    #SDI <- mean(m$SDI)
    list(size=size, SDI=SDI, n=n)
}
m2 <- m[,collapse(.SD),by=c('mouse_ID','group')]
m2[group=='Primary']
tst <- dunn_test_ES(m2, SDI ~ group)
#stat.test <- mywilcox2(m2, SDI ~ group, paired=F)
m2[size < 1, size:=1]

p <- ggplot(m2, aes(x=group, y=SDI))  +
    scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,by=0.25)) + 
    geom_point(position=position_jitter(width=0.15, height=0, seed=42), aes(size=size, color=group), pch=16) + 
    geom_boxplot(fill=NA, outlier.shape=NA, width=0.5) + 
    theme_ang(base_size=12) + 
    scale_size(range=c(3,12)) +
    #scale_size_area(max_size=10, breaks=seq(1,13,by=3)) +
    scale_color_manual(values=group_cols, name='Tissue type') +
    stat_compare_means(method = "kruskal.test", label.y = 1.05, geom = "label") + 
    stat_pvalue_manual(tst, label='label', y.position=c(0.95,1.0,1.05), tip.length=0) 
ggsave(here('figures/LEGO_boxplot_with_ROI_size_avg.pdf'),width=8, height=6)








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scrap? continue from here
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tiles <- pd[,c('id','n_pt','n_dm','n_td'),with=F]
tiles$frac_pt <- tiles$n_pt / max(tiles$n_pt)
tiles$frac_dm <- tiles$n_dm / max(tiles$n_dm)
tiles$frac_td <- tiles$n_td / max(tiles$n_td)
tiles_n <- data.table::melt(tiles[,c('id','n_pt','n_dm','n_td'),with=F], id.var='id')
tiles_n[,variable:=gsub('n_','',variable)]
tiles_frac <- data.table::melt(tiles[,c('id','frac_pt','frac_dm','frac_td'),with=F], id.var='id')
tiles_frac[,variable:=gsub('frac_','',variable)]
tiles <- merge(tiles_n, tiles_frac, by=c('id','variable'), all=T)
names(tiles) <- c('id','type','n','norm')
tiles[type=='pt',type:='Primary']
tiles[type=='td',type:='Tumor deposit']
tiles[type=='dm',type:='Distant']

p2 <- ggplot(tiles, aes(x=id, y=type)) +
    geom_tile(aes(fill=norm)) +
    scale_fill_gradient(low='white',high='steelblue', name='Fraction of row-max') + 
    theme_ang(base_size=12) +
    #theme(axis.text.x=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position='none') +
    labs(x='Patient', y='Number of samples') +


fit <- lm(obs ~ n_td + n_dm + n_pt, data=pd)
coefs <- coef(summary(fit))
ps <- summary(fit)$coef[2:4,4]
est <- coefs[, 1]
upr <- coefs[, 1] + 1.96 * coefs[, 2]
lwr <- coefs[, 1] - 1.96 * coefs[, 2]
info <- as.data.table(cbind(est, lwr, upr, coefs[, 4]))
info <- cbind(type=rownames(coefs), info)
colnames(info) <- c("type","est", "lwr.95ci", "upr.95ci","pval")
info <- info[-1,]
info[type=='n_pt',type:='Primary']
info[type=='n_td',type:='Tumor deposit']
info[type=='n_dm',type:='Distant']

p3 <- ggplot(info, aes(x=type, y=est)) +
    geom_point() +
    geom_errorbar(aes(min=lwr.95ci,max=upr.95ci),width=0.25) +
    geom_hline(yintercept=0, color='red', linetype='dashed', size=0.5) +
    coord_flip() +
    theme_ang(base_size=12) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.line.y=element_blank()) +
    scale_y_continuous(position='right') + labs(x=NULL,y='Effect size') 

p4 <- ggplot() + theme_ang() + theme_nothing()
p <- plot_grid(p1, p4, p2, p3, align='hv', axis='lrtb',ncol=2, nrow=2, rel_heights=c(3,1.5), rel_widths=c(3,1.5))
ggsave(here('figures/fig_5b_td_dm_distinct_origin_metlevel_withIDs.pdf'),width=20,height=11)




