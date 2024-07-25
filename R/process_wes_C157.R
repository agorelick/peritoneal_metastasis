



source(here::here('R/func.R'))
additional_packages <- c('igraph','rjson','reticulate','pheatmap','ggimage','MutationalPatterns','BSgenome','BSgenome.Hsapiens.UCSC.hg19')
suppressMessages(trash <- lapply(additional_packages, require, character.only = TRUE))
ref_genome="BSgenome.Hsapiens.UCSC.hg19"

d <- fread(here('original_data/wes/C157_merged_filtered.ccf.oncokb.maf.gz'))
d[,ShortVariantID:=paste0(Chromosome,':',Reference_Allele,Start_Position,Tumor_Seq_Allele2)]
d[is.na(clonality), clonality:='INDETERMINATE']
length(unique(d$ShortVariantID))
if(!dir.exists(here('figures/wes'))) dir.create(here('figures/wes'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add in substitution type + trinucleotide context
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

convert_to_gr <- function(dat) {
    require(MutationalPatterns)
    dat <- dat[,c('ShortVariantID','Chromosome','Start_Position','Reference_Allele','Tumor_Seq_Allele2'),with=F]
    names(dat) <- c('ShortVariantID','chr','pos','REF','ALT')
    dat <- dat[!duplicated(ShortVariantID),]
    dat$paramRangeID <- NA
    dat <- as.data.frame(dat)
    rownames(dat) <- paste0(dat$chr,':',dat$REF,dat$pos,dat$ALT)
    dat$chr <- paste0('chr',dat$chr)
    g <- makeGRangesFromDataFrame(dat, ignore.strand=T, seqnames.field='chr', start.field='pos', end.field='pos', keep.extra.columns=T)
    GenomeInfoDb::genome(g) = 'hg19'
    g
}

d$pos <- 1:nrow(d)
g <- convert_to_gr(d[!duplicated(ShortVariantID) & !Variant_Type %in% c('INS','DEL')])
tc <- type_context(g, ref_genome, extension=1)
types <- data.table(ShortVariantID=g$ShortVariantID, type=tc$type, context=tc$context)
types[,trinuc:=paste0(substr(context,1,1),'[',type,']',substr(context,3,3))]
types <- types[,c('ShortVariantID','type','trinuc'),with=F]
names(types) <- c('ShortVariantID','muttype6','muttype96')
d <- merge(d, types, by='ShortVariantID', all.x=T)
d <- d[order(pos),]
length(unique(d$ShortVariantID))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first flag any mutations with insufficient detection power
# - after reviewing them, safe to exclude them from all samples 
# (all but 3 are not oncogenic, 3 are possible oncogenic but never clonal)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## estimate the probability of detecting a mutation at that position given 20% CCF, MCN=1, and TCN/purity observed
d[Chromosome %in% 1:22,total_copies:=purity * tcn + (1-purity) * 2]
d[Chromosome %in% c('X','Y'),total_copies:=purity * tcn + (1-purity) * 1] # C157 is male
d[,prob_fragment_ccf0.2_MCN1 := 1 * purity * 0.2 / total_copies] 
d[,power_ccf0.2_MCN1:=1-dbinom(x=0, size=d$t_depth, prob=prob_fragment_ccf0.2_MCN1)]  ## prob of detecting 1+ alt reads
d[power_ccf0.2_MCN1 < 0.9 & t_alt_count==0, t_alt_count:=NA] ## ambiguous sites

## there are 3 potentially oncogenic mutations which have power detectability in a few samples. Are they ever clonal? How many samples are affected?
check_if_ever_clonal <- d[is.na(t_alt_count),c('ShortVariantID','Hugo_Symbol','HGVSp_Short','Variant_Classification','Tumor_Sample_Barcode','ONCOGENIC','MUTATION_EFFECT','GENE_IN_ONCOKB'),with=F]
check_if_ever_clonal[ONCOGENIC=='Likely Oncogenic' | c(GENE_IN_ONCOKB==T & MUTATION_EFFECT %in% c('Loss-of-function','Likely Loss-of-function'))]

## conclude: they are never clonal, they affect a splice region which is low evidence for oncogenicity. Safe to exclude.
d[ShortVariantID=='2:C48040367T',c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','t_alt_count','ccf_1copy')]
d[ShortVariantID=='X:C76876004T',c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','t_alt_count','ccf_1copy')]
d[ShortVariantID=='2:A141032173-',c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','t_alt_count','ccf_1copy')]

bad_variants <- check_if_ever_clonal$ShortVariantID
length(unique(bad_variants))
d <- d[!ShortVariantID %in% bad_variants]
d[Variant_Classification=='Splice_Site', Protein_position:=gsub('p[.]X','',gsub('_splice','',HGVSp_Short))]

## add ABSENT class to clonality
d[t_alt_count==0, clonality:='ABSENT']

## clarify INDETERMINATE clonality
d[clonality=='INDETERMINATE' & t_alt_count > 0 & ccf_expected_copies_upper < 0.9, clonality:='SUBCLONAL']
d[clonality=='INDETERMINATE' & t_alt_count > 0 & ccf_expected_copies_upper >= 0.9, clonality:='CLONAL']


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# next, remove potential FFPE artifacts.
# Here, we find all mutations that are BOTH of:
# 1. C>T
# 2. only detected (subclonally) in 1 sample
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

collapse_mutations <- function(d) {
    affected_clonal <- length(unique(d[t_alt_count >= 1 & clonality=='CLONAL',(Tumor_Sample_Barcode)]))
    affected_notclonal <- length(unique(d[t_alt_count >= 1 & clonality!='CLONAL',(Tumor_Sample_Barcode)]))
    list(clonal=affected_clonal, notclonal=affected_notclonal)
}
ever_clonal <- d[,collapse_mutations(.SD), by=c('ShortVariantID','muttype6','ONCOGENIC','MUTATION_EFFECT')]
bad_variants <- ever_clonal[clonal==0 & notclonal==1 & muttype6=='C>T',(ShortVariantID)]
length(unique(bad_variants))  # 57
d2 <- d[!ShortVariantID %in% bad_variants]
length(unique(d2$ShortVariantID))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subset for valid mutation types 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

valid_classes <- c('Missense_Mutation','Nonsense_Mutation','Splice_Site','Translation_Start_Site','Silent','In_Frame_Ins','In_Frame_Del','Frame_Shift_Ins','Frame_Shift_Del')
length(unique(d2[!Variant_Classification %in% valid_classes,(ShortVariantID)]))
d3 <- d2[Variant_Classification %in% valid_classes,]
length(unique(d3$ShortVariantID))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add clonality classification
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parse <- function(aapos) {
    x <- strsplit(aapos,'[/]')[[1]][1]
    strsplit(x,'[-]')[[1]][1]
}
d3$Amino_acid_position <- as.integer(sapply(d3$Protein_position, parse))
d3[Variant_Classification %in% valid_classes, tm:=paste(Hugo_Symbol, gsub('p[.]','',HGVSp_Short))]
d3[tm %in% c('PHF1 R58=','SACS L3700='), tm:=paste0(tm,' (',HGVSc,')')] ## a few mutations occur in separate cDNA changes

## add in sample names/groups
map <- fread(here('original_data/wes/sample_map.txt'))
map$Sample_ID <- gsub('C157','',map$Sample_ID)
d3 <- merge(d3, map, by.x='Tumor_Sample_Barcode', by.y='Sample_ID', all.x=T)

## these mutations have NA copy number data in some sample
bad_variants <- unique(d3[Real_Sample_ID!='LN2b' & (is.na(tcn) | tcn==0 | is.na(lcn)),(ShortVariantID)])
length(bad_variants) 

qc <- d3[ShortVariantID %in% bad_variants & (ONCOGENIC %in% c('Likely Oncogenic','Oncogenic') | c(GENE_IN_ONCOKB==T & MUTATION_EFFECT %in% c('Loss-of-function','Likely Loss-of-function')))]
maf <- d3[!ShortVariantID %in% bad_variants,]

## standardize data when 0 ALT reads or TCN==0
maf[t_alt_count==0 | tcn==0, c('ccf_expected_copies','clonality'):=list(0,'ABSENT')]

## make classifications
maf[clonality=='ABSENT', clonality_3cati:=0]
maf[clonality=='SUBCLONAL', clonality_3cati:=1]
maf[clonality=='CLONAL', clonality_3cati:=2]

## save CCF distance matrix
dat <- data.table::dcast(ShortVariantID ~ Real_Sample_ID, value.var='ccf_expected_copies', data=maf)
dat$Normal1 <- 0
dat <- d2m(dat)
dm <- as.matrix(dist(t(dat), method='euclidean'))
write_distance_matrix(dm, here('processed_data/wes/ccf_distance_matrix_allsamples.txt'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make CCF heatmap and NJ tree
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## heatmap showing the clonal status of mutations that are EVER clonal in any sample
## note: we remove sample LN2b as it has too many mutations and is likely artifactual
summarize <- function(maf) {
    clonal=length(unique(maf$ShortVariantID[maf$clonality=='CLONAL']))
    subclonal=length(unique(maf$ShortVariantID[maf$clonality=='SUBCLONAL']))
    list(clonal=clonal, subclonal=subclonal)
}
n <- maf[clonality!='ABSENT',summarize(.SD), by=Real_Sample_ID]
n[,total:=clonal+subclonal]
n <- n[order(total,decreasing=T),]
n$x_median_total <- n$total / median(n$total) 
n$x_median_clonal <- n$clonal / median(n$clonal)
n$x_median_subclonal <- n$subclonal / median(n$subclonal)

maf <- maf[Real_Sample_ID!='LN2b', ]
tmp <- maf[,c('Real_Sample_ID','tm','ShortVariantID','Hugo_Symbol','HGVSp_Short','Variant_Classification','ONCOGENIC','ccf_expected_copies','ccf_expected_copies_lower','clonality')]
tmp$ShortVariantID %>% unique %>% length
ever_clonal <- unique(tmp[clonality=='CLONAL',(tm)])
tmp <- tmp[tm %in% ever_clonal,]
tmp$group <- 'Other'
tmp[grepl('Per',Real_Sample_ID), group:='Peritoneum']
ever_clonal <- unique(tmp[clonality=='CLONAL',(tm)])
tmp <- tmp[tm %in% ever_clonal,]
ccf_matrix <- dcast( tm ~ Real_Sample_ID, value.var='ccf_expected_copies', data=tmp)
ccf_matrix <- d2m(ccf_matrix)
valid_mutations <- rownames(ccf_matrix)
write_distance_matrix(ccf_matrix, here('processed_data/wes/sample_mutation_ccf_matrix.txt'))


prep_data_for_pyclone <- function(maf, ccf_matrix) {
    out <- maf[,c('Real_Sample_ID','tm','t_ref_count','t_alt_count','Chromosome','tcn','lcn','purity'),with=F]
    out[Chromosome %in% 1:22, normal_cn:=2]
    out[Chromosome %in% c('X','Y'), normal_cn:=1]

    ## merge the CCF data to the output from the MAF
    pd <- as.data.table(reshape2::melt(ccf_matrix))
    names(pd) <- c('tm','Real_Sample_ID','ccf')
    out <- merge(pd, out, by=c('tm','Real_Sample_ID'), all.x=T)
    out <- out[Real_Sample_ID!='Normal1',]
    out[,major_cn:=tcn - lcn]
    out[,minor_cn:=lcn]
    out[,c('ccf','Chromosome','tcn','lcn'):=NULL]
    out <- out[,c('tm','Real_Sample_ID','t_ref_count','t_alt_count','normal_cn','major_cn','minor_cn','purity'),with=F]
    names(out) <- c('mutation_id','sample_id','ref_counts','alt_counts','normal_cn','major_cn','minor_cn','tumour_content')
    out
}


plot_pyclone_clusters <- function(maf, file_pyclone_output, title) { 
    ## load output from pyclone to visualize each cluster's CCF in each sample
    pyclone_output <- fread(file_pyclone_output)
    pyclone_output[,ccf:=100*cellular_prevalence]
    pyclone_output[,cluster_id:=cluster_id+1]
    setnames(pyclone_output,'cluster_id','cluster')
    clusters <- sort(unique(pyclone_output$cluster))
    cols <- brewer.pal(length(clusters), 'Paired'); names(cols) <- clusters

    ## similar plot but with the raw CCFs
    pyclone_output <- merge(pyclone_output, maf[,c('Real_Sample_ID','tm','ccf_expected_copies'),with=F], by.x=c('sample_id','mutation_id'), by.y=c('Real_Sample_ID','tm'), all.x=T)
    pyclone_output$cluster <- factor(pyclone_output$cluster, levels=sort(unique(pyclone_output$cluster)))
    pyclone_output$garbage <- 'No'
    pyclone_output$garbage[pyclone_output$mutation_id %in% garbage_muts] <- 'Yes'

    ggplot(pyclone_output, aes(x=cluster, y=100*ccf_expected_copies)) +
        geom_point(data=pyclone_output, aes(fill=cluster), color='white', size=2.5, pch=21, stroke=0.25, position=position_jitter(width=0.333,height=0, seed=42)) +
        geom_boxplot(aes(fill=cluster), color='black', alpha=0.3, width=0.8, linewidth=0.25, outlier.shape=NA) +
        facet_wrap(facets=~sample_id) + 
        theme_fit(base_size=12) +
        theme(panel.grid.major=element_line(linewidth=0.25, color='grey92')) +
        scale_fill_manual(values=cols, name='Cluster') +
        labs(x='Cluster', y='Cancer cell fraction (%)', title=title)
}


prep_data_for_pairtree <- function(maf, ccf_matrix, file_from_pyclone=NA) { 
    dat <- prep_data_for_pyclone(maf, ccf_matrix)

    ## prepare data for pairtree
    mutation_prob <- maf[,c('tm','Chromosome','Start_Position','Real_Sample_ID','tcn','expected_alt_copies'),with=F]
    mutation_prob[,prob:=round(expected_alt_copies / tcn,6)]
    mutation_prob[prob > 1, prob:=1]
    mutation_prob[,c('tcn','expected_alt_copies'):=NULL]

    ## load data submitted to pyclone and annotate with mutation prob
    dat <- merge(dat, mutation_prob, by.x=c('mutation_id','sample_id'), by.y=c('tm','Real_Sample_ID'), all.x=T)
    dat$sample_id <- factor(dat$sample_id)
    dat[,dp:=ref_counts + alt_counts]

    ## order the mutations genomically
    dat$Chromosome <- factor(dat$Chromosome, levels=c(1:22,'X','Y'))
    dat <- dat[order(Chromosome, Start_Position, sample_id),]
    dat$mutation_id <- factor(dat$mutation_id, levels=unique(dat$mutation_id))

    ## generate the .ssm file for pairtree
    collapse <- function(field, dat) {
        dd <- data.table::dcast(mutation_id ~ sample_id, value.var=field, data=dat)
        md <- d2m(dd)
        value <- apply(md, 1, paste, collapse=',')
        out <- data.table(name=unique(dat$mutation_id), value=value)
        names(out)[2] <- field
        out
    }
    x_alt <- collapse('alt_counts', dat)
    x_dp <- collapse('dp', dat)
    x_prob <- collapse('prob', dat)
    ssm <- merge(x_alt, x_dp, by='name', all=T)
    ssm <- merge(ssm, x_prob, by='name', all=T)
    ssm <- cbind(id=paste0('s',0:(nrow(ssm)-1)), ssm)
    names(ssm) <- c('id','name','var_reads','total_reads','var_read_prob')

    if(is.na(file_from_pyclone)) {
        samples <- paste0('{"samples": [',paste(paste0('"',levels(dat$sample_id),'"'), collapse=', '),']')
        params <- paste0(samples, ', "clusters": [], "garbage": []}')

    } else {
        ## generate the params.json (fake json)
        pyclone_output <- fread(file_from_pyclone)
        pyclone_output <- pyclone_output[!duplicated(mutation_id),c('mutation_id','cluster_id'),with=F]
        pyclone_output[,cluster_id:=cluster_id+1]
        tmp <- ssm[,c('id','name'),with=F]
        tmp$pos <- 1:nrow(tmp)
        tmp <- merge(tmp, pyclone_output, by.x='name', by.y='mutation_id', all.x=T)
        tmp <- tmp[order(pos),]
        tmp$cluster_id <- factor(tmp$cluster_id, levels=sort(unique(tmp$cluster_id)))
        clusterdat <- copy(tmp)
        extract_variants_per_cluster <- function(tmp) {
            clusters <- paste0('[',paste(paste0('"',tmp$id,'"'), collapse=', '),']')
            list(ids=clusters)
        }
        clusters <- tmp[,extract_variants_per_cluster(.SD), by=cluster_id]
        clusters <- clusters[order(cluster_id),]
        clusters <- paste0(' "clusters": [',paste(clusters$ids, collapse=', '),']')
        samples <- paste0('{"samples": [',paste(paste0('"',levels(dat$sample_id),'"'), collapse=', '),']')
        params <- paste0(samples, ',', clusters, ', "garbage": []}')
    }
    
    list(ssm=ssm, params=params)
}


## before any clustering, use pairtree's removegarbage utility to flag bad mutations
pt <- prep_data_for_pairtree(maf, ccf_matrix)
write_tsv(pt$ssm,here('processed_data/wes/data_for_pairtree.ssm'))
cat(pt$params, file=here('processed_data/wes/data_for_pairtree.json'))

#pairtree_output_garbage_flagged.json

#### flag garbage
# pairtree/bin/removegarbage processed_data/wes/data_for_pairtree.ssm processed_data/wes/data_for_pairtree.json processed_data/wes/pairtree_output_garbage_flagged.json --seed 123


## create data for pyclone-vi clustering AFTER removing the garbage mutations
j <- rjson::fromJSON(file=here('processed_data/wes/pairtree_output_garbage_flagged.json'))
garbage <- j$garbage
tmp <- pt$ssm[,c('name','id'),with=F]
garbage_muts <- as.character(tmp$name[tmp$id %in% garbage])
maf_nogarbage <- maf[!tm %in% garbage_muts,]
ccf_matrix_nogarbage <- ccf_matrix[!rownames(ccf_matrix) %in% garbage_muts,]
data_for_pyclone <- prep_data_for_pyclone(maf_nogarbage, ccf_matrix_nogarbage)
write_tsv(data_for_pyclone, here('processed_data/wes/data_for_pyclone.tsv'))


## run pyclone
# pyclone-vi fit -i processed_data/wes/data_for_pyclone.tsv -o processed_data/wes/pyclone_output.h5 -c 20 -d beta-binomial -r 100 --seed 123
# pyclone-vi write-results-file -i processed_data/wes/pyclone_output.h5 -o processed_data/wes/pyclone_output.tsv


## plot the pyclone clustering results (boxplot of mutation raw CCFs per cluster, sample)
file_pyclone_output <- here('processed_data/wes/pyclone_output.tsv')
p <- plot_pyclone_clusters(maf_nogarbage, file_pyclone_output, title='Pyclone-vi clusters (clustering after garbage mutations are removed)')
ggsave(here('figures/wes/pyclone_raw_ccfs.pdf'),width=11, height=8)


## for each cluster, show the corresponding CCF heatmap
x <- fread(file_pyclone_output)
x <- merge(x, maf[,c('Real_Sample_ID','tm','ccf_expected_copies'),with=F], by.x=c('sample_id','mutation_id'), by.y=c('Real_Sample_ID','tm'), all.x=T)
setnames(x,'cluster_id','cluster')
x[,cluster:=cluster+1]
x$cluster <- factor(x$cluster, levels=sort(unique(x$cluster)))
cluster_heatmap <- function(k, x) {
    qc <- data.table::dcast(mutation_id ~ sample_id, value.var='ccf_expected_copies', data=x[cluster==k,])
    qc <- d2m(qc)
    pdf(here(paste0('figures/wes/cluster',k,'_heatmap.pdf')), width=8, height=8)
    pheatmap(t(qc), cluster_rows=F)
    dev.off()
}
lapply(1:7, cluster_heatmap, x)


## manually edit the clusters from pyclone to split up the high-variance clusters
clusters <- fread(file_pyclone_output)
clusters <- clusters[!mutation_id %in% c('CCDC108 M633*','NCOR2 I976Sfs*87'),] # remove cluster 2 entirely (garbage based on heatmap)
clusters <- clusters[!mutation_id %in% c('PYHIN1 R380=','TEC X550_splice'),] # remove these mutations which don't fit into cluster 7 and are n=1 clusters

## create this new cluster from cluster 4
clusters$cluster_id[clusters$mutation_id %in% c('C1orf74 Q131*','TIGD5 G339=','MDM1 P84R','LRRC39 Y309C')] <- max(clusters$cluster_id) + 1

## create this new cluster from cluster 4
clusters$cluster_id[clusters$mutation_id %in% c('ANKRD30A S406R','AK7 A159=','LRIG1 R574H')] <- max(clusters$cluster_id) + 1

## create this new cluster from cluster 5
clusters$cluster_id[clusters$mutation_id %in% c('MAEL V32I','ADAMTSL4 S68R')] <- max(clusters$cluster_id) + 1

## create this new cluster from cluster 7
clusters$cluster_id[clusters$mutation_id %in% c('LTBP1 P943=','DLGAP2 R669Q','OR51S1 G96S')] <- max(clusters$cluster_id) + 1

clusters$cluster_id <- factor(clusters$cluster_id, levels=unique(sort(clusters$cluster_id)))
clusters$cluster_id <- as.integer(clusters$cluster_id)
clusters$cluster_id <- clusters$cluster_id - 1
clusters <- clusters[order(cluster_id, mutation_id, sample_id),]
write_tsv(clusters, here('processed_data/wes/pyclone_output_reclustered.tsv'))


## create heatmaps for each cluster again after using the new clusters
x <- merge(clusters, maf[,c('Real_Sample_ID','tm','ccf_expected_copies'),with=F], by.x=c('sample_id','mutation_id'), by.y=c('Real_Sample_ID','tm'), all.x=T)
setnames(x,'cluster_id','cluster')
x[,cluster:=cluster+1]
x$cluster <- factor(x$cluster, levels=sort(unique(x$cluster)))

## for each cluster, show the corresponding CCF heatmap
cluster_heatmap <- function(k, x) {
    qc <- data.table::dcast(mutation_id ~ sample_id, value.var='ccf_expected_copies', data=x[cluster==k,])
    qc <- d2m(qc)
    pdf(here(paste0('figures/wes/reclustered_',k,'_heatmap.pdf')), width=8, height=8)
    pheatmap(t(qc), cluster_rows=F)
    dev.off()
}
lapply(1:10, cluster_heatmap, x)


## plot CCF distributions for the new clusters
p <- plot_pyclone_clusters(maf, here('processed_data/wes/pyclone_output_reclustered.tsv'), title='C157 mutation clusters after manual review') 
ggsave(here('figures/wes/pyclone_raw_ccfs_reclustered.pdf'),width=11, height=8)


## create data for orchard/pairtree using the new clusters
maf_reclustered <- maf[tm %in% clusters$mutation_id,]
ccf_matrix_reclustered <- ccf_matrix[unique(clusters$mutation_id),]
write_distance_matrix(ccf_matrix_reclustered, here('processed_data/wes/sample_mutation_ccf_matrix_reclustered.txt'))
pt <- prep_data_for_pairtree(maf_reclustered, ccf_matrix_reclustered, file_from_pyclone=here('processed_data/wes/pyclone_output_reclustered.tsv'))
write_tsv(pt$ssm,here('processed_data/wes/data_for_pairtree_reclustered.ssm'))
cat(pt$params, file=here('processed_data/wes/data_for_pairtree_reclustered.json'))


# run orchard
# python3 orchard/bin/orchard processed_data/wes/data_for_pairtree_reclustered.ssm processed_data/wes/data_for_pairtree_reclustered.json processed_data/wes/orchard_output.npz -p --seed 123 

# run pairtree to plot tree
# pairtree/bin/plottree --runid C157 processed_data/wes/data_for_pairtree_reclustered.ssm processed_data/wes/data_for_pairtree_reclustered.json processed_data/wes/orchard_output.npz processed_data/wes/orchard_output_results.html --tree-json processed_data/wes/orchard_output_clonetree.json --omit-plots pairwise_mle,pairwise_separate


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot clone tree and extract the pairtree data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## make a final heatmap using the set of variants retained in clusters and used in the clone tree
res <- fread(here('processed_data/wes/pyclone_output_reclustered.tsv'))
res[,cluster_id:=cluster_id+1]
res <- res[,c('sample_id','mutation_id','cluster_id','cellular_prevalence'),with=F]
setnames(res,'mutation_id','name')
ssm <- fread(here('processed_data/wes/data_for_pairtree_reclustered.ssm'))
ssm <- ssm[,c('id','name'),with=F]
clusters <- merge(res, ssm, by='name', all.y=T)

tmp <- fread(here('processed_data/wes/data_for_pairtree.ssm'))
retained_variants <- unique(tmp$name)
ccf <- fread(here('processed_data/wes/C157_merged_filtered.ccf.oncokb.cleaned.maf'),select=c('tm','Real_Sample_ID','ccf_expected_copies_lower','ShortVariantID'))
ccf <- ccf[Real_Sample_ID!='LN2b',]
ccf <- ccf[tm %in% retained_variants,]
ccf[ShortVariantID=='2:G54093907A', tm:='PSME4 L1792= (1)']
ccf[ShortVariantID=='2:T54093905C', tm:='PSME4 L1792= (2)']
clusters <- merge(clusters, ccf, by.x=c('sample_id','name'), by.y=c('Real_Sample_ID','tm'), all.y=T)
clusters <- clusters[!is.na(cluster_id)] # 4 mutations that were dropped from the clustered data manually are still here, only 1 entry/mutation, with cluster=NA
ccf_matrix <- data.table::dcast(cluster_id + name ~ sample_id, value.var='ccf_expected_copies_lower', data=clusters)
ccf_matrix$cluster_id <- factor(ccf_matrix$cluster_id, levels=sort(unique(ccf_matrix$cluster_id)))

## counts the number of variants per cluster
cluster_muts <- table(ccf_matrix$cluster_id)

myclusters <- sort(unique(clusters$cluster_id))
nc <- length(myclusters)
cluster_cols <- c('black','darkgrey',brewer.pal(nc, 'Paired')); names(cluster_cols) <- c(-1,0,1:nc)
my_colors <- list(Cluster=cluster_cols)
mat <- t(as.matrix(ccf_matrix[,(3:ncol(ccf_matrix)),with=F]))
colnames(mat) <- ccf_matrix$name
cols <- data.frame(Cluster=ccf_matrix$cluster_id)
row.names(cols) <- colnames(mat)
heatmap_info <- pheatmap(mat, annotation_col = cols, annotation_colors = my_colors, filename=here('figures/ed_fig_1a_center.pdf'), width=20, height=10, main='ED Fig 1a, center')
sample_levels <- heatmap_info$tree_row$labels[heatmap_info$tree_row$order]

tmp <- rbind(mat, Normal1=0)
dm <- dist(tmp, method='euclidean')
tree <- nj(dm)
tree <- phytools::reroot(tree, which(tree$tip.label=='Normal1'))

## show the final clone tree, annotate the edges with the number of mutations
j <- rjson::fromJSON(file=here('processed_data/wes/orchard_output_clonetree.json'))
nc <- length(j$eta)-1
cols <- c('darkgrey',brewer.pal(nc, 'Paired')); names(cols) <- c(0,1:nc)
parents <- j$parents
edges <- data.table(parent=parents, child=1:length(parents))
graph <- data.tree::FromDataFrameNetwork(edges, c("parent", "child"))
ig <- data.tree::as.igraph.Node(graph, "p", c("level", "isLeaf"), directed=T, direction='climb')
md <- data.table(cluster=names(V(ig)))
md$pos <- 1:nrow(md)
md$color <- cols[md$cluster]
V(ig)$color <- md$color

pdf(here('figures/ed_fig_1b.pdf'))
plot(ig, layout = layout.reingold.tilford(ig, root=1), main='ED Fig 1b')
dev.off()

np <- reticulate::import("numpy")
npz2 <- np$load(here("processed_data/wes/orchard_output.npz"))
eta <- npz2$f[['eta']]
eta <- eta[1,,]
samplenames <- as.character(npz2$f[['sampnames.json']])
samplenames <- gsub('"','',as.character(samplenames))
samplenames <- trimws(samplenames, which = "both", whitespace = "[\\[\\]\n]")
samplenames <- strsplit(samplenames, ', ')[[1]]
colnames(eta) <- samplenames
rownames(eta) <- (1:nrow(eta))-1
eta <- eta[-1,]
for(i in 1:ncol(eta)) eta[,i] <- eta[,i] / sum(eta[,i])

d_eta <- as.data.table(reshape2::melt(eta))
names(d_eta) <- c('cluster','sample','prop')
d_eta$sample <- factor(d_eta$sample, levels=rev(sample_levels))
d_eta$cluster <- factor(d_eta$cluster, levels=1:nc)
write_tsv(d_eta, here('processed_data/wes/clone_proportions.tsv'))
p <- ggplot(data=d_eta, aes(x=sample, y=prop)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0,1,by=0.25)) + 
    geom_bar(stat='identity', aes(fill=cluster)) + 
    scale_fill_manual(values=cluster_cols, name='Clone') +
    coord_flip() +
    theme_ang(base_size=12) +
    ggtitle('ED Fig 1a, right')
ggsave(here('figures/ed_fig_1a_right.pdf'),width=4,height=10)


## add in sample names/groups
group_cols <- c("#000000","#008C45","#EB5B2B","#FAB31D","#FAB31D","#4C86C6","#4C86C6","#bfbfbf")
names(group_cols) <- c('Normal','Primary','Locoregional','Peritoneum','Lung','Liver','Distant (other)','Other')

groups <- map[,c('Real_Sample_ID','group'),with=F]
setnames(groups,'Real_Sample_ID','label')
p0 <- ggtree(tree, layout='rect') + theme_tree2()
p0 <- p0 %<+% groups
p0 <- p0 + scale_color_manual(values=group_cols,name='Tissue type') 
p0 <- p0 + labs(title='ED Fig 1c')
p0 <- p0 + theme(legend.position='none')

popfreq <- t(eta)
nodes <- as.data.table(p0$data)
nodes <- nodes[isTip==T,c('node','label')]
stats <- cbind(label=rownames(popfreq), as.data.table(popfreq))
stats <- merge(nodes, stats, by='label', all.x=T)
stats <- stats[label!='Normal1',]
stats[,label:=NULL]
pies <- nodepie(stats, cols = 2:ncol(stats))
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cluster_cols))

p <- p0 + geom_inset(pies, width = .1, height = .1)
p <- p + geom_tiplab(aes(color=group), angle=0, hjust=-1.5) 
ggsave(here('figures/ed_fig_1c.pdf'),width=10,height=6)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test tree similarity between WES (original CCF), Poly-G, and SCNA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2024-07-24
dm_wes <- read_distance_matrix(here('processed_data/wes/ccf_distance_matrix_allsamples.txt'))
dm_scna <- read_distance_matrix(here('processed_data/copynumber/C157/C157_cnv_distance_matrix.txt'))
dm_polyg <- read_distance_matrix(here('processed_data/angular_distance_matrices/C157.txt'))

## plot 3 C157 trees unrooted
groups <- fread(here('original_data/wes/sample_map.txt'), select=c('Real_Sample_ID','group'))
names(groups)[1] <- 'label'
group_cols <- c("#000000","#008C45","#EB5B2B","#FAB31D","#FAB31D","#4C86C6","#4C86C6","#bfbfbf")
names(group_cols) <- c('Normal','Primary','Locoregional','Peritoneum','Lung','Liver','Distant (other)','Other')

tree0 <- nj(dm_wes)
tree0 <- root(tree0, outgroup='Normal1')
tree0 <- ape::rotate(tree0, 18)
tree0 <- ape::rotate(tree0, 22)

tree0$edge.length[1] <- 6 ## scale the LN2b branch
p0 <- ggtree(tree0, layout='ape')
p0 <- p0 %<+% groups
p0 <- p0 + scale_color_manual(values=group_cols,name='Tissue type') 
p0 <- p0 + labs(subtitle='WES', title='ED Fig 1d')
p0 <- p0 + theme(legend.position='none')
p0 <- p0 + geom_tiplab(aes(color=group), angle=0, size=3.5) 

tree1 <- nj(dm_polyg)
tree1 <- phytools::reroot(tree1, which(tree1$tip.label=='Normal1'))
p1 <- ggtree(tree1, layout='ape')
p1 <- p1 %<+% groups
p1 <- p1 + scale_color_manual(values=group_cols,name='Tissue type') 
p1 <- p1 + labs(subtitle='Poly-G', title='')
p1 <- p1 + theme(legend.position='none')
p1 <- p1 + geom_tiplab(aes(color=group), angle=1, size=3.5) 

tree2 <- nj(dm_scna)
tree2 <- phytools::reroot(tree2, which(tree2$tip.label=='Normal1'))
tree2 <- ape::rotate(tree2, 26)
tree2 <- ape::rotate(tree2, 19)

p2 <- ggtree(tree2, layout='ape')
p2 <- p2 %<+% groups
p2 <- p2 + scale_color_manual(values=group_cols,name='Tissue type') 
p2 <- p2 + labs(subtitle='SCNA', title='')
p2 <- p2 + theme(legend.position='none')
p2 <- p2 + geom_tiplab(aes(color=group), angle=1, size=3.5) 

p <- plot_grid(p0, p1, p2, nrow=1)
ggsave(here('figures/ed_fig_1d.pdf'),width=9, height=6)



