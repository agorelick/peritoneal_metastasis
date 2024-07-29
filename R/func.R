
ncpus = 4
options(repr.plot.width=7, repr.plot.height=7)


# ~~~~~~~~~~~~~~~~~~~~~~~~~
# - load necessary libraries
# - create expected directories for output
# ~~~~~~~~~~~~~~~~~~~~~~~~~

message('Loading required libraries ...')
packages <- c('ACE','ape','binom','caper','clipr','coin','colorspace','cowplot','data.table','dendextend','DescTools','dplyr','ggplot2','ggpubr','ggrepel','ggsignif','gmp','gplots','here','pals','parallel','patchwork','phangorn','phytools','QDNAseq','Quartet','RColorBrewer','readxl','reshape2','rstatix','tidyverse','TreeDist','TreeTools','viridis', 'ACE','adephylo','Biobase','GenomicRanges','ggtree','phyloseq', 'ggbeeswarm', 'rds')
suppressMessages(trash <- lapply(packages, require, character.only = TRUE))

f=function(directory) {
    if(!dir.exists(directory)) {
        message('Creating directory ',directory)
        dir.create(directory, recursive=T)
    }
}
dirs <- c(here('processed_data/copynumber'),
          here('figures_and_tables/copynumber/tree_comparisons'),
          here('figures_and_tables/copynumber/distance_matrix_comparisons'),
          here('figures_and_tables/copynumber/bootstrapped_trees'),
          here('figures_and_tables/copynumber/heatmaps'))
trash <- lapply(dirs, f)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define functions used throughout
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

theme_ang <- function (base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_bw(base_size = base_size, base_line_size = base_line_size,
        base_rect_size = base_rect_size) %+replace% theme(line = element_line(colour = "black",
        linewidth = base_line_size, linetype = 1, lineend = "round"),
        text = element_text(colour = "black", size = base_size,
            lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0,
            margin = margin(), debug = F), axis.text = element_text(colour = "black",
            size = rel(0.8)), axis.ticks = element_line(colour = "black",
            linewidth = rel(1)), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",
            linewidth = rel(1)), legend.key = element_blank(), strip.background = element_blank())
}


theme_fit <- function (base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    theme_bw(base_size = base_size, base_line_size = base_line_size,
        base_rect_size = base_rect_size) %+replace% theme(line = element_line(colour = "black", linewidth = base_line_size, linetype = 1, lineend = "round"),
        text = element_text(colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = F), 
        axis.text = element_text(colour = "black", size = rel(0.8)), 
        axis.ticks = element_line(colour = "black", linewidth = rel(1)), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        legend.key = element_blank()) 
} 


read_distance_matrix <- function (file, return.as.matrix = T) {
    distance_matrix <- fread(file)
    rows <- distance_matrix[[1]]
    distance_matrix <- distance_matrix[, (2:ncol(distance_matrix)), with = F]
    m <- as.matrix(distance_matrix)
    rownames(m) <- rows
    if (return.as.matrix == F) {
        as.dist(m, diag = T)
    }
    else {
        m
    }
}


write_distance_matrix <- function (dm_df, filepath) {
    write.table(dm_df, file = filepath, sep = "\t", quote = FALSE, col.names = NA)
}


write_tsv <- function (d, file, sep = "\t", quote = F, row.names = F, ...) { 
    write.table(d, file = file, sep = sep, quote = quote, row.names = row.names, ...)
}


table_freq <- function (value) {
    if (is.null(value) | length(value) == 0) {
        tbl <- data.table(value = NA, N = NA)
    }
    else {
        tbl <- as.data.table(table(value))
        tbl <- tbl[order(tbl$N, decreasing = T), ]
    }
    tbl
}


d2m <- function(dt) {
    ## assumes first column should be rownames of a matrix made from columns 2:ncol
    rows <- dt[[1]]
    if('data.table' %in% class(dt)) {
        dt <- dt[,c(2:ncol(dt)),with=F]
    } else if(class(dt)=='data.frame') {
        ## assume this is a data.frame
        dt <- dt[,c(2:ncol(dt))]
    } else {
        stop('enter a data.table or a data.frame')
    }
    m <- as.matrix(dt)
    rownames(m) <- rows
    m
}


parse_barcode <- function(barcodes) {
    ## extract the main tissue type, lesion, and sample number from each sample's barcode
    .parse_barcode <- function(barcode) {
        str <- strsplit(gsub("([A-Za-z]*)([0-9]*)([A-Za-z]*)", "\\1 \\2 \\3", barcode), " ")[[1]]
        list(barcode=barcode,type=str[1],lesion=str[2],sample=str[3])
    }
    s <- rbindlist(lapply(barcodes, .parse_barcode))
    s[is.na(sample),sample:='']
    s[grepl('-A$',barcode),autopsy:=T]
    s[!grepl('-A$',barcode),autopsy:=F]
    s$sample <- gsub('-A$','',s$sample)
    s
}


group_samples <- function(input,lun=T,liv=T,per=T,primary_autopsy_is_distant=T,highlight_peritoneum_when_multi=T,color=F) {
    if(any(c('data.frame','matrix') %in% class(input))) {
        barcodes <- rownames(input)
    } else if('character' %in% class(input)) {
        barcodes <- input
    }
    info <- parse_barcode(barcodes)    
    info[grepl('^N[0-9]',barcode) | grepl('^Normal[0-9]',barcode),group:='Normal']
    info[grepl('^P[0-9]',barcode) | grepl('^PT[0-9]',barcode),group:='Primary']
    info[grepl('^L[0-9]',barcode) | grepl('^LN[0-9]',barcode) | grepl('^TD[0-9]',barcode),group:='Locoregional']
    info[grepl('^Lun[0-9]',barcode),group:='Lung']
    info[grepl('^Liv[0-9]',barcode),group:='Liver']
    info[grepl('^Per[0-9]',barcode) | grepl('^Di[0-9]',barcode) | grepl('^Om[0-9]',barcode) | 
         grepl('^PerOv[0-9]',barcode),group:='Peritoneum']
    info[is.na(group),group:='Distant (other)']

    if(lun==F) info[group=='Lung',group:='Distant (other)']
    if(liv==F) info[group=='Liver',group:='Distant (other)']
    if(per==F) info[group=='Peritoneum',group:='Distant (other)']
    if(primary_autopsy_is_distant==T) info[group=='Primary' & autopsy==T,group:='Distant (other)']
    out <- info[,c('barcode','group'),with=F]
   
    if(highlight_peritoneum_when_multi) {
        ## default coloring showing all major types, but highlighting peritoneum
        out[group=='Normal',color:='black']
        out[group=='Primary',color:='#008c45']
        out[group=='Locoregional',color:='#eb5b2b']
        out[group=='Liver',color:='#4c86c6']
        out[group=='Lung',color:='#ea6a8c']
        out[group=='Peritoneum',color:='#fab31d']
        out[group=='Distant (other)',color:='#534797']    
    } else {
        ## default coloring showing all major types, but highlighting lung
        out[group=='Normal',color:='black']
        out[group=='Primary',color:='#008c45']
        out[group=='Locoregional',color:='#eb5b2b']
        out[group=='Liver',color:='#4c86c6']
        out[group=='Peritoneum',color:='#ea6a8c'] ## find alternative?
        out[group=='Lung',color:='#fab31d']
        out[group=='Distant (other)',color:='#534797']    
    }

    ## depending on which were included in input argumens, highlight single type
    if(lun==T & liv==F & per==F) {
        out[group=='Lung',color:='#fab31d']
        out[group %in% 'Distant (other)',color:='#4c86c6']

    } else if(lun==F & liv==T & per==F) {
        out[group=='Liver',color:='#fab31d']
        out[group %in% 'Distant (other)',color:='#4c86c6']

    } else if(lun==F & liv==F & per==T) {
        out[group=='Peritoneum',color:='#fab31d']
        out[group %in% 'Distant (other)',color:='#4c86c6']
    }
    if(color==F) out[,color:=NULL]
    out 
}


extract_tree <- function(tmp,type) {
    if(!is.null(tmp) & !is.null(tmp$plot)) {
        p <- tmp$plot
        p <- p + ggplot2::ggtitle(paste(tmp$patient, type))
        p
    } else {
        NULL
    }
}   


extract_data <- function(tmp) {
    patient <- tmp$patient
    dat <- tmp$data
    dat$patient <- patient
    dat
}   


extract_gglegend <- function (p) {
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) > 0)
        leg <- tmp$grobs[[leg]]
    else leg <- NULL
    leg
    legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
    plot <- p + theme(legend.position = "none")
    list(plot = plot, legend = legend)
}

get_met_specific_distances <- function(si, ad_table, comparison, distance, return_tree, met_color='blue') { 
    if(!comparison %in% c('normal','primary','pairwise','rds','all')) 
        stop("comparison should be one of 'normal', 'primary', 'pairwise', 'rds', or 'all'")
    if(!distance %in% c('node','angular','rds','scna','patristic')) 
        stop("distance should be one of 'node', 'angular', 'rds', 'scna', 'patristic'")
    if(comparison=='normal' & distance %in% c('patristic','angular')) 
        stop("angular distance/patritic node distabce should not be used for distance to normal")

    get_mean_distance_to_root <- function(distmatrix) {
        root <- grep('^N',rownames(distmatrix))
        distances <- as.numeric(distmatrix[root,])
        mean(distances[distances > 0])
    }

    get_mean_overall_node_distance <- function(distmatrix) {
        distmatrix[lower.tri(distmatrix,diag=T)] <- NA
        distances <- as.data.table(reshape2::melt(distmatrix))
        distances <- distances[!is.na(value),]
        mu <- mean(distances$value, na.rm=T)
        mu
    }
    
    met_specific_rds <- function(tm, normal_sample, primary_samples, met_samples) {
        ## rename all mets to 'Met1-X' so that we can combine samples across a group for RDS (e.g. PerOv/Per)
        ## note: as usual, this requires sample_info to be subset for only Normal, PT, and a single Met group
        map <- data.table(samplename=rownames(tm), pos=1:nrow(tm))
        map$samplename[map$samplename %in% met_samples] <- paste0('Met',1:length(met_samples))
        map <- map[order(pos)]
        rownames(tm) <- map$samplename; colnames(tm) <- map$samplename
        tree <- nj(tm)
        tree <- root(tree,outgroup=normal_sample,resolve.root=TRUE)
        klm_dat <- rds(tree)
        if(nrow(klm_dat)==0) klm_dat <- NULL
        klm_dat
    }

    distance_to_normal <- function(tm, normal_sample, primary_samples, met_samples) {
        ## given a distmatrix, return a data.table with the distance to normal for both the mets and the primaries 
        dm <- as.data.table(reshape2::melt(tm))
        dm <- dm[Var1!=normal_sample & Var2==normal_sample]
        dm[Var1 %in% primary_samples, group1:='Primary']
        dm[Var1 %in% met_samples, group1:='Metastasis']
        dm$group2 <- 'Normal'
        setnames(dm,c('Var1','Var2','value'),c('sample1','sample2','distance'))
        dm
    }

    distance_to_primary <- function(tm, normal_sample, primary_samples, met_samples) {
        ## given a distmatrix, return a data.table with the distance to normal for each 
        dm <- as.data.table(reshape2::melt(tm))
        dm <- dm[Var1 %in% met_samples & Var2 %in% primary_samples]
        dm$group1 <- 'Metastasis'
        dm$group2 <- 'Primary'
        setnames(dm,c('Var1','Var2','value'),c('sample1','sample2','distance'))
        dm
    }

    pairwise_distances <- function(tm, normal_sample=NULL, primary_samples=NULL, met_samples=NULL) {
        ## given a distmatrix, return a data.table with the distance to normal for each 

        ## first remove duplicate comparisons in distmatrix
        tm[lower.tri(tm)] <- NA

        ## format each unique comparison in a long table
        dm <- as.data.table(reshape2::melt(tm))
        dm$Var1 <- as.character(dm$Var1); dm$Var2 <- as.character(dm$Var2); 
        toadd <- t(apply(dm[,c('Var1','Var2'),with=F], 1, sort))
        dm$sample1 <- toadd[,1]
        dm$sample2 <- toadd[,2]
        dm <- dm[Var1!=Var2,]
        dm <- dm[!is.na(value),]
        dm[,c('Var1','Var2'):=NULL]
        dm <- dm[,c('sample1','sample2','value'),with=F]
        setnames(dm,'value','distance')

        if(!is.null(primary_samples) & !is.null(met_samples)) {
            ## get pairwise distances for the primary samples
            pt <- dm[sample1 %in% primary_samples & sample2 %in% primary_samples,]
            pt <- pt[order(sample1, sample2),]
            pt$group1 <- 'Primary'; pt$group2 <- 'Primary'

            ## get pairwise distances for the met samples
            mt <- dm[sample1 %in% met_samples & sample2 %in% met_samples,]
            mt <- mt[order(sample1, sample2),]
            mt$group1 <- 'Metastasis'; mt$group2 <- 'Metastasis'

            dm <- rbind(pt, mt)
            #factor <- mean(pt$distance)
            #dm$distance <- dm$distance / factor
        }
        dm
    }


    get_tree <- function(tree, normal_sample, primary_samples, met_samples, met_color, label='') {
        if(length(met_samples) > 0 & length(normal_sample)==1 & length(primary_samples) > 0) { 
            groups <- rbind(data.table(label=normal_sample,group='Normal',color='#000000'),
                            data.table(label=primary_samples,group='Primary',color='#008C45'),
                            data.table(label=met_samples,group='Metastasis',color=met_color))
            tmp <- groups[!duplicated(group),]
            mycols <- tmp$color; names(mycols) <- tmp$group
            p <- ggtree(tree,layout='ape') %<+% groups 
            p <- p + geom_tiplab(angle=0,aes(color=group)) + scale_color_manual(values=mycols,name='Sample type')
            p$data$label <- gsub('Normal','N',p$data$label)
            p <- p + guides(color='none')
        } else {
            p <- NULL
        }
        p
    }


    met_specific_distance <- function(mat, comparison, distance, return_tree=F, met_color='blue') { 
        ## default usage: specify a distance matrix that assumes 3 types: N, PT, and a single met group of interest (everything else).
        ## if comparison='normal': 
        ##  - returns the met-specific node distance to normal for those mets, and also PT samples in the same patient.
        ## if comparison='primary':
        ##  - returns the met-specific node distances between all pairs of the specific met-samples and PT samples.
        ## if comaprison='pairwise':
        ##  - returns all pairwise node distances between the met-samples, and also between the PT samples.

        mat <- as.matrix(mat)
        normal_sample <- rownames(mat)[str_detect(rownames(mat),'^N[0-9]|^Normal[0-9]')]
        primary_samples <- rownames(mat)[str_detect(rownames(mat),'^PT[0-9]|^PT[a-z][0-9]|^P[0-9]')]
        primary_samples <- primary_samples[grepl('-A',primary_samples)==F]
        met_samples <- rownames(mat)[!rownames(mat) %in% c(normal_sample, primary_samples)]
        valid_samples <- c(normal_sample, primary_samples, met_samples) 

        if(length(normal_sample) > 1) {
            stop('More than one normal sample!')
            out <- NULL

        } else if(length(normal_sample) < 1 | length(primary_samples) == 0){
            #cat(' (Insufficient samples)')
            out <- NULL
        } else if(length(valid_samples) < 3) { 
            #cat(' (Insufficient samples to construct tree)')
            out <- NULL
        } else if(length(met_samples) < 1 & comparison %in% c('primary','rds')) {
            #cat(' (Insufficient mets for comparison)')
            out <- NULL
        } else { 
            smat <- mat[valid_samples, valid_samples]    
            tree <- ape::nj(smat)
            tree <- root(tree,outgroup=normal_sample,resolve.root=TRUE)

            if(distance %in% c('node','scna')) {
                distmatrix <- as.matrix(adephylo::distTips(tree,tips="all",method="nNodes"))
            } else if(distance %in% c('patristic')) {
                distmatrix <- as.matrix(adephylo::distTips(tree,tips="all",method="patristic"))
            } else if(distance=='angular') {
                distmatrix <- smat
            } else {
                stop('distance must be among "node", "angular", "scna", "patristic"!')
            }

            ## determine how to normalize the data based on the comparison and distance combination
            if(comparison=='normal' & distance!='angular') scaling_factor <- get_mean_distance_to_root(distmatrix)          
            if(comparison=='normal' & distance=='angular') scaling_factor <- 1
            if(comparison=='primary' & distance!='angular') scaling_factor <- get_mean_overall_node_distance(distmatrix)          
            if(comparison=='primary' & distance=='angular') scaling_factor <- 1
            if(comparison=='pairwise' & distance!='angular') scaling_factor <- get_mean_overall_node_distance(distmatrix)          
            if(comparison=='pairwise' & distance=='angular') scaling_factor <- 1
            if(comparison=='all' & distance=='angular') scaling_factor <- 1
            if(comparison=='rds') scaling_factor <- 1

            if(comparison=='normal') {
                result <- distance_to_normal(distmatrix, normal_sample, primary_samples, met_samples)
                result$distance <- result$distance / scaling_factor
            } else if(comparison=='primary') {
                result <- distance_to_primary(distmatrix, normal_sample, primary_samples, met_samples) 
                result$distance <- result$distance / scaling_factor
            } else if(comparison=='pairwise') {
                result <- pairwise_distances(distmatrix, normal_sample, primary_samples, met_samples)
                result$distance <- result$distance / scaling_factor
            } else if(comparison=='all') {
                result <- pairwise_distances(distmatrix)
                result$distance <- result$distance / scaling_factor
                result[sample1 %in% normal_sample, group1:='Normal']
                result[sample2 %in% normal_sample, group2:='Normal']
                result[sample1 %in% primary_samples, group1:='Primary']
                result[sample2 %in% primary_samples, group2:='Primary']
                result[sample1 %in% met_samples, group1:='Metastasis']
                result[sample2 %in% met_samples, group2:='Metastasis']
            } else if(comparison=='rds') {
                if(length(met_samples)==0) stop('No mets!')
                result <- met_specific_rds(distmatrix, normal_sample, primary_samples, met_samples)
            } 

            if(comparison!='rds') {
                result$sample1 <- as.character(result$sample1)
                result$sample2 <- as.character(result$sample2)
                result$group1 <- as.character(result$group1)
                result$group2 <- as.character(result$group2)
            }

            ## return results, either a list of data+plot, or just the data
            if(return_tree) {
                p <- get_tree(tree, normal_sample, primary_samples, met_samples, met_color) 
                out <- list(data=result, plot=p)
            } else {
                out <-  result
            }
        }
        #cat('\n')
        out
    }


    run_met_specific_distance_for_patient <- function(patient, si, ad_table, comparison, distance, return_tree, met_color) { 
        #browser()
        ## usage:
        ## - adjust sample_info to include only Normal, Primary, and a desired type of Met (for example, synchronous peritoneal mets)
        ## - sample_info can be left uncollapsed (but this should only be for angular distance, not node distance)
        ## - this function will generate the request met-specific node distances using the three groups

        #cat(patient)
        mat <- ad_table[[patient]]

        ## subset the matrix for the samples in the provided sample_info. 
        info <- si[Patient_ID==patient,]
        valid_samples <- intersect(info$Real_Sample_ID, rownames(mat))
        mat <- mat[valid_samples, valid_samples]
        info <- info[Real_Sample_ID %in% valid_samples,]

        dat <- met_specific_distance(mat, comparison=comparison, distance=distance, 
                                     return_tree=return_tree, met_color=met_color)
        if(!is.null(dat))  dat$patient <- patient
        dat
    }

    ## run for all patients in the si table also in the ad_table
    all_patients <- unique(si$Patient_ID)
    valid_ad <- names(which(sapply(ad_table, is.null)==F))
    patients <- intersect(all_patients, valid_ad)
    l <- lapply(patients, run_met_specific_distance_for_patient, si, ad_table, 
                comparison=comparison, distance=distance, return_tree=return_tree, met_color=met_color)

    ## either return data merged into a data.table or a list with the data and tree for each
    if(return_tree==F) {
        rbindlist(l)
    } else {
        l
    }
}


mywilcox2 <- function(data, formula, paired, facet_field=NULL, include_n=T) {
    ## use this to get 'stat.test', which can be added to ggplot via e.g.
    ## stat_pvalue_manual(stat.test, label = "label", tip.length = 0.02) 
    require(rstatix)
    require(coin)

    toSN <- function(f) {
        ifelse(f < 0.001, formatC(f, format = "e", digits = 2), as.character(round(f, 3)))
    }
    wilcox_test2 <- function(...) {
        result <- wilcox_test(...)
        es <- wilcox_effsize(...)
        result$effsize <- es$effsize
        result
    }
    if(!is.null(facet_field)) {
        data$.group <- data[[facet_field]]
        stat.test <- data %>% group_by(.group) %>% wilcox_test2(formula, paired=paired)
        setnames(stat.test,'.group',facet_field)
    } else {
        stat.test <- data %>% wilcox_test2(formula, paired=paired)
    }
    stat.test <- stat.test %>% add_y_position()
    stat.test$label <- paste0('p=',toSN(stat.test$p),', ES=',round(stat.test$effsize,3)) 
  
    ## show the N for each group (or single N for paired)
    if(paired==T & include_n==T) {
       stat.test$label <- paste0(stat.test$label,', n=',stat.test$n1)
    } else if(include_n==T) {
       stat.test$label <- paste0(stat.test$label,', n=',stat.test$n1,';',stat.test$n2) 
    }
    stat.test
}


## generate minimal trees for node distances and RDS using Kim et al. to supplement liver patients.
## These are only used for node-distance analyses.

load_kim_node_distances <- function() {
    CRC1 <- ape::read.tree(text='(((((PT3, PT1), PT4), PT2), ((Liv1, Liv2), Liv3)), Normal1);')
    CRC1 <- root(CRC1,outgroup='Normal1',resolve.root=TRUE)
    CRC1 <- ape::rotate(CRC1,node=10)
    CRC1 <- ape::rotate(CRC1,node=14)
    CRC1 <- ape::rotate(CRC1,node=11)
    CRC1 <- ape::rotate(CRC1,node=12)
    CRC1 <- ape::rotate(CRC1,node=13)
    CRC2 <- ape::read.tree(text='(((((PT3, PT5), PT4), PT2), ((Liv1, Liv2), PT1)), Normal1);')
    CRC2 <- root(CRC2,outgroup='Normal1',resolve.root=TRUE)
    CRC2 <- ape::rotate(CRC2,node=10)
    CRC2 <- ape::rotate(CRC2,node=15)
    CRC2 <- ape::rotate(CRC2,node=13)
    CRC3 <- ape::read.tree(text='((((PT1, PT2), PT5), ((PT3, PT4), (((Liv2, Liv6), Liv1), ((Liv3, Liv4), Liv5)))), Normal1);')
    CRC3 <- root(CRC3,outgroup='Normal1',resolve.root=TRUE)
    CRC3 <- ape::rotate(CRC3,node=14)
    CRC3 <- ape::rotate(CRC3,node=17)
    CRC3 <- ape::rotate(CRC3,node=18)
    CRC3 <- ape::rotate(CRC3,node=19)
    CRC3 <- ape::rotate(CRC3,node=23)
    CRC4 <- ape::read.tree(text='(((((Liv2, Liv3), Liv4), Liv1), (PT1, PT2)), Normal1);')
    CRC4 <- root(CRC4,outgroup='Normal1',resolve.root=TRUE)
    CRC4 <- ape::rotate(CRC4,node=13)
    CRC4 <- ape::rotate(CRC4,node=12)
    CRC5 <- ape::read.tree(text='(((Liv1, Liv2), (PT1, PT2)), Normal1);')
    CRC5 <- root(CRC5,outgroup='Normal1',resolve.root=TRUE)
    CRC5 <- ape::rotate(CRC5,node=8)
    kim_trees <- list(CRC1=CRC1, CRC2=CRC2, CRC3=CRC3, CRC4=CRC4, CRC5=CRC5)
    f=function(tree) as.matrix(adephylo::distTips(tree, method='patristic'))
    kim_node_distances <- lapply(kim_trees, f)
    kim_node_distances
}


## rescale extremely-long terminal-branches (only use for plotting trees)
scale_branch <- function(tips, tree, sf=2) {
    newlength <- median(tree$edge.length)*sf
    selected_tip <- which(tree$tip.label %in% tips)
    for(sel in selected_tip) {
        selected_branch <- which(tree$edge[,2] == sel)
        tree$edge.length[selected_branch] <- newlength
    }
    tree
}


get_lmer <- function(res, title) {
    require(lme4)
    require(lmerTest)

    get_xpos <- function(pos) {
        pos <- pos / length(pos)
        offset <- mean(pos)-0.5
        pos - offset
    }

    get_patient_pos <- function(res) {
        res <- res[order(distance),]
        res$patient_pos <- 1:nrow(res)
        res$patient_pos <- get_xpos(res$patient_pos)
        res
    }

    patient_order <- function(res) {
        mid <- median(res$distance,na.rm=T)
        list(mid=mid)
    }

    mid_per <- res[class=='Per:Per',patient_order(.SD),by=c('patient')]
    mid_per <- mid_per[order(mid,decreasing=F),]
    res_per <- res[class=='Per:Per',]
    res_per$patient <- factor(res_per$patient, levels=mid_per$patient)
    res_per <- res_per[order(patient,distance,decreasing=F),]
    mid_liv <- res[class=='Liv:Liv',patient_order(.SD),by=c('patient')]
    mid_liv <- mid_liv[order(mid,decreasing=F),]
    res_liv <- res[class=='Liv:Liv',]
    res_liv$patient <- factor(res_liv$patient, levels=mid_liv$patient)
    res_liv <- res_liv[order(patient,distance,decreasing=F),]
    res2 <- rbind(res_per, res_liv)
    res2[,group:=paste0(patient,': ',group1)]
    res2$group <- factor(res2$group, levels=unique(res2$group))
    res2$group1 <- factor(res2$group1, levels=c('Liver','Peritoneum'))
    res2 <- res2[,get_patient_pos(.SD),by=c('group')]
    res2$global_pos <- (as.integer(res2$group)-1) + res2$patient_pos

    get_label_pos <- function(res2) {
        xpos <- max(res2$global_pos)
        xmid <- median(res2$global_pos)
        xstart <- xmid-0.25
        xend <- xmid+0.25
        ypos <- max(res2$distance)
        mid <- median(res2$distance)
        list(xpos=xpos, ypos=ypos, mid=mid, xstart=xstart, xend=xend)
    } 
    info <- res2[,get_label_pos(.SD),by=c('group','patient','group1')]

    #browser()
    tst <- wilcox_test(distance ~ group1, data=res2)
    es <- wilcox_effsize(distance ~ group1, data=res2)
    effect <- round(as.numeric(es$effsize),3)
    pval <- round(as.numeric(tst$p),3)
    stats1 <- paste0('Wilcox test (w.r.t liver): ES=',effect,'; P=',pval)

    require(geepack)
    #x <- res2[order(group, distance),]
    #mod <- geeglm(distance ~ group1, id=x$group, data=res2); summary(mod)

    x <- res2[order(patient, group1, distance),]
    #mod <- geeglm(distance ~ group1, id=x$patient, data=x, corstr='unstructured', scale.fix=T); summary(mod)
    mod <- geeglm(distance ~ group1, id=x$patient, data=x); summary(mod)

    coef <- summary(mod)$coef
    effect <- round(as.numeric(coef[2,1]),3)
    pval <- round(as.numeric(coef[2,4]),3)
    stats2 <- paste0('GEE-GLM (controlling for repeated measurements per patient): ES=',effect,'; P=',pval)

    mod <- lmer(distance ~ group1 + (1 | patient), data=res2, REML=F); summary(mod)
    coef <- summary(mod)$coef
    effect <- round(as.numeric(coef[2,1]),3)
    pval <- round(as.numeric(coef[2,5]),3)
    stats3 <- paste0('Mixed-effects linear model (each patient has own intercept): ES=',effect,'; P=',pval)

    p1 <- ggplot(res2, aes(x=global_pos, y=distance, group=group)) + 
        scale_y_continuous(limits=c(0,(max(res2$distance)+0.25)), breaks=seq(0,(max(res2$distance)+0.25),by=0.25)) + 
        geom_line(aes(color=group1)) +
        geom_point(pch=16,size=1,aes(color=group1)) +
        geom_segment(data=info, aes(y=mid, yend=mid,x=xstart, xend=xend), color='black',linewidth=0.5) +
        geom_text(data=info, aes(x=xpos, y=ypos, color=group1, label=patient), angle=0, hjust=1, vjust=-0.5) +
        scale_color_manual(values=group_cols) + 
        guides(color='none') +
        theme_ang(base_size=10) +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()) +
        labs(x=NULL,y='Distance',title=title,
             subtitle=paste(stats1,stats2,stats3,sep='\n'))

    p2 <- ggplot(res2, aes(y=distance)) + 
        scale_y_continuous(limits=c(0,(max(res2$distance)+0.25)), breaks=seq(0,(max(res2$distance)+0.25),by=0.25)) + 
        geom_density(aes(fill=group1),alpha=0.3,linewidth=0.25) + 
        scale_fill_manual(values=group_cols) + 
        guides(fill='none') +
        theme_ang(base_size=10) +
        theme(
              axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
              axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
        labs(y=NULL)

    p <- plot_grid(p1, p2, align='h', rel_widths=c(3,0.5), ncol=2, axis='tblr')
    p


}


dunn_test_ES <- function(data, formula, p.adjust.method='holm') {
    ## dunn test but add wilcox-test effect sizes to each pairwise comparison
    s <- as.character(formula)
    yvar <- s[2]
    xvar <- s[3]
    tst <- dunn_test(data, formula, p.adjust.method=p.adjust.method)
    tst$effsize <- as.numeric(NA) 
    for(i in 1:nrow(tst)) {
        x1 <- data[[yvar]][ data[[xvar]]==as.character(tst[i,2]) ]
        x2 <- data[[yvar]][ data[[xvar]]==as.character(tst[i,3]) ]
        tmp <- rbind(data.table(value=x1, group=1), data.table(value=x2, group=2))
        es <- as.numeric(wilcox_effsize(formula = value ~ group, data=tmp)$effsize)
        tst$effsize[i] <- es
    }
    tst$label <- paste0('p=',prettyNum(tst$p.adj, digits=1),', ES=',round(tst$effsize,3))
    tst
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define functions used for ACE lpWGS copy number pipeline
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## get 1Mb bins with adjusted-bin-copies and adjusted-segment-means
get_copynumber_adjusted_bins_for_sample <- function(sample, obj_list, fits, sex, exclude_mt=T) { 

    if(is.list(obj_list) & sample %in% names(obj_list) & sample %in% fits$barcode) {
        object <- obj_list[[sample]]    
        template <- objectsampletotemplate(object, index=1)
    } else if(grep('QDNAseq',class(obj_list)) & sample %in% sampleNames(obj_list) & sample %in% fits$barcode) {
        index <- grep(sample, sampleNames(obj_list)) 
        template <- objectsampletotemplate(object, index=index)
    } else {
        stop('obj_list is not set up properly!')
    }
    purity <- fits$purity[fits$barcode==sample]
    ploidy <- fits$ploidy[fits$barcode==sample]
    if(exclude_mt) template <- template[!template$chr %in% c('MT','M','chrM'),]

    ## get single-copy chromosomes
    if(sex=='female') {
        sgc <- c()
    } else if(sex=='male') {
        sgc <- c('X','Y')
    }

    ## get adjusted segments
    bindata <- as.data.table(template)
    bindata$chr <- gsub('chr','',as.character(bindata$chr))
    standard_chromosomes <- median(template$segments[template$chr %in% c(1:22)],na.rm=T) ## chromosomes used for standard
    segmentdf <- getadjustedsegments(template, ploidy=ploidy, cellularity = purity, standard=standard_chromosomes, sgc=sgc) 
    segmentdf <- as.data.table(segmentdf)
    setnames(segmentdf,c('Chromosome','Start','End'),c('chr','seg_start','seg_end'))
    segmentdf[,seg_id:=paste0(chr,':',seg_start/1e6,'-',seg_end/1e6)]
    segmentdf$seg_num <- as.integer(factor(segmentdf$seg_id, levels=unique(segmentdf$seg_id)))
    bindata$chr <- factor(bindata$chr, levels=unique(bindata$chr))
    segmentdf$chr <- factor(segmentdf$chr, levels=unique(segmentdf$chr))

    ## define the genome copies
    gc <- rep(2, nrow(template))
    gc[template$chr %in% sgc] <- 1

    ## adjust the bin's copy number with the same formula as for getting adjustedsegments
    bindata$adjustedcopynumbers <- bindata$copynumbers * (ploidy + 2/purity - 2)/standard_chromosomes - gc/purity + gc
    bindata <- cbind(sample=sample, bindata)

    setkey(bindata,'chr','start','end')
    setkey(segmentdf,'chr','seg_start','seg_end')
    bindata2 <- foverlaps(bindata, segmentdf, type='any')
    bindata2 <- bindata2[!is.na(adjustedcopynumbers),]

    segmentdf <- cbind(sample=sample, segmentdf)
    list(bins=bindata2, segs=segmentdf)
}


## test original ACE segments for subclonality 
test_segment <- function(bins, max_copies_tested) {
    copies <- unique(bins$Copies)
    mu <- mean(bins$adjustedcopynumbers)
    s <- sd(bins$adjustedcopynumbers)
    n_bins <- nrow(bins)
    if(mu >= copies & n_bins > 1 & copies %in% 0:(max_copies_tested-1)) {
        p_copies <- pnorm(q=copies, mean=mu, sd=s, lower.tail=T)
        p_altcopies <- pnorm(q=copies+1, mean=mu, sd=s, lower.tail=F)
        direction <- 'greater'
    } else if(mu < copies & n_bins > 1 & copies %in% 1:max_copies_tested) {
        p_copies <- pnorm(q=copies, mean=mu, sd=s, lower.tail=F)
        p_altcopies <- pnorm(q=copies-1, mean=mu, sd=s, lower.tail=T)
        direction <- 'less'
    } else {
        p_copies <- as.numeric(NA)
        p_altcopies <- as.numeric(NA)
        direction <- as.character(NA)
    }
    p_vals <- c(p_copies, p_altcopies)
    q_vals <- p.adjust(p_vals, method='BH')
    subclonal <- all(q_vals < 0.1)
    list(mu=mu, s=s, n_bins=n_bins, subclonal=subclonal, q_copies=q_vals[1], q_altcopies=q_vals[2], direction=direction)
}


## main script to extract/format data from the QDNAseq/ACE copy number object
process_scna_data_for_patient <- function(samples, obj_list, fits, sex, max_copies_tested=5, phylo_chromosomes=c(1:22)) {
    if(is.list(obj_list)) {
        valid <- all(samples %in% names(obj_list)) & all(samples %in% fits$barcode)
    } else if(grep('QDNAseq',class(obj_list))) {
        valid <- all(samples %in% sampleNames(obj_list)) & all(samples %in% fits$barcode)
    }
    if(valid==F) stop('Specified samples are not included in the obj_list and fits!')

    l <- lapply(samples, get_copynumber_adjusted_bins_for_sample, obj_list, fits, sex)
    get_bins <- function(l) l$bins
    bins <- rbindlist(lapply(l, get_bins))

    ## get unique segments in each sample. test for subclonal segments. annotate the bins accordingly
    get_segs <- function(l) l$segs
    segs <- rbindlist(lapply(l, get_segs))
    sig <- bins[chr %in% phylo_chromosomes,test_segment(.SD, max_copies_tested), by=c('sample','seg_id')]
    segs <- merge(segs, sig, by=c('sample','seg_id'), all.x=T)
    segs$Copies_heatmap <- as.numeric(segs$Copies)
    segs[!is.na(subclonal) & subclonal==T & mu > Copies, Copies_heatmap:=Copies + 0.5]
    segs[!is.na(subclonal) & subclonal==T & mu < Copies, Copies_heatmap:=Copies - 0.5]
    segs[Copies_heatmap < 0, Copies_heatmap:=0]
    segs[Copies_heatmap > max_copies_tested, Copies_heatmap:=max_copies_tested]

    ## convert bins to a distance matrix
    d2 <- data.table::dcast( chr + start + end + bin ~ sample, value.var='Segment_Mean2', data=bins)
    d2 <- d2[chr %in% phylo_chromosomes]
    m <- as.matrix(d2[,(samples),with=F])
    m <- t(m)
    dm <- as.matrix(dist(m, method='euclidean'))

    list(bins=bins, segs=segs, mat=d2, dm=dm)
}


## make a heatmap of each samples copy number segments ordered by euclidean distance NJ tree
scna_segment_heatmap <- function(dm, segs, chr, sex, groups, group_cols) {
    ## use the copy-number euclidean distance tree for ordering the samples
    tree <- nj(dm)
    tree <- root(tree,outgroup=grep(paste0('^N'),tree$tip.label),resolve.root=TRUE) 
    p_tree <- ggtree(tree,layout='rect',linewidth=0.5)
    dat <- as.data.frame(p_tree$data)
    dat <- dat[dat$isTip==T,]
    dat <- dat[order(dat$y,decreasing=F),]
    sample_levels <- dat$label[dat$isTip==T]

    ## for each sample, merge the segments into the full chromosome lengths
    new_chr <- copy(chr)
    if(sex=='male') {
        valid_chromosomes <- c(1:22,'X','Y')
    } else {
        valid_chromosomes <- c(1:22,'X')
    }
    segs <- segs[chr %in% valid_chromosomes]
    segs$chr <- factor(segs$chr, levels=valid_chromosomes)
    new_chr <- new_chr[chr %in% valid_chromosomes]
    new_chr[,chr_start:=1]
    new_chr[,chr_end:=length]
    new_chr$global_start <- as.numeric(NA)
    new_chr$global_end <- as.numeric(NA)
    new_chr[chr==1, global_start:=chr_start]
    new_chr[chr==1, global_end:=chr_end]
    new_chr$prev_regions <- as.numeric(NA)
    new_chr$prev_regions[1] <- 0
    for(i in 2:nrow(new_chr)) {
        new_chr$global_start[i] <- new_chr$global_end[i-1] + 1
        new_chr$global_end[i] <- new_chr$global_start[i] + new_chr$length[i]
        new_chr$prev_regions[i] <- new_chr$global_end[i-1]
    }

    setnames(new_chr,'length','chr_length')
    setkey(new_chr,'chr','chr_start','chr_end')
    new_chr$chr <- factor(new_chr$chr, valid_chromosomes)
    new_chr <- new_chr[order(chr)]
    gr_chr <- makeGRangesFromDataFrame(new_chr,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='chr_start', end.field='chr_end')

    message('Expanding segments to include NA regions in each chromosome ...')
    expand_segments_to_complete_chromosome_for_sample <- function(this.sample, segs, gr_chr) {
        message(this.sample)
        mat_sample <- segs[sample==this.sample,]
        mat_sample$chr <- factor(mat_sample$chr, levels(seqnames(gr_chr)))
        mat_sample <- mat_sample[order(chr)]
        gr_mat_sample <- makeGRangesFromDataFrame(mat_sample,keep.extra.columns=T,ignore.strand=T,seqnames='chr',start.field='seg_start', end.field='seg_end')
        NA_regions <- BiocGenerics::setdiff(gr_chr, gr_mat_sample)
        complete_regions <- as.data.table(sort(c(gr_mat_sample, NA_regions)))
        complete_regions$sample <- this.sample
        complete_regions[,segment:=paste0(seqnames,':',start,'-',end)]
        setnames(complete_regions,'seqnames','chr')
        complete_regions <- complete_regions[order(chr, start, end),]
        ## collapse regions with no difference in copy number
        complete_regions$sample <- this.sample
        complete_regions
    }
    sample_list <- lapply(samples, expand_segments_to_complete_chromosome_for_sample, segs, gr_chr)
    complete_segs <- rbindlist(sample_list)
    complete_segs$i <- 1:nrow(complete_segs)
    complete_segs$sample <- factor(complete_segs$sample, levels=sample_levels)
    complete_segs$sample_num <- as.integer(complete_segs$sample) - 0.5
    complete_segs <- merge(complete_segs, new_chr[,c('chr','prev_regions'),with=F], by='chr', all.x=T)
    complete_segs[, segment_start:=start + prev_regions]
    complete_segs[, segment_end:=end + prev_regions]
    complete_segs$global_midpoint <- complete_segs$prev_regions + (complete_segs$start + complete_segs$end) / 2
    new_chr$global_midpoint <- (new_chr$global_start + new_chr$global_end) / 2
 
    ## get colors for the y-axis labels
    tmp <- data.table(label=sample_levels)
    tmp$pos <- 1:nrow(tmp)
    tmp <- merge(tmp, groups, by='label', all.x=T)
    tmp <- tmp[order(pos),]
    y_cols <- group_cols[tmp$group]

    blues <- brewer.rdbu(9)[9:6]
    reds <-  brewer.rdbu(13)[6:1]
    cols <- c(blues, 'white', reds)
    names(cols) <- seq(0,5,by=0.5)
    
    ## recode the copy-names
    vals <- data.table(Copies_heatmap=seq(0,5,by=0.5))
    vals[,lwr:=floor(Copies_heatmap)]
    vals[,upr:=ceiling(Copies_heatmap)]
    vals$Copies_heatmap_val <- as.character(vals$Copies_heatmap)
    vals[lwr!=upr,Copies_heatmap_val:=paste0(lwr,'-',upr)]
    names(cols) <- vals$Copies_heatmap_val
    complete_segs <- merge(complete_segs, vals[,c('Copies_heatmap','Copies_heatmap_val'),with=F], by='Copies_heatmap', all.x=T)
    complete_segs <- complete_segs[order(i), ]
    complete_segs$Copies_heatmap_val <- factor(complete_segs$Copies_heatmap_val, levels=vals$Copies_heatmap_val)

    p_heatmap <- ggplot(complete_segs) + 
        scale_x_continuous(expand=c(0,0), breaks=new_chr$global_midpoint, labels=new_chr$chr) + 
        scale_y_continuous(expand=c(0,0), breaks=seq(0.5,by=1,length.out=length(sample_levels)), labels=sample_levels, position='right') +
        geom_rect(aes(xmin=segment_start, xmax=segment_end, ymin=sample_num-0.5, ymax=sample_num+0.5, fill=Copies_heatmap_val)) +
        geom_point(data=complete_segs[subclonal==T & !is.na(subclonal)], size=0.5, pch=16, color='black', aes(x=global_midpoint,y=sample_num)) +
        geom_vline(xintercept=c(0,new_chr$global_end), color='black', linewidth=0.5) +
        scale_fill_manual(values=cols,na.value='#bfbfbf',name='Copies') +
        theme_ang(base_size=10) +
        labs(x=NULL,y=NULL) +
        theme(axis.line.x=element_blank(), axis.line.y=element_blank(), axis.text.y = element_text(colour = y_cols))

    plot_grid(p_tree, p_heatmap, align='h', rel_widths=c(1,5))
}


# use a permutation test to test for non-random similarity between two trees (via quartet distance)
test_tree_similarity <- function(test_mat,ref_mat,nperm,title=NULL) { 
    ## subset both matrices to the common samples 
    common_samples <- intersect(rownames(test_mat),rownames(ref_mat))
    test_mat <- test_mat[common_samples,common_samples]
    ref_mat <- ref_mat[common_samples,common_samples]

    permute_mat <- function(i, mat) {
        original_samples <- rownames(mat)
        if(i==0) {
            permuted_samples <- copy(original_samples)
        } else {
            permuted_samples <- sample(original_samples,replace=F)
        }
        perm_mat <- mat[permuted_samples,permuted_samples]    
        rownames(perm_mat) <- original_samples; colnames(perm_mat) <- original_samples
        nj(perm_mat)        
    }   

    permute_test_and_get_shared_info <- function(i, test_mat, ref_mat) {
        test_tree_perm <- permute_mat(i, test_mat) ## permute the test tree, keep ref tree unchanged
        ref_tree <- nj(ref_mat)
        info <- as.data.frame(QuartetStatus(test_tree_perm, cf=ref_tree))
        info$perm <- i
        info
    }
    l <- lapply(0:nperm, permute_test_and_get_shared_info, test_mat, ref_mat)
    res <- rbindlist(l)
    res$similarity <- res$s / res$Q
    obs <- res[perm==0,]
    exp <- res[perm > 0,]
    numerator <- sum(exp$similarity >= obs$similarity) + 1
    denominator <- nrow(exp) + 1    
    pval <- numerator / denominator
    pval <- prettyNum(pval,digits=2)
    median_exp <- prettyNum(median(exp$similarity),digits=2)
    confint_exp <- prettyNum(quantile(exp$similarity,c(0.025,0.975)),digits=2)
    label <- paste0('Exp: ',median_exp,', [',confint_exp[1],',',confint_exp[2],']; Obs: ',prettyNum(obs$similarity,digits=3),'; P=',pval,'; bs=',nperm)
    p <- ggplot(exp,aes(x=similarity)) +
        geom_histogram(bins=50,fill='#a6a6a6',linewidth=1,color='white') +
        theme_ang(base_size=12) +
        geom_segment(x=obs$similarity,xend=obs$similarity,y=0,yend=Inf,color='red')
    p <- p + labs(x='Tree similarity',y='N permutations',subtitle=label)   
    p <- p + scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25),expand=c(0,0))
    list(data=res, plot=p)
}


## performs distance matrix similarity test and plots the matrices side-by-side including summary of test results
compare_matrices <- function(cnv_distance_matrix, polyg_distance_matrix, this.patient, R) { 
    ## align CNV and polyG distance matrices
    common_samples <- intersect(rownames(cnv_distance_matrix),rownames(polyg_distance_matrix))
    common_samples <- common_samples[!grepl('^Normal',common_samples)]
    sample_levels <- sort(common_samples, decreasing=T)
    d_subset <- cnv_distance_matrix[common_samples,common_samples]
    g_subset <- polyg_distance_matrix[common_samples,common_samples]
    tst <- dist_similarity(test_dist=d_subset, ref_dist=g_subset, nperm=R, return_only_pval=F, method='spearman')

    d_subset <- as.data.table(reshape2::melt(tst$test_dist))
    d_subset$data <- 'CNV'
    d_subset$Var1 <- factor(d_subset$Var1, levels=sample_levels)
    d_subset$Var2 <- factor(d_subset$Var2, levels=sample_levels)
    g_subset <- as.data.table(reshape2::melt(tst$ref_dist))
    g_subset$data <- 'poly-G'
    g_subset$Var1 <- factor(g_subset$Var1, levels=sample_levels)
    g_subset$Var2 <- factor(g_subset$Var2, levels=sample_levels)

    p1 <- ggplot(d_subset, aes(x=Var1,y=Var2)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        theme_ang(base_size=10) +
        geom_tile(aes(fill=value)) +
        theme(
              axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(), 
              legend.position='right') +
        scale_fill_gradient(low='white',high='steelblue',name='Euclidean\ndistance') +
        labs(x=NULL,y=NULL,subtitle=paste(this.patient,'SCNA distance matrix'))

    p2 <- ggplot(g_subset, aes(x=Var1,y=Var2)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        theme_ang(base_size=10) +
        geom_tile(aes(fill=value)) +
        theme(
              axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(), 
              legend.position='right') +
        scale_fill_gradient(low='white',high='steelblue',name='Angular\ndistance') +
        labs(x=NULL,y=NULL,subtitle=paste(this.patient,'Polyguanine distance matrix'))
    p <- plot_grid(p1, p2, tst$plot, ncol=3)
    p
}


## based on if patient has male sex, decide which chromosomes to include and their germline ploidy
sex_chr_info <- function(sex) {
    get_sgc <- function(sex) {
        if (sex == "male") {
            c("X", "Y")
        }
        else if (sex == "female") {
            c()
        }
    }
    sgc <- get_sgc(sex)
    if ("Y" %in% sgc) {
        chrsubset = 1:24
        excluded_chr = c("X", "Y", "MT")
    }
    else {
        chrsubset = 1:23
        excluded_chr = c("Y", "MT")
    }
    list(sgc=sgc, plot_chr_included=chrsubset, model_chr_excluded=excluded_chr)
}


## main function to input the QDNAseq copy number object and apply ACE to see fit for some purity/ploidy combo.
## - if save=F, the a plot is shown to the user. if save=T, it is saved to the outpur directory.
## - if purity=NA, the local minima are printed as possible solutions via ACE.
## - all chromosomes are refit with the purity/ploidy and retained in the output; 
## - however, only chr1:22 are used in ACE to fit the purity/ploidy.
## - this is because QDNAseq corrected for GC bias based on estimates obtained only from autosomes
refit <- function(object, sex, purity=NA, ploidy, samplename=NA, sampleindex=1, save=F, bottom=-1, cap=10, output_dir='.') {
    get_sgc <- function(sex) {
        if(sex=='male') {
            c('X','Y')
        } else if(sex=='female') {
            c()
        }
    }
    sgc <- get_sgc(sex)
    if('Y' %in% sgc) {
        ## 1X, 1Y
        chrsubset=1:24
        excluded_chr=c('X','Y','MT')
    } else {
        ## 2X, 0Y
        chrsubset=1:23
        excluded_chr=c('Y','MT')
    }

    if(is.na(samplename)) {
        pd <- Biobase::pData(object)
        samplename <- gsub("_.*", "",pd$name)
    }

    if(is.na(purity)) {
        model <- singlemodel(object, QDNAseqobjectsample=sampleindex, ploidy=ploidy, exclude=excluded_chr); model$minima
        return(model$minima)
    }

    title = paste0(samplename,': ploidy=',ploidy,', purity=',purity)
    template <- objectsampletotemplate(object, index=1)
    standard_chromosomes <- median(template$segments[template$chr %in% c(1:22)],na.rm=T) ## only use these chrs for calculating the standard, but show all chromosomes
    p <- singleplot(object, QDNAseqobjectsample=sampleindex, cellularity=purity, standard=standard_chromosomes, ploidy=ploidy, bottom=bottom, cap=cap, sgc=sgc, onlyautosomes=F, chrsubset=chrsubset, title=title)

    message(title)
    if(save==F) {
        ## return the results
        p
    } else {
        refits_dir=file.path(output_dir,'fits')
        if(!dir.exists(refits_dir)) dir.create(refits_dir,recursive=T)

        ## save the plot for this fit
        plot_file <- paste0(refits_dir,'/',samplename,'_N=',ploidy,'_cellularity=',purity,'.pdf')
        #message('Saving plot: ',plot_file)
        ggsave(p,dev=cairo_pdf,filename=plot_file,width=8,height=5)

        ## save data (append it to existing file along with the date of the entry)
        fit_file <- paste0(refits_dir,'/purity_ploidy.txt')
        if(!file.exists(fit_file)) cat('barcode\tpurity\tploidy\tdate\n',file=fit_file,append=F) 
        when <- as.character(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
        #message('Saving data: ',fit_file)
        cat(paste(samplename,purity,ploidy,when,sep='\t'),'\n',file=fit_file,append=T)
    }
}


distance_matrix_correlation <- function(test_dist_full, ref_dist_full, method, return_only_r=F) {
    ref_dist <- ref_dist_full[order(rownames(ref_dist_full)),order(colnames(ref_dist_full))]    
    test_dist <- test_dist_full[order(rownames(test_dist_full)),order(colnames(test_dist_full))]    

    ## subset the distances for common samples 
    common_samples <- sort(intersect(rownames(test_dist), rownames(ref_dist)))
    test_dist <- test_dist[common_samples,common_samples]
    ref_dist <- ref_dist[common_samples,common_samples]
    
    ## get a long table of sample:sample distances with columns for test and ref data
    ref_dist_long <- ref_dist; test_dist_long <- test_dist; 
    ref_dist_long[upper.tri(ref_dist_long)] <- NA
    test_dist_long[upper.tri(test_dist_long)] <- NA
    test_dist_long <- as.data.table(reshape2::melt(test_dist_long))
    ref_dist_long <- as.data.table(reshape2::melt(ref_dist_long))
    merged <- merge(test_dist_long,ref_dist_long,by=c('Var1','Var2'))
    merged <- merged[Var1!=Var2 & !is.na(value.x) & !is.na(value.y)] 
    merged[,c('Var1','Var2'):=NULL]
    names(merged) <- c('test_dist','ref_dist')

    ## get the correlation
    r <- cor(merged$test_dist,merged$ref_dist,method=method)    

    if(return_only_r==T) {
        out <- r
    } else {
        out <- list(r=r, merged=merged, test_dist=test_dist, ref_dist=ref_dist, test_dist_full=test_dist_full, ref_dist_full=ref_dist_full)
    }
    out
}


dist_similarity <- function(test_dist, ref_dist, nperm, return_only_pval, method) {

    shuffled_correlation <- function(i, test_dist, ref_dist, return_only_r, method) { 
        if(i > 0) {
            ## shuffle among all the samples with data (not only the common ones)
            orig_samples <- rownames(test_dist)
            shuffled_samples <- sample(orig_samples,replace=F)
            test_dist <- test_dist[shuffled_samples,shuffled_samples]
            rownames(test_dist) <- orig_samples
            colnames(test_dist) <- orig_samples
        }
        info <- distance_matrix_correlation(test_dist, ref_dist, return_only_r=F, method=method)
        r <- info$r
        merged <- info$merged
        test_dist <- info$test_dist
        ref_dist <- info$ref_dist
        test_dist_full <- info$test_dist_full
        ref_dist_full <- info$ref_dist_full
        ## return either just the distance or the complete data
        if(return_only_r) { 
            r
        } else {
            list(test_dist=test_dist, ref_dist=ref_dist, test_dist_full=test_dist_full, ref_dist_full=ref_dist_full, merged=merged, method=method, r=r)
        }
    }

    observed_info <- shuffled_correlation(0, test_dist, ref_dist, return_only_r=F, method=method)
    ref_dist <- observed_info$ref_dist
    test_dist <- observed_info$test_dist
    observed <- observed_info$r
    merged <- observed_info$merged

    permuted <- unlist(lapply(1:nperm, shuffled_correlation, test_dist, ref_dist, return_only_r=T, method=method))
    perms_as_correlated <- sum(permuted >= observed)
    p.value <- (perms_as_correlated + 1) / (nperm + 1) ## add pseudo-count to avoid p=0

    if(return_only_pval==T) {
        p.value
    } else {
        p.value <- prettyNum(p.value,digits=2)
        median_exp <- prettyNum(median(permuted),digits=2)
        confint_exp <- prettyNum(quantile(permuted,c(0.025,0.975)),digits=2)
        label <- paste0('Exp: ',median_exp,' [',confint_exp[1],',',confint_exp[2],']; Obs: ',prettyNum(observed,digits=3),'; P=',p.value,'; bs=',nperm)
        tmp <- data.table(permuted=permuted)
        xmin <- min(c(tmp$permuted, observed))
        xmin <- xmin - 0.1*xmin        
        xmax <- max(c(tmp$permuted, observed))
        xmax <- xmax + 0.1*xmax        
        xrange <- c(xmin, xmax)
        plot <- ggplot(tmp,aes(x=permuted)) +
            scale_x_continuous(limits=xrange) + 
            geom_histogram(bins=50,fill='#a6a6a6',linewidth=1,color='white') +
            theme_ang(base_size=12) +
            geom_segment(x=observed,xend=observed,y=0,yend=Inf,color='red') +
            labs(x='Spearman correlation',y='N permutations',subtitle=label)
        list(p.value=p.value, observed=observed, permuted=permuted, merged=merged, ref_dist=ref_dist, test_dist=test_dist, plot=plot)
    }
}


break_axis <- function(y, maxlower, minupper=NA, lowerticksize, upperticksize, ratio_lower_to_upper) {
    if(is.na(minupper)) {
        breakpos <- maxlower
        lowerticklabels <- seq(0,breakpos,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(breakpos+upperticksize,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- breakpos + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > breakpos
        newy[ind] <- breakpos + uppertickspacing*((newy[ind]-breakpos) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    } else {
        lowerticklabels <- seq(0,maxlower,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(minupper,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- maxlower + 0.5*lowerticksize + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > maxlower
        newy[ind] <- maxlower + 0.5*lowerticksize + 1*uppertickspacing + uppertickspacing*((newy[ind]-minupper) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    }
}


extract_gglegend <- function(p){
    require(ggplot2)
    require(cowplot)

    ## extract the legend from a ggplot object
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if(length(leg) > 0) leg <- tmp$grobs[[leg]]
    else leg <- NULL
    leg

    ## return the legend as a ggplot object
    legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
    plot <- p + theme(legend.position='none')
    list(plot=plot,legend=legend)
}


plot_bootstrapped_tree_unrooted <- function(tree, groups, title, size=4, bs_size=3) {
    ## add the bootstrap values to the tree
    p <- ggtree(tree,layout='ape') 
    p <- p %<+% groups
    suppressWarnings(p$data$boot <- as.numeric(p$data$label))
    maxx <- max(p$data$x)*1.1
    #p <- p + geom_nodelab(angle=0,color='blue')
    if(bs_size > 0) p <- p + geom_text_repel(aes(label=boot),color='blue',min.segment.length=0.15,seed=42,size=bs_size,max.overlaps=100)
    p <- p + geom_tiplab(fontface=1,size=size,hjust=0,angle=0,aes(color=group)) 
    p <- p + scale_color_manual(values=group_cols,name='Tissue') 
    p <- p + theme(legend.position='bottom') + guides(color='none')
    p <- p + labs(title=title) 
    #p <- p + xlim(c(-0.2,maxx))
    p
}


plot_bootstrapped_tree <- function(tree, groups, title) {
    ## add the bootstrap values to the tree
    p <- ggtree(tree,layout='rect') 
    p <- p %<+% groups
    suppressWarnings(p$data$boot <- as.numeric(p$data$label))
    maxx <- max(p$data$x)*1.1
    p <- p + geom_label(aes(label=boot,fill=boot), color='black',
                        label.padding=unit(0.1, "lines"),
                        label.r = unit(0, "lines"))
    p <- p + scale_fill_gradient(low='white',high='steelblue',name='Bootstrap value',limits=c(0,100))
    p <- p + geom_tiplab(fontface=1,size=4,hjust=0,angle=F,aes(color=group)) 
    p <- p + scale_color_manual(values=group_cols,name='Tissue') 
    p <- p + theme(legend.position='bottom') + guides(color='none')
    p <- p + labs(title=title) 
    p <- p + xlim(c(0,maxx))
    p
}


get_bootstrapped_scna_tree <- function(bins, dm, patient, phylo_chromosomes=c(1:22), R=1000) {
    bins <- bins[chr %in% phylo_chromosomes,]
    bins <- data.table::dcast(chr + bin + start + end ~ sample, value.var='Segment_Mean2', data=bins)
    get_resampled_bin_tree <- function(i, bins, samples) {
        chrs_with_replacement <- sample(phylo_chromosomes,replace=T)
        get_bins_for_chromosome <- function(this.chr, bins) {
            bins <- bins[chr %in% this.chr] 
            bins  
        }
        resampled_bins <- lapply(chrs_with_replacement, get_bins_for_chromosome, bins)
        resampled_bins <- rbindlist(resampled_bins)
        resampled_dm <- as.matrix(dist(t(as.matrix(resampled_bins[,(samples),with=F])), method='euclidean'))
        tree <- nj(resampled_dm)
        tree
    }
    bstrees <- lapply(1:R, get_resampled_bin_tree, bins, samples)
    bstrees <- TreeTools::as.multiPhylo(bstrees)
    tree <- nj(dm)
    tree <- root(tree,outgroup=grep(paste0('^N'),tree$tip.label),resolve.root=TRUE) ## root before adding boostrap values!
    tree <- addConfidences(tree, bstrees) 
    tree$node.label <- round(100*tree$node.label)
    p <- plot_bootstrapped_tree(tree, groups, paste0(patient,' SCNA tree (1000 bootstrap replicates)'))
    list(tree=tree, bstrees=bstrees, plot=p)
}


get_distance_from_bootstrapped_trees_to_original_tree <- function(tree, bstrees) { 
    class(bstrees) <- 'list'

    tree_quartet_similarity <- function(bs_tree, orig_tree) {
        # Q - The total number of quartets for _n_ leaves.
        # s - The number of quartets that are resolved identically in both trees.
        info <- as.data.frame(QuartetStatus(bs_tree, cf=orig_tree))
        info
    }

    message('Calculating quartet similarity for each resampled tree to original tree ...')
    l <- rbindlist(lapply(bstrees, tree_quartet_similarity, tree))
    l$similarity <- l$s / l$Q
    l
}


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


make_tree <- function(this.patient, sample_info, collapsed, outdir, show.depth=F, show.timing=F, show.bsvals=F, min.bsval.shown=50, tree.layout='ape', overwrite=F) { 
    set.seed(42)
    if(this.patient!='C38') {
        ad_file <- here(paste0('processed_data/angular_distance_matrices/',this.patient,'.txt'))
    } else {
        ad_file <- here(paste0('processed_data/polyG/science/results/sample_exclusion_0.3_rep_cut_0.11/results_angular_distance_representativeReplicates/angular_dist_matrix_w_root_usedmarkers/C38_angular_dist_matrix_w_root_usedmarkers_repreReplicate_newnames.txt'))
    }
    ad <- read_distance_matrix(ad_file)
    outfile <- file.path(outdir,paste0(this.patient,'.pdf'))
    if(file.exists(outfile) & overwrite==F) {
        message(this.patient,': file already exists.')
        return(NULL)
    }

    if(!dir.exists(outdir)) dir.create(outdir, recursive=T)
    message(this.patient)
    if(this.patient!='E3') {
        si <- sample_info[Patient_ID %in% this.patient]
    } else {
        si <- sample_info[grepl('E[a-c]3',Patient_ID),]
        si <- si[!duplicated(Sample_ID),]   
    }
    if(collapsed) si <- si[in_collapsed==T,]
    if(this.patient!='C38') {
        valid_samples <- intersect(rownames(ad), si$Real_Sample_ID)
        ad <- ad[valid_samples, valid_samples]
    } else {
        valid_samples <- rownames(ad)
    }
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

    if(tree.layout=='ape' & show.bsvals==T) tree$node.label[tree$node.label < min.bsval.shown] <- NA
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

    tree_title <- paste0(this.patient,'. BS values ',min.bsval.shown,'%+ shown')
    if(truncated_normal) tree_title <- paste0(tree_title,'. Normal branch truncated.')
    p <- p + theme(legend.position='bottom') + ggtitle(tree_title)

    if(tree.layout=='ape') 
        ggsave(outfile, plot=p, width=10, height=8)
    else
        ggsave(outfile, plot=p, width=8, height=10)
}


run_chemo_simulation <- function(sim, cells_start, cells_end, frac_surviving_chemo, b, d) { 
    message('Simulation: ',sim)
    
    # Define relevant parameters
    clones <- c(1:10)   # create some categorical clones
    PMclones <- 5
    Lclones <- 5
    number_lesions <- 4

    diversity_PM <- sample(clones,PMclones,replace=F) # choose clones that can populate PM
    diversity_L <- sample(clones,Lclones,replace=F) # choose clones that can populate L

    ## for PMs, clone prob is uniform (or 0 for non-PM-seeding clones)
    PMprobs <- rep(0,length(clones))
    PMprobs[clones %in% diversity_PM] <- rep(1,PMclones)/PMclones

    ## for Ls, clone prob is skewed (or 0 for non-PM-seeding clones)
    Lprobs <- rep(0,length(clones))
    lowprob <- 0.05 / (Lclones-1)
    highprob <- 0.95
    Lprobs[clones %in% diversity_L] <- c(highprob, rep(lowprob, Lclones-1))

    # create pre-treatment PM and L lesions. These are tables of the number of cells in each lesion corresponding to each clone
    PMs <- t(rmultinom(number_lesions, cells_start, prob = PMprobs))
    PMs <- apply(PMs, 2, as.numeric)
    Ls <- t(rmultinom(number_lesions, cells_start, prob = Lprobs))
    Ls <- apply(Ls, 2, as.numeric)
    rownames(PMs) <- paste0('lesion',1:4)
    colnames(PMs) <- paste0('clone',1:10)
    rownames(Ls) <- paste0('lesion',1:4)
    colnames(Ls) <- paste0('clone',1:10)

    # name rows/cols for diagnosis
    PMs_pc <- matrix(nrow=nrow(PMs), ncol=ncol(PMs))
    Ls_pc <- matrix(nrow=nrow(Ls), ncol=ncol(Ls))
    PMs_regrow <- matrix(nrow=nrow(PMs), ncol=ncol(PMs))
    Ls_regrow <- matrix(nrow=nrow(Ls), ncol=ncol(Ls))
    names(PMprobs) <- colnames(PMs)
    names(Lprobs) <- colnames(Ls)

    # kill 2/3 of cells in each PM lesion
    for(i in 1:nrow(PMs)) {
        cell_bucket_prechemo <- rep(clones, PMs[i,])
        total_cells <- sum(PMs[i,])
        cells_surviving <- round(total_cells*frac_surviving_chemo)
        cell_bucket_postchemo <- sample(cell_bucket_prechemo, size=cells_surviving, replace=F)
        cell_bucket_postchemo <- factor(cell_bucket_postchemo, levels=clones)
        tbl <- table(cell_bucket_postchemo) 
        PMs_pc[i,] <- as.integer(tbl)
    }

    # kill 2/3 of cells in each L lesion
    for(i in 1:nrow(Ls)) {
        cell_bucket_prechemo <- rep(clones, Ls[i,])
        total_cells <- sum(Ls[i,])
        cells_surviving <- round(total_cells*frac_surviving_chemo)
        cell_bucket_postchemo <- sample(cell_bucket_prechemo, size=cells_surviving, replace=F)
        cell_bucket_postchemo <- factor(cell_bucket_postchemo, levels=clones)
        tbl <- table(cell_bucket_postchemo) 
        Ls_pc[i,] <- as.integer(tbl)
    }

    # my original attempt to regrow lesions by killing off/dividing large chunks of cells at once
    PMs_regrow <- regrow_lesions(PMs_pc, prob_death=d/(b+d), max_cells=cells_end)
    Ls_regrow <- regrow_lesions(Ls_pc, prob_death=d/(b+d), max_cells=cells_end)

    inter_lesion_diversity  <- function(mat) {
        dm <- dist(mat,method="euclidian")
        median(dm) 
    }

    # pre-treatment intra-lesion diversity
    pm_intra_prechemo <- median(diversity(PMs,MARGIN=1))
    l_intra_prechemo <- median(diversity(Ls,MARGIN=1))

    # pre-treatment, post-chemo intra-lesion diversity
    pm_intra_postchemo_preregrowth <- median(diversity(PMs_pc,MARGIN=1))
    l_intra_postchemo_preregrowth <- median(diversity(Ls_pc,MARGIN=1))

    # post-regrowth intra-lesion diversity
    pm_intra_postregrowth <- median(diversity(PMs_regrow,MARGIN=1))
    l_intra_postregrowth <- median(diversity(Ls_regrow,MARGIN=1))

    # pre-treatment inter-lesion diversity
    pm_inter_prechemo <- inter_lesion_diversity(PMs)
    l_inter_prechemo <- inter_lesion_diversity(Ls)

    # pre-treatment, post-chemo inter-lesion diversity
    pm_inter_postchemo_preregrowth <- inter_lesion_diversity(PMs_pc)
    l_inter_postchemo_preregrowth <- inter_lesion_diversity(Ls_pc)

    # post-regrowth inter-lesion diversity
    pm_inter_postregrowth <- inter_lesion_diversity(PMs_regrow)
    l_inter_postregrowth <- inter_lesion_diversity(Ls_regrow)

    intra_data <- rbind(
                        data.frame(group='PM', setting='Pre-chemo', het_type='Intra-lesion', value=pm_intra_prechemo),
                        data.frame(group='PM', setting='Post-chemo, pre-regrowth', het_type='Intra-lesion', value=pm_intra_postchemo_preregrowth),
                        data.frame(group='PM', setting='Post-chemo, post-regrowth', het_type='Intra-lesion', value=pm_intra_postregrowth),
                        data.frame(group='L', setting='Pre-chemo', het_type='Intra-lesion', value=l_intra_prechemo),
                        data.frame(group='L', setting='Post-chemo, pre-regrowth', het_type='Intra-lesion', value=l_intra_postchemo_preregrowth),
                        data.frame(group='L', setting='Post-chemo, post-regrowth', het_type='Intra-lesion', value=l_intra_postregrowth))

    inter_data <- rbind(
                        data.frame(group='PM', setting='Pre-chemo', het_type='Inter-lesion', value=pm_inter_prechemo),
                        data.frame(group='PM', setting='Post-chemo, pre-regrowth', het_type='Inter-lesion', value=pm_inter_postchemo_preregrowth),
                        data.frame(group='PM', setting='Post-chemo, post-regrowth', het_type='Inter-lesion', value=pm_inter_postregrowth),
                        data.frame(group='L', setting='Pre-chemo', het_type='Inter-lesion', value=l_inter_prechemo),
                        data.frame(group='L', setting='Post-chemo, pre-regrowth', het_type='Inter-lesion', value=l_inter_postchemo_preregrowth),
                        data.frame(group='L', setting='Post-chemo, post-regrowth', het_type='Inter-lesion', value=l_inter_postregrowth))

    out <- rbind(intra_data, inter_data)
    out$simulation <- sim
    out
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reload polyG and kim et al data for following
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sample_info <- fread(here('processed_data/sample_info.txt'))
sample_info <- sample_info[cohort!='lung',]
#sample_info$met_treated_type <- ''
#sample_info[met_treated %in% c('hipec','hipec after untreated'), met_treated_type:='hipec']
#sample_info[met_treated %in% c('systemic chemo','systemic chemo after untreated'), met_treated_type:='systemic chemo']
#sample_info[met_treated %in% c('systemic chemo','hipec'), met_treated:='treated']
#sample_info[met_treated %in% c('systemic chemo after untreated','hipec after untreated'), met_treated:='treated after untreated']

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


# ~~~~~~~~~~~~~~~~~~
# color schemes 
# ~~~~~~~~~~~~~~~~~~

group_cols <- c("#000000", "#008C45", "#EB5B2B", "#FAB31D", "#4C86C6", "#4C86C6", "#4C86C6")
names(group_cols) <- c('Normal', 'Primary', 'Locoregional', 'Peritoneum', 'Lung', 'Liver', 'Distant (other)')

## tissue color scheme
cols <- c("#008C45", "#EB5B2B", "#EB5B2B", "#EB5B2B", "#FAB31D", "#FAB31D", "#4C86C6", "#4C86C6", "#4C86C6", "#4C86C6","#4C86C6")
names(cols) <- c('PT','LN','TD','LR','Per','PerOv','Liv','Lun','OvH','dMet','PT-A')

## get colors for depth
cols_depth <- cols
cols_depth_luminal <- lighten(cols, 0.15)
names(cols_depth_luminal) <- paste0(names(cols_depth),'_luminal')
cols_depth_deep <- darken(cols, 0.15)
names(cols_depth_deep) <- paste0(names(cols_depth),'_deep')
cols_depth <- c(cols_depth, cols_depth_deep, cols_depth_luminal)

## timing colors
cols_timing <- c('white','black')
names(cols_timing) <- c('synchronous','metachronous')

cols_depth2 <- c('#BA1822','#2F3080')
names(cols_depth2) <- c('deep','mucosal/luminal')



