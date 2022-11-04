#Normalization functions:

noob <- function(rgSet){
  minfi::preprocessNoob(rgSet)
}

funn <- function(rgSet){
  minfi::preprocessFunnorm(rgSet)
}
noob_pq <- function(rgSet){
  minfi::preprocessNoob(rgSet)%>% minfi::preprocessQuantile()
}
pq <- function(rgSet){
  minfi::preprocessQuantile(rgSet)
}

Em2 <- function(rgSet, arraytype = NULL){
  pd <- rgSet@colData
  qc <- ENmix::QCinfo(rgSet,distplot = F)
  mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
  mSetSqn <- ENmix::mpreprocess(rgSet = rgSet,impute = T)
  if(is.null(arraytype)){
    arraytype <- ifelse(500000 < nrow(mSetSqn), "EPIC", "450K" )
  }   
  
  if(arraytype=="EPIC"){
    mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(
      mat = mSetSqn,
      array = "IlluminaHumanMethylationEPIC", 
      annotation = "ilm10b4.hg19"
    )
  }else mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(mSetSqn)
  
  mSetSqn@colData <- pd
  return(mSetSqn)
}

Em <- function(rgSet, arraytype = NULL){
  pd <- rgSet@colData
  qc <- ENmix::QCinfo(rgSet,distplot = F)
  mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
  mSetSqn <- minfi::mapToGenome(mdat)
  mSetSqn@colData <- pd
  return(mSetSqn)
}

filter<-function(targets, rgSet,sampGroups=NULL,sampNames="Sample_Name",frac=0.1,pval=0.01,remove_sex=TRUE,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
  # requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  # requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
  n <- ncol(rgSet)
  #o <- rev(order(sampNames))
  #rgSet <- rgSet[, o]
  #sampNames <- sampNames[o]
  if (is.null(sampGroups)) g <- rep(1, n) else g <- targets[[sampGroups]]
  if (is.null(g)) g<-1:n
  g <- factor(g)
  
  pal =  c(
    "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
    "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
    "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
    "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
  
  
  # Check quality of the combined probe signal (testing against negative controls)
  # 1. Remove low quality samples. Defaults more than 10% we throw samples away.
  detP <- minfi::detectionP(rgSet, type = "m+u")
  grDevices::pdf(file = paste0(qc_folder,"mean_detection_pvalues.pdf"),width = 7,height = 7)
  ylabels<-colnames(detP)
  par(mar=c(max(4.1,max(nchar(ylabels))/2.2) ,4.1 , 4.1, 2.1))
  
  barplot(colMeans(detP), col=pal[g], las=2,
          cex.names=0.8, ylim=c(0,max(0.002,max(colMeans(detP))*2)), main ="Mean detection p-values")
  graphics::abline(h=0.05,col="red")
  graphics::legend("topleft", legend=levels(g), fill=pal[1:length(levels(g))],
                   bg="white")
  grDevices::dev.off()
  
  # 2. Removing low-quality samples (with fraction of probes not passing pval)
  bad_samples <- colnames(detP)[colSums(detP >=pval)/nrow(detP) > frac]
  if(length(bad_samples)>0){
    warning("The following samples will be discarded since they fail to pass the p-value filter ( ",
            frac*100,"% of the probes with p-val >", pval, "): \n ", paste(bad_samples,collapse = ", " ))
    rgSet <- rgSet[,setdiff(colnames(detP),bad_samples)]
  }else{
    cat("All samples passed detection P-value filter")
  }
  # 3. Removing low-quality probes (with p-value below pval)
  bad_probes<-which(rowSums(detP < pval) < ncol(rgSet)*(1-frac))
  rgSet <- rgSet[-c(bad_probes),]
  if(length(bad_samples)>0){
    warning("The following probes will be discarded since more than", frac*100,
            "% of the samples have detection p-values > ", pval, "): \n ", paste(bad_samples,collapse = ", " ))
  }else{
    cat("All samples passed detection P-value filter")
  }
  return(rgSet)
}



prep<-function(mSetSqn,remove_sex=TRUE,pval=0.01,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
  
  
  # 4. Removing probes with known SNPs at CpG site
  mSetSqn <-  minfi::mapToGenome(mSetSqn)
  mSetSqn <- minfi::dropLociWithSnps(mSetSqn)
  
  # 5. Removing cross reactive probes
  mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)
  
  # 6. Sex. prediction & removal
  mSetSqn$predictedSex <- minfi::getSex(mSetSqn, cutoff = -2)$predictedSex
  if(remove_sex){
    if(!is.null(arraytype)){anno<-get_anno(arraytype)
    }else{
      anno<-minfi::getAnnotation(mSetSqn)
      anno<-anno[!(anno$chr %in% c("chrX","chrY")),]
    }
  }
  return(mSetSqn)
}


# Add info:

# Purity:
#cnv.methyl::purify()

# Celltype:
#1. cellCounts <- FlowSorted.Blood.450k::estimateCellCounts(rgSet)
#2. FlowSorted.Blood

# Copy Number Variation:
#cnv.methyl::Kc_get(ss )

# Sex:
# minfi::get_sex()

# Age:
# 
Enmix<-function(rgSet2){
  qcE<-ENmix::QCinfo(rgSet2)
  mdat<-ENmix::preprocessENmix(rgSet2, bgParaEst="oob", dyeCorr="RELIC",
                               QCinfo=qc, nCores=6)
}

# # obtaining the beta values
# beta_values <- getBeta(gmSet)
# colnames(beta_values) <- metadata$sample
# 
# saveRDS(beta_values, file = "results/beta_values.rds")
# 
# # PRINCIPAL COMPONENT ANALYSIS
# 
# # selecting the top 100 most variable CpG sites
top_beta <- function(beta_values, n=1000){
  sdv <- apply(beta_values, 1, sd)
  top100 <- names(head(sort(sdv,decreasing=T), n))
  beta_top100 <- beta_values[top100,]
  return(beta_top100)
}
pca_res <- function(beta_top100,scale=T, center=T){
  prcomp(t(beta_top100), scale=scale, center=center)
}



corpca <- function(beta_top100,metadata,vars=NULL){
  requireNamespace("PCAtools")
  p<-PCAtools::pca(beta_top100,metadata = metadata, removeVar = 0.1)
  if(is.null(vars)){
    vars<-names(p$metadata)[sapply(p$metadata,function(x){
      !any(is.na(x))& length(unique(x))>1
    })
    ]
    vars[sapply(p$metadata,function(x)length(unique(x)))>1]
  }
  PCAtools::eigencorplot(p,
               components = PCAtools::getComponents(p, 1:6),
               metavars = vars,
               col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
               cexCorval = 0.7,
               colCorval = 'white',
               fontCorval = 2,
               posLab = 'bottomleft',
               rotLabX = 45,
               posColKey = 'top',
               cexLabColKey = 1.5,
               scale = TRUE,
               main = 'PC1-6 clinical correlations',
               colFrame = 'white',
               plotRsquared = FALSE)
  
}


bplot<-function(pca,ss,colgroup,s,combs=NULL, tit= NULL,folder = "analysis/pca/bplots/"){
  library(ggplot2)
  library(gplots)
  library(ggrepel)
  library(ggfortify)
  ss<-droplevels.data.frame(ss)
  pal =  c(
    "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
    "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
    "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
    "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
  lapply(colgroup, function(f) {
    f_len <- length(unique(with(ss,get(f))))
    cols <- pal[1:f_len]
    names(cols)<- unique(with(ss,get(f)))
    if(is.null(combs))combs<-combn(4,2)
    ss<-as.data.frame(ss)
    #rownames(ss)<-ss$Sample_Name
    for(i in 1:dim(combs)[2]){
      if(is.null(tit))tit <- paste0("colored.by.",f, "_shape.", s)
      ap<-ggplot2::autoplot(pca, x=combs[1,i], y=combs[2,i], data = ss, colour=f,shape=s,alpha=0.7,size=1)+
        geom_text_repel(aes(label = Sample_Name, color = with(ss,get(f))),
                        show.legend = FALSE, size = 1.5,max.overlaps = Inf,segment.size=0.2,min.segment.length = 0.8,point.size = 0.5)+
        #scale_color_brewer(palette = "Paired")+
        scale_color_manual(values = cols)+
        #geom_point(aes(size=0.2))+
        labs(colour=f)+ 
        
        ggtitle(tit)+
        theme_bw(base_size = 7)+
        theme(legend.key=element_blank(), legend.key.size=unit(1,"point"))
      
      
      dir.create(folder)
      #plot(ap)
      ggsave(paste0(folder,tit,i,".png"),plot=ap,width = 5.56, height = 2.80,units="in" )
    }
  })
  return(folder)
}

# Surrogate analysis:
surrogate<-function(grset,pheno,condition){
  mval<- getM(grset)
  pheno<-pData(grset)
  mod <- model.matrix(~as.factor(condition), data=pheno)
  mod0 <- model.matrix(~1, data=pheno)
  sva.results <- sva::sva(mval, mod, mod0)
}


# # PCA on the top 100 sites
# pca_res <- prcomp(t(beta_top100), scale=T, center=T)
# pca_all <- prcomp(t(beta_values), scale=T, center=T)
# 
# 
# 
# 
# 
# ## plotting PC1 and PC2 by condition
# autoplot(pca_res, x=1, y=2, data=metadata, colour="vascular_type", shape="type")+
#   geom_text_repel(aes(label=sample, color=vascular_type),hjust=-0.2, vjust=0, show.legend=F, size=3.5)+
#   labs(colour="Tissue", shape="Type")+
#   xlim(c(-0.5,0.3))+
#   theme_bw()+
#   ggtitle("PCA by tissue")
# 

#' Generate models
#' @title construct models and contrasts with limma
#' @param object your object containing beta values
#' @param group_var the variable used as independent variable
#' @param covs the set of variables to use as confounders
#' @param metadata the metadata or sample sheet
#' @param set a boolean vector to subset the observations
#' @param gr the group
#' @return fit2 ebayes model 
#' @author izar de Villasante
#' @export
#'
mod <- function(object, group_var, covs=NULL, metadata,set = TRUE,gr=NULL,pairwise = T,
                singular=F){
  data.table::setDT(as.data.frame(metadata))
  cont_sing=cont_pair=gr_cont_sing=gr_cont_pair=NULL
  metadata<-subset(metadata,set)
  metadata<-droplevels(metadata)
  object <- object[,metadata$barcode]
  covs_formula<-NULL
  if (!is.null(covs)&length(covs)>0)covs_formula<-paste0("+",paste0(covs,collapse=" + ",sep=""))
  design <- model.matrix( 
    formula(
      paste("~ 0 +" , paste0(group_var),covs_formula,sep= " " )
    ),
    data = metadata
  )
  fit <- limma::lmFit(object,design)
  cols <- with(metadata,paste0(group_var, unique(get(group_var))))
  if(pairwise == T){
    cont_pair <- apply(combn(cols,2),2,function(x) paste(x,collapse = "-"))
    
    if(!is.null(gr)) { 
      gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
      gr_cont_pair <- apply(combn(gr_cols,2),2,function(x) paste(x,collapse = "-"))
    }
  }
  
  
  if(singular == T){
    cont_sing<-apply(combn(cols,length(cols)-1),2,function(x){
      var <- setdiff(cols,x)
      group <- contgroup(group_var,levels=x)
      contrast <- paste0(var,"-", group)
      return(contrast)
    } )
    if(!is.null(gr)) { 
      gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
      gr_cont_sing <- apply(combn(gr_cols,length(gr_cols)-1),2,function(x){
        var <- setdiff(gr_cols,x)
        group <- contgroup(group_var,levels=x)
        contrast <- paste0(var,"-", group)
        return(contrast)
      } )
    }
  }
  
  cont <- c(cont_sing,cont_pair)
  gr_cont <- c(gr_cont_sing,gr_cont_pair)
  contMatrix <- limma::makeContrasts(
    contrasts=c(cont,gr_cont),
    levels=colnames(design)
  )
  # rename contrasts:
  # GR: 
  if(!is.null(gr)) colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), gr_cols, gr, vectorize_all = FALSE)
  
  # Singular 1vsMean.
  large <- colnames(contMatrix) %in% c(cont_sing,gr_cont_sing)
  colnames(contMatrix)[large] <- sapply(
    colnames(contMatrix)[large], function(x) paste0("sing_",strsplit(x,"-")[[1]][1]))
  
  # remove group_var prefix:
  colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), group_var, "", vectorize_all = FALSE)
  fit2 <- limma::contrasts.fit(fit, contMatrix)
  fit2 <- limma::eBayes(fit2)
  return(fit2)
}


contgroup<-function(name,levels){
  cols<-levels[grepl(name, levels, fixed = TRUE)]
  l<-length(cols)
  paste0("(",paste(cols,collapse = "+"),")/",l)
}

plotDMP <- function(DMPann,names){
  library(ggplot2)
  DMPresults <- data.frame(table(DMPann[ ,c("Contrast","Type")]))
  # plot DMPs (hypo/hyper)
  g1<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
    geom_bar(position="dodge", stat= "identity")+
    theme_bw()+
    scale_fill_manual(values=c("red", "skyblue")) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(x = "", y = "count", fill='Methylation')+
    ggtitle('Differently methylated probes')
  
  # ggplot2::ggsave(g1,paste0("analysis/DMPplots/",names,".png"))
  
  # plot with facets
  g2<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
    geom_bar(position="dodge", stat= "identity")+
    facet_wrap(.~Type, scales = "free_x") +
    theme_bw()+
    scale_fill_manual(values=c("red", "skyblue")) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(x = "", y = "count", fill='Methylation')+
    ggtitle('Differently methylated probes')
  # ggplot2::ggsave("analysis/DNA_methylation/DMP_p0.01_m0.3.facet.png")
  
  # plot proportion of DMPs in CGI
  DMP_annCGI <- data.frame(DMPann[ ,c("Contrast","Type", "Relation_to_Island")])
  g3<-ggplot2::ggplot(DMP_annCGI, aes(Contrast, fill = Relation_to_Island)) +
    facet_wrap(.~Type, scales = "free_x") +
    geom_bar(position ="fill", width = 0.8) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    ylab("DMPs") +
    xlab("")
  # ggplot2::ggsave("analysis/DNA_methylation/DMP_annCGI.png")
  
  #table(DMPann[ ,c("Relation_to_Island","Type", "Contrast")])
  
  # plot proportion of DMPs in genomic elements
  DMPann$UCSC_RefGene_Group[which(DMPann$UCSC_RefGene_Group == "")] <- "."
  DMPann$UCSC_RefGene_Group_short <- unlist(lapply(strsplit(DMPann$UCSC_RefGene_Group, ";"),'[[', 1))
  
  DMP_annGenomic <- data.frame(DMPann[ ,c("Contrast","Type", "UCSC_RefGene_Group_short")])
  g4<-ggplot2::ggplot(DMP_annGenomic, aes(Contrast, fill = UCSC_RefGene_Group_short)) +
    facet_wrap(.~Type, scales = "free_x") +
    geom_bar(position = "fill", width = 0.8) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(fill = "RefGene") +
    ylab("DMPs") +
    xlab("")
  # ggplot2::ggsave("analysis/DNA_methylation/DMP_annGenomic.png")
  
  # table(DMPann[ ,c("UCSC_RefGene_Group_short","Type", "Contrast")])
  
  return(list(g1,g2,g3,g4))
  
}

find_dmrs<-function(object, model, pcutoff = 0.01, betacutoff = 0.3, min.cpg=5){
  contrasts <- colnames(model$contrasts)
  results<-lapply(contrasts,function(x){
    myAnnotation <- DMRcate::cpg.annotate(
      object = object,                    
      datatype = "array", 
      what = "Beta", 
      analysis.type = "differential", 
      design = model$design, 
      contrasts = TRUE,
      cont.matrix = model$contrasts, 
      coef = x, 
      arraytype = "EPIC"
    )
    
    DMRs <- DMRcate::dmrcate(myAnnotation, 
                             pcutoff = pcutoff,
                             betacutoff= betacutoff,
                             min.cpgs =min.cpg)
    results.ranges <- DMRcate::extractRanges(DMRs)
    results.ranges$Contrast = x
    # data.frame(Contrast = x, results.ranges)
    return(results.ranges)
  }
  )
  
  # return(data.table::rbindlist(results) )
  res<-do.call("c",results)
  return(res)
}

options(future.globals.maxSize= 1891289600)
#gopath(dmrs_ANA,all.cpg=rownames(betas[,!is.na(ss_clean$ANA_dom)]),n=10,ann=ann)->pat
gopath <- function(ranges,all.cpg=NULL,n=20,ann=NULL){
  #require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  requireNamespace(c("S4Vectors","future","future.apply","missMethyl","data.table"))
  future::plan("multisession")
  cont<-unique(ranges$Contrast)
  pathways <- future.apply::future_lapply(cont, function(x){
    
    results.ranges<- ranges[ranges$Contrast == x,]
    gst_go <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
                                                    all.cpg = all.cpg,
                                                    collection = "GO", 
                                                    array.type = "EPIC",
                                                    anno = ann,
                                                    sig.genes = T)
    g<-missMethyl::topGSA(gst, n=n)
    g$Contrast <- x
    g$method <- "GO"
    g
    },
    error = function(e) {
      data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO")
    }
    )
    gst_prom <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
                                                      all.cpg = all.cpg,
                                                      collection = "GO", 
                                                      array.type = "EPIC",
                                                      genomic.features=c("TSS200","TSS1500","1stExon"),
                                                      anno = ann,
                                                      sig.genes = T)
    g<-missMethyl::topGSA(gst, n=n)
    g$Contrast <- x
    g$method <- "GO_prom"
    g
    },
    
    error = function(e) {
      data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO_prom")
    }
    )
    gst_kegg <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
                                                      all.cpg = all.cpg,
                                                      collection = "KEGG", 
                                                      array.type = "EPIC",
                                                      anno = ann,
                                                      sig.genes = T)
    g<-missMethyl::topGSA(gst, n=n)
    g$Contrast <- x
    g$method <- "KEGG"
    g$ONTOLOGY <- NA
    g$TERM<-g$Description
    g$Description<-NULL
    g
    },
    
    error = function(e) {
      data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="KEGG")
    }
    )
    result <- data.table::rbindlist(list(gst_go,gst_prom,gst_kegg),fill=T)
    return(result)
  },future.packages = c("S4Vectors","GenomicRanges"))
  
  
  return(data.table::rbindlist(pathways))
  
}

# cpgs <- GenomicRanges::GRanges(seqnames = anno$chr, 
#                                ranges = IRanges::IRanges(start = anno$pos, 
#                                                          end = anno$pos),
#                                strand = anno$strand,
#                                name = anno$Name)
# 
# overlaps <- GenomicRanges::findOverlaps(cpgs,regions)
# sig.cpg <- cpgs$name[overlaps@from]
# .getGO <- function(){
#   if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
#     stop("org.Hs.eg.db package required but not installed.")
#   egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
#   GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
#   d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
#   GeneID.PathID <- GeneID.PathID[d, ]
#   GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
#                                                       keys=unique(GeneID.PathID$go_id), 
#                                                       columns=c("GOID","ONTOLOGY","TERM"), 
#                                                       keytype="GOID"))
#   go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
#   
#   list(idList=go, idTable=GOID.TERM)
# }
# go <- .getGO()
# result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=go$idList, 
#                   array.type=array.type, plot.bias=plot.bias, 
#                   prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
#                   fract.counts=fract.counts, 
#                   genomic.features = genomic.features,
#                   sig.genes = sig.genes)
# result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
# rownames(result) <- result$GOID

# subset(ss,!is.na(ss$ANA_dom1))
# technical_features <- c("Sample_Name", "organism", "Basename", "barcode")
# 
# 
#  mval<- getM(npq[,!is.na(ss$ANA_dom1)])
#  pheno<-ss[!is.na(ANA_dom1)]
#  mod <- model.matrix(
#    formula(paste(" ~ ANA_dom1 +", paste(batch,sep="+",collapse="+")))
#    ,data = pheno
#  )
#  mod0 <- model.matrix(
#    formula(paste(" ~", paste(batch,sep="+",collapse="+")))
#    ,data = pheno
#  )
#  sva.results <- sva(mval, mod, mod0)
# # 
# design <- model.matrix(
#   formula(paste(" ~", paste(covs,sep="+",collapse="+")))
#   ,data = ss
#   )
# n.sv = sva::num.sv(betas,design,method="leek")
# 
# svobj = sva(rgSet,mod,mod0,n.sv=n.sv)
# 

#