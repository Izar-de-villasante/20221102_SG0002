# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)

library(tarchetypes) # Load other packages as needed. # nolint





# Set target options:
tar_option_set(
  packages = c("tibble","foreach","S4Vectors"), # packages that your targets need to run
  #imports = "cnv.methyl",
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore"# multiprocess or multicore, LSF, SGE, Slurm etc.
)

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)

# Load the R scripts with your custom functions:
#lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)
# source("other_functions.R") # Source other scripts as needed. # nolint
source("R/functions.R")

# produce_data <- function() {
#   expand.grid(samplesheet = c("a", "b"), model = c("c", "d"), normalization = c(1, 2, 3))
# }
# list(
#   tar_group_by(data, produce_data(), samplesheet, model),
#   tar_target(group, data, pattern = map(data))
# )

samplesheets <- c(
  sample_Sheet1="data/ss_H2009.rds",
  sample_Sheet2="data/ss_H358.rds"
)

values <- tibble::tibble( # Use all possible combinations of input settings.
  method_function = rlang::syms(c("noob")),#, "pq","funn","noob_pq","Em")),
  data_paths = samplesheets,
  data_names = c("H2009","H358")
  )



targets <- tar_map(
  values = values,
  names = data_names, #"data_source", # Select columns from `values` for target names.

  tar_target(samplesheet_path, data_paths, format = "file"),
  tar_target(samplesheet, readRDS(samplesheet_path)),
  tar_target(ss,samplesheet[-1,]),
  tar_target(category,samplesheet[1,]),
  tar_target(rgSet, cnv.methyl::read.metharray.exp.par(ss,folder="ana",extended=T)),
  # Qc report:
  tar_target(QC_plots, cnv.methyl::qc(rgSet,sampGroups = "condition",qc_folder = paste0("analysis/intermediate/QC/",data_names,"/"))),

  # Calculate purity:
  tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
  tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="Condition")),

  tar_target(normalize, method_function(filtered)),
  tar_target(clean, prep(normalize)),
  tar_target(ss_clean, droplevels.data.frame( cbind(clean@colData,purity))),
  tar_target(ann, minfi::getAnnotation(clean)),
  tar_target(betas, minfi::getBeta(clean)),
  tar_target(top,top_beta(betas,n=1000)),
  tar_target(pca, pca_res(top)),
  tar_target(pca_corrplot,corpca(top,metadata=ss_clean)),
  tar_target(bplotvars, c("predictedSex",colnames(as.matrix(category[1,]))[as.matrix(category[1,])%in%c("batch","covs")])),
  tar_target(bplots, bplot(pca,ss=ss_clean,col=bplotvars,s="",
                           folder = paste0("analysis/pca/bplots/",data_names,"/"))
             ),

  tar_target(model, mod(object = betas, group_var = "Sample_Group",
                          metadata = ss_clean,gr=c("Control","Case","Treated","Untreated"))
             ),

  tar_target(dmps_mod1, cnv.methyl::DMPextr(fit = model,
                                           ContrastsDM = colnames(model$contrasts),
                                           beta_normalized = betas,
                                           p.value = 0.01,
                                           mDiff = 0.15,
                                           ann = ann,
                                           writeOut = F,
                                           w
                                            )),
  tar_target(dmps_summary_ana,{
    data.table::setDT(dmps_mod1)
    dmps_mod1[,list(Hyper=sum(Type=="Hyper"),Hypo=sum(Type=="Hypo")),by=c("Contrast")]
    }),


  tar_target(dmpplot_mod1, plotDMP(dmps_mod1)),

  tar_target(dmrs, find_dmrs(betas,model,pcutoff = 0.1, betacutoff = 0.2, min.cpg=3)),
  # tar_target(gopath, gopath(dmrs,all.cpg=rownames(betas),n=40,ann=ann)),
  # 
  # tar_target(save_gopath, writexl::write_xlsx(gopath[FDR<0.1,.SD,by="Contrast",
  #                                                            .SDcols=c("TERM","N","DE","FDR","Contrast","SigGenesInSet","method")],
  #                                                 paste0("results/gopath","_noob",".xlsx"))),
  NULL
)
# list(
#   tar_target(samplesheet_path, "data/ss_H358.rds", format = "file"),
#   tar_target(samplesheet, readRDS(samplesheet_path)),
#   # tar_target(ss,samplesheet),
#   tar_target(ss,samplesheet[-1,]),
#   tar_target(category,samplesheet[1,]),
#   tar_target(rgSet, cnv.methyl::read.metharray.exp.par(ss,folder="ana",extended=T)),
#   # Qc report:
#   tar_target(QC_plots, cnv.methyl::qc(rgSet,sampGroups = "condition")),
#
#   # Calculate purity:
#   tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
#   tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="condition")),
#   targets,
#   # tar_target( alldmps,
#   #             list(apply(tidyr::expand_grid(c("dmps_ANA","dmps_SANDRA","dmps_WT"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
#   # )
#
#   #tar_combine(combined_gopath_ANA,
#   #list(apply(tidyr::expand_grid(c("gopath_ANA"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
#   NULL
#   )
#
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
# # subset ss
# BiocManager::install(c("Biobase", "conumee", "GenomicRanges", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "impute", "IRanges", "limma", "maxprobes", "minfi", "SummarizedExperiment"))
