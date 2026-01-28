## subsampling-PhyloThin: R-script for subsampling large phylogenetic trees and applying PhyloThin

# by Hannah Götsch

# Compile this code using:
# Rscript phylothin_subsampling.r path_to_folder input_tree (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)

#### TUNING PARAMETER #############################################################################################

# proportion on how often a sample has to be in a oversampling-clade such that it gets classified as oversampled
alpha_5 <- 0.9 

#### INPUT ########################################################################################################
suppressPackageStartupMessages(
library(ape) # package for the analysis of phylogenetics and evolution
)

# extract command-line arguments:
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 7) {
  stop("Wrong command line input\n
       Usage is \"Rscript phylothin_subsampling.r path_to_folder input_tree 
       (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)\"", call.=FALSE)
}
default <- 1
runs <- NULL
subsamplesize <- NULL
basepath <- NULL
input_tree_file <- NULL
pathd8 <- 1
i <- 1
while (i <= length(args)) {
  if (args[i] == "-r") {
    runs <- as.numeric(args[i + 1]) # number of subsamples
    default <- 0
    i <- i + 2
  } else if (args[i] == "-s") {
    subsamplesize <- as.numeric(args[i + 1]) # size of subsamples
    default <- 0
    i <- i + 2
  } else if (args[i] == "no_PATHd8") {
    pathd8 <- 0 # skip PATHd8
    i <- i + 1
  } else {
    if (is.null(basepath)) { # first positional argument assumed to be basepath
      basepath <- args[i]
      if (substr(basepath, nchar(basepath), nchar(basepath)) == "/") { # path should NOT end with /
        basepath <- substr(basepath, 1, nchar(basepath)-1)
      }
    } else if (is.null(input_tree_file)){ # second positional argument assumed to be input tree
      input_tree_file <- args[i]
      if (substr(input_tree_file, nchar(input_tree_file)-3, nchar(input_tree_file)) == ".nwk") {
        tree_name <- substr(input_tree_file, 1, nchar(input_tree_file)-4)
      } else {
        tree_name <- input_tree_file
      }
    }
    i <- i + 1
  }
}

# sanity check:
if (is.null(basepath)) {
  stop("You need to specify the path to the data folder\n
       Usage is \"Rscript phylothin_subsampling.r path_to_folder input_tree 
       (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)\"", call.=FALSE)
} 
if (is.null(input_tree_file)) {
  stop("Missing argument input_tree\n
       Usage is \"Rscript phylothin_subsampling.r path_to_folder input_tree 
       (-r number_of_subsamples) (-s subsample_size) (no_PATHd8)\"", call.=FALSE)
} 

# load (ultrametric) tree:
input_tree <- read.tree(file.path(basepath, input_tree_file))
save_full_tree <- input_tree
num_sample_full <- length(save_full_tree$tip.label) # sample size

# sanity check:
if (input_tree$Nnode < 2){ # "cintervals" can not be computed
  stop("Your tree has only one node. PhyloThin can not be used.")
}

if (!file.exists(file.path(basepath, "phylothinoutput"))){ # output-folder for PhyloThin
  dir.create(file.path(basepath, "phylothinoutput"))
}

###################################################################################################################
#### STEP 1: PHYLOTHIN ON FULL TREE ###############################################################################
###################################################################################################################

print("start PhyloThin on full input tree.")

Sys.setenv(phylothin = paste(basepath, "/phylothin.r", sep = ""))
Sys.setenv(path_to_folder = basepath)
Sys.setenv(input_tree = input_tree_file)

if (pathd8==0){ # skip PATHd8 since requested by input
  system('Rscript $phylothin $path_to_folder $input_tree no_PATHd8 no_clade')
} else{ # make ultrametric tree with PATHd8 (Britton et al 2007)
  system('Rscript $phylothin $path_to_folder $input_tree no_clade')
}

# load reduced (ultrametric) tree:
input_tree <- read.tree(paste(basepath, "/phylothinoutput/reduced_tree_", 
                              tree_name, ".nwk", sep = ""))
save_red_tree <- input_tree
num_sample_red <- length(save_red_tree$tip.label) # reduced sample size

print(paste("PhyloThin on full input tree removed", num_sample_full-num_sample_red, "samples."))

###################################################################################################################
#### STEP 2: SUBSAMPLING ##########################################################################################
###################################################################################################################

print("start subsampling on (by phylothin) reduced input tree.")

# sanity check:
if (input_tree$Nnode < 2){ # "cintervals" can not be computed
  stop("Your by PhyloThin reduced tree has only one node. Subsampling can not be used.")
}

#### DEFINE PARAMETERS ############################################################################################

if (pathd8==0){ # skip PATHd8 since requested by input
  print("Making the reduced tree ultrametric with PATHd8 is skipped.")
  um_tree <- input_tree
} else{ # make ultrametric tree with PATHd8 (Britton et al 2007)
  print("start calculating ultrametric tree for subsampling with PATHd8.")
  # define paths needed in bash-script
  Sys.setenv(PATHd8 = paste(basepath, "/PATHd8", sep = ""))
  Sys.setenv(tree_file = paste(basepath, "/phylothinoutput/reduced_tree_", tree_name, ".nwk", sep = ""))
  Sys.setenv(PATHd8_tree_file = paste(basepath, "/phylothinoutput/pathd8_reduced_tree_", input_tree_file, sep = ""))
  Sys.setenv(um_tree_file = paste(basepath, "/phylothinoutput/um_reduced_tree_", input_tree_file, sep = ""))
  # bash-script for PATHd8
  system('chmod u+x $tree_file') # to get permission
  system('$PATHd8 $tree_file $PATHd8_tree_file')
  system('sed -n "/d8.*;/p" $PATHd8_tree_file > $um_tree_file')
  # ultrametric tree
  um_tree <- read.tree(paste(basepath, "/phylothinoutput/um_reduced_tree_", input_tree_file, sep = ""))
  if(class(um_tree)=="multiPhylo"){um_tree <- um_tree$`d8tree:`}
}

# compute default setting
if (default == 1) {
  # compute default subsamplesize (on by PhyloThin reduced tree)
  if (is.null(subsamplesize)) {
    external_br_index <- which(um_tree$edge[,2] <= num_sample_red) # find external branches
    external_br_length <- um_tree$edge.length[external_br_index] # length of external branches
    if(length(which(external_br_length == 0))>0){
      external_br_length <- external_br_length[-which(external_br_length == 0)] # external branch length > 0
    }
    subsamplesize <- round(sqrt(coalescent.intervals(um_tree)$total.depth/min(external_br_length)))
    # sanity check:
    if (subsamplesize > num_sample_full){
      stop(paste("Error: the default subsample size", subsamplesize, "is larger than the (reduced) sample size", 
                 num_sample_red, ". Skip the subsampling procedure of PhyloThin (maybe not needed for your tree) 
                or define your own subsample size: 'Rscript phylothin_subsampling.r path_to_folder input_tree 
                (-r number_of_subsamples) -s subsample_size (no_PATHd8)'."))
    }
  }
  # compute default number of runs
  if (is.null(runs)) {
    runs <- round(100*num_sample_red/subsamplesize) # on average every tip will be sampled 100 times
  }
}
# sanity check:
if (subsamplesize > num_sample_red){
  stop(paste("Error: subsample size", subsamplesize, "has to be smaller than (reduced) sample size", num_sample_red, "."))
}
print(paste("The following parameter setting is used for subsampling:", 
            runs, "drawings with a subsample size of", subsamplesize, "."))

if (!file.exists(file.path(basepath, "phylothinoutput/subsampling"))){ # folder for sub-sampling output
  dir.create(file.path(basepath, "phylothinoutput/subsampling"))
}

#### SUBSAMPLING ##################################################################################################

subsample_table <- data.frame(sample = um_tree$tip.label, 
                              num_sampled = rep(0, num_sample_red), 
                              num_in_cluster = rep(0, num_sample_red))
write.csv(subsample_table, file = paste(basepath, "/phylothinoutput/subsample_table_", tree_name, ".csv", sep = "" ),
          quote = F, row.names = F)

for (i in 1:runs) {
  print(paste("%%%%%%% START Subsampling", i, "%%%%%%%"))
  set.seed(i)
  nosample_tips <- um_tree$tip.label[sample(1:num_sample_red, num_sample_red-subsamplesize)]
  sample_tree <- drop.tip(um_tree, nosample_tips)
  write.tree(sample_tree, file = paste(basepath, "/phylothinoutput/subsampling/um_", tree_name, "_", i, ".nwk", sep = ""))

###################################################################################################################
#### STEP 3: PHYLOTHIN ON SUBSAMPLES ##############################################################################
###################################################################################################################
  
  Sys.setenv(phylothin = paste(basepath, "/phylothin.r", sep = ""))
  Sys.setenv(path_to_folder = paste(basepath, "/phylothinoutput/subsampling", sep = ""))
  Sys.setenv(input_tree = paste("um_", tree_name, "_", i, ".nwk", sep = ""))
  system('Rscript $phylothin $path_to_folder $input_tree no_PATHd8')
  
  # # remove unecessary output
  # Sys.setenv(remove_file = paste(basepath, "/phylothinoutput/subsampling/phylothinoutput/k*", sep = ""))
  # system('rm $remove_file')
  # Sys.setenv(remove_file = paste(basepath, "/phylothinoutput/subsampling/phylothinoutput/r*", sep = ""))
  # system('rm $remove_file')
  # Sys.setenv(remove_file = paste(basepath, "/phylothinoutput/subsampling/phylothinoutput/t*", sep = ""))
  # system('rm $remove_file')
  # Sys.setenv(remove_file = paste(basepath, "/phylothinoutput/subsampling/phylothinoutput/check/*", sep = ""))
  # system('rm $remove_file')

###################################################################################################################
#### OUTPUT #######################################################################################################
###################################################################################################################

  subsample_table <- read.csv(paste(basepath, "/phylothinoutput/subsample_table_", tree_name, ".csv", sep = "" ), header=T)
  for (subsample_sample in setdiff(um_tree$tip.label, nosample_tips)) {
    subsample_table[subsample_table$sample == subsample_sample,]$num_sampled <- 1 + 
      subsample_table[subsample_table$sample == subsample_sample,]$num_sampled
  }
  subtree_name <- paste("um_", tree_name, "_", i, sep = "")
  if (file.exists(paste(basepath, "/phylothinoutput/subsampling/phylothinoutput/clades_", subtree_name, ".csv", sep = "" ))) {
    clades_file <- read.csv(paste(basepath, "/phylothinoutput/subsampling/phylothinoutput/clades_", subtree_name, ".csv", sep = "" ),
                            header=T)
    clades_tips <- clades_file[!is.na(clades_file$clade),]$samples
    for (clades_sample in clades_tips) {
     subsample_table[subsample_table$sample == clades_sample,]$num_in_cluster <- 1 + 
        subsample_table[subsample_table$sample == clades_sample,]$num_in_cluster
    }
  }
  write.csv(subsample_table, file = paste(basepath, "/phylothinoutput/subsample_table_", tree_name, ".csv", sep = "" ),
            quote = F, row.names = F)
}

# Sys.setenv(remove_path = paste(basepath, "/phylothinoutput/subsampling/phylothinoutput/check", sep = ""))
# system('rmdir $remove_path')

subsample_table <- read.csv(paste(basepath, "/phylothinoutput/subsample_table_", tree_name, ".csv", sep = "" ), header=T)
subsample_table$oversampled <- subsample_table$num_in_cluster/subsample_table$num_sampled
write.csv(subsample_table, file = paste(basepath, "/phylothinoutput/subsample_table_", tree_name, ".csv", sep = "" ),
          quote = F, row.names = F)

# samples classified as oversampled
if (file.exists(paste(basepath, "/phylothinoutput/removed_ids_", tree_name, ".txt", sep = ""))){
  removed_labels_full <- read.table(paste(basepath, "/phylothinoutput/removed_ids_", 
                                          tree_name, ".txt", sep = ""))$V1  # removed tips by phylothin on full tree
} else {
  removed_labels_full <- c()
}
removed_labels_sub <- subsample_table[subsample_table$oversampled >= alpha_5,]$sample # removed tips after subsampling
removed_labels <- c(removed_labels_full,removed_labels_sub)
num_removed <- length(removed_labels)

print(paste("subsampling removed another", num_removed-(num_sample_full-num_sample_red), "tips. In total", 
            num_removed, "tips have been removed."))

# removed sample-ids
if(num_removed>0){
  write.table(removed_labels, 
              file = paste(basepath, "/phylothinoutput/removed_ids_", tree_name, "_subsampling.txt" ,sep = "" ),
              quote = F, row.names = F, col.names = F)
}

# kept sample-ids
write.table(setdiff(save_full_tree$tip.label,removed_labels), 
            file = paste(basepath, "/phylothinoutput/kept_ids_", tree_name, "_subsampling.txt" ,sep = "" ),
            quote = F, row.names = F, col.names = F)

# tree with only kept samples
dropped_tree <- drop.tip(save_full_tree, tip = removed_labels) # tree where tips removed
write.tree(dropped_tree, 
           file = paste(basepath, "/phylothinoutput/reduced_tree_", tree_name, "_subsampling.nwk", sep = ""))

# comparison of tree before and after removing samples
if (num_sample_full - num_removed == 2){
  warning = "too many tips removed"
}else{
  warning = ""
}
# input-tree:
pdf(paste(basepath, "/phylothinoutput/treecomparison_", tree_name, "_subsampling.pdf", sep = ""))
par(mfrow = c(1,3))
plot(save_full_tree, show.tip.label = F, main = warning, sub = tree_name)
removed_ones <- which(is.element(save_full_tree$tip.label, removed_labels_full))
tiplabels(tip = removed_ones, col = "red" , pch = 4)
plot(save_red_tree, show.tip.label = F, sub = "after applying PhyloThin once", col.sub = "red")
removed_ones <- which(is.element(save_red_tree$tip.label, removed_labels_sub))
tiplabels(tip = removed_ones, col = "blue" , pch = 4)
plot(dropped_tree, show.tip.label = F, 
     main = paste(num_removed, "of", num_sample_full, "removed" ),
     sub = "after subsampling", col.sub = "blue")
invisible(dev.off())

###################################################################################################################