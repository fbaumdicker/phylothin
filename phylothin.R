## PhyloThin: R-script for removing oversampled genomes from phylogenetic tree

# by Franz Baumdicker and Hannah Götsch

# Compile this code using:
# Rscript phylothin.r path_to_folder input_tree (priority_list) (no_PATHd8) (no_clade) (-m number_variable_sites)

###################################################################################################################

## download data-set from e.g. NCBI reference sequence database ###################################################
## use e.g. panX (Ding, Baumdicker, Neher 2018) to make a coalescence/phylogenetic tree ###########################

#### TUNING PARAMETERS ############################################################################################
# The parameters have been optimized based on simulations.
alpha_0 <- 0.5
alpha_1 <- 0.1 
alpha_2 <- 0.25*alpha_1 
alpha_3 <- 0.95
alpha_4 <- 0.9

#### INPUT ########################################################################################################

suppressPackageStartupMessages({
  library(ape) # package for the analysis of phylogenetics and evolution
  library(phangorn) # for nnls.tree() and cophenetic.phylo()
  library(dplyr) # for near()
})

# extract command-line arguments:
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 7) {
  stop("Wrong command line input\n
       Usage is \"Rscript phylothin.r path_to_folder input_tree 
       (priority_list) (no_PATHd8) (no_clade) (-m number_variable_sites)\"", call.=FALSE)
}
basepath <- NULL
input_tree_file <- NULL
prio <- F
pathd8 <- T
no_clade <- F
mutation_sensitive <- F
i <- 1
while (i <= length(args)) {
  if (args[i] == "no_PATHd8") {
    pathd8 <- F # skip PATHd8
    i <- i + 1
  } else if (args[i] == "no_clade") {
    no_clade <- T # no clade-output
    i <- i + 1
  } else if (args[i] == "-m") {
    num_snp <- as.numeric(args[i + 1]) # number of variable sites
    mutation_sensitive <- T
    i <- i + 2
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
    } else { # third positional argument assumed to be priority list
      prio <- T # priority list given
      # load (optional) priority list
      priority_list_file <- args[i]
      priority_list <- read.csv(file.path(basepath, priority_list_file), header=F)
      priority_list[,2] <- as.numeric(priority_list[,2])
      # check if priority list has the right format
      priority_list_check <- priority_list[-1, ]
      if (sum(colSums(!is.na(priority_list_check)) > 0) > 2 || # no more than two non-empty columns
          !is.numeric(priority_list_check[!is.na(priority_list_check[2]),2])){ # second column has numeric entries
        stop("The given priority list is not in the right format.")
      }
    }
    i <- i + 1
  }
}

# sanity check:
if (is.null(basepath)) {
  stop("You need to specify the path to the data folder\n
  Usage is \"Rscript phylothin.r path_to_folder input_tree 
  (priority_list) (no_PATHd8) (no_clade) (-m number_variable_sites)\"", call.=FALSE)
} 
if (is.null(input_tree_file)) {
  stop("Missing argument input_tree\n
  Usage is \"Rscript phylothin.r path_to_folder input_tree 
  (priority_list) (no_PATHd8) (no_clade) (-m number_variable_sites)\"", call.=FALSE)
} 

# load (ultrametric) tree:
input_tree <- read.tree(file.path(basepath, input_tree_file))

if (input_tree$Nnode < 2){ # "cintervals" can not be computed
  stop("Your tree has only one node. PhyloThin can not be used.")
}

if (length(input_tree$tip.label) > 10000) { # clade-output not possible for big trees (>10,000)
  no_clade <- T # no clade-output
  print("Your tree is large (> 10,000 tips), clade-output skipped.")
}

if (pathd8){ # make ultrametric tree with PATHd8 (Britton et al 2007)
  print("start calculating ultrametric tree with PATHd8.")
  # define paths needed in bash-script
  Sys.setenv(PATHd8 = paste(basepath, "/PATHd8", sep = ""))
  Sys.setenv(tree_file = file.path(basepath, input_tree_file))
  Sys.setenv(PATHd8_tree_file = paste(basepath, "/pathd8_", input_tree_file, sep = ""))
  Sys.setenv(um_tree_file = paste(basepath, "/um_", input_tree_file, sep = ""))
  # bash-script for PATHd8
  system('chmod u+x $tree_file') # to get permission
  system('$PATHd8 $tree_file $PATHd8_tree_file')
  system('sed -n "/d8.*;/p" $PATHd8_tree_file > $um_tree_file')
  # ultrametric tree
  um_tree <- read.tree(paste(basepath, "/um_", input_tree_file, sep = ""))
  if(class(um_tree)=="multiPhylo"){um_tree <- um_tree$`d8tree:`}
} else { # skip PATHd8 since requested by input
  print("Making the given tree ultrametric with PATHd8 is skipped.")
  um_tree <- input_tree
}

# CAUTION: lot of memory needed for the following (can be skipped if no priority list is given and no clade-output is needed)
if (is.ultrametric(um_tree) == F){ # make the ultrametric tree really ultrametric (precision/rounding errors)
  um_tree <- nnls.tree(cophenetic(um_tree), um_tree, method = "ultrametric", 
                    rooted = is.rooted(um_tree), trace = 0)
}

# create output folders (if not already exist)
if (!file.exists(file.path(basepath, "phylothinoutput"))){
  dir.create(file.path(basepath, "phylothinoutput"))
  print("output-folder phylothinoutput/ has been created")
}
if (!file.exists(file.path(basepath, "phylothinoutput/check"))){
  dir.create(file.path(basepath, "phylothinoutput/check"))
}

###################################################################################################################

# start PhyloThin ...
if (!no_clade) {
  clades <- data.frame(samples = um_tree$tip.label, clade = rep(NA, length(um_tree$tip.label))) # tip-labels
}
save_tree <- um_tree # save the initial (ultrametric) tree 

###################################################################################################################
#### SCALING ######################################################################################################
###################################################################################################################

## scale the tree to coalescent based on T_i for i in a subset of {2, ..., k}

if (is.binary(um_tree) == F){ # test if the tree is binary
  um_tree <- multi2di(um_tree) # transform tree into binary by adding branches of length zero
}
num_start <- um_tree$Nnode+1 # number of nodes + 1 (m+1; binary tree = number of samples n)
cintervals <- coalescent.intervals(um_tree) # information about coalescent intervals
                                           # (number of lineages, interval lengths, interval count, total depth)
                                           # from the ultrametric phylogenetic tree

if (length(c(which(cintervals$interval.length < 0))) > 0){
  print("Warning: some coalescent intervals are negative? check your tree!")
}

# assume that approx. (at least) k (< m+1) samples don't contribute to oversampling

cs1 <- cumsum(choose(2:um_tree$Nnode,2) * cintervals$interval.length[um_tree$Nnode:2])
f_i <- cs1/(cs1 + choose(3:num_start,2) *
              rev(cumsum(cintervals$interval.length[1:(um_tree$Nnode-1)]))) *
    (2:um_tree$Nnode)/(1:(um_tree$Nnode-1))

kmin <- min(6, ceiling(0.1*num_start+1)) # min(5, 0.1n)
k <- max(kmin, which(f_i < alpha_3)+1) # k > kmin = min(5, 0.1n)
l_from_fi <- max(2, which(f_i < alpha_4)+1) # f_i for i=2,...,n-1 -> which(...+1)
  
# scale_index = which coal. times we use for scaling the tree to coalescent time scale
  
rev_scale_index <- which(cintervals$interval.length > 0) # start with times > 0
# only look at the upper k-1 coal. intervals of the tree (to avoid to "include oversampling"):
rev_scale_index <- rev_scale_index[rev_scale_index >= um_tree$Nnode - k + 2]
scale_index <- um_tree$Nnode - rev_scale_index + 2 # the corresponding indices of T_i
  
if (length(scale_index) == 0){
  scale_index <- 2:kmin
} else if (length(scale_index) < kmin-1){
  while (length(scale_index) < kmin-1){
    neg <- (2:k)[-(scale_index-1)]
    scale_index <- c(scale_index, min(neg))
  }
}

scale_index <- rev(sort(scale_index))
rev_scale_index <- um_tree$Nnode - scale_index + 2

###################################################################################################################

if (sum(choose(scale_index,2) * cintervals$interval.length[rev_scale_index]) == 0) {
  print("Error: scaling parameter is zero! scaling skipped")
  scaling <- 1
  scale_index <- c()
  rev_scale_index <- c()
} else {
  
  scaling <- length(scale_index) / sum(choose(scale_index,2) * cintervals$interval.length[rev_scale_index]) 
  # = 1/(average over chosen coal. times)
  um_tree$edge.length <- um_tree$edge.length * scaling # scale edge length
  
  if (length(scale_index) <= kmin-1) {
    print("Warning: scaling maybe not good enough")
    print(paste(c("average over",length(scale_index),
                  "coalescent times is used for the scaling parameter", 1/scaling,
                  "(T_i for i in {",rev(scale_index),"}). sample size: ", num_start), collapse=" "))
  }
}

if (mutation_sensitive) { # compute mutation rate
  theta <- num_snp/sum(um_tree$edge.length)
  print(paste("mutation rate", theta, 
              "(number of variable sites divided by total branch length of scaled ultrametric tree) computed."))
  print("PhyloThin uses mutation sensitive thinning of the tree.")
}

###### CHECK ######################################################################################################

## plot "coalescent times"
pdf(paste(basepath, "/phylothinoutput/check/treecutting_", tree_name, ".pdf", sep = "")) 
  ctimes <- choose((um_tree$Nnode+1):2,2) * coalescent.intervals(um_tree)$interval.length # scaled coal. times
  plot(ctimes, ylab = "choose(i,2) T_i", xlab = "counting starts from tips")
  cs_ctimes <- cumsum(ctimes) # cumulative sums starting from tips (under neutrality: should follow gamma distr.)
  plot(cs_ctimes, ylab = "cumsum(choose(i,2) T_i)", xlab = "summation starts from tips")
  plot(pgamma(cs_ctimes,1:length(cs_ctimes)), ylab = "pgamma(cumsum(choose(i,2) T_i))", 
       xlab = "summation starts from tips") # gamma distribution functions / cumulative probabilities (vector)
  revcs_ctimes <- cumsum(rev(ctimes)) # cumulative sums for reversed ctimes (starting from root)
  plot(pgamma(revcs_ctimes,1:length(revcs_ctimes)), ylab = "pgamma(cumsum(choose(i,2) T_i))", 
       xlab = "summation starts from root") 
invisible(dev.off())

###################################################################################################################
#### REMOVE SAMPLES ###############################################################################################
###################################################################################################################

removed <- T # control
num_removed <- 0 # number of removed tips/samples
removed_labels <- c() # ids of removed samples

# how many cumulative sums are used for testing:
# cumulative sums (starting from root) of scaled coalescent times (reversed, revcs_ctimes from above)
test_rev <- cumsum(rev(choose((um_tree$Nnode+1):2,2) * coalescent.intervals(um_tree)$interval.length))
testprob_rev <- pgamma(test_rev,1:length(test_rev)) # cumulative prob.
local_length_rev <- max(which(testprob_rev > alpha_0),0)
local_length <- length(test_rev)-local_length_rev

if(length(test_rev) > l_from_fi){
  if(local_length == 0){
    local_length <- length(test_rev)+1-l_from_fi
  }else{
    local_length <- min(local_length, length(test_rev)+1-l_from_fi)
  }
}

print(paste("At maximum", local_length, "coalescent intervals are used for testing."))
#write.table(local_length, file = paste(basepath, "/phylothinoutput/testing_bound_", tree_name, ".txt" ,sep = "" ), quote = F, row.names = F, col.names = F)

if (prio){ # priority list given
  # id of tips don't want to keep (priority 0)
  no_wantkeep_all <- priority_list[priority_list$NCBI_refseq == 0,]$tip_label
  # id of tips with no priority
  nopriority_all <- union(priority_list[is.na(priority_list$priority)]$tip_label,
                          setdiff(um_tree$tip.label, priority_list$tip_label))
}

if (!no_clade) {clade_num <- 1}

# removing algorithm ##################

cutting_step <- function() {
    external_br_index <<- which(um_tree$edge[,2] <= length(um_tree$tip.label)) # find external branches
    tip_index <<- um_tree$edge[,2][external_br_index] # tip indices (of external branches)
    
    # positions of all (usually two) smallest edge length of external branch 
    # -> one of these tips should be removed; they are all in the same clade
    all_min_branch <<- which(near(um_tree$edge.length[external_br_index], min(um_tree$edge.length[external_br_index])))
    
    # adjust clades-dataframe for this clade
    clade_tips <<- um_tree$tip.label[tip_index[all_min_branch]] # ids of clade-tips
    if (!no_clade) { 
      clade_samples <<- clades$samples %in% clade_tips
      clade_clades <<- unique(clades[clade_samples,]$clade)
      if (length(clade_clades) == 1){
        if (is.na(clade_clades)){ # it is a new clade
          clades[clade_samples,]$clade <<- clade_num
          clade_num <<- clade_num +1
        } # otherwise: one already defined clade, nothing to do
      } else if (length(clade_clades[!is.na(clade_clades)]) == 1){ # additional tips are joining an existing clade
        clades[clade_samples,]$clade <<- clade_clades[!is.na(clade_clades)]
      } else { # two or more clades are merging: one new big clade
        clades[clades$clade %in% clade_clades,]$clade <<- clade_num
        clade_num <<- clade_num +1
      }
    }
    
    # priority list given
    if (prio){
      # candidates to remove which have priority 0
      no_wantkeep <<- intersect(clade_tips, no_wantkeep_all) 
      # candidates to remove which have no priority
      nopriority <<- intersect(clade_tips, nopriority_all) 
      if(length(no_wantkeep) > 0){
        min_branch_id <<- no_wantkeep[1]
      } else if (length(nopriority) > 0){ # there are samples with no priority
        min_branch_id <<- nopriority[1] # id of one sample with no priority
      } else { # sample with "highest" priority (gets removed):
        match_row <<- which.max(priority_list[priority_list$tip_label %in% clade_tips,]$priority) # row in priority list
        min_branch_id <<- priority_list[priority_list$tip_label %in% clade_tips,]$tip_label[match_row] # id of sample with "highest" priority
      }
    } else { # no priority list: remove one tip
      min_branch_id <<- clade_tips[1]
    }
    
    removed_labels <<- c(removed_labels, min_branch_id) # add id of to remove tip
    um_tree <<- drop.tip(um_tree, min_branch_id) # remove tip from tree
    
    num_removed <<- num_removed + 1 
    local_length <<- local_length - 1
}
#######################################

while(removed && num_start - num_removed > 2 && local_length > 0){ # ensure to have at least two tips/samples
  
  if (!mutation_sensitive) { # basic PhyloThin thinning
    
    # cumulative sums (starting from tips) of scaled coalescent times (cs_ctimes from above)
    test <- cumsum(choose((um_tree$Nnode+1):2,2) * coalescent.intervals(um_tree)$interval.length)
    # (negativ coalescent intervals don't matter here (zero and negative have both prob zero))
    
    myprob <- pgamma(test,1:length(test))[1:local_length] # cumulative prob. of interest
    
    # remove tip if sum(T_i[(n-l):n]) (Gamma[l+1,1], l < l_max) or T_n (Exp[1]) are much smaller then expected:
    if (min(myprob) < alpha_2 || myprob[1] < alpha_1){ 
      cutting_step()
    } else {
      removed <- F # desired distr. reached; no further tip should be removed, stop the while loop
    }
    
  } else { # mutation sensitive thinning
    
    # cumulative sums (starting from tips) of coalescent times
    test <- cumsum(coalescent.intervals(um_tree)$interval.length)
    test[test < 0] <- 0 # set negativ (due to rounding errors?) coalescent intervals to zero

    geom_p <- 1/(1+2*theta/(((um_tree$Nnode+1):2)*(((um_tree$Nnode+1):2)-1))) # parameters of geometric distribution (starting from the tips)
    # geom_p smaller for higher theta & closer to root (2)
    myprob <- pgeom(test,geom_p)[1:local_length]
    # mass of geom distr closer to zero for higher geom_p
    
    # geom_max: until when geom. distr. should be computed for convolution (larger = better)
    geom_p_min <- 1/(1+2*theta/(((um_tree$Nnode+1):2)*(((um_tree$Nnode+1):2)-1)))[local_length] # smallest used parameter
    geom_max <- 0
    while (dgeom(geom_max,geom_p_min) > 1e-20) { # (smaller = better)
      geom_max <- geom_max +1
    }
    
    if (myprob[1] < alpha_1) { # last coal interval too small
      print("Test on last coalescent interval has been successful.")
      cutting_step() # cutting
    } else { # cumsum probably too small
      testing_sum <- T # control
      j <- 1
      geom_conv <- dgeom(0:geom_max,geom_p[j])
      while (testing_sum) {
        j <- j+1
        # Fourier transform to compute convolution
        geom_conv <- convolve(geom_conv, rev(dgeom(0:geom_max,geom_p[j])), type = "open") 
        geom_conv <- geom_conv[1:(geom_max+1)] # reduce length back to geom_max
        # cumulative probabilities (cdf)
        conv_cdf <- cumsum(geom_conv)
        if (conv_cdf[floor(test[j])+1] < alpha_2) {
          print(paste("Test on sum of", j, "last coalescent interval has been successful."))
          cutting_step() # cutting
          testing_sum <- F
        } 
        if (testing_sum && j == local_length) {
          testing_sum <- F
          removed <- F # desired distr. reached; no further tip should be removed, stop the while loop
        }
      }
    }
  }
}

print(paste("PhyloThin removed", num_removed, "of", num_start, "samples."))

if (num_start - num_removed == 2){
  print(paste("Warning: only two tips left. sample size: ", num_start))
}

# check and format clades-dataframe
if (!no_clade) { 
for (i in unique(clades[!is.na(clades$clade),]$clade)){
  clade_ids <- subset(clades$samples, clades$clade == i)
  clade_keeper <- setdiff(clade_ids, removed_labels) # keepers in the current clade
  if (length(clade_keeper) == 0){ # no keeper in the clade
    print("Error: Something went wrong in computing the clades. (there exists a clade with no keeper)")
  } else if(length(clade_keeper) > 1){ # more than one keeper in the clade
    clade_tree <- drop.tip(save_tree, tip = setdiff(save_tree$tip.label,clade_ids)) # subtree of the current clade
    distmatrix <- cophenetic.phylo(clade_tree) # distance-matrix; CAUTION: does not work for big trees (>10,000)
    distmatrix_red <- distmatrix[clade_keeper,setdiff(clade_ids,clade_keeper)] # relevant distances
    clade_removed_more_keeper <- c()
    if (length(setdiff(clade_ids,clade_keeper)) == 1){ # only one removed tip in clade
      if (length(unique(distmatrix_red)) > 1){ # remove "unnecessary" keeper from clade
        clades[clades$samples %in% names(distmatrix_red[distmatrix_red != min(distmatrix_red)]),]$clade <- NA
        #print("Error: Something went wrong in computing the clades.")
      }
    } else {
    if (any(apply(distmatrix_red, 2, function(x) sum(x == min(x)) > 1))){
      clade_removed_more_keeper <- names(which(apply(distmatrix_red, 2, function(x) sum(x == min(x)) > 1)))
      distmatrix_red_real <- distmatrix_red
      distmatrix_red <- distmatrix_red[,-which(apply(distmatrix_red, 2, function(x) sum(x == min(x)) > 1))]
    }
    if (!is.matrix(distmatrix_red)) { # gives in which row(s) the min distance is (if only one column left)
      min_dist <- which(distmatrix_red == min(distmatrix_red))
    } else { # gives in which row(s) the min distance is
      min_dist <- apply(distmatrix_red, 2, FUN = function(x) which(x == min(x)))
    }
    for (j in clade_keeper){
      # removed tips with shortest distance to this keeper
      clade_removed <- names(which(min_dist == which(rownames(distmatrix_red) == j)))
      if (length(clade_removed) == 0){ # keeper is currently alone in clade
        clades[clades$samples == j,]$clade <- "assign"
      } else {
        clades[clades$samples == j,]$clade <- clade_num # assign keeper to new clade
        clades[clades$samples %in% clade_removed,]$clade <- clade_num # assign removed tips with shortest dist to same clade
        clade_num <- clade_num +1
      }
    }
    if (length(clade_removed_more_keeper) > 0){ # removed tips with more than one keeper in clade
      for (l1 in clade_removed_more_keeper){
        l2 <- which(colnames(distmatrix_red_real) == l1)
        more_keeper <- names(which(distmatrix_red_real[,l2] == min(distmatrix_red_real[,l2]))) # closest keepers for removed one
        more_keeper_clades <- unique(clades[clades$samples %in% more_keeper,]$clade)
        if ("assign" %in% more_keeper_clades){
          if (length(more_keeper_clades) == 1){ # no clade assigned to keepers -> start new clade
            clades[clades$samples == l1,]$clade <- clade_num
            clades[clades$samples %in% more_keeper,]$clade <- clade_num
            clade_num <- clade_num +1
          } else { # assign clade to unassigned keepers and removed tip
            clades[clades$samples == l1,]$clade <- more_keeper_clades[more_keeper_clades != "assign"][1]
            clades[clades$samples %in% more_keeper & clades$clade == "assign",]$clade <- more_keeper_clades[more_keeper_clades != "assign"][1]
          }
          more_keeper_clades <- unique(clades[clades$samples %in% more_keeper,]$clade)
        }
        if (length(more_keeper_clades) == 1){ # keepers are already in the same clade
          clades[clades$samples == l1,]$clade <- more_keeper_clades[1] # assign removed tip to clade
        } else { # join clades and add removed tip
          clades[clades$clade %in% more_keeper_clades[-1],]$clade <- more_keeper_clades[1]
          clades[clades$samples == l1,]$clade <- more_keeper_clades[1]
        }
      }
    }
    if (length(which(clades$clade == "assign")) > 0){
      print("Error: Something went wrong in computing the clades.")
      clades[clades$clade %in% c("assign"),]$clade <- NA # assign keeper to no clade (should not happen)
    }
  }}
}

# "rename" clade numbers to lowest possible
clade_names <- sort(as.numeric(unique(clades[!is.na(clades$clade),]$clade)))
clades$clade <- match(clades$clade, clade_names)

# check if every removed sample is assigned to a clade
if (length(setdiff(removed_labels, clades[!is.na(clades$clade),]$samples)) > 0){
  print("Warning: Something went wrong with the clade-ouput! Not all removed samples are assigned to a clade.")
}
# check if in every clade is exactly one keeper
clade_keeper_all <- setdiff(clades[!is.na(clades$clade),]$samples, removed_labels)
if (length(clade_keeper_all) < length(unique(subset(clades$clade, clades$samples %in% clade_keeper_all)))){
  print("Warning: Maybe something went wrong with the clade-ouput. There are less keepers in the clades than different clades.")
} else if (length(clade_keeper_all) > length(unique(subset(clades$clade, clades$samples %in% clade_keeper_all)))){
  print("Warning: Maybe something went wrong with the clade-ouput. There are more keepers in the clades than different clades.")
}
}

###### CHECK ######################################################################################################

## visualization of k, scale_index, scaling and num_removed
pdf(paste(basepath, "/phylothinoutput/check/scaling_", tree_name, ".pdf", sep = "")) # export
plot(2:(cintervals$interval.count+1), 
     choose(2:(cintervals$interval.count+1),2) * cintervals$interval.length[(cintervals$interval.count):1],
     pch = 1, col = "black", ylab = "choose(i,2) T_i", xlab = "index i") 
if((cintervals$interval.count-k+1)>0){
  fi_alpha3_points <- which(cintervals$interval.length[(cintervals$interval.count-k+1):1]>0)+k
  points(fi_alpha3_points, 
         choose(fi_alpha3_points,2) * 
           (cintervals$interval.length[(cintervals$interval.count-k+1):1])[fi_alpha3_points-k],
         pch = 16, col = "black")
}
if(length(scale_index)>0){
points(scale_index, 
       choose(scale_index,2) * cintervals$interval.length[rev_scale_index],
       pch = 16, col = "blue") }
abline(h = 1/scaling, col = "blue", lty=1)
abline(v = local_length_rev+1.5, col = "darkgreen", lty=1)
abline(v = l_from_fi+0.5, col = "green", lty=1)
abline(v = num_start-num_removed+0.5, col = "red", lty=1)
if (num_removed > 0) {
  points((num_start-num_removed+1):num_start, 
         choose((num_start-num_removed+1):num_start,2) * 
           cintervals$interval.length[num_removed:1],
         pch = 4, col = "red")
}
legend( x="topright", 
        legend=c("T_i=0","T_i used for scaling",paste("scaling factor:",1/scaling),
                 "bound defined through f_i","bound defined through pgamma", "removed"), 
        col=c("black","blue","blue","green","darkgreen","red"), lwd=1,
        lty=c(0,0,1,1,1,1), pch=c(1,16,NA,NA,NA,4), merge=FALSE )
invisible(dev.off()) # close the plotting device

###################################################################################################################
#### OUTPUT #######################################################################################################
###################################################################################################################

# oversampled clades
if (!no_clade) { 
if(length(removed_labels)>0){
  write.csv(clades, file = paste(basepath, "/phylothinoutput/clades_", tree_name, ".csv" ,sep = "" ))
}
}

# removed sample-ids
if(length(removed_labels)>0){
  write.table(removed_labels, file = paste(basepath, "/phylothinoutput/removed_ids_", tree_name, ".txt" ,sep = "" ),
              quote = F, row.names = F, col.names = F)
}

# kept sample-ids
write.table(um_tree$tip.label, file = paste(basepath, "/phylothinoutput/kept_ids_", tree_name, ".txt" ,sep = "" ),
            quote = F, row.names = F, col.names = F)

# tree with only kept samples
dropped_tree <- drop.tip(input_tree, tip = removed_labels) # tree where tips removed
write.tree(dropped_tree, file = paste(basepath, "/phylothinoutput/reduced_tree_", tree_name, ".nwk", sep = ""))

# comparison of tree before and after removing samples
removed_ones <- which(is.element(save_tree$tip.label, removed_labels))
if (num_start - num_removed == 2){
  warning = "too many tips removed"
}else{
  warning = ""
}
text1 = paste(num_removed, "of", length(save_tree$tip.label), "removed (", round(100*num_removed/length(save_tree$tip.label)), "%).")
text2 = paste(length(save_tree$tip.label)-num_removed, "remaining tips.")

# ultrametric tree:
pdf(paste(basepath, "/phylothinoutput/check/um_treecomparison_", tree_name, ".pdf", sep = ""))
  par(mfrow = c(1,2)) # create a 1 x 2 plotting matrix
  plot(save_tree, show.tip.label = F, main = warning, sub = tree_name) # plot initial tree
  tiplabels(tip = removed_ones, col = "red" , pch = 4) # add red X at tips (initial tree) which have been removed
  plot(drop.tip(save_tree, tip = removed_labels), show.tip.label = F,  # plot final tree
     main = text1, sub = text2)
invisible(dev.off())
# input-tree:
pdf(paste(basepath, "/phylothinoutput/treecomparison_", tree_name, ".pdf", sep = ""))
  par(mfrow = c(1,2))
  plot(input_tree, show.tip.label = F, main = warning, sub = tree_name)
  tiplabels(tip = removed_ones, col = "red" , pch = 4)
  #tiplabels(tip = which(is.element(save_tree$tip.label, clades[!is.na(clades$clade),]$samples)), col = "blue" , pch = 1) # tips which belong to a oversampled clade
  plot(dropped_tree, show.tip.label = F, main = text1, sub = text2)
invisible(dev.off())

# full input-tree
pdf(paste(basepath, "/phylothinoutput/tree_", tree_name, ".pdf", sep = ""))
  plot(input_tree, show.tip.label = T, 
       main = paste(num_removed, "of", length(save_tree$tip.label), "removed. ", warning),
       sub = tree_name)
  tiplabels(tip = removed_ones, col = "red" , pch = 4)
invisible(dev.off())

print("PhyloThin done.")

###################################################################################################################

