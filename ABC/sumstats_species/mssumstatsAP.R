#' Calculates several summary statistics from microsats
#'
#' @param data microsatellites, whereby every locus is represented by two adjacent columns. Or a gtypes object from strataG
#' @param by_pop name of population variable. If specified, all summary statistics will be calculated within populations
#' @param start_geno integer, specifying the first column with genotypes. If now specified, all column are expected to be genotypes
#' @param mratio defaults to "strict". if "loose", the mratio will be calculated differently, see ?m_ratio
#' @param rarefaction if TRUE, calculates mean and sd number of alleles as mean of n_samp individuals and n_loc loci over nresamp bootstraps
#' @param nsamp number of samples to subssample
#' @param nloc number of loci to subsample
#' @param nboot number of bootstraps
#'
#' @examples
#'
#'
#' data(fur_seal)
#' out <- mssumstats(fur_seal, by_pop = "pop", start_geno = 4, mratio = "loose", rarefaction = TRUE,
#' nresamp = 10, nind = 20, nloc = 5)
#'
#' @export
#'
#'

mssumstatsAP <- function(data, by_pop = NULL, start_geno = NULL,
                       mratio = c("strict", "loose"), rarefaction = FALSE, nresamp = 1000, nind = NULL, nloc = NULL) {
  
  # define how to calculate MRatio
  if (length(mratio) == 2) mratio  <- mratio[1]
  # start of actual genotypes
  if (is.null(start_geno)) start_geno <- 1
  
  genotypes <- data
  
  # function to calculate summary statistics from a genotype data.frame, where
  # the genotypes start with the 1st column
  calculate_sumstats <- function(genotypes){
    
    # transform to gyptes object
    g_types_geno <- strataG::df2gtypes(as.data.frame(genotypes),
                                       ploidy = 2, id.col = NULL, strata.col = NULL, loc.col = 1)
    
    # calc summary statistics
    # num_alleles, allel_richness, prop_unique_alleles, expt_het, obs_het
    # mean and sd
    num_alleles <- strataG::numAlleles(g_types_geno)
    num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
    num_alleles_sd <- sd(num_alleles, na.rm = TRUE)
    # exp_het
    exp_het <- strataG::exptdHet(g_types_geno)
    exp_het_mean <- mean(exp_het, na.rm = TRUE)
    exp_het_sd <- sd(exp_het, na.rm = TRUE)
    # obs_het
    obs_het <- strataG::obsvdHet(g_types_geno)
    obs_het_mean <- mean(obs_het, na.rm = TRUE)
    obs_het_sd <- sd(obs_het, na.rm = TRUE)
    
    # allele frequencies
    afs <- strataG::alleleFreqs(g_types_geno, by.strata = FALSE)
    # in case anything is NA
    if (any(is.na(names(afs)))) afs <- afs[-which(is.na(names(afs)))]
    
    # prop low frequency alleles
    prop_low_af <- function(afs){
      # low_afs <- (afs[, "freq"] / sum(afs[, "freq"])) < 0.05
      low_afs <- afs[, "prop"] <= 0.05
      prop_low <- sum(low_afs) / length(low_afs)
    }
    # and mean/sd for all
    prop_low_afs <- unlist(lapply(afs, prop_low_af))
    prop_low_afs_mean <- mean(prop_low_afs, na.rm = TRUE)
    prop_low_afs_sd <- stats::sd(prop_low_afs, na.rm = TRUE)
    
    # for empirical data (with repeat sizes being not 1 as in the simulated microsimr data)
    # find the most common repeat size per locus
    rpt_size <- 2:8
    freqs <- strataG::alleleFreqs(g_types_geno)
    
    # what is the (most likely) repeat size?
    allele_sizes <- lapply(freqs, function(x) sort(as.numeric(row.names(x))))
    allele_size_diffs <- lapply(allele_sizes, function(x) diff(x))
    
    find_repeat_size <- function(size_diff, rpt_size){
      
      all_possible <- matrix(data = NA, nrow = length(rpt_size), ncol = length(size_diff))
      count <- 1
      for (r in rpt_size) {
        all_possible[count, ] <- size_diff%%r == 0
        count <- count + 1
      }
      # instead of checking whether all alleles stick to a certain repeat size take the most
      # common repeat size
      r <- rpt_size[which.max(rowSums(all_possible))]
    }
    repeat_size_per_locus <- as.numeric(lapply(allele_size_diffs, find_repeat_size, rpt_size))
    # allele range
    allele_range <- unlist(lapply(afs, function(x) diff(range(as.numeric(row.names(x))), na.rm = TRUE)))
    # delete loci where range is 0
    to_keep <- !(allele_range < 2)
    allele_range <- allele_range[ to_keep ]
    rep_per_loc <- repeat_size_per_locus[ to_keep ]
    mean_allele_range <- mean(allele_range /   rep_per_loc, na.rm = TRUE)
    sd_allele_range <- sd(allele_range /   rep_per_loc, na.rm = TRUE)
    
    # allele size variance and kurtosis
    # create vector of all alleles per locus
    all_alleles <- function(afs_element){
      alleles <- as.numeric(rep(row.names(afs_element), as.numeric(afs_element[, "freq"])))
      size_sd <- stats::sd(alleles)
      size_kurtosis <- moments::kurtosis(alleles, na.rm = TRUE)
      out <- data.frame(size_sd = size_sd, size_kurtosis = size_kurtosis)
    }
    all_allele_size_ss <- do.call(rbind, lapply(afs, all_alleles))
    
    mean_allele_size_sd <- mean(all_allele_size_ss$size_sd / repeat_size_per_locus, na.rm = TRUE)
    sd_allele_size_sd <- sd(all_allele_size_ss$size_sd / repeat_size_per_locus, na.rm = TRUE)
    
    mean_allele_size_kurtosis <- mean(all_allele_size_ss$size_kurtosis, na.rm = TRUE)
    sd_allele_size_kurtosis <- sd(all_allele_size_ss$size_kurtosis, na.rm = TRUE)
    
    
    # measure similar to mRatio
    if (mratio == "strict") {
      rpt_size <- 8:2
      mratio_all <- strataG::mRatio(g_types_geno, by.strata = FALSE, rpt.size = rpt_size)
      mratio_mean <- mean(mratio_all, na.rm = TRUE)
      mratio_sd <- stats::sd(mratio_all, na.rm = TRUE)
    } else if (mratio == "loose") {
      mratio_all <- m_ratio(g_types_geno)
      # mratio <- mratio[mratio != 1] # not sure if makes sense
      mratio_mean <- mean(mratio_all, na.rm = TRUE)
      mratio_sd <- NA
      if (!is.na(mratio_mean))  mratio_sd <- stats::sd(mratio_all, na.rm = TRUE)
      
    } else {
      stop("specify whether mratio is calculated strict or loose ")
    }
    # mratio might be larger than 1 in "loose"
    if (!is.na(mratio_mean) & mratio_mean > 1) mratio_mean <- 1
    
    out <- data.frame(
      num_alleles_mean, num_alleles_sd,
      exp_het_mean, exp_het_sd,
      obs_het_mean, obs_het_sd,
      mean_allele_size_sd, sd_allele_size_sd,
      mean_allele_size_kurtosis, sd_allele_size_kurtosis,
      mean_allele_range, sd_allele_range,
      mratio_mean, mratio_sd,
      prop_low_afs_mean,  prop_low_afs_sd)
  }
  
  
  
  # functions to draw nresamp times nind individuals and nloc loci from the genotypes and recalculate summary stats
  single_rarefaction <- function(resamp_num, geno, nind, nloc){
    # defaults to all individuals and all loci in case nothing is specified
    if (is.null(nind)) nind <- nrow(geno)
    if (is.null(nloc)) nloc <- ncol(geno) / 2
    # in case that nind or nloc is greater than the actual dataset, set nind or nloc down
    if(nind > nrow(geno)) nind <- nrow(geno)
    if(nloc > ncol(geno)/2) nloc <- ncol(geno/2)
    # subsamples nind individuals from the data
    sample_inds <- sample(1:nrow(geno), nind, replace = FALSE)
    # subsamples ncol loci from the data
    locus_vec <- seq(from = 1, to = ncol(geno), by = 2)
    sample_locs <- sample(locus_vec, nloc, replace = FALSE)
    full_sample_locs <- sort(c(sample_locs, sample_locs + 1))
    # calc ss
    out <- calculate_sumstats(geno[sample_inds, full_sample_locs])
  }
  
  # CI function
  CI <- 0.95
  calc_CI <- function(x) {
    out <- stats::quantile(x, c((1 - CI)/2, 1 - (1 - CI)/2), na.rm = TRUE)
    names(out) <- c("CIlow", "CIhigh")
    out <- data.frame(t(out))
  }
  
  # draw nresamp times nind individuals and nloc loci from the genotypes and recalculate summary stats
  calc_ss_with_rarefaction <- function(geno, nresamp, nind, nloc){
    all_rarefac_sumstats <- lapply(1:nresamp, single_rarefaction, geno, nind, nloc)
    all_rarefac_df <- do.call(rbind, all_rarefac_sumstats)
    out <- as.data.frame(lapply(all_rarefac_df, function(x) out <- data.frame("PE" = mean(x, na.rm = TRUE), calc_CI(x))))
    # out <- c("PE" = mean(x), calc_CI(x))))
    names(out) <- gsub(".PE", "", names(out))
    #out
    out2<-all_rarefac_df
    mylist <- list(out, out2)
    mylist
  }
  
  
  # summary statistics within populations
  if (!is.null(by_pop)){
    # make sure pop is character
    genotypes[[by_pop]] <- as.character(genotypes[[by_pop]])
    # list all populations (or clusters or whatever)
    all_pops <- names(table(genotypes[[by_pop]]))
    # delete populations with just one individual from the list
    # discard clusters with few individuals
    one_ind_pop <- any(as.numeric(table(genotypes[[by_pop]])) <= 10)
    if (one_ind_pop) {
      all_pops <- all_pops[-(which(as.numeric(table(genotypes[[by_pop]])) == 1))]
    }
    
    if (rarefaction == TRUE){
      all_sumstats <- lapply(all_pops, function(x)  calc_ss_with_rarefaction(genotypes[genotypes[[by_pop]] == x,
                                                                                       start_geno:ncol(genotypes)], nresamp, nind, nloc))
      out <- as.data.frame(do.call(rbind, all_sumstats))
    } else if (rarefaction == FALSE){
      all_sumstats <- lapply(all_pops, function(x)  calculate_sumstats(genotypes[genotypes[[by_pop]] == x, start_geno:ncol(genotypes)]))
      out <- as.data.frame(do.call(rbind, all_sumstats))
    }
    
    row.names(out) <- all_pops
  }
  # summary statistics for the full dataset
  if (is.null(by_pop)){
    
    if (rarefaction == TRUE){
      out <- calc_ss_with_rarefaction(genotypes[, start_geno:ncol(genotypes)], nresamp, nind, nloc)
    } else if (rarefaction == FALSE) {
      out <- calculate_sumstats(genotypes[, start_geno:ncol(genotypes)])
    }
    
  }
  
  out
  
}