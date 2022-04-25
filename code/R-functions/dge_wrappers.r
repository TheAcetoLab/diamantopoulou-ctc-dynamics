# Check for pairing
# Use only selected samples

# Function run_edgeR_comp
##' @title Run differential using a DGEList or SingleCellExperiment as input
##' @description Wrapper to run edgeR from a DGEList or SCEset using different configuration options.
##' @param object either a DGEList or a SingleCellExperiment
##' @param use_samples character or logical array of samples to analyse.
##'
##' Desing configuration for differential expression
##' @param group_var a character with name of the variable from sample annotation tested for differential expression.
##' @param group_sample a character string with the name of evaluated level from group_var.
##' @param group_ref a character string with the name of referece level from group_var.
##' @param numeric_covar a character string with a name of a numerical variable to block the differential expression analysis.
##' @param batch_vars a character string with the name of a categorical variable to block the differential expression analysis.
##' @param design_formula a character string with the design formular, e.g.: "~ 0 + batch_vars + group_var".
##' @param contrast a character string with the contrast passed to glmQLFTest or glmLRT. Contrast should be a numeric vector or
##' matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero. If NULL and coef==NULL,
##' contrast will be created as 'group_var''group_sample'-'group_var''group_ref'
##' @param coef integer or character index vector passed to glmQLFTest or glmLRT indicating which coefficients of the linear model are to be tested
##' equal to zero. Ignored if contrast is not NULL. If 'last', the last column of design matrix will be used
##'
##' Conversion from SingleCellExperiment to DGEList
##' @param spike_normalization normalize data using spike-in size factors for SingleCellExperiment object.
##' @param assay_to_DGEList a character indicating the type of assay matrix from SingleCellExperiment put into the DGEList object.
##' Default "counts".
##' @param assay_to_row_filter a character indicating the type of assay matrix used for filtering genes. Default, "tpm".
##' If "object" is DGEList, only "counts" is available at the moment.
##' @param use_colData a numeric array with column indexes from colData slot of the SingleCellExperiment object to add to DGEList.
##' @param use_rowData a numeric array with column indexes from rowData slot of the SingleCellExperiment object to add to DGEList.
##'
##' Feature filtering parameters
##' @param use_filterByExpr a logical, wether to runedgeR::filterByExpr to filter features by expression level. Only available
##' if group_var is not numeric.
##' @param min_counts a integer value of the minimum expression level for a feature to be considered. Ignored if
##' use_filterByExpr is not NULL.
##' @param min_present_prop a numeric value of the minimum proportion samples for a feature to be selected.
##' Expression determined by  >= min.count level for each compared group. Ignored if use_filterByExpr is not NULL.
##'
##' EdgeR workflow configuration
##' @param run_calcNormFactors character string or logical indicating wether to run internal edgeR normalization
##' using calcNormFactors. If TRUE, TMM method will be run. Alternatively, a string indicating the normalization
##' method to be used "TMM","TMMwsp","RLE","upperquartile","none". If NULL and object is SingleCellExperiment class
##' sizeFactors stored from the original object will be used for normalization.
##' @param estimateDisp_robust logical indicating if the estimation of prior.df be robustified against outliers.
##' See edgeR::estimateDisp.
##' @param estimateDisp_trend.method character string specifying the trend method used in estimateDisp function. S
##' ee edgeR::estimateDisp.
##' @param glmQLFit_robust logical, whether to estimate the prior QL dispersion distribution robustly. Only used
##' if glm_approach is QL. See edgeR::glmQLFit.
##' @param glm_approach character string specifying the approach used for edgeR. Possible values are "LRT" for
##' the likelihood ratio tests or "QLF" for the  quasi-likelihood F-tests.
##'
##' Output configuration
##' @param adjust_method character string passed to p.adjust method. Possible values are "none", "BH", "fdr"
##' (equivalent to "BH"), "BY" and "holm". Default 'BH'.
##' @param assays_from_SingleCellExperiment a character array indicating the names of the assays from the original
##' SingleCellExperiment to return, but not used for the differential expression analysis. Only used if input object
##' is of SingleCellExperiment class
##' @template roxygen-template
##' @details
##' This function is a wrapper to run edgeR differential expression wirkflow using a DGEList or SingleCellExperiment as
##' input and multiple parameters.
##' @return object of class DGEcomp.
edgeR_dge <- function(
  object,
  use_samples = NULL,
  # Desing configuration for differential expression
  group_var = NULL,
  group_sample = NULL,
  group_ref = NULL,
  numeric_covar = NULL,
  batch_vars = NULL,
  design_formula = NULL,
  contrast = NULL,
  coef = NULL,
  # Conversion from SingleCellExperiment to DGEList
  spike_normalization = FALSE,
  assay_to_DGEList = 'counts',
  assay_to_row_filter = "tpm",
  use_colData = NULL,
  use_rowData = NULL,
  # Feature filtering parameters
  use_filterByExpr = TRUE,
  min_counts = 1,
  min_present_prop = 0.40,
  # EdgeR workflow configuration
  run_calcNormFactors = 'TMM',
  estimateDisp_robust = FALSE,
  estimateDisp_trend.method = "locfit",
  glmQLFit_robust = FALSE,
  glm_approach = "QLF",
  # Output configuration
  adjust_method = 'BH',
  assays_from_SingleCellExperiment = NULL

  )
{

  require(edgeR)
  require(biomaRt)

  # Subset of samples used in the comparison
  if(!is.null(use_samples)){
    de_obj <- object[, use_samples]
  } else {
    de_obj <- object
  }


  # Processing SingleCellExperiment object
  if(class(de_obj) == 'SingleCellExperiment'){
    require(scran)

    # Spike-in normalization
    if(spike_normalization){
      de_obj <- computeSpikeFactors(de_obj, general.use=TRUE)
      de_obj <- normalize(de_obj)
    }

    # Get assay data from SingleCellExperiment for filtering genes
    if(assay_to_row_filter != assay_to_DGEList)
      exprs_mat <- assay(de_obj, assay_to_row_filter)

    # Creating DGEList from SingleCellExperiment
    # In convertTo edgeR, size factors are turned into normalization factors. This conversion is performed by :
    # source : getMethod("convertTo","SingleCellExperiment")
    # sf <- colData(de_obj)$size_factor
    # dds <- convertTo(de_obj, type = "edgeR")
    # nf <- log(sf/dds$samples$lib.size)
    # nf <- exp(nf - mean(nf))

    if(is.null(use_colData))
      use_colData <- 1:ncol(colData(de_obj))
    if(is.null(use_rowData))
      use_rowData <- 1:ncol(rowData(de_obj))

    de_obj <- convertTo(
      de_obj,
      type = "edgeR",
      assay.type = assay_to_DGEList
      )

  }


  # Extract info from DGEList
  sample_annot <- de_obj$sample
  counts_mat <- de_obj$counts
  if(assay_to_row_filter == assay_to_DGEList)
    exprs_mat <- counts_mat

  # Comparison name
  comparison_name <- ifelse(
    !is.numeric(sample_annot[, group_var]),
    paste0(group_var, '_', group_sample, "--over--", group_ref),
    group_var
  )


  # Create variables to design matrix
  if(!is.numeric(sample_annot[, group_var])) {
    assign(group_var, factor(sample_annot[, group_var]) %>% relevel(., group_ref))
  } else {
    assign(group_var, sample_annot[, group_var])
  }
  if(!is.null(batch_vars))
    for(i in batch_vars)
      assign(i, factor(sample_annot[, i]))
  if(!is.null(numeric_covar))
    for( i in numeric_covar)
      assign(i, sample_annot[, i])

  # Create desing variable and contrast
  design_mat <- model.matrix(as.formula(design_formula))

  if(coef == 'last')
    coef <- ncol(design_mat)
  if(is.null(contrast) & is.null(coef))
      contrast <- paste(
        paste0(group_var, group_sample),
        paste0(group_var, group_ref),
        sep = "-")

  # Add samples compared to the annotation
  if(!is.numeric(sample_annot[, group_var])) {
    de_obj$sample$use_sample <- sample_annot[, group_var] %in% c(group_ref, group_sample)
  } else {
    de_obj$sample$use_sample <- TRUE
  }


  # Filter genes by expression level using filterByExpr : The filterByExprfunction keeps rows that have worthwhile counts
  # in a minumum number of samples. The function accesses the group factor contained in y in order to compute the minimum
  # group size, but the filtering is performed independently of which sample belongs to which group so that no bias is introduced.
  # The group factor or the experimental design matrix can also be given directly to the filterByExpr function if not already
  # set in the DGEList object. It is also recommended to recalculate the library sizes of the DGEListobject after the filtering,
  # although the downstream analysis is robust to whether this is done or not. The actual filtering uses CPM values rather than counts
  # in order to avoid giving preference to samples with large library sizes.

  if(use_filterByExpr & !is.numeric(sample_annot[, group_var])) {
    use_groups <- sample_annot[, group_var]
    if(!is.factor(use_groups))
      use_groups <- factor(use_groups)
    keep_features <- filterByExpr(de_obj, group = sample_annot[, group_var] %>% factor, min.count = min_counts, min.prop = min_present_prop)
  }

  # Filter genes by expression with manually defined parameters
  features_prop_matrix <- NULL
  if(!use_filterByExpr | is.numeric(sample_annot[, group_var])) {
    is_present <- exprs_mat >= min_counts
    if(!is.numeric(sample_annot[, group_var])) {
      is_ref <- sample_annot[, group_var] == group_ref
      is_sample <- sample_annot[, group_var] == group_sample

      prop_present_all <- apply(is_present, 1, mean)
      is_present_all <- prop_present_all >= min_present_prop
      n_present_all <- apply(is_present, 1, sum)

      prop_present_ref <- apply(is_present[,is_ref, drop = FALSE], 1, mean)
      is_present_ref <- prop_present_ref >= min_present_prop
      n_present_ref <- apply(is_present[,is_ref, drop = FALSE], 1, sum)

      prop_present_sample <- apply(is_present[,is_sample, drop = FALSE], 1, mean)
      is_present_sample <- prop_present_sample >= min_present_prop
      n_present_sample <- apply(is_present[,is_sample, drop = FALSE], 1, sum)

      keep_features <- is_present_all | is_present_ref | is_present_sample
      features_prop_matrix <- cbind(
        p.all = prop_present_all,
        p.ref = prop_present_ref,
        p.sample = prop_present_sample,
        f.all = n_present_all,
        f.ref = n_present_ref,
        f.sample = n_present_sample
        ) %>%
        data.frame %>%
        rownames_to_column('feature')
    } else {
      n_present_all <- rowSums2(is_present)
      prop_present_all <- apply(is_present, 1, mean)
      keep_features  <- prop_present_all >= min_present_prop
      use_samples <- rep(TRUE, nrow(sample_annot))
      features_prop_matrix <- cbind(p.all = prop_present_all)  %>%
        data.frame %>%
        rownames_to_column('feature')
    }
  }

  # Filter features
  de_obj$keep_features <- keep_features
  keep_features <- de_obj$keep_features # named logical array
  de_obj <- de_obj[de_obj$keep_features,,keep.lib.sizes=FALSE]

  # TMM normalization
  if(!is.null(run_calcNormFactors)){
    if(is.logical(run_calcNormFactors)){
      de_obj <- calcNormFactors(de_obj, method = 'TMM')
    } else {
      de_obj <- calcNormFactors(de_obj, method = run_calcNormFactors)
    }
  }

  # Estimate the negative binomial (NB) dispersions
  de_obj <- estimateDisp (
    de_obj,
    design = design_mat,
    robust = estimateDisp_robust,
    trend.method = estimateDisp_trend.method
  )

  # QLF approach
  if(glm_approach == 'QLF') {
    fit <- glmQLFit(
      de_obj,
      design = design_mat,
      robust = glmQLFit_robust
      )

    if(!is.null(contrast)){
      contrast_mat <- makeContrasts(contrasts=contrast, levels=design_mat)
      dge_lrt <- glmQLFTest(fit, contrast=contrast_mat)
    } else {
      dge_lrt <- glmQLFTest(fit, coef=coef)
      }
  }

  # LRT approach
  if(glm_approach == 'LRT') {
    # Fit negative bionomial GLM
    fit <- glmFit(
      de_obj,
      design = design_mat
      )
    # Carry out Likelihood ratio test
    if(!is.null(contrast)){
      contrast_mat <- makeContrasts(contrasts=contrast, levels=design_mat)
      dge_lrt <- glmLRT(fit, contrast=contrast_mat)
    } else {
      dge_lrt <- glmLRT(fit, coef=coef)
    }
  }


  # P-value adjustment
  dge_lrt$table$FDR <- p.adjust(dge_lrt$table$PValue, method = adjust_method)

  # Create table with differential expression results
  de_table <- cbind(dge_lrt$genes, dge_lrt$table) %>%
    data.frame %>%
    rownames_to_column('feature') %>%
    dplyr::arrange(FDR, PValue)
  if(!is.null(features_prop_matrix))
    de_table <- left_join(de_table, features_prop_matrix, by = 'feature')

  de_table <-  de_table %>%
    dplyr::mutate(rownames = feature) %>%
    column_to_rownames('rownames')

  # Collect arguments
  run_args <-  as.list(match.call())

  # Remove object from match.call to reduce size of the output object
  run_args$object <- NULL

  # Prepare output list
  res <- list(
    run_args = run_args,
    results = de_table,
    dge_obj = list(
      DGEList = de_obj,
      DGEGLM = fit,
      DGELRT = dge_lrt
    ),
    comparison = comparison_name
  )

  # Add additional expression matrices from input DGEList
  if(class(object) == 'DGEList') {
    res$assays <- list(
      logcpm = edgeR::cpm(de_obj, log=TRUE, normalized.lib.sizes = FALSE),
      norm_logcpm = edgeR::cpm(de_obj, log=TRUE, normalized.lib.sizes = TRUE)
    )
  }

  # Add additional expression matrices from input SingleCellExperiment
  if(class(object) == 'SingleCellExperiment'){
    if(!is.null(assays_from_SingleCellExperiment)){
      assays_from_SingleCellExperiment <- intersect (
        assays_from_SingleCellExperiment,
        assayNames(object)
        )
      if(length(assays_from_SingleCellExperiment) > 0) {
        res$assays <- NULL
        for (i in assays_from_SingleCellExperiment) {
          res[[i]] <- assay(object, i)[de_table$feature, ]
        }
      }
    }
  }

  class(res) <- append(class(res), "DGEcomp")
  return(res)
}






