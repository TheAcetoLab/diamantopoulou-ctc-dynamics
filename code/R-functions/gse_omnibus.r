# Function gse_omnibus
# TDO : Harmonize column names (p.value, fdr ....)
##' @title Run multiple tools for gene and pathway analysis from differential expression results
##' @description Wrapper to run multiple gense set and pathway enrichment tools. Input must be a data frame containing
##' log fold-changes, P values and adjusted P values. The function runs different type of analyses:
##'
##' Hypergeometric test : topGO, goseq, G:Profiler, clusterProfiler enricher (enricher, enrichGO, enrichKEGG, enrichPathway)
##'
##' GSEA-like : clusterProfiler gse (GSEA, gseGO, gseKEGG, gsePathway)
##'
##' Topological : ROntoTools
##'
##'
##' @param p numeric, array of p values (either unadjusted or adjusted). Must have the same length and order than 'fc' and 'feature_names'
##' @param fc numeric, array of log fold-changes. Must have the same length and order than 'p' and 'feature_names'
##' @param feature_names character, names of features. Must have the same length and order than 'p' and 'fc'
##' @param feature_type character indicating the feature type id :
##' @param p_thrs
##' @param fc_thrs
##' @param annot_db
##' @param organism
##' @param kegg_organism
##' @param reactome_organism
gse_omnibus <- function(
  p = NULL,
  fc = NULL,
  feature_names = NULL,
  feature_type = 'SYMBOL',
  p_thrs = 0.05,
  fc_thrs = 0.5,
  return_sets = c('abs', 'up', 'down'),

  annot_db = 'org.Hs.eg.db',
  organism = 'hsapiens',
  kegg_organism = 'hsa',
  reactome_organism = 'human',
  go_ontologies = c('BP', 'MF', 'CC'),
  gmt_files = NULL,
  save_intermediates = NULL,

  run_all_ora = FALSE,
  run_all_gsea = FALSE,

  run_topGO = FALSE,
  run_goseq = FALSE,
  run_gprofiler = FALSE,

  run_enricher = FALSE,
  run_enrichGO = FALSE,
  run_enrichKEGG = FALSE,
  run_enrichMKEGG = FALSE,
  run_enrichReactome = FALSE,

  run_GSEA = FALSE,
  run_gseGO = FALSE,
  run_gseKEGG = FALSE,
  run_gseMKEGG = FALSE,
  run_gseReactome = FALSE,
  run_gprofiler_gsea = FALSE,

  run_ROntoTools = FALSE,
  run_ROntoTools_all = FALSE,

  run_GSEAPreranked = FALSE,

  args_topGO = list(algorithm = "elim", statistic = "fisher", nodeSize = 10, topNodes = 100, pvalueCutoff = 1),
  args_goseq = list(genome = "hg38", id = "knownGene", method = "Wallenius"),
  args_gprofiler = list(
    correction_method = "g_SCS"
  ),
  args_enricher = list(qvalueCutoff = 0.25, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1), #common arguments for enricher, enrichGO, enrichKEGG, enrichMKEGG
  args_gse = list(minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1),
  args_ROntoTools = list(nboot = 5000, seed = 12345, verbose = FALSE),
  args_GSEAPreranked = list(
    msigdb_collection = c('h.all.v7.0', 'c2.cp.v7.0'),
    gsea_cli_path = '~/apps/gsea/4.0.3/gsea-cli.sh',
    run_script = TRUE,
    rnd_seed = 1029384756)
)
{

  # Sanity checks
  if(!length(feature_names) == length(p) & length(p) == length(fc))
    stop("different length of feature_names, p and fc")

  # Generate directory for intermediate files
  if(!is.null(save_intermediates)){
    if(!file.exists(save_intermediates)) {
      dir.create(save_intermediates, recursive = TRUE)
    }
  }

  # List of results to return
  gse_res <- list()

  # Collect arguments
  gse_res$run_args <-  as.list(match.call())
  for (i in c('feature_type', 'p_thrs', 'fc_thrs')) {
    gse_res$run_args[[i]] <- get(i)
  }

  # Table with number of selected features
  gse_res$ora_thrs_tab <- table(
    p <= p_thrs,
    abs(fc) >= fc_thrs
  ) %>% data.frame %>% set_names(c('p_pass', 'fc_pass', 'freq'))

  # Generate data.frame with input data
  object <- data.frame(fc, p, feature_names, stringsAsFactors = FALSE) %>%
    set_names(c('fc', 'p', feature_type))
  object <- object[!duplicated(object[[feature_type]]),]

  # Generate entrez id
  if(feature_type != "ENTREZID") {
    eg <- clusterProfiler::bitr(object[[feature_type]], fromType=feature_type, toType="ENTREZID", OrgDb=annot_db)
    eg <- eg[!duplicated(eg[[feature_type]]),] # Remove duplicated mappings
    eg <- eg[!duplicated(eg[["ENTREZID"]]),]
    object <- left_join(object, eg)
  }

  # obtain list of features and statistics required for the different GSE methods
  feature_lists_de_fdr_fc <- generate_de_gene_list(object,
                                                   pvalue = p_thrs,
                                                   fc = fc_thrs,
                                                   feature_var = feature_type,
                                                   fc_var = 'fc',
                                                   pvalue_var = 'p',
                                                   return_sets = return_sets)
  feature_lists_de_fdr_fc_names <- map(feature_lists_de_fdr_fc, function(x) names(x)[x==1]) %>% unlist %>% unique
  feature_stats_de_fdr_fc <- generate_stats_list(object,
                                                 features = feature_lists_de_fdr_fc_names,
                                                 feature_var = feature_type,
                                                 fc_var = 'fc',
                                                 pvalue_var = 'p')
  feature_stats_de_fdr <- generate_stats_list(object,
                                              pvalue = p_thrs,
                                              feature_var = feature_type,
                                              pvalue_var = 'p')
  feature_fc_list <- generate_fc_list(object, feature_var = feature_type, fc_var = 'fc')
  feature_stats <- generate_stats_list(object,
                                       pvalue = 1,
                                       fc = 0,
                                       feature_var = feature_type,
                                       fc_var = 'fc',
                                       pvalue_var = 'p')

  # obtain list of features (entrez id) and statistics required for the different GSE methods
  if(feature_type != "ENTREZID") {
    eg_lists_de_fdr_fc <- generate_de_gene_list(object,
                                                pvalue = p_thrs,
                                                fc = fc_thrs,
                                                feature_var = 'ENTREZID',
                                                fc_var = 'fc',
                                                pvalue_var = 'p',
                                                return_sets = return_sets)
    eg_lists_de_fdr_fc_names <- map(eg_lists_de_fdr_fc, function(x) names(x)[x==1]) %>% unlist %>% unique
    eg_stats_de_fdr_fc <- generate_stats_list(object,
                                              features = eg_lists_de_fdr_fc_names,
                                              feature_var = 'ENTREZID',
                                              fc_var = 'fc',
                                              pvalue_var = 'p')
    eg_stats_de_fdr <- generate_stats_list(object,
                                           pvalue = p_thrs,
                                           feature_var = 'ENTREZID',
                                           pvalue_var = 'p')
    eg_fc_list <- generate_fc_list(object, feature_var = 'ENTREZID', fc_var = 'fc')
    eg_stats <- generate_stats_list(object,
                                    pvalue = 1,
                                    fc = 0,
                                    feature_var = 'ENTREZID',
                                    fc_var = 'fc',
                                    pvalue_var = 'p')
  } else {
    eg_lists_de_fdr_fc <- feature_lists_de_fdr_fc
    eg_stats_de_fdr_fc <- feature_stats_de_fdr_fc
    eg_stats_de_fdr <- feature_stats_de_fdr
    eg_fc_list <- feature_fc_list
    eg_stats <- feature_stats
  }

  # Return gene lists
  gse_res$gene_list <- list(
    feature_lists_de_fdr_fc = feature_lists_de_fdr_fc,
    feature_lists_de_fdr_fc_names = feature_lists_de_fdr_fc_names,
    feature_stats_de_fdr_fc = feature_stats_de_fdr_fc,
    feature_stats_de_fdr = feature_stats_de_fdr,
    feature_fc_list = feature_fc_list,
    feature_stats = feature_stats,
    eg_lists_de_fdr_fc = eg_lists_de_fdr_fc,
    eg_stats_de_fdr_fc = eg_stats_de_fdr_fc,
    eg_stats_de_fdr = eg_stats_de_fdr,
    eg_fc_list = eg_fc_list,
    eg_stats = eg_stats
  )


  # define GO ontologies to use
  names(go_ontologies) <- go_ontologies

  # define other gene-sets from input gmt files
  if(!is.null(gmt_files)) {
    gmt_gsets <- lapply(gmt_files, clusterProfiler::read.gmt)
    if(feature_type != "ENTREZID") {
      gmt_gsets_eg <- lapply(gmt_gsets, function(gg) {
        eg <- clusterProfiler::bitr(gg$gene, fromType=feature_type, toType="ENTREZID", OrgDb=annot_db)
        eg <- eg[!duplicated(eg[[feature_type]]),] # Remove duplicated mappings
        eg <- eg[!duplicated(eg[["ENTREZID"]]),]
        gg <- left_join(gg, eg, by = c( 'gene' = feature_type))
        gg <- na.omit(gg)
        gg$gene <- gg$ENTREZID
        gg$ENTREZID <- NULL
        gg
      })
    } else{
      gmt_gsets_eg <- gmt_gsets
    }
  }

  # Run topGO
  if(run_topGO | run_all_ora) {
    message("Running topGO")
    gse_res$topGO <- map(feature_lists_de_fdr_fc, function(x) {
      res <- foreach(ontology = go_ontologies) %do% {
        do.call('wrapper_topGO', c(
          list(
            x,
            ID = feature_type,
            ontology = ontology,
            annot_db = annot_db
          ),
          args_topGO
        ))
      }
      names(res) <- go_ontologies
      return(res)
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$topGO, file = file.path(save_intermediates, 'topGO.rds'))
  }

  # Run goseq GO terms
  if(run_goseq | run_all_ora) {
    require(goseq)
    message("Running goseq for GO terms")
    gse_res$goseq <- map(eg_lists_de_fdr_fc, function(x) {
      pwf <- do.call('nullp', c(list(x, plot.fit = FALSE), args_goseq[c('genome', 'id')]))

      res <- foreach(ontology = go_ontologies) %do% {
        enr.res <- list('pwf' = pwf)
        enr.res$result <- do.call('goseq', c(
          list(
            pwf,
            test.cats  = paste0("GO:", ontology)
          ),
          args_goseq
          )
        )
        enr.res$result <- enr.res$result %>%
          dplyr::rename("Count" = "numDEInCat", "setSize" = "numInCat") %>%
            mutate(
              GeneProp = Count / setSize,
              pvalue = p.adjust(over_represented_pvalue, method = 'BH')
              )
        return(enr.res)
      }
      names(res) <- go_ontologies
      return(res)
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$goseq, file = file.path(save_intermediates, 'goseq.rds'))
  }

  # Run goseq using gene sets from gmt_files
  if((run_goseq | run_all_ora) & !is.null(gmt_files))  {
    require(goseq)
    message("Running goseq from GMT files")
    gse_res$goseq_gmt <- map(eg_lists_de_fdr_fc, function(x) {
      pwf <- do.call('nullp', c(list(x, plot.fit = FALSE), args_goseq[c('genome', 'id')]))
      res <- foreach(use_gsets = gmt_gsets_eg) %do% {
        enr.res <- list('pwf' = pwf)
        enr.res$result <- do.call('goseq', c(
          list(
            pwf,
            gene2cat = use_gsets
          ),
          args_goseq
        )
        )
        enr.res$result <- enr.res$result %>%
          dplyr::rename("Count" = "numDEInCat", "setSize" = "numInCat") %>%
          mutate(
            GeneProp = Count / setSize,
            pvalue = p.adjust(over_represented_pvalue, method = 'BH')
          )
        return(enr.res)
      }
      names(res) <- names(gmt_gsets_eg)
      return(res)
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$goseq_gmt, file = file.path(save_intermediates, 'goseq_gmt.rds'))
  }

  # Run g:Profiler over enrichment
  # https://biit.cs.ut.ee/gprofiler/page/docs
  if(run_gprofiler | run_all_ora) {
    require(gprofiler2)
    message("Running g:Profiler")
    gse_res$gprofiler <- map(feature_lists_de_fdr_fc, function(x) {
      do.call('gost', c(
        list(
          names(x)[x == '1'],
          custom_bg = names(x),
          organism = organism
        ),
        args_gprofiler
      )
      )
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$gprofiler, file = file.path(save_intermediates, 'gprofiler.rds'))
  }

  # Run g:Profiler ogsea
  # https://biit.cs.ut.ee/gprofiler/page/docs
  if(run_gprofiler_gsea | run_all_gsea) {
    require(gprofiler2)
    message("Running g:Profiler")
    gse_res$gprofiler_gsea <- do.call('gost', c(
      list(
        names(feature_fc_list),
        custom_bg = names(feature_fc_list),
        organism = organism
      ),
      args_gprofiler
    )
    )
    if(!is.null(save_intermediates))
      saveRDS(gse_res$gprofiler_gsea, file = file.path(save_intermediates, 'gprofiler_gsea.rds'))
  }


  # Run clusterProfiler::enricher
  if((run_enricher | run_all_ora) & !is.null(gmt_files) ) {
    require(clusterProfiler)
    message("Running clusterProfiler::enricher")
    gse_res$enricher <- map(feature_lists_de_fdr_fc, function(x) {
      res <- foreach(use_gsets = gmt_gsets) %do% {
        enr.i <- do.call('enricher', c(
          list(
            names(x)[x == '1'],
            TERM2GENE=use_gsets,
            universe = names(x)
          ),
          args_enricher
        ))
        if(!is.null(enr.i)) {
          enr.res <- enr.i@result
          enr.res$setSize <- as.numeric(gsub("/.*", "", enr.res$BgRatio))
          enr.res$GeneProp <- enr.res$Count / enr.res$setSize
          enr.i@result <- enr.res

        }
        enr.i
      }
      names(res) <- names(gmt_gsets)
      return(res)
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$enricher, file = file.path(save_intermediates, 'enricher.rds'))
  }

  # Run clusterProfiler::enrichGO
  if(run_enrichGO | run_all_ora) {
    require(clusterProfiler)
    gse_res$enrichGO <- map(feature_lists_de_fdr_fc, function(x) {
      res <- foreach(ontology = go_ontologies) %do% {
        message("Running clusterProfiler::enrichGO ", ontology)
        # enr.i <- enrichGO(names(x)[x == '1'], universe = names(x), ont = ontology, OrgDb = org.Hs.eg.db,
        # keyType = "SYMBOL", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1)
        enr.i <- do.call('enrichGO', c(
          list(
            names(x)[x == '1'],
            universe = names(x),
            ont = ontology,
            OrgDb = get(annot_db),
            keyType = feature_type
          ),
          args_enricher
        ))
        if(!is.null(enr.i)) {
          enr.i.sim <- simplify(enr.i, cutoff=0.7, by="p.adjust", select_fun=min)

          enr.res <- enr.i@result
          enr.res$setSize <- as.numeric(gsub("/.*", "", enr.res$BgRatio))
          enr.res$GeneProp <- enr.res$Count / enr.res$setSize
          enr.res$simplify <- enr.res$ID %in% enr.i.sim@result$ID
          enr.i@result <- enr.res
        }
        enr.i
      }
      names(res) <- go_ontologies
      return(res)
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$enrichGO, file = file.path(save_intermediates, 'enrichGO.rds'))
  }

  # Run clusterProfiler::enrichKEGG
  if(run_enrichKEGG | run_all_ora) {
    require(clusterProfiler)
    message("Running clusterProfiler::enrichKEGG")
    gse_res$enrichKEGG <- map(eg_lists_de_fdr_fc, function(x) {
      enr.i <- do.call('enrichKEGG', c(
        list(
          names(x)[x == '1'],
          universe = names(x),
          organism = kegg_organism
        ),
        args_enricher
      ))
      if(!is.null(enr.i)) {
        enr.i <- setReadable(enr.i, OrgDb = annot_db, keyType="ENTREZID")
        enr.res <- enr.i@result
        enr.res$setSize <- as.numeric(gsub("/.*", "", enr.res$BgRatio))
        enr.res$GeneProp <- enr.res$Count / enr.res$setSize
        enr.i@result <- enr.res
      }
      enr.i
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$enrichKEGG, file = file.path(save_intermediates, 'enrichKEGG.rds'))
  }

  # Run clusterProfiler::enrichMKEGG
  if(run_enrichMKEGG | run_all_ora) {
    require(clusterProfiler)
    message("Running clusterProfiler::enrichMKEGG")
    gse_res$enrichMKEGG <- map(eg_lists_de_fdr_fc, function(x) {
      enr.i <- do.call('enrichMKEGG', c(
        list(
          names(x)[x == '1'],
          universe = names(x),
          organism = kegg_organism
        ),
        args_enricher
      ))
      if(!is.null(enr.i)) {
        enr.i <- setReadable(enr.i, OrgDb = annot_db, keyType="ENTREZID")
        enr.res <- enr.i@result
        enr.res$setSize <- as.numeric(gsub("/.*", "", enr.res$BgRatio))
        enr.res$GeneProp <- enr.res$Count / enr.res$setSize
        enr.i@result <- enr.res
      }
      enr.i
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$enrichMKEGG, file = file.path(save_intermediates, 'enrichMKEGG.rds'))
  }

  # Run ReactomePA::enrichPathway
  if(run_enrichReactome | run_all_ora) {
    require(clusterProfiler)
    require(ReactomePA)
    message("Running ReactomePA::enrichPathway")
    gse_res$enrichReactome <- map(eg_lists_de_fdr_fc, function(x) {
      enr.i <- do.call('enrichPathway', c(
        list(
          names(x)[x == '1'],
          universe = names(x),
          organism = reactome_organism
        ),
        args_enricher
      ))
      if(!is.null(enr.i)) {
        enr.i <- setReadable(enr.i, OrgDb = annot_db, keyType="ENTREZID")
        enr.res <- enr.i@result
        enr.res$setSize <- as.numeric(gsub("/.*", "", enr.res$BgRatio))
        enr.res$GeneProp <- enr.res$Count / enr.res$setSize
        enr.i@result <- enr.res
      }
      enr.i
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$enrichReactome, file = file.path(save_intermediates, 'enrichReactome.rds'))
  }

  # Run clusterProfiler::GSEA
  if((run_GSEA | run_all_gsea) & !is.null(gmt_files) ) {
    require(clusterProfiler)
    message("Running clusterProfiler::GSEA")
    gse_res$GSEA <- map(gmt_gsets, function(use_gsets) {
      enr.i <- do.call('GSEA', c(
        list(
          feature_fc_list,
          TERM2GENE=use_gsets
        ),
        args_gse
      ))
      enr.res <- enr.i@result
      enr.res$Count <- strsplit(enr.res$core_enrichment, "\\/") %>% lapply(., length) %>% unlist
      enr.res$.sign <- ifelse(enr.res$NES < 0, 'supressed', 'activated')
      enr.res$GeneRatio <- enr.res$Count / enr.res$setSize
      enr.i@result <- enr.res

      return(enr.i)
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$GSEA, file = file.path(save_intermediates, 'GSEA.rds'))
  }

  # Run clusterProfiler::gseGO
  if(run_gseGO | run_all_gsea){
    require(clusterProfiler)
    gse_res$gseGO <- map(go_ontologies, function(ontology) {
      message("Running clusterProfiler::gseGO ", ontology)
      enr.i <- do.call('gseGO', c(
        list(
          feature_fc_list,
          ont = ontology,
          OrgDb = get(annot_db),
          keyType = feature_type
        ),
        args_gse
      ))
      enr.i.sim <- simplify(enr.i, cutoff=0.7, by="p.adjust", select_fun=min)

      enr.res <- enr.i@result
      enr.res$simplify <- enr.res$ID %in% enr.i.sim@result$ID
      enr.res$Count <- strsplit(enr.res$core_enrichment, "\\/") %>% lapply(., length) %>% unlist
      enr.res$.sign <- ifelse(enr.res$NES < 0, 'supressed', 'activated')
      enr.res$GeneRatio <- enr.res$Count / enr.res$setSize
      enr.i@result <- enr.res
      return(enr.i)
    })
    if(!is.null(save_intermediates))
      saveRDS(gse_res$gseGO, file = file.path(save_intermediates, 'gseGO.rds'))
  }

  # Run clusterProfiler::gseKEGG
  if(run_gseKEGG | run_all_gsea){
    require(clusterProfiler)
    message("Running clusterProfiler::gseKEGG")
    enr.i <-do.call('gseKEGG', c(
      list(
        eg_fc_list,
        organism = kegg_organism
      ),
      args_gse
    ))
    if(!is.null(enr.i)){
      enr.i <- setReadable(enr.i, OrgDb = annot_db, keyType="ENTREZID")
      enr.res <- enr.i@result
      enr.res$Count <- strsplit(enr.res$core_enrichment, "\\/") %>% lapply(., length) %>% unlist
      enr.res$.sign <- ifelse(enr.res$NES < 0, 'supressed', 'activated')
      enr.res$GeneRatio <- enr.res$Count / enr.res$setSize
      enr.i@result <- enr.res
      gse_res$gseKEGG <- enr.i
      if(!is.null(save_intermediates))
        saveRDS(gse_res$gseKEGG, file = file.path(save_intermediates, 'gseKEGG.rds'))
    }
  }

  # Run clusterProfiler::gseMKEGG
  if(run_gseMKEGG | run_all_gsea){
    require(clusterProfiler)
    message("Running clusterProfiler::gseMKEGG")
    enr.i <-do.call('gseMKEGG', c(
      list(
        eg_fc_list,
        organism = kegg_organism
      ),
      args_gse
    ))
    if(!is.null(enr.i)){
      enr.i <- setReadable(enr.i, OrgDb = annot_db, keyType="ENTREZID")
      enr.res <- enr.i@result
      enr.res$Count <- strsplit(enr.res$core_enrichment, "\\/") %>% lapply(., length) %>% unlist
      enr.res$.sign <- ifelse(enr.res$NES < 0, 'supressed', 'activated')
      enr.res$GeneRatio <- enr.res$Count / enr.res$setSize
      enr.i@result <- enr.res
      gse_res$gseMKEGG <- enr.i
      if(!is.null(save_intermediates))
        saveRDS(gse_res$gseMKEGG, file = file.path(save_intermediates, 'gseMKEGG.rds'))
    }
  }

  # Run ReactomePA::gsePathway
  if(run_gseReactome | run_all_gsea){
    require(clusterProfiler)
    require(ReactomePA)
    message("Running ReactomePA::gsePathway")
    enr.i <-do.call('gsePathway', c(
      list(
        eg_fc_list,
        organism = reactome_organism
      ),
      args_gse
    ))
    if(!is.null(enr.i)){
      enr.i <- setReadable(enr.i, OrgDb = annot_db, keyType="ENTREZID")
      enr.res <- enr.i@result
      enr.res$Count <- strsplit(enr.res$core_enrichment, "\\/") %>% lapply(., length) %>% unlist
      enr.res$.sign <- ifelse(enr.res$NES < 0, 'supressed', 'activated')
      enr.res$GeneRatio <- enr.res$Count / enr.res$setSize
      enr.i@result <- enr.res
      gse_res$gseReactome <- enr.i
      if(!is.null(save_intermediates))
        saveRDS(gse_res$gseReactome, file = file.path(save_intermediates, 'gseReactome.rds'))
    }
  }

  # Run ROntoTools with significant genes
  if(run_ROntoTools) {
    message("Running ROntoTools")
    gse_res$ROntoTools <- do.call('wrapper_ROntoTools', c(
      list(
        fc = eg_stats_de_fdr_fc$fc,
        pv = eg_stats_de_fdr_fc$padj,
        kegg_organism = kegg_organism,
        ref = names(eg_stats$fc)
      ),
      args_ROntoTools
    ))
    if(!is.null(save_intermediates))
      saveRDS(gse_res$ROntoTools, file = file.path(save_intermediates, 'ROntoTools.rds'))
  }

  # Run ROntoTools with all genes
  if(run_ROntoTools_all) {
    message("Running ROntoTools with all genes")
    gse_res$ROntoTools_all <- do.call('wrapper_ROntoTools', c(
      list(
        fc = eg_stats$fc,
        pv = eg_stats$padj,
        kegg_organism = kegg_organism
      ),
      args_ROntoTools
    ))
    if(!is.null(save_intermediates))
      saveRDS(gse_res$ROntoTools_all, file = file.path(save_intermediates, 'ROntoTools_all.rds'))
  }

  # Run GSEA pre-ranked
  if(run_GSEAPreranked) {
    message("Running GSEA Preranked")
    gse_res$GSEAPreranked <- do.call('wrapper_GSEAPreranked', c(
      list(
        feature_fc_list,
        feature_type = feature_type
      ),
      args_GSEAPreranked
    ))
    if(!is.null(save_intermediates))
      saveRDS(gse_res$GSEAPreranked, file = file.path(save_intermediates, 'GSEAPreranked.rds'))
  }

  # Return all data
  return(gse_res)
}


# Function gse_edger_omnibus
# Run gene-se enrichment analysis for edgeR objects DGEList and DGELRT
gse_edger_omnibus <- function(
  object,
  design = NULL,
  feature_type = 'SYMBOL',
  organism = 'Hs',
  annot_db = 'org.Hs.eg.db',
  go_ontologies = c('BP', 'MF', 'CC'),
  gmt_files = NULL,
  run_all = FALSE,
  run_fry = FALSE,
  run_camera = FALSE,
  run_romer = FALSE,
  run_goana = FALSE,
  run_kegga = FALSE,
  minGSSize = 10,
  maxGSSize = 500,
  ora_arg = list(FDR = 0.05, restrict.universe = TRUE)
) {

  # https://github.com/davismcc/fibroblast-clonality/blob/master/src/Rmd/DE_pathways_FTv62_callset_baseclone_vs_others.Rmd
  # https://www.biorxiv.org/content/biorxiv/early/2018/09/10/413047.full.pdf
  # https://f1000research.com/articles/5-1438/v2

  # I haven't provided this functionality in camera() because camera() is designed to work with directional pathways and F-statistics are not directional. camera() is not designed to work with GO either, because of its redundancy and lack of directionality. (https://support.bioconductor.org/p/87794/)

  # Self-contained versus Competitive gene set test
  ## A competitive test compares differential expression of the gene set to a standard defined by the complement of that gene set. A self-contained test, in contrast, compares the gene set to a fixed standard that does not depend on the measurements of genes outside the gene set. The competitive test is most popular. (https://academic.oup.com/bioinformatics/article/23/8/980/198511#1980970)
  ## Self-contained : Using the terminology of Goeman and Buhlmann (1), ‘self-contained’ gene set tests examine a set of genes in their own right without reference to other genes in the genome (3–8), whereas ‘competitive’ gene set tests compare genes in the test set relative to all other genes. Self-contained tests are of interest for assessing the relevance of an individual biological process to the experiment at hand (8), whereas the competitive tests focus more on distinguishing the most important biological processes from those that are less important. Competitive tests are overwhelmingly more commonly used in the genomic literature (9). (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/#!po=78.5714)
  # Self-contained
  ## fry() : function performs ROAST gene set tests [40]. It is a self-contained gene set test. Given a gene set, it tests whether the majority of the genes in the set are DE across the comparison of interest.
  # Competitive gene set test
  # camera() : The camera() function performs a competitive gene set test accounting for inter-gene correlation. It tests whether a set of genes is highly ranked relative to other genes in terms of differential expression
  # romer(): The romer() function performs a gene set enrichment analysis. It implements a GSEA approach [38] based on rotation instead of permutation

  # Directionality
  ## camera() is designed to work with directional pathways and F-statistics are not directional. camera() is not designed to work with GO either, because of its redundancy and lack of directionality.

  require(clusterProfiler)
  require(GSA)

  # List of results to return
  gse_res <- list()

  # Collect arguments
  gse_res$run_args <- as.list(match.call())
  gse_res$run_args$object <- NULL

  # Clean arguments
  if(is.null(minGSSize))
    minGSSize <- -Inf
  if(is.null(minGSSize))
    maxGSSize <- Inf

  # Generate entrez id
  if(feature_type != "ENTREZID") {
    eg <- clusterProfiler::bitr(rownames(object$DGELRT), fromType=feature_type, toType="ENTREZID", OrgDb=annot_db)
    eg <- eg[!duplicated(eg[[feature_type]]),] # Remove duplicated mappings
    eg <- eg[!duplicated(eg[["ENTREZID"]]),]

    # use only genes with entrez id
    use_genes <- rownames(object$DGELRT) %in% eg[[feature_type]]
    object$DGELRT <- object$DGELRT[use_genes,]
    object$DGEList <- object$DGEList[use_genes,]

    # change rownames to entrez id
    rownames(object$DGELRT) <- eg$ENTREZID
    rownames(object$DGEList) <- eg$ENTREZID
  }

  # Load GMT files
  if(!is.null(gmt_files)){
    gmt_gsets <- lapply(gmt_files, function(x){
      invisible(capture.output( res <- GSA::GSA.read.gmt(x)))
      names(res$genesets) <- res$geneset.names
      geneset_length <- lapply(res$genesets, length)
      geneset_use <- geneset_length >= minGSSize & geneset_length <= maxGSSize
      return(res$genesets[geneset_use])
    })
    gmt_gsets_indices <- lapply(gmt_gsets, ids2indices, identifiers = rownames(object$DGEList))
  }

  # Run fry on DGEList
  if((run_fry | run_all) & !is.null(gmt_files))
    gse_res$fry <- map(gmt_gsets_indices, function(x){
      fry(y = object$DGEList, index = x, design = design) %>%
        rownames_to_column('Term')
    })

  # Run camera on DGEList
  if((run_camera | run_all) & !is.null(gmt_files))
    gse_res$camera <- map(gmt_gsets_indices, function(x){
      camera(y = object$DGEList, index = x, design = design) %>%
        rownames_to_column('Term')
    })

  # # Run romer on DGEList
  # if((run_romer | run_all) & !is.null(gmt_files))
  #   gse_res$romer <- map(gmt_gsets_indices, function(x){
  #     romer(y = object$DGEList, index = x) %>%
  #       rownames_to_column('Term')
  #   })

  # Run goana
  if(run_goana | run_all){
    gse_res$goana <- do.call('goana', c(
      list(
        object$DGELRT,
        species = organism
      ),
      ora_arg
    )) %>%
      rownames_to_column('ID')
    pmin_ord <- apply(gse_res$goana, 1, function(x) min(x['P.Up'], x['P.Down'])) %>% as.numeric %>% order
    gse_res$goana <- gse_res$goana[pmin_ord,]
    gse_res$goana$Direction <- ifelse(gse_res$goana$Down > gse_res$goana$Up, 'Down', 'Up')
    gse_res$goana$GeneRatio.up <- gse_res$goana$Up / gse_res$goana$N
    gse_res$goana$GeneRatio.down <- gse_res$goana$Down / gse_res$goana$N
  }

  # Run kegga
  if(run_kegga | run_all){
    gse_res$kegga <- do.call('kegga', c(
      list(
        object$DGELRT,
        species = organism
      ),
      ora_arg
    )) %>%
      rownames_to_column('ID')
    pmin_ord <- apply(gse_res$kegga, 1, function(x) min(x['P.Up'], x['P.Down'])) %>% as.numeric %>% order
    gse_res$kegga <- gse_res$kegga[pmin_ord,]
    gse_res$kegga$Direction <- ifelse(gse_res$kegga$Down > gse_res$kegga$Up, 'Down', 'Up')
    gse_res$kegga$GeneRatio.up <- gse_res$kegga$Up / gse_res$kegga$N
    gse_res$kegga$GeneRatio.down <- gse_res$kegga$Down / gse_res$kegga$N
  }

  # Return all data
  return(gse_res)

}


# Function generate_de_gene_list
##' @title Obtain lists of genes differentially expressed.
##' @description From a differential expression results dataframe, obtain list of features that are differentially
##' expressed, upregulated and downregulated according to filtering parameters
##' @return a list of 0/1 arrays with features as names indicating genes that are differentially expressed (abs slot),
##'  upregulated (up slot) and downregulated (down slot) according to filtering parameters
generate_de_gene_list <- function(
  object,
  fc = 0,
  pvalue = 1,
  feature_var = 'SYMBOL',
  fc_var = 'fc',
  pvalue_var = 'PValue',
  return_sets = c('abs', 'up', 'down')
)
{

  is_significant <- object[[pvalue_var]] <= pvalue

  res <- list(
    abs = abs(object[[fc_var]]) >= fc & is_significant,
    up = object[[fc_var]] >= fc & is_significant,
    down = object[[fc_var]] <= fc & is_significant
  )

  res <- map(res, function(x) x %>% as.integer %>% as.factor %>% set_names(object[[feature_var]]))
  res <- map(res, function(x) x[!is.na(names(x))])


  return(res[return_sets])
}


# Function generate_stats_list
##' @title Obtain lists of with statistics from an input list of features or for differentially expressed features.
##' @description From a differential expression results dataframe, obtain list of with statistics differential expression statistics :
##' fold-change, p-value, p-adj. The list will be obtained from an input list of features or from a list of features differentially expressed
##' according to inut parameters.
##' @details the list of features indicated in the parameters has priority over pvalue, fdr and pvalue
##' @return a list of numerical arrays with features as names showing the log2 fold-change (fc slot), P value (pvalue slot) or Adjusted P value (padj slot)
generate_stats_list <- function(
  object,
  fc = 0,
  pvalue = 1,
  features = NULL,
  feature_var = 'SYMBOL',
  fc_var = 'fc',
  pvalue_var = 'PValue'
)
{

  if(!is.null(features)){
    use_features <- object[[feature_var]] %in% features
    y <- object[use_features,]
  }
  if(is.null(features)){
    is_significant <- object[[pvalue_var]] <= pvalue & abs(object[[fc_var]]) >= fc
    y <- object[is_significant,]
  }
  res <- list(
    fc = y[[fc_var]],
    pvalue = y[[pvalue_var]],
    padj = y[[pvalue_var]]
  )

  res <- map(res, function(x) x  %>% set_names(y[[feature_var]]))
  res <- map(res, function(x) x[!is.na(names(x))])

  return(res)
}


# Function generate_fc_list
##' @title Obtain an array of fold-change from an input list of features or for differentially expressed features.
##' @description From a differential expression results dataframe, obtain an array of fold-changes names with feature names.
generate_fc_list <- function(
  object,
  feature_var = 'feature',
  fc_var = 'logFC'
)
{
  res <- object[[fc_var]] %>% set_names(object[[feature_var]]) %>% sort(., decreasing = T)
  res <- res[!is.na(names(res))]
  return(res)
}



# Function wrapper_topGO
##' @title Run topGO and return a lsit with results
##' @description Perform topGO analysis from an array of 0/1 named with feature name (with ID type) and return original data and s
##' ummarised table of results with gene names for each go.
##' @param ID Entrez,  GenBank,  Alias,  Ensembl,  GeneSymbol, GeneName and UniGene.
##' @param ... Other parameters for topGO::runTest.
wrapper_topGO <- function(
  gene_list,
  ID = "SYMBOL",
  ontology = 'BP',
  annot_db = 'org.Hs.eg.db',
  nodeSize = 10,
  topNodes = 100,
  pvalueCutoff = 1,
  ...
){
  require(topGO)
  require(annot_db, character.only = TRUE)
  genes_GOdata <- new( "topGOdata", ontology = ontology,
                       allGenes = gene_list, nodeSize = nodeSize,
                       annot = annFUN.org,
                       mapping = annot_db,
                       ID=ID)
  # genes_GOresult <- runTest( genes_GOdata, algorithm = "elim", statistic = "fisher" )
  genes_GOresult <- runTest( genes_GOdata, ... )
  genes_GOtab <- GenTable( genes_GOdata, elimFisher = genes_GOresult, topNodes = topNodes ) %>% dplyr::filter(elimFisher < pvalueCutoff)
  genes_GOtab$GeneRatio <- genes_GOtab$Significant / genes_GOtab$Annotated
  GO_to_genes <- genesInTerm(genes_GOdata, genes_GOtab$GO.ID) %>%
    map(., function(gn) names(gene_list[gn][gene_list[gn] == 1]) %>% paste(collapse = ';')) %>%
    unlist
  genes_GOtab$features <-  GO_to_genes[genes_GOtab$GO.ID]
  genes_GOtab$elimFisher <- as.numeric(genes_GOtab$elimFisher)

  if(ID == 'Ensembl') {
    genes_GOtab$SYMBOL <- genes_GOtab$features %>%
      map(function(x){
        x %>% strsplit(";") %>% unlist() %>%
          mapIds(org.Hs.eg.db, keys = ., keytype = "ENSEMBL", column="SYMBOL") %>% paste(collapse = ";")}
      )
  }

  return(
    list(
      GOtab = genes_GOtab,
      GOresult = genes_GOresult,
      GOdata = genes_GOdata
    )
  )

}


# Function run wrapper_ROntoTools
##' @title Run ROntoTools topological bses gene set enrichment analysis.
##' @description Perform ROntoTools Pathway-Express (PE) and Primary dis-regulatio (PDis) analysis from a given arrays
##' of fold-change and p-values. The input list of genes could be prefiltered list of genes by significance or all genes.
##' This perform the first tool to perform analysis of signaling pathways using important biological factors like all the
##' interactions between the genes, the type of interactionbetween them and the position  and  magnitude of expression change
##' for all the differentiallyexpressed genes. It uses KEGG pathway databases.
##' @param fc numeric array of fold-changes named with entrezId.
##' @param pv numeric array of pvalues named with entrezId.
##' @param kegg_organism character, name of KEGG organism id. Default kegg_organism="hsa".
##' @param ref character, array of reference genes with entrezId.
##' @param ... Other parameters for topGO::pe and topGO::pDis.
##' @return list, with the following elements:
##' pe : object of class peRes, returned from ROntoTools:pe.
##' pe_summary : data.frame with peRes results.
##' pdis : object of class pDisRes, returned from ROntoTools:pDis.
##' pdis_summary : data.frame with pDis results.
##' kegg_version : KEGG version used.
wrapper_ROntoTools <- function(fc, pv, kegg_organism, ref = NULL, ...){
  require(graph)
  require(ROntoTools)
  res <- list()

  # Generate graph for pathway database
  res$kegg_version <- keggInfo("pathway") %>% str_split(., "\n") %>% unlist %>% grep('path ', ., value = TRUE) %>% gsub("path *(.*)", "\\1", .)
  kpg <- keggPathwayGraphs(kegg_organism, verbose = TRUE)
  kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
                        edgeWeightByType = list(activation = 1, inhibition = -1, expression = 1, repression = -1),
                        defaultWeight = 0)
  kpn <- keggPathwayNames(kegg_organism)

  # add organism prefix to entrez id names
  fc <- fc %>% set_names(paste('hsa', names(.), sep = ':'))
  pv <- pv %>% set_names(paste('hsa', names(.), sep = ':'))

  # setting the node weight
  kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)

  # The reference contains all the genes measured in the analysis
  if(!is.null(ref)){
    ref <- paste('hsa', ref, sep = ':')
  } else {
    ref <- names(fc)
  }

  # Pathway-Express analysis
  res$pe <- pe(x = fc, graphs = kpg, ref = ref, ...)

  res$pe_summary <- Summary(res$pe,
                            pathNames = kpn,
                            order.by = "pPert")
  # Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE, pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert")

  # Pathway-Express analysis
  res$pdis <- pDis(x = fc, graphs = kpg, ref = ref, ...)

  res$pdis_summary <- Summary(res$pdis,
                              pathNames = kpn,
                              order.by = "pPert")

  # res$kpn <- kpn


  return(res)


  # p <- res$pe@pathways[["path:hsa04064"]]
  # g <- layoutGraph(p@map, layoutType = "dot")
  # graphRenderInfo(g) <- list(fixedsize = FALSE)
  # edgeRenderInfo(g) <- peEdgeRenderInfo(p)
  # nodeRenderInfo(g) <- peNodeRenderInfo(p)
  # renderGraph(g)
}



# Function run wrapper_GSEAPreranked
##' @title Run GSEA pre-ranked desktop
##' http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
##' GSEAPreranked calculated p-value based on gene permutation. Gene permutation breaks the structure of intergene correlations
##' within a gene set and, in fact, reflects the unrealistic assumption that genes are independent. This makes pre-ranked GSEA
##' highly prone to false positive results. On the other hand, the original GSEA method permutes samples (rather than genes),
##' thereby pre- serving intergene correlations within each gene set . Therefore, this method is likely less sensitive to the length bias.
##' However, sample permutation is only effec- tive for datasets with a large number ofreplicate samples, whereas small datasets have too
##' few samples to support their robust permutation
##'
wrapper_GSEAPreranked <- function (
  fc = NULL,
  feature_type = 'SYMBOL',
  path_dir = '.',
  msigdb_collection = c('h.all.v7.0', 'c2.cp.v7.0'),
  gsea_cli_path = '~/apps/gsea/4.0.3/gsea-cli.sh',
  gsea_ftp = 'ftp.broadinstitute.org://pub/gsea/gene_sets',
  minGSSize = 10,
  maxGSSize = 500,
  rnd_seed = 1029384756,
  plot_top_x = 20,
  nperm = 1000,
  run_script = TRUE,
  analysis_prefix = NULL,
  run_slurm = FALSE
)
{
  if(!file.exists(path_dir))
    dir.create(path_dir, recursive = TRUE)

  # RNK: Ranked list file format (*.rnk)
  rnk_data <- data.frame(gene=names(fc), FC=fc)
  if(!is.null(analysis_prefix)) {
    rnk_file <- file.path(path_dir, paste0(analysis_prefix, '.rnk'))
  } else {
    rnk_file <- tempfile(pattern = "file", tmpdir = path_dir, fileext = ".rnk")
  }
  rnk_prefix <- basename(rnk_file) %>% gsub("\\.rnk", "", .)
  write_tsv(rnk_data, path = rnk_file, col_names = FALSE)

  # Write bash script
  res_list <- list()

  for (use_collection in msigdb_collection) {

    # msigdb configuration
    if(feature_type == 'SYMBOL')
      gmt_file <- paste0(use_collection, '.symbols.gmt')
    if(feature_type == 'ENTREZID')
      gmt_file <- paste0(use_collection, '.entrez.gmt')

    # output paths
    out_path <- gsub("\\.rnk", "_gseaPreranked", rnk_file)
    # command file
    cmd_file <- gsub("\\.rnk", paste0("_", use_collection, "_gseaPreranked.sh"), rnk_file)

    sink(cmd_file)

    cat('module purge; module load Java\n')
    cat(
      gsea_cli_path, 'GSEAPreranked',
      '-rnk', rnk_file,
      '-gmx', file.path(gsea_ftp, gmt_file),
      '-collapse No_Collapse',
      '-mode Max_probe',
      '-norm meandiv',
      paste('-nperm', nperm),
      '-scoring_scheme weighted',
      '-rpt_label', use_collection,
      'order descending',
      '-create_svgs false',
      '-include_only_symbols true',
      '-make_sets true',
      paste('-plot_top_x', plot_top_x),
      paste('-rnd_seed', rnd_seed),
      paste('-set_min', minGSSize),
      paste('-set_max', maxGSSize),
      '-zip_report false',
      '-out', out_path,
      '\n'
    )
    sink()

    if(run_slurm) {
      system(paste('sbatch', cmd_file))
    } else if (run_script) {
      system(paste('bash', cmd_file))
      out_collection <- list.files(out_path, pattern = paste0(use_collection, '.GseaPreranked'), full.names = TRUE)

      res <- list()
      res$result_neg <- list.files(out_collection, pattern ='gsea_report_for_na_neg.*xls', full.names = TRUE) %>%
        read_tsv()
      res$result_pos <- list.files(out_collection, pattern ='gsea_report_for_na_pos.*xls', full.names = TRUE) %>%
        read_tsv()
      res$params_rpt <- list.files(out_collection, pattern ='.*rpt', full.names = TRUE) %>%
        read_tsv(col_names = c('param', 'value'))
      res$index_file <- list.files(out_collection, pattern ='index.html', full.names = TRUE)
      res_list[[use_collection]] <- res
    } else {
      res_list[[use_collection]] <- cmd_file
    }
  }
  return(res_list)
}
