#' Run GSEA on multiple list using fgsea
#'
#' @param gene_list named list object with named numeric vectors of gene symbols and logFC
#' @param gmt_ls list object with gene ontology data
#' @param nperm number of permutations for P-value calculations. Default is 1000
#' @param group Which Geneset database you want to run GSEA on
#' @param subgroup Specify subgroup of Geneset DB you want to run GSEA on \cr
#' \tabular{rrrrr}{
#'  \strong{Group} \tab    \strong{Sub-group}\cr
#'  C1  \tab      \cr
#'  C2  \tab      CGP\cr
#'  C2  \tab      CP\cr
#'  C2  \tab      CP:BIOCARTA\cr
#'  C2  \tab      CP:KEGG\cr
#'  C2  \tab      CP:PID\cr
#'  C2  \tab      CP:REACTOME\cr
#'  C2  \tab      CP:WIKIPATHWAYS\cr
#'  C5  \tab      GO:BP\cr
#'  C5  \tab      GO:CC\cr
#'  C5  \tab      GO:MF\cr
#'  C5  \tab      HPO\cr
#'  C6  \tab     \cr
#'  H   \tab      \cr
#'  }
#' @param S Specify specie for Geneset database default is human
#'
#' @return Output of a list of list with GSEA results
#' @export
#'
#' @examples
GSEA_run <- function(gene_list, gmt_ls=NULL, nperm=1000,group = NULL, subgroup = NULL, S = "human"){
  pathway <- pval <- padj <- ES <- NES <- size <- leadingEdge <- fgsea.FC  <- NULL
  #Blank list to hold results
  all.results <- list()

  #### Data ####
  #Load gene ontology
  if(!is.null(group) & !is.null(subgroup)){
    #myGO <- fgsea::gmtPathways(gmt_file)
    gene_sets = msigdbr::msigdbr(species = S, category = group, subcategory = subgroup)
    myGO <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  } else if(!is.null(group) & is.null(subgroup)){
    gene_sets = msigdbr::msigdbr(species = S, category = group)
    myGO <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
  }else if(!is.null(gmt_ls)){
    myGO <- gmt_ls
  } else {
    stop("Please provide gene set data as file for list object.")
  }

  #### Loop ####
  #Loop through each list in the gene_list object
  for(genes in names(gene_list)){
    message(genes)
    #Extract 1 gene list
    genes.temp <- gene_list[[genes]]
    #Order by fold change
    genes.temp <- sort(genes.temp, decreasing = TRUE)

    #### FGSEA ####
    #Set score type based on fold change
    if(min(genes.temp) < 0 & max(genes.temp) > 0){
      scoreType <- "std"
    } else if(max(genes.temp) <= 0){
      scoreType <- "neg"
    } else if(min(genes.temp) >= 0){
      scoreType <- "pos"
    } else{
      stop("Could not determine score type from fold changes.")
    }

    #Run GSEA with fgsea
    fg.result <- fgsea::fgseaSimple(pathways = myGO,
                                    stats = genes.temp,
                                    nperm=nperm,
                                    #eps=0,
                                    scoreType=scoreType) %>%
      as.data.frame()




    #### Combine results ####
    gsea.result <- fg.result %>%
      dplyr::select(pathway, pval, padj, ES, NES, size, leadingEdge) %>%
      dplyr::mutate(fgsea.FC = ifelse(NES < 0, "down","up")) %>%
      #format leading edge list
      tidyr::unnest(cols = c(leadingEdge)) %>%
      dplyr::group_by(pathway,pval,padj,ES,
               NES,size,fgsea.FC) %>%
      dplyr::summarise(leadingEdge=paste(unique(leadingEdge), collapse=";"),
                       .groups="drop")

    #### Save ####
    all.results[[genes]] <- gsea.result
  }

  #### Format output ####
  #Unlist results into 1 df



  if(group!= "c6_GSEA.result"){
    all.results.df <- do.call(rbind.data.frame, all.results) %>%
      tibble::rownames_to_column("group") %>%
      dplyr::mutate(group = gsub("[.][0-9]{0,4}","",group)) %>% dplyr::mutate(pathway = sub("[A-Z]*_","",pathway)) %>% dplyr::mutate(pathway = gsub("_"," ",pathway))
  } else {
    all.results.df <- do.call(rbind.data.frame, all.results) %>%
      tibble::rownames_to_column("group") %>%
      dplyr::mutate(group = gsub("[.][0-9]{0,4}","",group))
  }

  return(all.results.df)
}
