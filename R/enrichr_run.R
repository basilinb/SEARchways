



##' onLoad hook to setup package options
##'
##' onLoad hook to setup package options and to check connection to website
##' @title onLoad hook to setup package options
##' @param libname (Required). Library name
##' @param pkgname (Required). Package name
##' @return NULL
##' @author Wajid Jawaid \email{wj241@alumni.cam.ac.uk}
.onAttach <- function(libname, pkgname) {
  options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  options(enrichR.live = TRUE)
  packageStartupMessage("Welcome to enrichR\nChecking connection ... ", appendLF = TRUE)
  options(modEnrichR.use = TRUE)
  options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr", "OxEnrichr"))
  if (getOption("modEnrichR.use")) {
    listEnrichrSites()
  } else {
    getEnrichr(url=paste0(getOption("enrichR.base.address"), "datasetStatistics"))
    packageStartupMessage("Enrichr ... ", appendLF = FALSE)
    if (getOption("enrichR.live")) packageStartupMessage("Connection is Live!")
  }
}


##' Helper function
##'
##' Helper function for GET
##' @title Helper function for GET
##' @param url (Required). URL address requested
##' @param ... (Optional). Additional parameters to pass to GET
##' @return same as GET
##' @author Wajid Jawaid \email{wj241@alumni.cam.ac.uk}
##' @author I-Hsuan Lin \email{i-hsuan.lin@manchester.ac.uk}
##' @importFrom httr GET
##' @importFrom httr status_code
##' @importFrom httr http_status
getEnrichr <- function(url, ...) {
  options(enrichR.live = FALSE)
  tryCatch({
    x <- GET(url = url, ...)
    code <- status_code(x)
    if(code != 200) {
      # Error with status code
      message(http_status(code)$message)
    } else {
      # OK/success
      options(enrichR.live = TRUE)
      invisible(x)
    }
  },
  # Warning message
  warning = function(warn) {
    message(warn); message("") # force newline
  },
  # Error without status code
  error = function(err) {
    message(err); message("") # force newline
  },
  finally = function() {
    invisible(x)
  })
}

##' List modEnrichr Websites
##'
##' List Enrichr Websites
##' @title List Enrichr Websites
##' @return print Enrichr Website status
##' @author Alexander Blume
##' @param ... (Optional  Additional parameters)
##' @export
listEnrichrSites <- function(...) {
  for (site in getOption("enrichR.sites")) {
    getEnrichr(url = paste0(getOption("enrichR.sites.base.address"), site, "/", "datasetStatistics"))
    packageStartupMessage(paste0(site, " ... "), appendLF = FALSE)
    if (paste0(getOption("enrichR.sites.base.address"), site, "/")  == getOption("enrichR.base.address")) {
      if (getOption("enrichR.live")) packageStartupMessage("Connection is Live!")
    } else
      if (getOption("enrichR.live")) packageStartupMessage("Connection is available!")

  }
}










#' Run EnrichR function on multiple Geneset databases on either gene list or modules of intertest
#'
#' @param gene_mod_list Here input either Genelist in HGNC or Ensembl ID format or input module object output from module function
#' @param type Type can be either gene or module
#' @param gene_id gene id is either HGNC or Ensembl
#' @param dbs Enter the geneset database you want to run Enrichr on the default is MSigDB Hallmark 2020 to check available database run enrichR::listEnrichrDb()
#' @param mod_of_interest Input your modules of interest by default will run modules

#' @return list with enrichr result for each geneset database
#'
#' @export
#'
#' @examples
enrichr_run <- function(gene_mod_list, type = "gene",gene_id = "HGNC",dbs = c("MSigDB_Hallmark_2020"),mod_of_interest = "all") {
  #Set variables to Null
  gene_biotype <- module <- ensembl_gene_id <- hgnc_symbol <- NULL
  #Download Ensembl gene list to get HGNC symbols
  ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  all_genes <- biomaRt::getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol'), mart = ensembl) %>%
    dplyr::filter(gene_biotype == "protein_coding")
  #Setting the website for Enrichr to run
  enrichR::setEnrichrSite("Enrichr")
  #Defining two empty lists
  enrich_list <- list()
  gene_list <- list()

  #Checking if user wants to run on gene list or modules and formatting the data
  if (type!="gene"){
    if(mod_of_interest == "all"){
      mod_of_interest = gene_mod_list$mods %>% dplyr::pull(module) %>% unique() %>% sort()
    }
  } else {
    if (gene_id == "HGNC") {
      gene_list <- gene_mod_list
    }else{
      # Check to see if its HGNC symbol or Ensembl ID and formatting the data
      gene_list <- all_genes %>% dplyr::filter(ensembl_gene_id %in% gene_mod_list) %>% dplyr::pull(hgnc_symbol)
    }
  }

  #Loop through the databases provided by user
  for (d in dbs){
    if (type == "gene"){
      #Running Enrichr
      enriched <- enrichR::enrichr(gene_list,d)
      # Append to our main list we want to return, So if we provide multiple databases this will append them to the main list.
      enrich_list <-  append(enrich_list,enriched)
    }
    else {
      #Loop through the modules and genes within the modules
      for (m in mod_of_interest) {
        #Getting gene list from the module
        gene_list <- gene_mod_list$mods %>% subset(module == m) %>% dplyr::pull(hgnc_symbol)
        #Running Enrichr
        enriched <- enrichR::enrichr(gene_list,d)
        #Append to our main list we want to return, So if we provide multiple databases this will append them to the main list.
        enrich_list[[paste("module", m,sep="_")]] <- enriched
      }
    }
  }
  return(enrich_list)
}
