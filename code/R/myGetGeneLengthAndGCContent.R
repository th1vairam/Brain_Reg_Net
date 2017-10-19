myGetGeneLengthAndGCContent <- function(id){
  
  org = 'hsa'
  id.type = 'ensembl_gene_id'
  host = 'dec2016.archive.ensembl.org'
  
  message("Connecting to BioMart ...")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = host)
  ds <- listDatasets(ensembl)[, "dataset"]
  ds <- grep(paste0("^", org), ds, value = TRUE)
  
  if (length(ds) == 0){
    stop(paste("Mart not found for:", org))
  } else if (length(ds) > 1) {
    message("Found several marts")
    sapply(ds, function(d) message(paste(which(ds == d), d, sep = ": ")))
    n <- readline(paste0("Choose mart (1-", length(ds), ") : "))
    ds <- ds[as.integer(n)]
  }
  
  ensembl <- useDataset(ds, mart = ensembl)
  message(paste0("Downloading sequence", ifelse(length(id) > 1, "s", ""), " ..."))
  if (length(id) > 100)
    message("This may take a few minutes ...")
  attrs <- c(id.type, "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end")
  coords <- getBM(filters = id.type, attributes = attrs, values = id, mart = ensembl)
  id <- unique(coords[, id.type])
  coords <- sapply(id, function(i) {
    i.coords <- coords[coords[, 1] == i, 3:5]
    g <- GRanges(i.coords[, 1], IRanges(i.coords[, 2], i.coords[, 3]))
    return(g)
  })
  coords <- lapply(coords[id], reduce)
  len <- plyr::ldply(coords, function(x) sum(IRanges::width(x)), .id = 'ensembl_gene_id') %>%
    dplyr::rename(gne.length = V1)
  
  gc.cont <- getBM(filters = id.type, attributes = c(id.type, 'hgnc_symbol', 'percentage_gc_content'), values = id, mart = ensembl)
  
  res <- dplyr::full_join(gc.cont, len)
  return(res)
}
