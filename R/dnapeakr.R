# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# argument to znf_peaks w postaci tabelki (wczytane uprzednio przez read.table)
mapping <- function(znf_peaks, epsilon=5000) {
  ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
  result <- c()
  for (i in 1:dim(znf_peaks)[1]){
    # chr id: as numeric
    chromosome_id<-ifelse(
      stringi::stri_sub(znf_peaks[i, 1], 4)=="X",
      "X",
      as.numeric(stringi::stri_sub(znf_peaks[i, 1], 4)))
    st <- znf_peaks[i, 2] - epsilon
    en <- znf_peaks[i, 3] + epsilon
    gens<- biomaRt::getBM(attributes = c('ensembl_gene_id','chromosome_name','start_position', 'end_position'),
                 filters = c('chromosome_name','start','end'),
                 values = list(chromosome_id,st,en),
                 mart = ensembl)
    result<-rbind(result,
                 if (nrow(gens) == 0) {
                   data.frame(ensembl_gene_id=NA, chromosome_name=chromosome_id, start_position=st, end_position=en)
                   } else gens[1, ])
  }
  return (result)
}

