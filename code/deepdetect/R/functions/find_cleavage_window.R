find.cleavageWindow = function(peptides, proteins, terminal.size = 7) {
  proteinIndexes = match(peptides$protein, proteins$accession)
  
  cleavages = do.call(rbind, lapply(1:length(proteinIndexes), function(i) {
    proteinSequence = proteins$sequence[proteinIndexes[i]]
    if ('start' %in% colnames(peptides)) {
      start = peptides$start[i]
    }
    else {
      start = gregexpr(proteinSequence, pattern = peptides$sequence[i], fixed = TRUE)[[1]][1]
    }
    if (start < 0) {
      stop(paste0(peptides$sequence[i], ' not found in ', proteins$accession[proteinIndexes[i]]))
    }
    nTerm = gsub(' ', '_', format(
      substring(proteinSequence, start - terminal.size, start - 1),
      width = terminal.size,
      justify = 'right'
    ))
    end = start + nchar(peptides$sequence[i]) - 1
    cTerm = gsub(' ', '_', format(
      substring(proteinSequence, end + 1, end + terminal.size),
      width = terminal.size,
      justify = 'left'
    ))
    cbind(
      start = start,
      end = end,
      nTerminal = nTerm,
      cTerminal = cTerm
    )
  }))
  cleavages = data.frame(cleavages, stringsAsFactors = FALSE)
  cleavages$start = as.integer(cleavages$start)
  cleavages$end = as.integer(cleavages$end)
  cleavages
}
