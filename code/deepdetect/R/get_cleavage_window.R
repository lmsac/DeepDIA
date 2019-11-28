library(readr)

negative.detectability.file = list.files(pattern = '_negative\\.detectability\\.csv$')[1]
detectability.file = list.files(pattern = sub(
  '_negative\\.detectability\\.csv$', '\\.detectability\\.csv$', negative.detectability.file
))[1]

proteins = read.fasta(list.files(pattern = sub('\\.detectability\\.csv$', '.fasta', detectability.file)))
proteins = data.frame(
  accession = sapply(proteins, function(x) x$name),
  sequence = sapply(proteins, function(x) x$sequence),
  description = sapply(proteins, function(x) sub('^>[^ ]* ', '', x$description)),
  stringsAsFactors = FALSE
)

lapply(c(detectability.file, negative.detectability.file), function(f) {
  peptides = read_csv(f)
  if (is.na(match(peptides$protein[1], proteins$accession)) &&
      grepl('^[A-Za-z0-9_]+\\|([A-Za-z0-9_]+)\\|.*', proteins$accession[1]) &&
      !grepl('^[A-Za-z0-9_]+\\|([A-Za-z0-9_]+)\\|.*', peptides$protein[1])) {
    proteins$accession = sub('^[A-Za-z0-9_]+\\|([A-Za-z0-9_]+)\\|.*', '\\1', proteins$accession)
  }
  
  cleavages = find.cleavageWindow(peptides, proteins)
  peptides = cbind(peptides, cleavages)
  
  peptides = peptides[which(
    !grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptides$sequence) &
      !grepl('[^_ACDEFGHIKLMNPQRSTVWY]', peptides$nTerminal) &
      !grepl('[^_ACDEFGHIKLMNPQRSTVWY]', peptides$cTerminal)
  ), ]
  
  peptides = peptides[nchar(peptides$sequence) <= 50, ]
  
  write.csv(
    peptides,
    f,
    row.names = FALSE
  )
})
