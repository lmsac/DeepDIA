library(readr)
digested.file = list.files(pattern = '_digested.*\\.txt$')[1]
digested = read_delim(digested.file, '\t', escape_double = FALSE, trim_ws = TRUE)

proteins = read.fasta(list.files(pattern = sub('_digested.*\\.txt$', '.fasta', digested.file)))

peptides = data.frame(
  protein = digested$Protein_Name,
  sequence = digested$Sequence,
  position = digested$Tryptic_Name,
  miss = as.integer(sub('^t[0-9]+\\.([0-9]+).*', '\\1', digested$Tryptic_Name)) - 1,
  stringsAsFactors = FALSE
)

proteins = data.frame(
  accession = sapply(proteins, function(x) x$name),
  sequence = sapply(proteins, function(x) x$sequence),
  description = sapply(proteins, function(x) sub('^>[^ ]* ', '', x$description)),
  stringsAsFactors = FALSE
)

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
  sub('_digested(.*)\\.txt$', '\\1.peptide.csv', digested.file),
  row.names = FALSE
)
