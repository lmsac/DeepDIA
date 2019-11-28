library(readr)

digested.files = list.files(pattern = '_digested.*\\.txt$')

peptides = do.call(rbind, lapply(digested.files, function(file) {
  read_delim(file, delim = '\t', trim_ws = TRUE)
}))

peptides = data.frame(
  protein = peptides$Protein_Name,
  sequence = peptides$Sequence,
  stringsAsFactors = FALSE
)

peptides = unique(peptides)

peptides = peptides[nchar(peptides$sequence) <= 50, ]
peptides = peptides[!grepl('[^ACDEFGHIKLMNPQRSTVWY]', peptides$sequence), ]

write.csv(
  peptides,
  sub('_digested.*\\.txt$', '.peptide.csv', digested.files[1]),
  row.names = FALSE
)
