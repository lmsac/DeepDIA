library(readr)
detectability.file = list.files(pattern = '\\.detectability\\.csv$')[1]
peptides = read_csv(detectability.file)

digested.file = list.files(pattern = sub('\\.detectability\\.csv$', '_digested.*\\.txt$', detectability.file))[1]
digestedPeptides = read_delim(
  digested.file, 
  '\t', escape_double = FALSE, trim_ws = TRUE
)

peptides.negative = local({
  digestedPeptides = digestedPeptides[!grepl('[^ACDEFGHIKLMNPQRSTVWY]', digestedPeptides$Sequence), ]
  digestedPeptides = digestedPeptides[nchar(digestedPeptides$Sequence) <= 50, ]
  
  digestedPeptides = digestedPeptides[!duplicated(digestedPeptides$Sequence), ]
  
  indexes = which(!(digestedPeptides$Sequence %in% peptides$sequence))
  
  data.frame(
    protein = digestedPeptides$Protein_Name[indexes],
    sequence = digestedPeptides$Sequence[indexes],
    detectability = 0,
    stringsAsFactors = FALSE
  )
})

peptides.negative$protein = sub('^[A-Za-z0-9]+\\|([A-Z0-9]+)\\|.*', '\\1', peptides.negative$protein)

write.csv(
  peptides.negative, 
  sub('\\.detectability\\.csv$', '_negative.detectability.csv', detectability.file),
  row.names = FALSE
)
