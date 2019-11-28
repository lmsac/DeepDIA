library(readr)
report.file = list.files(pattern = '\\.PeptideReport\\.csv$')[1]
peptideReport = read_csv(report.file)

peptideReport.HPRP = peptideReport[grepl('HPRP', peptideReport$R.FileName), ]
peptideReport = peptideReport[!grepl('HPRP', peptideReport$R.FileName), ]

peptideQuantity = local({
  runs = unique(peptideReport$R.FileName)
  peptideQuantity = Reduce(
    function(x, y) list(
      run = '',
      data = merge(
        x$data, y$data, 
        by = 'sequence', all = TRUE
      )
    ), 
    lapply(runs, function(run) {
      peptideReport = peptideReport[peptideReport$R.FileName == run, ]
      list(
        run = run,
        data = data.frame(
          sequence = peptideReport$PEP.StrippedSequence,
          intensity = peptideReport$`PEP.Label-Free Quant`,
          stringsAsFactors = FALSE
        )
      )
    })
  )$data
  colnames(peptideQuantity)[-1] = runs
  peptideQuantity
})

peptides = local({
  sequence = peptideQuantity$sequence
  intensity = apply(peptideQuantity[, -1], 1, function(x) median(x, na.rm = TRUE))

  rowIndexes = match(sequence, peptideReport$PEP.StrippedSequence)
  missCleavages = peptideReport$PEP.NrOfMissedCleavages[rowIndexes]
  qValue = peptideReport$PEP.QValue[rowIndexes]
  protein = sapply(strsplit(peptideReport$PG.ProteinAccessions[rowIndexes], ';'), function(x) x[1])
  start = sapply(strsplit(peptideReport$PEP.StartingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
  end = sapply(strsplit(peptideReport$PEP.EndingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))

  data.frame(
    protein = protein,
    sequence = sequence,
    start = start,
    end = end,
    missCleavages = missCleavages,
    intensity = intensity,
    stringsAsFactors = FALSE
  )
})

peptides = local({
  peptides = peptides[order(peptides$protein), ]
  protein.indexes = which(!duplicated(peptides$protein))
  peptides$relativeIntensity = do.call(c, lapply(1:length(protein.indexes), function(i) {
    if (i < length(protein.indexes)) {
      rowIndexes = protein.indexes[i]:(protein.indexes[i + 1] - 1)
    }
    else {
      rowIndexes = protein.indexes[i]:nrow(peptides)
    }
    
    intensity = peptides$intensity[rowIndexes]
    relativeIntensity = intensity / max(intensity)
  }))
  peptides
})
peptides$detectability = pmax((5 + log10(peptides$relativeIntensity)) / 5 * 0.5, 0) + 0.5


peptides.HPRP = local({
  peptideReport = peptideReport.HPRP
  
  sequence = setdiff(peptideReport.HPRP$PEP.StrippedSequence, peptides$sequence)
  
  rowIndexes = match(sequence, peptideReport$PEP.StrippedSequence)
  missCleavages = peptideReport$PEP.NrOfMissedCleavages[rowIndexes]
  qValue = peptideReport$PEP.QValue[rowIndexes]
  protein = sapply(strsplit(peptideReport$PG.ProteinAccessions[rowIndexes], ';'), function(x) x[1])
  start = sapply(strsplit(peptideReport$PEP.StartingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
  end = sapply(strsplit(peptideReport$PEP.EndingPositions[rowIndexes], ';'), function(x) as.integer(gsub('^\\(([0-9]+)\\).*', '\\1', x[1])))
  
  value = min(peptides$relativeIntensity) / 2
  detectability = pmax((5 + log10(value)) / 5 * 0.5, 0) + 0.5
  
  data.frame(
    protein = protein,
    sequence = sequence,
    start = start,
    end = end,
    missCleavages = missCleavages,
    intensity = NA,
    relativeIntensity = NA,
    detectability = detectability,
    stringsAsFactors = FALSE
  )
})


proteinReport = read_csv(sub('\\.PeptideReport\\.csv$', '.ProteinReport.csv', report.file))

proteinAccession = local({
  coverageThreshold = 25
  coverage = sapply(strsplit(proteinReport$PG.Coverage, ';'), function(x) as.numeric(sub('%', '', x[1])))
  unique(sapply(strsplit(proteinReport$PG.ProteinAccessions[which(
    proteinReport$PG.UniquePeptides > 1 &
      coverage >= coverageThreshold
  )], ';'), function(x) x[1]))
})

writeLines(
  proteinAccession, 
  paste0(sub('\\.PeptideReport\\.csv', '', report.file), '_excludeSingleHit_coverage25.proteinAccession.txt')
)


peptides = peptides[peptides$protein %in% proteinAccession, ]

peptides.HPRP = peptides.HPRP[peptides.HPRP$protein %in% proteinAccession, ]


peptides = rbind(peptides, peptides.HPRP)

peptides = local({
  proteinHit = table(peptides$protein)
  proteins = names(proteinHit)[proteinHit > 1]
  peptides[peptides$protein %in% proteins, ]
})

write.csv(
  peptides, 
  paste0(sub('\\.PeptideReport\\.csv', '', report.file), '_excludeSingleHit_coverage25.detectability.csv'),
  row.names = FALSE
)
