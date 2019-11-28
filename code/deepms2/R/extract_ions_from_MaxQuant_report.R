charge = 2:3

msms = do.call(rbind, lapply(
  list.files(pattern = 'msms.txt', full.names = T, recursive = TRUE), 
  function(msmsFile) {
    msms = read.delim(
      msmsFile, header = TRUE, sep = '\t', 
      fileEncoding = 'UTF-8', 
      stringsAsFactors = FALSE
    )
    msms = msms[which(
      msms$Reverse != '+' &
        msms$Modifications == 'Unmodified' &
        msms$Charge %in% charge &
        msms$Type %in% c('MSMS', 'MULTI-MSMS')
    ), ]
    msms = cbind(Result.file = msmsFile, msms)
  })
)


lapply(charge, function(ch) {
  indexs = which(msms$Charge == ch)
  indexs = indexs[order(msms$Score[indexs], decreasing = TRUE)]
  indexs = indexs[!duplicated(msms$Sequence[indexs])]
  
  ions = lapply(indexs, function(i) {
    peptide = msms[i, ]$Sequence
    charge = msms[i, ]$Charge
    score = msms[i, ]$Score
    assigned = msms[i, ]$Peak.coverage
    
    matches = strsplit(msms[i, ]$Matches, ';')[[1]]
    intensities = as.numeric(strsplit(msms[i, ]$Intensities, ';')[[1]])
    
    if (length(matches) == 0 || length(intensities) == 0 || all(intensities == 0))
      return(NULL)
    
    peptideLength = nchar(peptide)
    ionTypes = c(
      'b1', 'b2', 'bn1', 'bn2', 'bo1', 'bo2', 
      'y1', 'y2', 'yn1', 'yn2', 'yo1', 'yo2'
    )
    ions = lapply(ionTypes, function(x) {
      rep(0, nchar(peptide) - 1)
    })
    names(ions) = ionTypes
    
    lapply(1:length(matches), function(i) {
      label = matches[i]
      type = substring(label, 1, 1)
      label = substring(label, 2)
      charge = if (endsWith(label, '(2+)')) 2 else 1
      label = sub('\\(2\\+\\)$', '', label)
      nh3 = if (endsWith(label, '-NH3')) TRUE else FALSE
      label = sub('-NH3$', '', label)
      h2o = if (endsWith(label, '-H2O')) TRUE else FALSE
      label = sub('-H2O$', '', label)
      number = suppressWarnings(as.integer(label))
      
      if (!is.na(number)) {
        index = if (type %in% c('b')) 
          number
        else if (type %in% c('y')) 
          peptideLength - number
        else
          number
        
        label = paste0(
          type,
          if (nh3) 'n' else '',
          if (h2o) 'o' else '',
          charge
        )
        
        if (label %in% ionTypes) {
          ions[[label]][index] <<- intensities[i]
        }
      }
    })
    
    if (all(unlist(ions) == 0))
      return(NULL)
    
    list(
      peptide = peptide,
      charge = charge,
      score = score,
      assigned = assigned,
      ions = ions,
      resultFile = msms[i, ]$Result.file,
      rawFile = msms[i, ]$Raw.file,
      scan = msms[i, ]$Scan.number
    )
  })
  ions = ions[!sapply(ions, is.null)]
  
  writeLines(rjson::toJSON(ions), sub('XX', ch, 'chargeXX.ions.json'))
  
  message(paste0('Charge ', ch, ': ', length(ions), ' peptides.'))
})
