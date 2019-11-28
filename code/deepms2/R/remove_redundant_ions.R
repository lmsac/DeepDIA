lapply(2:3, function(charge) {
  files = list.files(pattern = paste0('_charge', charge, '\\.ions\\.RData$'))
  
  extracted.ions = do.call(c, lapply(files, function(file) {
    load(file)
    extracted.ions
  }))
  
  indexes = local({
    scores = sapply(extracted.ions, function(x) -x$qvalue)
    rank = order(scores, decreasing = TRUE)
    dup = duplicated(sapply(extracted.ions[rank], function(x) x$peptide))
    uni = rank[!dup]
    filtered = uni # [scores[uni] <= 0.01]
    indexes = filtered[order(filtered)]
  })
  
  extracted.ions = extracted.ions[indexes]

  writeLines(
    rjson::toJSON(extracted.ions, indent = 2), 
    sub('\\.RData$', '.json', files[1])
  )
})