read.fasta = function(file) {
  lines = readLines(file)
  start.indexes = grep('^>', lines)
  if (length(start.indexes) == 0) {
    return(list())
  }
  end.indexes = c(start.indexes[-1] - 1, length(lines))
  
  lapply(1:length(start.indexes), function(i) {
    first.line = lines[start.indexes[i]]
    s = strsplit(first.line, ' ')[[1]]
    name = sub('^>', '', s[1])
    description = first.line
    
    if (end.indexes[i] > start.indexes[i]) {
      sequence = paste0(trimws(lines[(start.indexes[i] + 1):(end.indexes[i])]), collapse = '')
    }
    else {
      sequence = ''
      warning(paste0('Line ', start.indexes[i], ': empty sequence.'))
    }
    list(
      name = name,
      sequence = sequence,
      description = description
    )
  })
}

# read.fasta = function(file) {
#   fasta = seqinr::read.fasta(file, seqtype = 'AA', as.string = TRUE)
#   lapply(fasta, function(entry) {
#     list(name = attr(entry, 'name'),
#          sequence = entry[1],
#          description = attr(entry, 'Annot'))
#   })
# }