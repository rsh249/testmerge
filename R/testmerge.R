## R function for merging paired end DNA sequence data when you *know* that all pairs should overlap. That is, when fragment length is smaller than 2x read length.


# broken, proably won't fix
mergeReads = function(R1, R2, minoverlap=10) {
  require(Biostrings, quietly=TRUE)
  require(stringr, quietly=TRUE)
  r1 = readQualityScaledDNAStringSet(path2R1)
  r2 = readQualityScaledDNAStringSet(path2R2)
  #complement R2
  r2 = reverseComplement(r2)
  reads = data.frame(
    l1 = 1:length(r1),
    l2 = 1:length(r2),
    lm = 1:length(r1),
    seq = 1:length(r1)
  )
  for (i in 1:length(r1)) {
    dist=vector()
    nr1 = nchar(r1[[i]])
    nr2 = nchar(r2[[i]])
    p=proc.time()
    #maxlen = 0
    for(z in minoverlap:length(r1[[i]])){
      forw = (str_sub(r1[[i]], -z))
      reve = (str_sub(r2[[i]], 1, z))
      dist[[z]] = adist(forw, reve)
      currlen = nchar(dist[[z]])
      #if (currlen > maxlen) { maxlen = currlen }
      #if (currlen < maxlen && currlen < minsincemax) { minsincemax = currlen }
      if (z > 50){
        if (dist[[z-30]] > dist[[z - 10]] && dist[[z]] > dist[[z-10]]) { optimstatus = 1; break }
      }
      
    }
    if (optimstatus == 0 ) { print("No overlap found"); next }
    
    minedit = which(dist/1:length(dist) == min(dist/1:length(dist), na.rm=T))
    for.unmatched = str_sub(r1[[i]], 1, nchar(r1[[i]]) - minedit)
    for.matched = str_sub(r1[[i]], -minedit)
    rev.matched = str_sub(r2[[i]], 1,minedit)
    rev.unmatched = str_sub(r2[[i]], -(nchar(r2[[i]]) - minedit))
    mcat = paste(for.unmatched, 
                 for.matched, 
                 rev.unmatched, sep='')
    proc.time() - p
    #print(mcat)
    #can get ties in read fuzzy alignment resulting in more than one mcat element
    #should consider what happens if no overlap is actually present (then min distance is probably at starting overlap)
    reads$l1[i] = nchar(r1[[i]])
    reads$l2[i] = nchar(r2[[i]])
    reads$lm[i] = nchar(mcat[1])
    reads$seq[i] = mcat[1]
  }
  return(reads)
}
