#' Test for motif enrichment significance
#'
#'@param target.promoters A dataframe with target promoters.
#'@param universe.promoters A dataframe with the total number of promotors to test against.
#'@param motifs An output from DNAStringSet()
#'@return A data frame of tested motifs with sorted p-values
#'@example motifEnrichment(my.targets, my.universe, all.counts = F, motifs=motifsSS)

motifEnrichment <- function(target.promoters,universe.promoters,all.counts=F,motifs=motifsSS) {
  
  #use vcountPDict to count the occurences of each motif in each promoter
  target.counts <- vcountPDict(motifs,target.promoters,fixed=F) + 
    vcountPDict(motifsSS,reverseComplement(target.promoters),fixed=F)
  universe.counts <- vcountPDict(motifs,universe.promoters,fixed=F) + 
    vcountPDict(motifsSS,reverseComplement(universe.promoters),fixed=F)

  if (all.counts) { 
    #count all occurences of a motif instead of the number of promoters that it occurs in
    target.counts.sum <- apply(target.counts,1,sum)
    universe.counts.sum <- apply(universe.counts,1,sum)
  } else {
    target.counts.sum <- apply(ifelse(target.counts > 0,1,0),1,sum)
    universe.counts.sum <- apply(ifelse(universe.counts > 0 , 1, 0),1,sum)
  }
  n.motifs <- length(target.counts.sum)
  results <- vector(mode="numeric",length=n.motifs)
  for (i in 1:n.motifs) {
    if (all.counts) { #the contigency tables are different depending on whether we are looking at promoters or overall occurences
      #test if ratio of occurences to promoters is the same in the target and the universe
      m <- matrix(c(
        target.counts.sum[i],                       #number of occurences within target
        dim(target.counts)[2],                      #number of promoters in target
        universe.counts.sum[i],                  #number of occurences within universe
        dim(universe.counts)[2]                  #number of promoters in universe
      ),ncol=2)
    } else { #looking at promoters with and without hits
      m <- matrix(c(
        target.counts.sum[i],                            #number of promoters in target with hit
        dim(target.counts)[2]-target.counts.sum[i],      #number of promoters in target with no hit
        universe.counts.sum[i],                          #number of promoters in universe with hit
        dim(universe.counts)[2]-universe.counts.sum[i]   #number of promoters in universe with no hit
      ),ncol=2)
    } #else
    results[i] <- fisher.test(m,alternative="greater")$p.value
  } #for loop
  results.table <- data.frame(
    motif=names(motifs),
    universe.percent = round(universe.counts.sum/dim(universe.counts)[2],3)*100,
    target.percent = round(target.counts.sum/dim(target.counts)[2],3)*100,
    p.value =  results)
  results.table <- results.table[order(results.table$p.value),]
  results.table
}