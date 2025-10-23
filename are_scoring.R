## An implementation of the algorithm to score ARE's as described here:
## http://arescore.dkfz.de/info.html
## http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002433
library(Biostrings)
library(IRanges)

## Returns an integer vector containing number of matches for each pattern (MIndex object)
# previously in Biostrings, but now deprecated?
elementLengths <- function(mi_index_obj) {
  
  return(length(mi_index_obj))
}


## basal=1;overlapping=1.5;d1.3=.75;d4.6=0.4;d7.9=0.2;within.AU=0.3; aub.min.length=20; aub.p.to.start=0.8; aub.p.to.end=0.55
AREscore <- function(x, basal=1.0, overlapping=1.5, d1.3=0.75, d4.6=0.4,
                     d7.9=0.2, within.AU=0.3,
                     aub.min.length=20, aub.p.to.start=0.8, aub.p.to.end=0.55, pentamer="ATTTA", overmer="ATTTATTTA") {
  #xtype <- match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))
  xtype <- ifelse(grepl("T", x), "DNA", "RNA")
  
  if (xtype == "RNA") {
    pentamer <- gsub("T", "U", "ATTTA")
    overmer <- gsub("T", "U", "ATTTATTTA")
}

  x <- as(x, sprintf("%sString", xtype))

  pmatches <- matchPattern(pentamer, x)
  omatches <- matchPattern(overmer, x)

  basal.score <- elementLengths(pmatches) * basal
  over.score <- elementLengths(omatches) * overlapping
  

  no.cluster <- data.frame(d1.3=0, d4.6=0, d7.9=0)
  
  if (length(pmatches) < 2) {
    clust <- no.cluster
  } else {
    wg <- width(gaps(pmatches))
    clust <- data.frame(d1.3=sum(wg <= 3), d4.6=sum(wg >= 4 & wg <= 6),
               d7.9=sum(wg >= 7 & wg <= 9))
  }
    

  dscores <- clust$d1.3 * d1.3 + clust$d4.6 * d4.6 + clust$d7.9 * d7.9

  au.blocks <- identifyAUBlocks(x, aub.min.length, aub.p.to.start, aub.p.to.end)
  aub.score <- sum(countOverlaps(pmatches, au.blocks) * within.AU)
  
  #print(basal.score)
  #print(over.score)
  #print(aub.score)
  score <- basal.score + over.score + dscores + aub.score
  #print(score)
  #print(data.frame(clust))
  #ans <- DataFrame(score=score, n.pentamer=elementLengths(pmatches),
  #                 n.overmer=elementLengths(omatches), au.blocks=list(au.blocks),
  #                 n.au.blocks=length(au.blocks))
  #returnval <- cbind(ans, DataFrame(clust))
 
  return(score)
}


##' In order to account for an AU-rich context, AREScore identifies AU-blocks as
##' regions that are generally rich in As and Us. It does so by sliding a window
##' whose size is defined by \code{min.length} NTs over the entire length
##' of the input sequence, one nucleotide at a time.
##'
##' For each window, the algorithm calculates percent AU, and if the AU
##' percentage is equal to or greater than the number specified in
##' \code{p.to.start}, it marks the position of the first nucleotide of that
##' window as the beginning of a new AU-block. The algorithm continues scanning
##' downstream until it discovers a window with AU percentage equal to or
##' smaller than \code{p.to.end} (their implementation looks like its only
##' smaller than!). The last nucleotide of this window is marked
##' as the last nucleotide of the AU-block.
identifyAUBlocks <- function(x, min.length=20, p.to.start=0.8, p.to.end=0.55) {
  #xtype <- match.arg(substr(class(x), 1, 3), c("DNA", "RNA"))
  xtype <- ifelse(grepl("T", x), "DNA", "RNA")
  
  stopifnot(isSingleNumber(min.length) && min.length >= 5 && min.length <= 50)
  stopifnot(isSingleNumber(p.to.start) && p.to.start >= 0.50 &&
            p.to.start <= 0.95)
  stopifnot(isSingleNumber(p.to.end) && p.to.end >= 0.20 && p.to.end <= 0.70)
  stopifnot(p.to.start > p.to.end)

  if (xtype == "DNA") {
    AU <- "AT"
  } else {
    AU <- "AU"
  }

  x <- as(x, sprintf("%sString", xtype))

  au.freq <- letterFrequencyInSlidingView(x, min.length, AU, as.prob=TRUE)
  
  starts <- c()
  stops <- c()
  
  curr_start <- 1
  curr_end <- 1
  
  while (curr_end <= nrow(au.freq)) {
    valid_start <- au.freq[curr_start] >= p.to.start 
    valid_end <- au.freq[curr_end] <= p.to.end
    
    if (valid_start) {
    if (valid_end) {
      starts <- c(starts, curr_start)
      stops <- c(stops, curr_end)
      curr_start <- curr_end+1
      curr_end <- curr_end+1
    } else {
      curr_end <- curr_end+1
    }
    } else {
    curr_start <- curr_start+1
    curr_end <- curr_start
    }
  }
  starts <- c(starts, curr_start)
  stops <- c(stops, curr_end)
  
  blocks <- IRanges(starts, stops+min.length-1)

  return(blocks)
}
