#' @title Calculate Attribute Scores
#' @description Calculates weighted average attribute scores for individual sensitivity attributes, with or without bootstrapping. Based on code from the Highly Migratory Species (HMS) CVA team.

#' @param attribute a data frame consisting of all expert scores for a single sensitivity attribute and species.
#' @param bootstrap binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores. Will impact outputs.
#' @param samples number of samples to run for bootstrapping. Default is 10,000.

#' @return weighted averaged scores for each attribute. If \code{bootstrap = TRUE}, this is a data frame with the number of rows equal to the number of samples, and the number of columns equal to the number of attributes. If \code{bootstrap = FALSE}, this is a vector equal to the length of the number of attributes.

calculate_attribute_score <- function(
  attribute,
  bootstrap = TRUE,
  samples = 10000
) {
  scores.only <- attribute[, (dim(attribute)[2] - 3):dim(attribute)[2]] #isolate scores

  if (bootstrap) {
    #if bootstrap is TRUE, perform bootstrap
    set.seed(2026)
    scoresB <- NULL

    if (sum(scores.only, na.rm = T) == 0) {
      #check for data (if all are NAs, sum will be equal to 0)
      scoresB <- rep(0, length = samples)
    } else {
      for (i in 1:samples) {
        ## from HMS CVA code
        #calculate number of tallies into L, M, H, V
        L <- sum(scores.only$Scoring.Rank1, na.rm = T)
        M <- sum(scores.only$Scoring.Rank2, na.rm = T)
        H <- sum(scores.only$Scoring.Rank3, na.rm = T)
        V <- sum(scores.only$Scoring.Rank4, na.rm = T)

        #create draw pile based on number of tallies in each category
        draw_pile <- c(rep(1, L), rep(2, M), rep(3, H), rep(4, V))

        #randomly sample with replacement 25 (5 tallies x 5 experts) from draw pile
        samp <- sample(
          x = draw_pile,
          size = 5 * nrow(scores.only),
          replace = TRUE
        ) #replaced 5 with 5 * nrow(scores.only)

        #table of 1,2,3,4 in samp
        samp.tab <- data.frame(table(samp))

        sampL <- ifelse(
          any(unique(samp) == 1),
          samp.tab$Freq[which(samp.tab$samp == 1)],
          0
        )
        sampM <- ifelse(
          any(unique(samp) == 2),
          samp.tab$Freq[which(samp.tab$samp == 2)],
          0
        )
        sampH <- ifelse(
          any(unique(samp) == 3),
          samp.tab$Freq[which(samp.tab$samp == 3)],
          0
        )
        sampV <- ifelse(
          any(unique(samp) == 4),
          samp.tab$Freq[which(samp.tab$samp == 4)],
          0
        )

        #calculate weighted average
        score <- ((sampL * 1) + (sampM * 2) + (sampH * 3) + (sampV * 4)) /
          (nrow(scores.only) * 5) #replaced num.experts with nrow(scores.only)
        scoresB <- c(scoresB, score)
      } #end loop
    } #end if

    return(scoresB)
  } else {
    #if bootstrap if FALSE, just get weighted mean
    score <- (sum(scores.only$Scoring.Rank1, na.rm = T) +
      sum(scores.only$Scoring.Rank2, na.rm = T) * 2 +
      sum(scores.only$Scoring.Rank3, na.rm = T) * 3 +
      sum(scores.only$Scoring.Rank4, na.rm = T) * 4) /
      sum(scores.only, na.rm = T)
    #names(score) <- attribute$Attribute.Name
    return(score)
  }
}
