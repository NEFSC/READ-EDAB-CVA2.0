#' @title Calculate Directionality Score
#' @description Calculates weighted average directionality scores for individual sensitivity attributes, with or without bootstrapping. Based on code from the Highly Migratory Species (HMS) CVA team.

#' @param species a data frame consisting of all expert directionality scores for a single species.
#' @param bootstrap binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores. Will impact outputs.
#' @param samples number of samples to run for bootstrapping. Default is 10,000.

#' @return weighted averaged directionality score. If \code{bootstrap = TRUE}, this is a vector containing bootstrapped scores, with a length equal to the number of samples. If \code{bootstrap = FALSE}, this is a single score.

calculate_directionality <- function(species, bootstrap = TRUE, samples = 10000){
  scores.only <- species[,(dim(species)[2]-2):
                             dim(species)[2]] #isolate scores

  if(bootstrap){ #if bootstrap is TRUE, perform bootstrap
    set.seed(2026)
    scoresB <- NULL

    if(sum(scores.only,na.rm = T) == 0){ #check for data (if all are NAs, sum will be equal to 0)
      scoresB <- rep(0, length=samples)
    }else{
      for(i in 1:samples){
        #calculate number of tallies into Neg(itive), Neu(trals), Pos(itive)
        num_negs<-sum(scores.only$Negative, na.rm = T)
        num_neu<-sum(scores.only$Neutral, na.rm = T)
        num_pos<-sum(scores.only$Positive, na.rm = T)

        #create draw pile based on number of tallies in each category
        draw_pile<-c(rep("neg",num_negs),rep("neu",num_neu),rep("pos",num_pos))

        #randomly sample with replacement 20 (4 tallies x 5 experts) from draw pile
        samp<-sample(x=draw_pile,size=4 * nrow(scores.only),replace=TRUE) #replaced 4 with 5 * nrow(scores.only)

        #table of 1,2,3,4 in samp
        samp.tab<-data.frame(table(samp))

        sampNEG<-ifelse(any(unique(samp)=="neg"),samp.tab$Freq[which(samp.tab$samp=="neg")],0)
        sampNEU<-ifelse(any(unique(samp)=="neu"),samp.tab$Freq[which(samp.tab$samp=="neu")],0)
        sampPOS<-ifelse(any(unique(samp)=="pos"),samp.tab$Freq[which(samp.tab$samp=="pos")],0)

        #calculate weighted average
        score <-((sampNEG*-1) + (sampNEU*0) + (sampPOS*1)) / (nrow(scores.only) * 4)

        #ranking
        rank<-ifelse(score<= -0.33,"negative",
                     ifelse(score>=0.33,"positive","neutral"))

        scoresB <- c(scoresB, rank)
      } #end loop
    } #end if

    return(scoresB)

  } else { #if bootstrap if FALSE, just get weighted mean
    #calculate number of tallies into Neg(itive), Neu(trals), Pos(itive)
    num_negs<-sum(scores.only$Negative, na.rm = T)
    num_neu<-sum(scores.only$Neutral, na.rm = T)
    num_pos<-sum(scores.only$Positive, na.rm = T)

    #calculate weighted average
    score <- ((num_negs*-1) + (num_neu*0) + (num_pos*1)) / (nrow(scores.only) * 4)

    #ranking
    rank<-ifelse(score<= -0.33,"negative",
                 ifelse(score>=0.33,"positive","neutral"))

    return(rank)
  }

}
