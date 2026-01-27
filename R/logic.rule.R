#' @title Calculate Total Sensitivity Score with the Logic Rule
#' @description Calculates total sensitivity from expert counts or bootstrapped values. Based on code from the Highly Migratory Species (HMS) CVA team. 

#' @param attribute.vec output from \code{attribute.score}
#' @param bootstrap binary TRUE/FALSE to turn on/off bootstrapping of sensitivity scores. Will impact outputs. 

#' @return returns the final sensitivity scores based on the logic rule. The output is similar to \code{attribute.score}, where if \code{bootstrap = TRUE}, the function returns a vector containing the final sensitivity score for each sample, and if \code{bootstrap = FALSE}, it returns a single final value.


logic.rule <- function(attribute.vec, bootstrap){
  if(bootstrap){ #if bootstrap = T, attribute.vec will be a MATRIX; apply the logic rule to each ROW 
    total.score <- vector(length = nrow(attribute.vec))
    for(x in 1:nrow(attribute.vec)){
      if(length(which(attribute.vec[x,]>=3.5))>=3){ # Very High 
        total.score[x] <- 4
      }else if(length(which(attribute.vec[x,]>=3))>=2){ # High
        total.score[x] <- 3
      }else if(length(which(attribute.vec[x,]>=2.5))>=2){ # Moderate
        total.score[x] <- 2
      }else{ # Low
        total.score[x] <- 1
      } #end if/elses
    } # end for 
  } else { #if bootstrap == F, attribute.vec will be a vector; apply logic rule to vector only
    total.score <- NULL
    if(length(which(attribute.vec>=3.5))>=3){ # Very High 
      total.score <- 4
    }else if(length(which(attribute.vec>=3))>=2){ # High
      total.score <- 3
    }else if(length(which(attribute.vec>=2.5))>=2){ # Moderate
      total.score <- 2
    }else{ # Low
      total.score <- 1
    } #end if/elses
  }
  
  return(total.score) #same thing to return either way :) 
}

