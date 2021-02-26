# read trait information from .docx
library(qdapTools)
raw_traits <- qdapTools::read_docx(
  "C:/Users/bruno/Documents/RESEARCH/PhD/EvolutionProject/taxon_traits.docx"
  )

# fill a data frame with the 360 morphological traits
df <- data.frame(matrix(NA, ncol = 361))

# first one is all NA anyway, but we need the name
df[1, 1] <- sub("\\ .*", "", raw_traits[[1]])
for (i in 2:length(raw_traits)) {
  # create a placeholder for this row
  df_sub <- data.frame(matrix(0, ncol = 361))
  
  # get the taxon name
  df_sub[1, 1] <- sub("\\ .*", "", raw_traits[[i]])
  
  # get the character data
  chars <- strsplit(sub(".*\\ ", "", raw_traits[[i]]), "")[[1]]
  
  # very confusing while loop necessary since 
  # a lot of the characters are uncertain (e.g. [12])
  
  # count (index of the df we're in)
  count <- 2
  
  # candidates for character (inside [])
  cands <- c()
  
  # whether we're inside a []
  running <- FALSE
  
  # loop count
  w <- 1
  
  # 362 since we have 360 characters + column 1
  while (count < 362) {
    # if the char is [, need to get in it
    if (chars[w] == "[") {
      running <- TRUE
      
      # increase loop counter
      w <- w + 1
      
      # continue to the next while iteration
      next
    }
    
    # if we are inside a []
    if (running) {
      # if we already checked all of the numbers inside
      if (chars[w] == "]") {
        # sample one out of them to be the trait value
        df_sub[1, count] <- sample(cands, 1)
        
        # reset all the necessary variables
        cands <- c()
        running <- FALSE
        
        # increase counters
        count <- count + 1
        w <- w + 1
        
        # continue to the next while iteration
        next
      } else {
        # if we're in the interior of a [], add the value
        # to the cands vector, including NA if it's a ?
        cands <- c(cands, ifelse(chars[w] == "?", NA, chars[w]))
        
        # increase loop counter
        w <- w + 1
        
        # continue to next while iteration
        next
      }
    } else {
      # if we are not inside [], just add value
      # to the matrix (remembering ? = NA)
      df_sub[1, count] <- ifelse(chars[w] == "?", NA, chars[w])
      
      # increase counters
      count <- count + 1
      w <- w + 1
    }
  }
  
  # add row to our data frame
  df <- rbind(df, df_sub)
}

# informative column names
colnames(df) <- c("species", paste0("trait_", 1:360))

# get the taxon list with time + geographic data (from Table 1)
taxon <- read.csv(
  "C:/Users/bruno/Documents/RESEARCH/PhD/EvolutionProject/taxon_list.txt",
                  sep = "\t", header = FALSE)

# select only the taxon list
taxon_list <- as.character(taxon[, 1])

# make genus species separation a _ instead of space
taxon_list <- unlist(lapply(taxon_list, function(s) sub("\\ ", "_", s)))

# create new data frame to hold only the caninae species
final_df <- data.frame()
for (i in 1:nrow(df)) {
  # if the species name is in the list, put it in the new df
  if (df[i, 1] %in% taxon_list) {
    final_df <- rbind(final_df, df[i, ])
  }
}

# expect 78 rows, got
print(nrow(final_df))

# so we know some are missing, but which?
for (i in 1:length(taxon_list)) {
  if (!(taxon_list[i] %in% final_df[, 1])) {
    print(taxon_list[i])
  }
}

# looking over Appendix S2 and the paper, seems like there were 
# some typos and/or name disagreements, so we correct that

# add the remaining speces
final_df <- rbind(final_df, df[df$species %in% 
                                 c("Lupullella_adusta", 
                                   "Lupullella_mesomelas", 
                                   "Xenocyon_lycanoides"), ])

# make their name equal to the paper (to agree with tree)
final_df[final_df$species %in% c("Lupullella_adusta", 
                                 "Lupullella_mesomelas", 
                                 "Xenocyon_lycanoides"), 1] <-
  c("Lupulella_adusta", "Lupulella_mesomelas", "Xenocyon_lycaonoides")

# finally, make all trait values numeric
for (i in 2:ncol(final_df)) {
  final_df[, i] <- as.numeric(final_df[, i])
}

# for how many traits do we have complete information?
info <- unlist(lapply(2:ncol(final_df), function(c) sum(!is.na(final_df[, c]))))
which(info == nrow(final_df))

# 193 (194 in df) is body size, so we have a continuous trait

# categorical traits for which we have every taxa are 
# really low variance, so we find a vector of variances
vars <- unlist(lapply(2:ncol(final_df), function(x) var(final_df[!is.na(final_df[, x]), x])))

# and find traits such that we have more than 0.1 variance
# and more than 95% information on taxa
which((vars > 0.1) & (info > 0.95*nrow(final_df)))

# 146 (147 in df) seems to work well