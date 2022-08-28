library(lubridate)

# Data manager

#
partition_eeg_data <- function(df, start, end)
{
   s_ind <- which(df$Time == start)
   e_ind <- which(df$Time == end)
   
   return (df[s_ind:e_ind, ])
}

lower_resolution <- function(df, steps)
{
  I <- c(1) # Include one since following computation skips it.
  for (n in 1:(nrow(df)/steps))
  {
    #print(paste(round((n * steps + 1) / nrow(df), 3), '%', sep = ''))
    I <- append(I, n * steps + 1)  
  }
  
  return (df[I, ])
}

remove_channels <- function(df, channels)
{
  for (c in channels)
  {
    print(c)
    df[c] <- NULL
  }
  return (df)
}


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#remove_repeated_anoms(subset(POINT_ANOMS, POINT_ANOMS$variate == 1))
