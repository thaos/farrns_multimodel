library(magrittr)

df <- read.table("GHCND/UK000056225.dly", skip = 291, stringsAsFactors = FALSE)
df <- subset(df, grepl(pattern = "(09|10|11)PRCP", x = df[, 1])) 
ichar <- sapply(df, is.character) 

df_pr <- df[, !ichar]
df_flag <- df[, ichar][, -1]
df_year <- substr(df[, ichar][, 1], 12, 15) %>% as.integer
df_month <- substr(df[, ichar][, 1], 16, 17) %>% as.integer
df_pr[df_pr < 0] <- NA

df_pr <- apply(df_pr, 1, max, na.rm = TRUE)
df_pr <- tapply(df_pr, df_year, max)

df_pr <- cbind(year = unique(df_year), pr = df_pr / 10)
# in mm
saveRDS(df_pr, file = "pr_ghcdn.rds")
