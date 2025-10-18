library(tidyverse)

df <- read.csv("subsystem_significance_summary_only_pval.csv", check.names = FALSE)

colnames(df) <- gsub("p \\?0.05_Pos", "p>0.05_Pos", colnames(df))
colnames(df) <- gsub("p\\?0.05_Neg", "p>0.05_Neg", colnames(df))

df <- df %>%
  mutate(
    `p<0.05_Pos` = as.numeric(`p<0.05_Pos`),
    `p>0.05_Pos` = as.numeric(`p>0.05_Pos`),
    `p<0.05_Neg` = as.numeric(`p<0.05_Neg`),
    `p>0.05_Neg` = as.numeric(`p>0.05_Neg`)
  ) %>%
  replace_na(list(
    `p<0.05_Pos` = 0,
    `p>0.05_Pos` = 0,
    `p<0.05_Neg` = 0,
    `p>0.05_Neg` = 0
  ))

chi_results <- df %>%
  rowwise() %>%
  mutate(
    a = ifelse(`p<0.05_Pos` == 0, 0.5, `p<0.05_Pos`),
    b = ifelse(`p>0.05_Pos` == 0, 0.5, `p>0.05_Pos`),
    c = ifelse(`p<0.05_Neg` == 0, 0.5, `p<0.05_Neg`),
    d = ifelse(`p>0.05_Neg` == 0, 0.5, `p>0.05_Neg`),
    
    chi_p = tryCatch({
      tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      suppressWarnings(chisq.test(tbl)$p.value)
    }, error = function(e) NA),
   
    fisher_p = tryCatch({
      tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      fisher.test(tbl)$p.value
    }, error = function(e) NA)
  ) %>%
  ungroup()

write.csv(chi_results, "Statistical test output p value only.csv", row.names = FALSE)

head(chi_results)