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

results <- df %>%
  rowwise() %>%
  mutate(
    # Laplace correction for zeros only
    a = ifelse(`p<0.05_Pos` == 0, 0.5, `p<0.05_Pos`),
    b = ifelse(`p>0.05_Pos` == 0, 0.5, `p>0.05_Pos`),
    c = ifelse(`p<0.05_Neg` == 0, 0.5, `p<0.05_Neg`),
    d = ifelse(`p>0.05_Neg` == 0, 0.5, `p>0.05_Neg`),
    
    # Chi-square test
    chi_p = tryCatch({
      tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      suppressWarnings(chisq.test(tbl)$p.value)
    }, error = function(e) NA_real_()),
    
    # Fisher's Exact Test p-value
    fisher_p = tryCatch({
      tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      fisher.test(tbl)$p.value
    }, error = function(e) NA_real_()),
    

    odds_ratio = tryCatch({
      (a * d) / (b * c)
    }, error = function(e) NA_real_()),
    
    # Wald Confidence Interval
    se_log_or = tryCatch({
      sqrt(1/a + 1/b + 1/c + 1/d)
    }, error = function(e) NA_real_()),
    
    log_or = log(odds_ratio),
    
    ci_low = exp(log_or - 1.96 * se_log_or),
    ci_high = exp(log_or + 1.96 * se_log_or)
  ) %>%
  ungroup() %>%
  select(-se_log_or, -log_or)

# Save output
write.csv(results, "1.csv", row.names = FALSE)

# Preview
head(results)
