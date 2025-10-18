# Load libraries
library(dplyr)

# Read your files
tpms <- read.csv("TPMs.csv")
mapping <- read.csv("mouse_human.csv")

# Merge by the mouse gene symbol column
merged_df <- tpms %>%
  left_join(
    mapping %>%
      select(mouse_Symbol,
             human_geneid,
             human_EnsemblID,
             human_Symbol,
             human_type_of_gene,
             human_description),
    by = c("mgi_symbol" = "mouse_Symbol")
  )

# Save merged file
write.csv(merged_df, "TPMs_with_human_info.csv", row.names = FALSE)
