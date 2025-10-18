library(dplyr)
library(tidyr)
library(pheatmap)
library(stringr)

# === Load files ===
tpms <- read.csv("TPMs_with_human_info.csv", check.names = FALSE)
top75 <- read.csv("top75_subsystems.csv", check.names = FALSE)

# === Step 1: Extract Purine catabolism genes ===
purine_genes <- top75 %>%
  filter(str_trim(tolower(subsystem)) == "purine catabolism") %>%
  separate_rows(associated_genes, sep = ";\\s*") %>%
  pull(associated_genes) %>%
  unique()

# === Step 2: Filter TPMs for these genes ===
purine_tpms <- tpms %>%
  filter(human_Symbol %in% purine_genes)

# === Step 3: Detect expression columns ===
expr_cols <- setdiff(names(purine_tpms),
                     c("ensembl_gene_id", "mgi_symbol", "ensembl_gene_id_version",
                       "start_position", "end_position", "size",
                       "human_geneid", "human_EnsemblID",
                       "human_Symbol", "human_type_of_gene", "human_description"))

expr_mat <- purine_tpms[, expr_cols]
rownames(expr_mat) <- purine_tpms$human_Symbol

# === Step 4: Log2 transform ===
expr_mat_log <- log2(expr_mat + 1)

# === Step 5: Save heatmap with square cells and legend ===
png("7777777.png", width = 2000, height = 1200, res = 200)

pheatmap(expr_mat_log,
         cluster_rows = FALSE,     
         cluster_cols = TRUE,      
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 8,
         cellwidth = 15,
         cellheight = 15,
         legend = TRUE,   # ensures color legend is drawn
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Purine Catabolism")

dev.off()
