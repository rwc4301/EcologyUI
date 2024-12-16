library(dplyr)

tax_table_1 <- readxl::read_xlsx("Genus master absolute.xlsx", sheet = "TAX.1")
tax_table_2 <- readxl::read_xlsx("Genus master absolute.xlsx", sheet = "TAX.2")
tax_table_3 <- readxl::read_xlsx("Genus master absolute.xlsx", sheet = "TAX.3")

tax_table_1 <- t(tax_table_1)
tax_table_2 <- t(tax_table_2)
tax_table_3 <- t(tax_table_3)

tax_table <- unique(rbind(tax_table_1, tax_table_2, tax_table_3))

abund_table <- readxl::read_xlsx("Genus master absolute.xlsx", sheet = "OTU.1")
abund_table <- abund_table %>% tibble::column_to_rownames(var = "Taxonomy")

abund_table_2 <- readxl::read_xlsx("Genus master absolute.xlsx", sheet = "OTU.2")
abund_table_2 <- abund_table_2 %>% tibble::column_to_rownames(var = "Taxonomy")

abund_table_3 <- readxl::read_xlsx("Genus master absolute.xlsx", sheet = "OTU.3")
abund_table_3 <- abund_table_3 %>% tibble::column_to_rownames(var = "Taxonomy")

abund_table <- as.data.frame(t(abund_table))
abund_table_2 <- as.data.frame(t(abund_table_2))
abund_table_3 <- as.data.frame(t(abund_table_3))

abund_table$row.names <- rownames(abund_table)
abund_table_2$row.names <- rownames(abund_table_2)
abund_table_3$row.names <- rownames(abund_table_3)

abund_table <- full_join(abund_table, abund_table_2, by = "row.names")
abund_table <- full_join(abund_table, abund_table_3, by = "row.names")

rownames(abund_table) <- abund_table$row.names

abund_table <- abund_table[, !names(abund_table) %in% "row.names"]
abund_table <- as.data.frame(t(abund_table))

abund_table[is.na(abund_table)] <- 0

write.csv(abund_table, file = "genus-master-absolute.csv")
write.csv(tax_table, file = "genus-tax-table.csv")