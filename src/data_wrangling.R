library(dplyr)

abund_table <- readxl::read_xlsx("../Phylum master absolute.xlsx", sheet = "Absolute.1")
abund_table <- abund_table %>% tibble::column_to_rownames(var = "Taxonomy")

abund_table_2 <- readxl::read_xlsx("../Phylum master absolute.xlsx", sheet = "Absolute.2")
abund_table_2 <- abund_table_2 %>% tibble::column_to_rownames(var = "Taxonomy")

abund_table_3 <- readxl::read_xlsx("../Phylum master absolute.xlsx", sheet = "Absolute.3")
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

write.csv(abund_table, file = "phylum-master-absolute.csv")