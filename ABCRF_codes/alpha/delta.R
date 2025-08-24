args <- commandArgs(trailingOnly = TRUE)
sim_number <- args[1]


output_files <- list.files(pattern = paste0("EMIBD_output_", sim_number, "_\\d+\\.txt"))
file_nums <- as.numeric(sub(paste0("EMIBD_output_", sim_number, "_(\\d+)\\.txt"), "\\1", output_files))
output_files <- output_files[order(file_nums)]

delta_list <- vector("list", 9)  # Delta1 to Delta9
for (i in 1:9) delta_list[[i]] <- list()
r_list <- list()


F_Mean_matrix <- matrix(0, nrow = 0, ncol = 15)  # Assuming 15 individuals
F_SD_matrix   <- matrix(0, nrow = 0, ncol = 15)

for (file in output_files) {
  lines <- readLines(file)

  start <- grep("IBD coefficient estimates", lines)
  if (length(start) == 0) next

  data_start <- start + 3
  data_lines <- lines[data_start:(data_start + 224)]  # 225 dyads
  data_con <- textConnection(data_lines)
  df <- read.table(data_con, header = FALSE)
  close(data_con)

  colnames(df) <- c("Dyad", "Indiv1", "Indiv2",
                    paste0("IIS", 1:9),
                    paste0("Delta", 1:9),
                    "r")

  for (i in 1:9) {
    delta_list[[i]][[length(delta_list[[i]]) + 1]] <- df[[paste0("Delta", i)]]
  }

  r_list[[length(r_list) + 1]] <- df$r

  start2 <- grep("Indiv genotypes at polymorphic loci & inbreeding", lines)
  if (length(start2) == 0) next

  data_start2 <- start2 + 3
  data_lines2 <- lines[data_start2:(data_start2 + 14)]  # 15 individuals
  data_con2 <- textConnection(data_lines2)
  df2 <- read.table(data_con2, header = FALSE)
  close(data_con2)

  F_Mean_matrix <- rbind(F_Mean_matrix, t(df2[[6]]))  # F_Mean column
  F_SD_matrix   <- rbind(F_SD_matrix, t(df2[[7]]))    # F_SD column
}

delta_matrices <- lapply(delta_list, function(x) do.call(rbind, x))  # 200 × 225
r_matrix <- do.call(rbind, r_list)  # 200 × 225

for (i in 1:9) {
  write.table(delta_matrices[[i]], file = paste0("Delta", i, "_", sim_number, ".txt"),
              row.names = FALSE, col.names = FALSE)
}

write.table(r_matrix, file = paste0("r_", sim_number, ".txt"),
            row.names = FALSE, col.names = FALSE)

write.table(F_Mean_matrix, file = paste0("F_Mean_", sim_number, ".txt"),
            row.names = FALSE, col.names = FALSE)

write.table(F_SD_matrix, file = paste0("F_SD_", sim_number, ".txt"),
            row.names = FALSE, col.names = FALSE)