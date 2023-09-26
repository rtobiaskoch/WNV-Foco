


rm(list = ls())
pacman::p_load(tidyverse, Biostrings, readr)
source("https://raw.githubusercontent.com/rtobiaskoch/TK_useful_functions/main/fasta_filter.R")


mdata = read_rds("data_output/mdata_c_cmplt.rds")

fasta = Biostrings::readDNAStringSet("data_output/fasta_ns.fasta")

larimer = mdata %>%
  filter(grepl("Larimer", division) & grepl("Larimer", trap))

csu = mdata %>%
  filter(date  %in% larimer$date)

t = rbind(csu, larimer)         

f = filter_fasta(fasta, t$strain)

Biostrings::writeXStringSet(f, "data_mid/larimer_fasta_test.fasta")
