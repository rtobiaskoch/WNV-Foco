rm(list = ls())
library(tidyverse)

n1 = read.csv("../../WNV-Foco-nextstrain1/data/headers.csv", header = F)
n2 = read.csv("../data/metadata.csv") %>% select(-X)

df_p = n2 %>%
  mutate(year = year(date)) %>%
  group_by(year, state) %>%
  count()
  
ggplot(df_p, aes(x = year, y = n, fill = state)) +
  geom_col()

n1_not_n2 = anti_join(n1, n2, by = c("V1" = "accession"))

n2_not_n1 = anti_join(n2, n1, by = c("accession" = "V1"))

n2_in_n1 = semi_join(n2, n1, by = c("accession" = "V1"))

n2_in_n1.2 = semi_join(n2, n1, by = c("accession" = "V1", "date" = "V2"))
