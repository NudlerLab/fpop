# parser of uniprot formatted fasta file with protein sequences into dataframe

# header is split in columns:
# 1. db
# 2. UniqueIdentifier
# 3. EntryName
# 4. ProteinName
# 5. OS=OrganismName
# 6. OX=OrganismIdentifier
# 7. GN=GeneName
# 8. PE=ProteinExistence
# 9. SV=SequenceVersion
# 10. sequence
# 11. length
# 12. mw

library(tidyverse)

calc_mw <- function(protein_sequence) {
  
  aa_weights <- c("G" = 57.02147,"A" = 71.03712,"S" = 87.03203,
                  "P" = 97.05277,"V" = 99.06842,"T" = 101.04768,
                  "C" = 103.00919,"U" = 168.96420,"I" = 113.08407,
                  "L" = 113.08407,"N" = 114.04293,"D" = 115.02695,
                  "Q" = 128.05858,"K" = 128.09497,"E" = 129.0426,
                  "M" = 131.04049,"H" = 137.05891,"F" = 147.06842,
                  "R" = 156.10112,"Y" = 163.06333,"W" = 186.07932)
  
  current_sequence <- str_split(string = protein_sequence, pattern = "")[[1]]
  current_mw <- map_dbl(current_sequence, function(x) aa_weights[x])
  current_mw[is.na(current_mw)] <- 110.0 # 100 Da for average amino acid
  
  return(sum(current_mw))
}

parse_uniprot_fasta <- function(fasta_file_name) {
  
  fasta <- read_file(file = fasta_file_name) # reading whole file in one string of chars
  
  fasta <- str_split(string = fasta, pattern = ">")[[1]] # splitting at the beginning of each header
  
  fasta <- as_tibble(fasta) # converting list of header_sequence strings into 1-col dataframe
  
  fasta <- fasta %>%
    separate(col = value, into = c("fasta_header", "protein_sequence"), sep = "\n", extra = "merge", fill = "right") %>%
    filter(!is.na(protein_sequence)) %>%
    mutate(protein_id = str_match(string = fasta_header, pattern = "\\|(.*?)\\|")[,2],
           protein = str_match(string = fasta_header, pattern = "\\|.*?\\|(.*?\\_.*?)\\s")[,2],
           gene = str_match(string = fasta_header, pattern = "GN=(.*?)\\s")[,2],
           description = str_match(string = fasta_header, pattern = "\\s(.*?)\\sOS")[,2]) %>%
    select(-fasta_header) %>%
    mutate(protein_sequence = str_replace_all(string = protein_sequence, pattern = "\n", replacement = "")) %>%
    mutate(len = str_length(string = protein_sequence),
           mw = map_dbl(protein_sequence, calc_mw))
 
  return(fasta)
}