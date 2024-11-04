setwd("~/PAth/to/WorkDir/")
library(xlsx)
library("tidyverse")

## select features from annotation table

lapply(1:length(list.files()), function(xl){
  annotation = read.delim(list.files()[xl][1], header = T)
  colnames(annotation) =  c("CHROM", "POS", "REF", "ALT", "EFFECT",
                          "AA_Change", "NUC_Change", "IMPACT",
                          "BIOTYPE", "QUAL", "PASS", "GENE_NAME", "GENE_ID" )

  annotation[annotation==""]<-NA

  annotation_select = annotation %>%
  dplyr::filter(!is.na(AA_Change))


  annotation_select_pass = annotation_select %>%
  dplyr::filter(PASS == "true")

 

  annotation_select_pass_effect = annotation_select_pass %>%
    dplyr::filter(EFFECT %in% c("chromosome_number_variation", "coding_sequence_variant", "conservative_inframe_deletion",
                                "conservative_inframe_insertion", "disruptive_inframe_deletion", "disruptive_inframe_insertion",
                                "exon_loss", "exon_loss_variant", "exon_variant", "frameshift_variant", "gene_variant",
                                "initiator_codon_variant", "missense_variant", "missense_variant&splice_region_variant",
                                "rare_amino_acid_variant", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
                                "start_lost", "start_retained", "stop_gained", "stop_lost", "stop_retained_variant", "synonymous_variant",
                                "transcript_ablation", "transcript_amplification", "transcript_variant", "frameshift_variant",
                                "rare_amino_acid_variant", "frameshift_variant&stop_gained", "frameshift_variant&start_lost",
                                "frameshift_variant&splice_region_variant", "frameshift_variant&stop_lost",
                                "stop_gained&disruptive_inframe_deletion", "frameshift_variant&stop_gained&splice_region_variant",
                                "splice_acceptor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant",
                                "start_lost&conservative_inframe_deletion", "stop_gained&conservative_inframe_insertion",
                                "missense_variant&splice_region_variant", "splice_region_variant&synonymous_variant",
                                "stop_lost&splice_region_variant", "start_lost&conservative_inframe_deletion", "frameshift_variant&stop_lost",
                                "frameshift_variant&splice_region_variant", "splice_acceptor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant",
                                "stop_gained&splice_region_variant", "splice_region_variant&synonymous_variant",
                                "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&stop_gained",
                                "stop_gained&splice_region_variant", "stop_lost&splice_region_variant",
                                "conservative_inframe_deletion&splice_region_variant", "start_lost&conservative_inframe_deletion",
                                "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant",
                                "splice_donor_variant&missense_variant&splice_region_variant&intron_variant",
                                "splice_donor_variant&splice_region_variant&synonymous_variant&intron_variant",
                                "stop_gained&conservative_inframe_insertion", "stop_gained&splice_region_variant",
                                "splice_donor_variant&splice_region_variant&synonymous_variant&intron_variant","stop_gained&disruptive_inframe_insertion",
                                "splice_acceptor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant",
                                "splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
                                "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant",
                                "splice_acceptor_variant&missense_variant&splice_region_variant&intron_variant",
                                "splice_donor_variant&splice_region_variant&synonymous_variant&intron_variant",
                                "start_lost&conservative_inframe_deletion", "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant",
                                "splice_acceptor_variant&missense_variant&splice_region_variant&intron_variant", "stop_gained&splice_region_variant",
                                "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant", "frameshift_variant&splice_region_variant"))


  annotation_select_pass_effect = annotation_select_pass_effect %>%
  add_count(GENE_NAME) %>%
  arrange(desc(n))

  annotation = annotation_select_pass_effect %>%
  group_by(GENE_NAME) %>%
  mutate(AA_Change_collapse = paste0(AA_Change, collapse = ";")) %>%
  mutate(EFFECT_collapse = paste0(EFFECT, collapse = ";")) %>%
  mutate(NUC_Change_collapse = paste0(NUC_Change, collapse = ";"))

  annotation_distinct_AA = annotation %>%
  distinct(AA_Change_collapse, .keep_all = TRUE)

  file_name = paste(list.files()[xl][1], "annotation_distinct_AA", sep = "_")
  
  file.path = paste0("~/PAth/to/WorkDir//", file_name,".csv" )
  file.path1 = paste0("~/PAth/to/WorkDir/", file_name,".xlsx" )

  write.csv(annotation_distinct_AA, file = file.path )
  writexl::write_xlsx(annotation_distinct_AA, path = file.path1)
  
  
  })


list_files = list.files( pattern = ".xlsx")

mut_list_dff = lapply(1:length(list_files), function(l){
  print(l)
  name = gsub("annotation_distinct_AA.xlsx", "", list_files[l])
  df = t(readxl::read_excel(list_files[l], sheet = 1) %>%
           as.data.frame() %>%
           dplyr::select(GENE_NAME, AA_Change_collapse)) %>%
    as.data.frame() 
  
  rownames(df)[2] = name
  
  colnames(df) = df %>%
    slice(1)
  
  df = df %>%
    slice(-1)
  
  print(head(df))

}
) %>%
  bind_rows()
class(mut_list_dff)
##write.csv (mut_list_dff, "Mutations_samplewise.csv")

dat.fr = mut_list_dff

dat.fr_mod = data.frame(lapply(dat.fr, function(x) { 
  if(is.factor(x)) {
    x <- as.character(x)
  }
  x[x!="0"]="1"
  x}))

dat.fr_mod[is.na(dat.fr_mod)] <- 0

rownames(dat.fr_mod) = rownames(mut_list_dff)
##write.csv (dat.fr_mod, "Mutations_samplewise_01.csv")


Effect_list_dff = lapply(1:length(list_files), function(l){
  print(l)
  name = gsub("annotation_distinct_AA.xlsx", "", list_files[l])
  df = t(readxl::read_excel(list_files[l], sheet = 1) %>%
           as.data.frame() %>%
           dplyr::select(GENE_NAME, EFFECT_collapse)) %>%
    as.data.frame() 
  
  rownames(df)[2] = name
  
  colnames(df) = df %>%
    slice(1)
  
  df = df %>%
    slice(-1)
  
  print(head(df))
  
}
) %>%
  bind_rows()


##write.csv (Effect_list_dff, "EFFECT_samplewise.csv")


## Calculate Stats
effect_stats = table(annotation_select_pass$EFFECT) %>%
  as.data.frame()
colnames(effect_stats) = c("effect", "frequency")
effect_stats = effect_stats %>%
  arrange(desc(frequency))

biotype_stats = table(annotation_select_pass$BIOTYPE) %>%
  as.data.frame()
colnames(biotype_stats) = c("biotype", "frequency")

biotype_stats = biotype_stats %>%
  arrange(desc(frequency))

GENE_NAME_stats = table(annotation_select_pass$GENE_NAME)%>%
  as.data.frame()
colnames(GENE_NAME_stats) = c("GENE_NAME", "frequency")

GENE_NAME_stats = GENE_NAME_stats %>%
  arrange(desc(frequency))

Impact_stats = table(annotation_select_pass$IMPACT)%>%
  as.data.frame()
colnames(Impact_stats) = c("IMPACT", "frequency")

annotation_list = list(variants = annotation_select_pass_effect,
     biotype = biotype_stats,
     effect = effect_stats,
     genes = GENE_NAME_stats,
     Impact = Impact_stats)

theList = annotation_list

wb <- createWorkbook()
sheet <- createSheet(wb,"T_N_PASS")

currRow <- 1
for(i in 1:length(theList)){
  
  cs <- CellStyle(wb) + Font(wb, isBold=TRUE) + Border(position=c("BOTTOM", "LEFT", "TOP", "RIGHT"))
  
  addDataFrame(theList[[i]],
               sheet=sheet,
               startRow=currRow,
               row.names=FALSE,
               colnamesStyle=cs)
  
  currRow <- currRow + nrow(theList[[i]]) + 2 
}

saveWorkbook(wb,file = "T_N_PASS.xlsx")
