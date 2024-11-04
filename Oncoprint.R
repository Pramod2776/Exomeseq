setwd("~/Path/to/Workdir/")
library(xlsx)
library("tidyverse")

## select the features from annotation files
lapply(1:length(list.files()), function(xl){

 
  annotation = read.delim(list.files()[xl][1], header = T)
  colnames(annotation) =  c("CHROM", "POS", "REF", "ALT", "AF", "DP", "ALT_Allele", "EFFECT",
                          "AA_Change",  "IMPACT",  "NUC_Change",
                          "BIOTYPE",  "PASS", "GENE_NAME", "GENE_ID" )

annotation[annotation==""]<-NA
annotation_select = annotation %>%
  dplyr::filter(!is.na(AA_Change))

annotation_select_pass = annotation_select %>%
  dplyr::filter(PASS == "true")


annotation_select_pass_effect = annotation_select_pass %>%
  dplyr::filter(grepl('missense|frameshift', EFFECT)) %>%
  dplyr::filter(DP >= 100) %>%
  dplyr::filter(BIOTYPE == "true")


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
file.path = paste0("~/Path/to/Workdir/", file_name,".csv" )
file.path1 = paste0("~/Path/to/Workdir/", file_name,".xlsx" )

write.csv(annotation_gb_distinct_AA, file = file.path )
writexl::write_xlsx(annotation_gb_distinct_AA, path = file.path1)
  
  })

## Arrange Mutations Sample wise
list_files = list.files( pattern = ".xlsx")

mut_list_dff = lapply(1:length(list_files), function(l){
  name = gsub("_annotation_distinct_AA.xlsx", "", list_files[l])
  df = t(readxl::read_excel(list_files[l], sheet = 1) %>%
           as.data.frame() %>%
           dplyr::select(GENE_NAME, AA_Change_collapse)) %>%
    as.data.frame() 
  
  rownames(df)[2] = name
  
  colnames(df) = df %>%
    slice(1)
  
  df = df %>%
    slice(-1)

}
) %>%
  bind_rows()

##write.csv (mut_list_dff, "Mutations_samplewise_missense.csv")

dat.fr = mut_list_dff

dat.fr_mod = data.frame(lapply(dat.fr, function(x) { 
  if(is.factor(x)) {
    x <- as.character(x)
  }
  x[x!="0"]="1"
  x}))

dat.fr_mod[is.na(dat.fr_mod)] <- 0

rownames(dat.fr_mod) = rownames(mut_list_dff)
##write.csv (dat.fr_mod, "Mutations_samplewise_01_missense.csv")

Effect_list_dff = lapply(1:length(list_files), function(l){

  name = gsub("_annotation_distinct_AA.xlsx", "", list_files[l])
  df = t(readxl::read_excel(list_files[l], sheet = 1) %>%
           as.data.frame() %>%
           dplyr::select(GENE_NAME, EFFECT_collapse)) %>%
    as.data.frame() 
  
  rownames(df)[2] = name
  
  colnames(df) = df %>%
    slice(1)
  
  df = df %>%
    slice(-1)
  
}
) %>%
  bind_rows()


##write.csv (Effect_list_dff, "EFFECT_samplewise_missense.csv")

### oncoprint heatmap

###mat = read.table(system.file("extdata", package = "ComplexHeatmap", 
                             "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
                 header = TRUE, stringsAsFactors = FALSE, sep = "\t")

matt1 = Effect_list_dff
matt2 = matt1[,1:28]

matt2 = t(as.matrix(matt2))
matt2[is.na(matt2)] = ""

##write.csv(matt2, "matt2.csv",)

col1 = c("missense_variant" = "blue", "frameshift_variant" = "red")

alter_fun1 = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  missense_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["missense_variant"], col = NA))
  },
  # big red
  frameshift_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "pt"), h-unit(0.5, "pt"), 
              gp = gpar(fill = col["frameshift_variant"], col = NA))
  }
   # ,
  # # small green
  # MUT = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
  #             gp = gpar(fill = col["missense_variant, frameshift_variant"], col = NA))
  # }
  
)


alter_fun1 = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  missense_variant = alter_graphic("rect", fill = col["missense_variant"]),
  frameshift_variant = alter_graphic("rect", fill = col["frameshift_variant"])
)

column_title = "OncoPrint"
heatmap_legend_param1 = list(title = "Alternations", at = c("missense_variant", "frameshift_variant"), 
                            labels = c("Missense", "Frame_shift"))
oncoPrint(as.matrix(matt2),
          alter_fun = alter_fun1, col = col1, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param1,
          alter_fun_is_vectorized = FALSE)


oncoPrint(matt2,
          alter_fun = alter_fun1, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          pct_side = "right", row_names_side = "left",
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)




























# amainoacid_change = annotation_select$AA_Change
# amainoacid_change = amainoacid_change %>%
#   as.data.frame()
# 
# 
# annotation_select_high = annotation_select %>%
#   dplyr::filter(IMPACT == "HIGH" | IMPACT == "MODERATE")


effect_stats = table(annotation_select_pass$EFFECT) %>%
  as.data.frame()

colnames(effect_stats) = c("effect", "frequency")

effect_stats = effect_stats %>%
  arrange(desc(frequency))

head(effect_stats)








biotype_stats = table(annotation_select_pass$BIOTYPE) %>%
  as.data.frame()
colnames(biotype_stats) = c("biotype", "frequency")

biotype_stats = biotype_stats %>%
  arrange(desc(frequency))

head(biotype_stats)

GENE_NAME_stats = table(annotation_select_pass$GENE_NAME)%>%
  as.data.frame()
colnames(GENE_NAME_stats) = c("GENE_NAME", "frequency")

head(GENE_NAME_stats)

GENE_NAME_stats = GENE_NAME_stats %>%
  arrange(desc(frequency))


Impact_stats = table(annotation_select_pass$IMPACT)%>%
  as.data.frame()

colnames(Impact_stats) = c("IMPACT", "frequency")

head(Impact_stats)






annotation_list = list(variants = annotation_select_pass_effect,
     biotype = biotype_stats,
     effect = effect_stats,
     genes = GENE_NAME_stats,
     Impact = Impact_stats)



theList = annotation_list

wb <- createWorkbook()
sheet <- createSheet(wb,"TS600_NS_PASS")

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

saveWorkbook(wb,file = "TS600_NS_PASS.xlsx")
