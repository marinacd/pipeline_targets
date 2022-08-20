get_data <- function(file) {
  read_xlsx(file)
}

subset <- function(data){
  data[503:703,]
}

# input_exp <-"~/Desktop/TFG/ProtegéPipeline/TEST_27may_Exposome/data/inputs/raw/metabolomic_associations.xlsx"
# input_phyto <-"~/Desktop/TFG/ProtegéPipeline/PhytoHub/phytohub_associations.xlsx"

#preprocessing subfunctions
split_data <- function(data,sep) {
  data %>%
    mutate(Intake = strsplit(Intake, sep)) %>%
    unnest(Intake) %>%
    filter(!is.na(Intake))
}

add_IDcol <- function(data) {
  data %>%
    mutate(ID = 1:length(data$Biomarker))
}


##PREPROCESSING##
preprocessing <- function(data, split) {
  
  if(split){
    data <- split_data(data,"/")
  }
  
  # if (!('ID' %in% colnames(data))){
    data <- add_IDcol(data)
  
  ids <- c("InChIKey", "ChEBI", "ChemSpider", "InChI", "KEGG", "PubChem", "HDMB" )

  data %>% 
    dplyr::rename(food = Intake) %>%
    dplyr::rename(biomarker = Biomarker) %>%
    mutate(biomarker = tolower(biomarker)) %>%
    select(ID, biomarker, food, which(colnames(data) %in% ids))
  # unique(by = c(biomarker, food)) #no repeated relationships plssssssss
}

#ids subfunctions
get_inchis <- function(data){
  cts_convert(c(data$biomarker), from = "Chemical Name", to = "InChIKey", match = "first")
}
get_hdmb <- function(data){
  cts_convert(c(data$biomarker), from = "Chemical Name", to = "human metabolome database", match = "first")
}
get_chemspider <- function(data){
  cts_convert(c(data$biomarker), from = "Chemical Name", to = "ChemSpider", match = "first")
}
get_chebi <- function(data){
  cts_convert(c(data$biomarker), from = "Chemical Name", to = "ChEBI", match = "first")
}
get_inchicode <- function(data){
  cts_convert(c(data$biomarker), from = "Chemical Name", to = "inchi code", match = "first")
}
get_kegg <- function(data){
  cts_convert(c(data$biomarker), from = "Chemical Name", to = "KEGG", match = "first")
}
get_pubchem <- function(data){
  cts_convert(c(data$biomarker), from = "Chemical Name", to = "PubChem CID", match = "first")
}

##IDs##
get_IDs <- function(data, InChIKey, HDMB, ChemSpider, ChEBI, InChI , KEGG, PubChem){
  if(InChIKey){
    data <- data %>%
      mutate(InChIKey = unlist(get_inchis(data)))    #maybe i could add this to preprocessing since it is needed for classification
  }
  if(HDMB){
    data <- data %>%
      mutate(HDMB = unlist(get_hdmb(data)))
  }
  if(ChemSpider){
    data <- data %>%
      mutate(ChemSpider = unlist(get_chemspider(data)))
  }
  if(ChEBI){
    data <- data %>%
      mutate(ChEBI = unlist(get_chebi(data)))
  }
  if(InChI){
    data <- data %>%
      mutate(InChI = unlist(get_inchicode(data)))
  }
  if(KEGG){
    data <- data %>%
      mutate(KEGG = unlist(get_kegg(data)))
  }
  if(PubChem){
    data <- data %>%
      mutate(`PubChem CID` = unlist(get_pubchem(data)))
  }
}

save_ids <- function(data){
  openxlsx::write.xlsx(data, "annotationsFooDB5.xlsx", overwrite = TRUE)
}


fobi_met <- function(){
  fobitools::parse_fobi(terms = "FOBI:01501", get = "des")   #"FOBI:01501" is the root of metabolites
  
}

#classification subfunctions

classifiable <- function(data){
  data <- data[which(unlist(lapply(data$InChIKey, is.inchikey))),] #can only handle one inchikey string
  data  %>%
    select(biomarker , `InChIKey`) %>%    
    filter(!is.na(InChIKey)) %>%
    filter(!duplicated(.))
}

# not_classifiable <- function(data){
#   data  %>%
#     select(biomarker , `InChIKey`) %>%
#     filter(is.na(InChIKey)) %>%
#     filter(!duplicated(.))
# }

classify <- function(data){
  df <- purrr::map(data$InChIKey, get_classification)
  
  data <- purrr::map(df, classification) %>%             
    set_names(data$InChIKey) %>%                
    enframe() %>%                                           
    unnest(cols = c(value)) %>%                               
    filter(Classification != "Organic compounds")   
}

get_classes <- function(data){
  classify(classifiable(data))
  
}

out_fobi <- function(data, met_fobi){
  #which of the CLASSES found are not in fobi (not metabolites) ADD TO SÜMMARYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
  
  unique(data$Classification[!data$Classification %in% met_fobi$name])  
}

###CLASSIFICATION##

class_df <- function(data, out_fobi){
  ls_class <- list() #create a list
  
  for (i in 1:nrow(data)) {                                    #for every row in the classification dataframe
    if(data$Classification[i] %in% out_fobi){                  #if the classification assigned is not in fobi
      ls_class[[i]] <- data.frame(parent = data$Classification[i-1], child = data$Classification[i])
    }
  }
  ls_class <- bind_rows(ls_class) %>%   #binds all df in the classification list by rows
    # the result will be a df with only the classifications we dont have in fobi bc the other ones are null
    #so nothing to bind
    filter(!duplicated(.))
}


# relationships subfunctions
annotate_food <- function(data){
  # if(duplicated(data$ID)){
  #   data <- data %>%
  #     mutate(ID = row_number())
  # }
  data %>%
    select(ID, food) %>%
    fobitools::annotate_foods(similarity = 0.85, reference = fobitools::foods)
  
}



associations <- function(data,foods){
  foods_ann <- foods$annotated %>%
    dplyr::rename(food = FOOD_NAME)
  
  ids <- data %>%
    select(biomarker,InChIKey)
  
  data <- data %>%
    right_join(foods_ann, by = "food") %>%
    select(FOBI_ID, FOBI_NAME, food, biomarker)
  
  data <- data %>%
    left_join(ids, by = "biomarker") %>%
    mutate(biomarker = tolower(biomarker))
  
  return(data)
}

ass_anno <- function(ass,mets){
  mets <- mets %>%
    select(1:4) %>%
    mutate(biomarker = tolower(name))   
  
  ass <- ass %>%
    left_join(mets, by = "biomarker") %>%
    filter(!duplicated(.)) %>%
    dplyr::rename(original_name = name)
}


#### relations to include (metabolites already present in FOBI)

relations_fobi <- function(ass_ann){
  ass_ann %>%
    filter(!is.na(original_name)) %>%
    filter(!duplicated(.))
}


save_fobi_rel <- function(fobi_rel){
  openxlsx::write.xlsx(fobi_rel, "relations_fobiFooDB5.xlsx", overwrite = TRUE)
}


#### relations to include (new FOBI metabolites)

relations_fobi_new <- function(ass_ann){
  ass_ann %>%
  filter(is.na(original_name)) %>%                  #retain rows that AREN'T in fobi
  filter(!is.na(InChIKey)) %>%                      #retain rows that have and inchikey
  filter(!duplicated(.)) %>%
  mutate(biomarker = gsub(":", "-", biomarker))
}

save_new_rel <- function(new_rel){
  openxlsx::write.xlsx(new_rel, "relations_fobi_newFooDB5.xlsx", overwrite = TRUE)
}



new_classes <- function(classes){
  classes %>%
    dplyr::rename(InChIKey = name) %>%
    group_by(InChIKey) %>%             #next operations will be done by Inchikey group
    dplyr::slice(n()) %>%              #index rows by location, location = size of group
    #aka will get the last row for each inchikey so the last classification
    filter(!duplicated(InChIKey))
} 

####Classification of new metabolites

class_new_met <- function(relations_fobi_new, df_class2){
  relations_fobi_new %>%
    left_join(df_class2, by = "InChIKey") %>%
    select(biomarker, Classification)   %>%
    mutate(biomarker = gsub(":", "-", biomarker))
}




summary <- function(data, foods, ann, fobi_rel, new_rel, class_data){
  
  info <- list()
  
  info$Input_Data <- paste0("Found ", nrow(data) , " relationships in input data, associating ", 
                            length(unique(data$biomarker)), " biomarkers with ", length(unique(data$food)), " food items.")
  
  info$Food_Annotation <-  paste0( "Able to annotate: ", nrow(foods$annotated), " foods. ", 
                                   "Unable to annotate: ", nrow(foods$unannotated), " foods.")
  
  info$Missing_classification <- paste0("Can't retrieve ", sum(is.na(ann$InChIKey)), " InChIKeys, unable to classify those metabolites.")
  
  a <- sum(!is.na(ann$InChIKey), na.rm = TRUE)
  b <- sum(!is.na(ann$HDMB), na.rm = TRUE)
  c <- sum(!is.na(ann$ChemSpider), na.rm = TRUE)
  d <- sum(!is.na(ann$ChEBI), na.rm = TRUE)
  e <- sum(!is.na(ann$InChI), na.rm = TRUE)
  f <- sum(!is.na(ann$KEGG), na.rm = TRUE)
  g <- sum(!is.na(ann$`PubChem CID`), na.rm = TRUE)
  
  info$Total_IDS <- tibble( "InChIKey" = a, "HDMB" = b, "ChemSpider" = c, "ChEBI" = d, InChI = e, "KEGG" = f , "PubChem CID" = g)
  
  info$Relationships_Found <- paste0("Found ",  nrow(fobi_rel), " new relationships from metabolites already present in fobi." )
  info$New_Relationships_Found <- paste0("Found " , nrow(new_rel), " new relationships from metabolites already present in fobi.")
  
  info$classes <- paste0("Found ", nrow(class_data), " new classes.")
  
  print(info)
  
  
}

add_fobi <- function(data){
  seq <- seq(from = 50340, to = 50340 + sum(is.na(data$FOBI_ID)) -1 )
  ids <- paste0("FOBI:0", seq)
  data$FOBI_ID[is.na(data$FOBI_ID)] <- ids
  return(data)
}

new_merge <- function(rel,ann,class){
  rel <- rel %>% select(1:5)
  ann <- ann %>% select(-ID) %>% filter(!duplicated(biomarker))
  b <- left_join(ann,rel)
  c <- left_join(b, class ) %>% distinct(.keep_all = TRUE)
  add_fobi(c)
  return(c)
}



save <- function(data) {
  openxlsx::write.xlsx(data, "new_mets_and_relsFooDB5.xlsx")
}





save_data <- function(data) {
  openxlsx::write.xlsx(data, "composition_data_processedFooDB5.xlsx")
}
save__new_classes <- function(data){
  openxlsx::write.xlsx(data, "fobi_new_classesFooDB5.xlsx")
}
save_new_classes <- function(data){
  openxlsx::write.xlsx(data, "fobi_new_met_classesFooDB5.xlsx")
}




