get_data <- function(file) {
  read_xlsx(file)
}

subset <- function(data){
  data[1:10,]
}

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
    split_data(data,"/")
  }
  if (!('ID' %in% colnames(data))){
    add_IDcol(data)}
  data %>% select(ID,Biomarker,Intake)%>%
    dplyr::rename(food = Intake) %>%
    dplyr::rename(biomarker = Biomarker) %>%
    mutate(biomarker = tolower(biomarker)) #%>%
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
get_IDs <- function(data, InChIKey, HDMB, ChemSpider, ChEBI, InChI, KEGG, PubChem){
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
  openxlsx::write.xlsx(data, "annotations.xlsx")
}


fobi_met <- function(){
  fobitools::parse_fobi(terms = "FOBI:01501", get = "des") %>%    #"FOBI:01501" is the root of metabolites
    mutate(biomarker = tolower(name))   
}

#classification subfunctions

classifiable <- function(data){
  data  %>%
    select(biomarker , `InChIKey`) %>%    
    filter(!is.na(InChIKey)) %>%
    filter(!duplicated(.))
}

not_classifiable <- function(data){
  data  %>%
    select(biomarker , `InChIKey`) %>%
    filter(is.na(InChIKey)) %>%
    filter(!duplicated(.))
}

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
  data %>%
    select(ID, food) %>%
    fobitools::annotate_foods(similarity = 0.85, reference = fobitools::foods)

}



associations_fobi_anno <- associations_fobi %>%
  left_join(metabolites_fobi_2, by = "biomarker") %>%
  filter(!duplicated(.)) %>%
  dplyr::rename(original_name = name.y)





















# 
# 
# summary <- function(data, foods,ids){
# 
#   n1 <- length(unique(data()$biomarker))
#   n2 <- length(unique(data$food))
#   n3 <- nrow(data)
#   rel <- paste0("Found ", n3 , " relationships in input data, associating ", n1, " biomarkers with ", n2, " food items.")
# 
#   info <- list()
#   info$relationships <- rel
#   info$foods <-  paste0( "Found ", nrow(foods$annotated), " foods that could be annotated and ", nrow(foods$unannotated), " that couldn't.")
#   
#   info$InChIKey <- paste0("Can't retrieve ", sum(is.na(ids$InChIKey)), " InChIKeys, unable to classify those metabolites.")
#   
#   
#   print(info)
#   
#   
#   
# }
# 

save_data <- function(data) {
  openxlsx::write.xlsx(data, "composition_data_processed.xlsx")
}
save__new_classes <- function(data){
  openxlsx::write.xlsx(data, "fobi_new_classes.xlsx")
}




