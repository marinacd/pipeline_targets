# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tibble", "tidyverse", "readxl", "dplyr", "classyfireR", "fobitools", "webchem"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Load the R scripts with your custom functions:
lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
tar_target(file, "~/Desktop/TFG/ProtegéPipeline/FooDB/data/associations.xlsx", format = "file"),
  tar_target(dataa, get_data(file)),
  tar_target(data, subset(dataa)),
  tar_target(prep, preprocessing(data, TRUE)),
  tar_target(ids, get_IDs(prep, FALSE, TRUE,TRUE,TRUE,FALSE,TRUE, TRUE)),
  
  tar_target(classes, get_classes(ids)),
  tar_target(foods, annotate_food(prep)),
  tar_target(fobi_mets, fobi_met()),
  
  tar_target(class_out, out_fobi(classes,fobi_mets)),
  tar_target(class_data, class_df(classes, class_out)), #works until here
  
  tar_target(ass, associations(ids, foods)),
  tar_target(ass_ann, ass_anno(ass, fobi_mets)),
  
  tar_target(fobi_rel, relations_fobi(ass_ann)),
  tar_target(new_rel, relations_fobi_new(ass_ann)),
  
  tar_target(df_class, new_classes(classes)),
  tar_target(new_met_class, class_new_met(new_rel, df_class)),

  tar_target(new, new_merge(new_rel, ids, new_met_class)),
  tar_target(fobiID, add_fobi(new)),

  tar_target(save_dat, save_data(prep)),
  tar_target(save_ann, save_ids(ids)),
  tar_target(save_rel, save_fobi_rel(fobi_rel)),
  tar_target(save_rel_new, save_new_rel(new_rel)),
  tar_target(save_class, save__new_classes(class_data)),
  tar_target(save_new_class,save_new_classes(new_met_class)),
  tar_target(save_new, save(new)),

  tar_target(sum, summary(prep,foods,ids, fobi_rel, new_rel, class_data))
  
)


#MAYBE I CAN CREATE A SCRIPT CALLED OPTIONS FOR STUFF SUCH AS THE SEPARATOR OR OTHER INPUTS?

