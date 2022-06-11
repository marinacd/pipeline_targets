# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tibble", "tidyverse", "readxl", "classyfireR", "fobitools", "webchem"), # packages that your targets need to run
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
  tar_target(file, "~/Desktop/TFG/ProtegéPipeline/TEST_27may_Exposome/data/inputs/raw/metabolomic_associations.xlsx", format = "file"),
  tar_target(dataaaa, get_data(file)),
  tar_target(data, subset(dataaaa)),
  tar_target(prep, preprocessing(data, FALSE)),
  tar_target(ids, get_IDs(prep, TRUE, TRUE,TRUE,TRUE,TRUE,TRUE, TRUE)),
  tar_target(save_ann, save_ids(ids)),
  tar_target(classes, get_classes(ids)),
  tar_target(foods, annotate_food(prep)),
  tar_target(fobi_mets, fobi_met()),
  tar_target(class_out, out_fobi(classes,fobi_mets)),
  tar_target(class_data, class_df(classes, class_out)),
  # tar_target(sum, summary(prep,foods,ids)),
  tar_target(save_dat, save_data(prep)),
  tar_target(save_class, save__new_classes(class_data))
  
  )


#MAYBE I CAN CREATE A SCRIPT CALLED OPTIONS FOR STUFF SUCH AS THE SEPARATOR OR OTHER INPUTS?

