library(curatedMetagenomicData)
library(dplyr)
library(tidyr)

# 1. All variable names
colnames(sampleMetadata)

# 2. Which variables are populated (at least one non-NA) in each study
vars_per_study <- sampleMetadata |>
  group_by(study_name) |>
  summarise(across(everything(), ~ any(!is.na(.x)))) |>
  pivot_longer(-study_name, names_to = "variable", values_to = "present") |>
  filter(present) |>
  select(-present)

# as a named list: study -> character vector of variables
vars_list <- split(vars_per_study$variable, vars_per_study$study_name)
vars_list[["AsnicarF_2017"]]

# 3. In how many studies is each variable available?
vars_per_study |>
  count(variable, sort = TRUE)

# 4. Study x variable presence matrix (TRUE/FALSE)
presence <- sampleMetadata |>
  group_by(study_name) |>
  summarise(across(everything(), ~ any(!is.na(.x))))
# 4. Blood pressure availability (systolic_p is the only BP variable)
bp <- sampleMetadata |> filter(!is.na(systolic_p))
n_distinct(bp$study_name, bp$subject_id)  # number of people

bp |>
  group_by(study_name) |>
  summarise(samples = n(), people = n_distinct(subject_id)) |>
  arrange(desc(people))
