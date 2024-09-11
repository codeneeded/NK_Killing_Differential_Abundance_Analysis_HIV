library(readxl)
library(dplyr)
library(lme4)
library(ggplot2)
library(ggeffects)
library(lmerTest)
library(tidyr)
library(reshape2)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(plotly)
library(broom.mixed)
library(extrafont)
loadfonts(device = "win")
########### Read in Data and Viral Data ############
setwd("C:/Users/ammas/Documents/NK_Killing_Differential_Abundance_Analysis_HIV")
data <- read_excel("NK Function MetaDATA_REVISED0724v2.xlsx", col_names = TRUE)
v_data <- read.csv("Viral_Titres.csv")
str(data)
########## Clean Column Names and Delete uneeded columns ##############

# Preliminary setup
names(data)[names(data) == "Donor ID"] <- "PID" # Rename column to PID

# Insert a new column "viral load" as the 11th column initialized with NA
data$viral_load <- NA
# This step involves creating a new column order vector
new_order <- c(names(data)[1:10], "viral_load", names(data)[11:(ncol(data)-1)])

# Reorder the dataframe columns according to the new order
data <- data[, new_order]

# Fix errored entries
data$`age (yrs)` <- as.numeric(data$`age (yrs)`)
data$gender <- as.factor(data$gender)
data$HIV <- as.factor(data$HIV)
data$Group <- as.factor(data$Group)
data$Treatment <- as.factor(data$Treatment)
data$Batch <- as.factor(data$Batch)

data <- data %>%
  mutate(
    Timepoint = if_else(grepl("CE037 V11", `Sample Name`), "12", as.character(Timepoint)),
    `age (yrs)` = if_else(grepl("CE037 V11", `Sample Name`), 0.85, `age (yrs)`)
  )
str(data)

data <- data %>%
  filter(!(`Sample Name` == "A1 CE023 V1 NK only.fcs" & Batch == "Batch 9"))

######### Add VIral Load Data ###########

# Step 1 & 2: Update "viral load" based on "Group" and "HIV"
data$viral_load[data$Group != "TARA"] <- NA
data$viral_load[data$HIV == "negative"] <- 0

# Step 3 & 4: Extract and set "viral load" for specific conditions
for (i in 1:nrow(data)) {
  if (data$Group[i] == "TARA" && data$HIV[i] == "positive") {
    if (data$Timepoint[i] == "12") {
      # Extract from v_data where PID matches and Age is 12
      matched <- v_data[v_data$PID == data$PID[i] & v_data$Age == 12,]
    } else if (data$Timepoint[i] == "Entry" && grepl(" V1", data$`Sample Name`[i])) {
      # Extract from v_data where PID matches and Age is 1
      matched <- v_data[v_data$PID == data$PID[i] & v_data$Age == 1,]
    } else if (data$Timepoint[i] == "Entry" && grepl(" V2", data$`Sample Name`[i])) {
      # Extract from v_data where PID matches and Age is 2
      matched <- v_data[v_data$PID == data$PID[i] & v_data$Age == 2,]
    } else {
      next
    }
    
    # Update the "viral load" if matched entry found
    if (nrow(matched) > 0) {
      data$viral_load[i] <- matched$Viral.Titre[1] # Assuming the first match is taken if multiple
    }
  }
}

# Restore original column names with spaces or special characters as necessary
names(data) <- gsub("viral_load", "viral load", names(data))

data_cleaned <- data %>% select(-c(1, 3, 4))

levels (as.factor(data_cleaned$PID))

######### Split Data based on MFI and Frequency ############

col_names <- names(data_cleaned)
col_names_excluding_killing <- col_names[-((ncol(data_cleaned)-2):ncol(data_cleaned))]

# Identify indices for the shared first 9 columns and the last 3 columns
shared_cols <- 1:11
killing_additional_cols <- (ncol(data_cleaned)-2):ncol(data_cleaned) # Last 3 columns

# Identify columns for MFI and Freq excluding the last three columns
mfi_cols <- grep("Median", col_names_excluding_killing)
freq_cols <- grep("\\(\\%\\)", col_names_excluding_killing, perl = TRUE)

# Map these columns back to the original dataframe's column indices
mfi_indices <- match(col_names_excluding_killing[mfi_cols], col_names)
freq_indices <- match(col_names_excluding_killing[freq_cols], col_names)

# Creating the dataframes
MFI <- data_cleaned[, unique(c(shared_cols,killing_additional_cols, mfi_indices))]
Freq <- data_cleaned[, unique(c(shared_cols, killing_additional_cols,freq_indices))]

#### Clean Frequency Dataset prior to Splitting, and convert % to raw counts ########
# Replace portions of column names
colnames(Freq) <- gsub('Lymphocytes/Single Cells/Live/CD3-Label-/Total NK', 'Total NK', colnames(Freq))
colnames(Freq) <- gsub('Lymphocytes/Single Cells/Live/CD3-Label-', 'P4', colnames(Freq))
colnames(Freq) <- gsub('Lymphocytes/Single Cells/Live', 'P3', colnames(Freq))
colnames(Freq) <- gsub('Lymphocytes/Single Cells', 'P2', colnames(Freq))
colnames(Freq) <- gsub('Lymphocytes', 'P1', colnames(Freq))

names(Freq) <- gsub(" \\| Freq\\. of Parent \\(%\\)", "", names(Freq))

# Convert percentages to actual values
Freq$P1 <- Freq$Count * (Freq$P1 / 100) # P1 as a percentage of Count
Freq$P2 <- Freq$P1 * (Freq$P2 / 100)    # P2 as a percentage of the new P1 value
Freq$P3 <- Freq$P2 * (Freq$P3 / 100)    # P3 as a percentage of the new P2 value
Freq$P4 <- Freq$P3 * (Freq$P4 / 100)    # P4 as a percentage of the new P3 value
Freq$`Total NK` <- Freq$P4 * (Freq$`Total NK` / 100)    # P4 as a percentage of the new P3 value

str(Freq)
# Step 1: Identify column groups
single_slash_cols <- grep("^Total NK/[^/]+$", names(Freq), value = TRUE)
double_slash_cols <- grep("^Total NK/.+/[^/]+$", names(Freq), value = TRUE)

##### Checking for non-numeric values
# Combine all columns of interest into a single vector
all_cols <- c("Total NK", single_slash_cols, double_slash_cols)

# Identify non-numeric columns from the list
non_numeric_cols <- all_cols[sapply(Freq[all_cols], function(x) !is.numeric(x))]
# Convert identified non-numeric columns to numeric
Freq[non_numeric_cols] <- lapply(Freq[non_numeric_cols], function(x) as.numeric(as.character(x)))

# Verify changes
str(Freq[non_numeric_cols])

######## Split based on cohort ###########

# Split the dataframe based on the Group column
split_Freq <- split(Freq, Freq$Group)

# Access the dataframe for "Florah"
florah_Freq <- split_Freq[["Florah"]]

# Access the dataframe for "PAVE KL"
pave_kl_Freq <- split_Freq[["PAVE KL"]]

# Access the dataframe for "TARA"
tara_Freq <- split_Freq[["TARA"]]
# Set all values in columns  that are less than 1% to 0
tara_Freq[, 15:ncol(tara_Freq)][tara_Freq[, 15:ncol(tara_Freq)] < 1] <- 0



#### Impute NA Values #########



# Columns to impute (example; adjust as needed)
columns_to_impute <- c("Total NK", single_slash_cols, double_slash_cols)

# Tara

# Impute NAs based on the median for the HIV status group
tara_Freq <- tara_Freq %>%
  mutate(across(all_of(columns_to_impute), ~if_else(is.na(.),
                                                    ave(., HIV, FUN = function(x) median(x, na.rm = TRUE)),
                                                    .),
                .names = "{.col}"))

# Florah
florah_Freq <- florah_Freq %>%
  mutate(across(all_of(columns_to_impute), ~if_else(is.na(.),
                                                    ave(., HIV, FUN = function(x) median(x, na.rm = TRUE)),
                                                    .),
                .names = "{.col}"))
# Pave Kl

pave_kl_Freq <- pave_kl_Freq %>%
  mutate(across(all_of(columns_to_impute), ~if_else(is.na(.),
                                                    ave(., HIV, FUN = function(x) median(x, na.rm = TRUE)),
                                                    .),
                .names = "{.col}"))



# Step 2: Calculate values for direct subsets of "Total NK"

# Tara
for(col in single_slash_cols) {
  tara_Freq[[col]] <- tara_Freq[["Total NK"]] * (tara_Freq[[col]] / 100)
}


# Srep 3 Double / cols
for(col in double_slash_cols) {
  # Extract the parent subset name from the column name
  parts <- strsplit(col, "/")[[1]]
  parent_subset_name <- paste(parts[1:length(parts)-1], collapse="/")
  
  # Ensure parent_subset_name is in single_slash_cols or double_slash_cols
  if(parent_subset_name %in% names(tara_Freq)) {
    parent_value <- tara_Freq[[parent_subset_name]]
    tara_Freq[[col]] <- parent_value * (tara_Freq[[col]] / 100)
  }
}

# Update the dataframe column names to reflect the calculated values
names(tara_Freq) <- gsub("Total NK/", "", names(tara_Freq))

# Florah
for(col in single_slash_cols) {
  florah_Freq[[col]] <- florah_Freq[["Total NK"]] * (florah_Freq[[col]] / 100)
}


# Srep 3 Double / cols
for(col in double_slash_cols) {
  # Extract the parent subset name from the column name
  parts <- strsplit(col, "/")[[1]]
  parent_subset_name <- paste(parts[1:length(parts)-1], collapse="/")
  
  # Ensure parent_subset_name is in single_slash_cols or double_slash_cols
  if(parent_subset_name %in% names(florah_Freq)) {
    parent_value <- florah_Freq[[parent_subset_name]]
    florah_Freq[[col]] <- parent_value * (florah_Freq[[col]] / 100)
  }
}

# Update the dataframe column names to reflect the calculated values
names(florah_Freq) <- gsub("Total NK/", "", names(florah_Freq))

# Pave KL
for(col in single_slash_cols) {
  pave_kl_Freq[[col]] <- pave_kl_Freq[["Total NK"]] * (pave_kl_Freq[[col]] / 100)
}


# Srep 3 Double / cols
for(col in double_slash_cols) {
  # Extract the parent subset name from the column name
  parts <- strsplit(col, "/")[[1]]
  parent_subset_name <- paste(parts[1:length(parts)-1], collapse="/")
  
  # Ensure parent_subset_name is in single_slash_cols or double_slash_cols
  if(parent_subset_name %in% names(pave_kl_Freq)) {
    parent_value <- pave_kl_Freq[[parent_subset_name]]
    pave_kl_Freq[[col]] <- parent_value * (pave_kl_Freq[[col]] / 100)
  }
}

# Update the dataframe column names to reflect the calculated values
names(pave_kl_Freq) <- gsub("Total NK/", "", names(pave_kl_Freq))

### Standardise all values to Total NK

# Step 1: Identify the `Total NK` column
total_nk_col <- tara_Freq$`Total NK`

# Step 2: Convert values after column 19 to percentages of `Total NK`
tara_Freq[, 20:ncol(tara_Freq)] <- round(tara_Freq[, 20:ncol(tara_Freq)] / total_nk_col * 100,2)

######
# Assuming your dataframe is named df and the column range is from column A to column B
cols_to_check <- 15:ncol(tara_Freq)

# Calculate the proportion of 0s in the specified columns
cols_to_keep <- colMeans(tara_Freq[, cols_to_check] == 0) <= 0.8

# Keep only the columns within the specified range that have 80% or fewer 0s
tara_Freq <- tara_Freq[, c(names(tara_Freq)[-cols_to_check], names(tara_Freq)[cols_to_check][cols_to_keep])]
########
###### NK KILLING Analysis #############

### Cleaning Dataframe ###


tara_Freq$Treatment <- droplevels(tara_Freq$Treatment)
levels(tara_Freq$Treatment)

### Mixed Model Functions ####
filter_top_effects <- function(data_frame, significance_threshold = 0.05) {
  # Step 1: Filter significant effects
  significant_effects <- data_frame %>%
    filter(P.Value <= significance_threshold)
  
  # Step 2: Calculate how many more effects are needed to reach 10
  remaining_needed <- 10 - nrow(significant_effects)
  
  # Step 3: If fewer than 10 significant effects, pick the remaining non-significant effects
  if (remaining_needed > 0) {
    # Filter non-significant effects
    non_significant_effects <- data_frame %>%
      filter(P.Value > significance_threshold)
    
    # Split non-significant effects into positive and negative
    positive_non_significant <- non_significant_effects %>%
      filter(Estimate > 0) %>%
      arrange(desc(Estimate)) # Sort positive by descending Estimate
    
    negative_non_significant <- non_significant_effects %>%
      filter(Estimate < 0) %>%
      arrange(Estimate) # Sort negative by ascending Estimate
    
    # Determine how many positive and negative we can pick
    max_to_pick <- ceiling(remaining_needed / 2)
    pos_to_pick <- min(nrow(positive_non_significant), max_to_pick)
    neg_to_pick <- min(nrow(negative_non_significant), remaining_needed - pos_to_pick)
    
    # If we couldn't get enough negatives, pick more positives
    if (pos_to_pick + neg_to_pick < remaining_needed) {
      pos_to_pick <- remaining_needed - neg_to_pick
    }
    
    # Select the rows
    selected_positives <- positive_non_significant %>% head(pos_to_pick)
    selected_negatives <- negative_non_significant %>% head(neg_to_pick)
    
    # Combine significant effects with the selected non-significant effects
    top_effects <- bind_rows(significant_effects, selected_positives, selected_negatives)
  } else {
    # If there are already 10 or more significant effects, just return the top 10
    top_effects <- significant_effects %>%
      arrange(desc(abs(Estimate))) %>%
      head(10)
  }
  
  return(top_effects)
}
filter_top_effects_fas <- function(data_frame, significance_threshold = 0.05) {
  # Filter significant effects
  significant_effects <- data_frame %>%
    filter(P.Value <= significance_threshold)
  
  # Ensure CD56dimCD16+/FasL is included if it's not in significant effects
  if (!"CD56dimCD16+/FasL" %in% significant_effects$Subset) {
    fasl_row <- data_frame %>%
      filter(Subset == "CD56dimCD16+/FasL")
    significant_effects <- bind_rows(significant_effects, fasl_row)
  }
  
  # Calculate how many more effects are needed
  remaining_needed <- 10 - nrow(significant_effects)
  
  # Pick remaining non-significant effects to fill up to 10
  if (remaining_needed > 0) {
    non_significant_effects <- data_frame %>%
      filter(P.Value > significance_threshold) %>%
      filter(Subset != "CD56dimCD16+/FasL") %>%
      arrange(desc(abs(Estimate))) %>%
      head(remaining_needed)
    
    top_effects <- bind_rows(significant_effects, non_significant_effects)
  } else {
    top_effects <- significant_effects %>%
      arrange(desc(abs(Estimate))) %>%
      head(10)
  }
  
  return(top_effects)
}
fit_lmer_models <- function(data, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05) {
  models_list <- list()
  
  # Loop through each subset, ensuring column names are correctly handled
  for (subset_name in flow_population_columns) {
    # Escape subset_name if it contains special characters
    escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
    
    # Construct the formula string with escaped column names
    formula_str <- paste(fixed_part_of_formula, "+", escaped_subset_name)
    
    # Convert the string to a formula
    current_formula <- as.formula(formula_str)
    
    # Fit the model using the current formula
    model <- tryCatch({
      lmer(current_formula, data = data)
    }, error = function(e) {
      cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
      return(NULL)  # Return NULL if there was an error fitting the model
    })
    
    # Store the model if successfully fitted
    if (!is.null(model)) {
      models_list[[subset_name]] <- model
    }
  }
  
  # Initialize an empty data frame to store the results
  results_df <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                           Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each model to extract information
  for (subset_name in names(models_list)) {
    model <- models_list[[subset_name]]
    summary_model <- summary(model)
    df <- as.data.frame(summary_model$coefficients)
    df$Subset <- subset_name
    df$Effect <- rownames(df)
    results_df <- rbind(results_df, df)
  }
  
  # Rename columns to make them consistent
  results_df <- results_df %>%
    rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
    select(Subset, Effect, Estimate, `Std.Error`, P.Value)
  
  # Escape subset names to match Effect names in the results
  escaped_subset_names <- paste0("`", gsub("/", "\\/", flow_population_columns), "`")
  
  # Filter results to include only significant effects
  filtered_results_df <- results_df %>%
    filter(Effect %in% escaped_subset_names)
  
  # Filter for significant effects
  significant_results <- filtered_results_df %>%
    #filter(P.Value < p_value_threshold) %>%
    filter(Std.Error < 1) %>%
    mutate(Significance = ifelse(P.Value < 0.05, "*", ""),
           Color = ifelse(Estimate > 0, "lightblue", "lightgreen"))
  
  return(significant_results)
}

tara_Freq$`Specific Killing` <- as.numeric(tara_Freq$`Specific Killing`)
tara_Freq_plot <- tara_Freq %>% drop_na(14)
tara_Freq_plot_filtered <- tara_Freq_plot %>%
  filter(!Treatment %in% "untreated", HIV != "HUU")
tara_Freq_plot_filtered <- tara_Freq_plot_filtered %>%
  mutate(HIV = recode(HIV, "positive" = "HEI", "negative" = "HEU"))
# Ensure Timepoint is a factor and reorder it so "Entry" comes before "12"
tara_Freq_plot_filtered <- tara_Freq_plot_filtered %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Entry", "12")))


######## Paired Plots #####
# Plot the updated data
library(ggpubr)
setwd("C:/Users/ammas/Documents/NK_Killing_Differential_Abundance_Analysis_HIV")

ggplot(tara_Freq_plot_filtered, aes(x = Timepoint, y = `Specific Killing`, group = PID, color = HIV)) +
  geom_line(aes(linetype = HIV), size = 1) +  # Set line size
  geom_point(size = 2) +  # Set point size
  facet_grid(HIV ~ Treatment, scales = "free_y") +  # Splitting by both HIV and Treatment
  theme_minimal(base_size = 15) +  # Increase base font size for readability
  labs(
    title = "Specific Killing Across Treatments, Timepoints, and HIV Status",
    y = "Specific Killing (%)", 
    x = "Timepoint"
  ) +
  scale_color_manual(values = c("HEI" = "red", "HEU" = "blue")) +  # Custom colors for HEI and HEU
  scale_linetype_manual(values = c("HEI" = "solid", "HEU" = "dashed")) +  # Custom linetypes for HEI and HEU
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    strip.text = element_text(size = 14),  # Adjust facet label size
    legend.title = element_blank(),  # Remove legend title
    legend.position = "top"  # Move the legend to the top
  )
ggsave("Specific_Killing_Paired_Plot.png", width = 10, height = 8, dpi = 300,bg='white')

##### HUT Mechanism ########

#### Clean and create the HUT dataframe ####
#### Seperate HUT Data

# Separate the HUT data
TARA_HUT78_Freq <- tara_Freq %>%
  filter(Treatment == "HUT78")

# Filter unneeded columna
TARA_HUT78_Freq$Group <- NULL
TARA_HUT78_Freq$`age group` <- NULL
TARA_HUT78_Freq$`age (yrs)` <- NULL
TARA_HUT78_Freq$Target <- NULL
TARA_HUT78_Freq$`Target/Dead`<- NULL
TARA_HUT78_Freq$Timepoint <- as.factor(TARA_HUT78_Freq$Timepoint)
TARA_HUT78_Freq$Timepoint <- factor(TARA_HUT78_Freq$Timepoint, levels = rev(levels(TARA_HUT78_Freq$Timepoint)))

TARA_HUT78_Freq <- TARA_HUT78_Freq %>%
  filter(HIV != "HUU")

# Clean

TARA_HUT78_Freq <- TARA_HUT78_Freq[!is.na(TARA_HUT78_Freq$`Specific Killing`), ]
TARA_HUT78_Freq <- droplevels(TARA_HUT78_Freq)
TARA_HUT78_Freq$`Specific Killing`<- as.numeric(TARA_HUT78_Freq$`Specific Killing`)
TARA_HUT78_Freq$`Specific Killing`<- round(TARA_HUT78_Freq$`Specific Killing`, 2)


# Standardize the 'viral load' variable because of high values

TARA_HUT78_Freq$`viral load` <- scale(TARA_HUT78_Freq$`viral load`)

# List of columns representing flow populations, starting from column 13 onwards
flow_population_columns <- colnames(TARA_HUT78_Freq)[15:ncol(TARA_HUT78_Freq)] 

##### FLOW HUT #######

setwd("C:/Users/ammas/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/Mixed_Model_Plots/HUT_Mechanism")
# Initialize lists to store results

fixed_part_of_formula <- "`Specific Killing` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_HUT78_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05) 
filter_top_effects(HUT78)
# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on HUT78 Specific Killing (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill
filter_top_effects(HUT78)
ggsave("Flow_Effects_on_Specific_Killing_HUT78_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')
library(ggridges)

ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, size = abs(Estimate), color = Estimate)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  scale_color_gradient(low = "cyan", high = "magenta") +
  geom_text(aes(label = Significance), colour = "darkred", hjust = .5, vjust = .75, size = 4) +  # Center significance star
  coord_flip() +
  labs(title = "Bubble Plot of HIV Effect on Specific Killing", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  theme(
    text = element_text(family = "Calibri"),  # Change font to Calibri
    axis.text.y = element_text(size = 12),  # Make row names larger and thicker
    axis.title.y = element_text(size = 14),  # Larger y-axis title
    axis.title.x = element_text(size = 14),  # Larger x-axis title
    plot.title = element_text(size = 16, face = "bold"),  # Larger and bold title
    legend.position = "right"  # Keep the color legend
  ) +
  guides(size = "none")  # Remove the legend for size

HUT78_filtered <- HUT78 %>%
  filter(nchar(Subset) <= 30)
HUT78_filtered <- HUT78_filtered %>%
  filter(Subset != "CD107a-Granzyme B-Perforin+")
HUT78_filtered <- HUT78_filtered %>%
  filter(P.Value <0.05)

ggplot(HUT78_filtered, aes(x = reorder(Subset, Estimate), y = Estimate, color = Estimate)) +
  geom_segment(aes(x = reorder(Subset, Estimate), xend = reorder(Subset, Estimate), y = 0, yend = Estimate), size = 1) +
  geom_point(size = 4) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on HUT78 Specific Killing", 
       x = "Subsets", 
       y = "Effect Size") +
  ylim(-1.5, 2) +  # Set the y-axis limits from -1.5 to 1.5
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(family = "Calibri"), # Set font to Calibri
        axis.text.x = element_text(size = 8, family = "Calibri"), 
        axis.text.y = element_text(size = 9, family = "Calibri"), 
        plot.title = element_text(hjust = 0.5, size = 15, family = "Calibri")) # Center th

ggsave("Flow_Effects_on_Specific_Killing_HUT78_TARA_Freq_significant_only_plot_lolipop.png", width = 7, height = 4, dpi = 300,bg='white')

##### IFNg HUT #######
fixed_part_of_formula <- "`CD56dimCD16+/IFNy` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_HUT78_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05) 

# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/IFNy` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Flow_Effects_on_IFNy_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')

##### CD56dimCD16+/MIP-1B HUT #######
# Initialize lists to store results

fixed_part_of_formula <- "`CD56dimCD16+/MIP-1B` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_HUT78_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05) 

# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/MIP-1B` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Flow_Effects_on_MIP-1B_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')



##### CD56dimCD16+/CD107A HUT #######

fixed_part_of_formula <- "`CD56dimCD16+/CD107a` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_HUT78_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05)

# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/CD107a` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Flow_Effects_on_CD107a_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')


#### K562 Mechanism ####


#### Clean and create the HUT dataframe ####
#### Seperate K562  Data

# Separate the K562 data
TARA_K562_Freq <- tara_Freq %>%
  filter(Treatment == "K562")

# Filter unneeded columna
TARA_K562_Freq$Group <- NULL
TARA_K562_Freq$`age group` <- NULL
TARA_K562_Freq$`age (yrs)` <- NULL
TARA_K562_Freq$Target <- NULL
TARA_K562_Freq$`Target/Dead`<- NULL
TARA_K562_Freq$Timepoint <- as.factor(TARA_K562_Freq$Timepoint)
TARA_K562_Freq$Timepoint <- factor(TARA_K562_Freq$Timepoint, levels = rev(levels(TARA_K562_Freq$Timepoint)))

TARA_K562_Freq <- TARA_K562_Freq %>%
  filter(HIV != "HUU")

# Clean

TARA_K562_Freq <- TARA_K562_Freq[!is.na(TARA_K562_Freq$`Specific Killing`), ]
TARA_K562_Freq <- droplevels(TARA_K562_Freq)
TARA_K562_Freq$`Specific Killing`<- as.numeric(TARA_K562_Freq$`Specific Killing`)
TARA_K562_Freq$`Specific Killing`<- round(TARA_K562_Freq$`Specific Killing`, 2)


# Standardize the 'viral load' variable because of high values

TARA_K562_Freq$`viral load` <- scale(TARA_K562_Freq$`viral load`)

# List of columns representing flow populations, starting from column 13 onwards
flow_population_columns <- colnames(TARA_K562_Freq)[15:ncol(TARA_K562_Freq)] 

##### FLOW K562 #######

setwd("C:/Users/ammas/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/Mixed_Model_Plots/K562_Mechanism")
# Initialize lists to store results

fixed_part_of_formula <- "`Specific Killing` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
K562 <- fit_lmer_models(TARA_K562_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05) 

# Create the plot with only significant effects
ggplot(filter_top_effects(K562), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on K562 Specific Killing (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

#K562_filtered <- K562 %>%
#  filter(nchar(Subset) <= 30)
K562_filtered <- K562 %>%
  filter(Subset != "CD56dimCD16+/CD107a total" & Subset != "CD107a-"
         & Subset != "CD107a total"& Subset != "CD107a+Granzyme B+Perforin+"
         & Subset != "CD107a-Granzyme B+Perforin+")


K562_filtered <- K562_filtered %>%
  filter(P.Value <0.05)

ggplot(K562_filtered, aes(x = reorder(Subset, Estimate), y = Estimate, color = Estimate)) +
  geom_segment(aes(x = reorder(Subset, Estimate), xend = reorder(Subset, Estimate), y = 0, yend = Estimate), size = 1) +
  geom_point(size = 4) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  ylim(-1.5, 2) +  # Set the y-axis limits from -1.5 to 1.5
  labs(title = "Effect of Flow Subsets on K562 Specific Killing", 
       x = "Subsets", 
       y = "Effect Size") +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(family = "Calibri"), # Set font to Calibri
        axis.text.x = element_text(size = 8, family = "Calibri"), 
        axis.text.y = element_text(size = 9, family = "Calibri"), 
        plot.title = element_text(hjust = 1, size = 15, family = "Calibri")) # Center th
ggsave("Flow_Effects_on_Specific_Killing_K562_TARA_Freq_significant_plot_lolipop.png", width = 7, height = 4, dpi = 300,bg='white')


##### IFNg K562 #######
fixed_part_of_formula <- "`CD56dimCD16+/IFNy` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
K562 <- fit_lmer_models(TARA_K562_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05) 

# Create the plot with only significant effects
ggplot(filter_top_effects(K562), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/IFNy` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Flow_Effects_on_IFNy_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')

##### CD56dimCD16+/MIP-1B K562 #######
# Initialize lists to store results

fixed_part_of_formula <- "`CD56dimCD16+/MIP-1B` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
K562 <- fit_lmer_models(TARA_K562_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05) 

# Create the plot with only significant effects
ggplot(filter_top_effects(K562), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/MIP-1B` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Flow_Effects_on_MIP-1B_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')



##### CD56dimCD16+/CD107A K562 #######

fixed_part_of_formula <- "`CD56dimCD16+/CD107a` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
K562 <- fit_lmer_models(TARA_K562_Freq, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05)

# Create the plot with only significant effects
ggplot(filter_top_effects(K562), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/CD107a` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Flow_Effects_on_CD107a_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')






### Combine untreated with Specific Killing ###


# Separate the untreated data
untreated_data <- tara_Freq %>%
  filter(Treatment == "untreated")

# Create a function to merge specific killing values for each treatment
merge_specific_killing <- function(data, treatment) {
  data %>%
    left_join(
      tara_Freq %>%
        filter(Treatment == treatment) %>%
        select(PID, Timepoint, !!paste0(treatment, "_Specific_Killing") := `Specific Killing`),
      by = c("PID", "Timepoint")
    )
}

# List of treatments to be added as columns
treatments <- levels(tara_Freq$Treatment)[levels(tara_Freq$Treatment) != "untreated"]

# Loop over each treatment and merge specific killing values
TARA_Killing <- untreated_data
for (treatment in treatments) {
  TARA_Killing <- merge_specific_killing(TARA_Killing, treatment)
}

# Removing NA specific Killings as well as redundant columns


TARA_Killing$Group <- NULL
TARA_Killing$`age group` <- NULL
TARA_Killing$`age (yrs)` <- NULL
TARA_Killing$Target <- NULL
TARA_Killing$`Target/Dead`<- NULL
TARA_Killing$Timepoint <- as.factor(TARA_Killing$Timepoint)
TARA_Killing$Timepoint <- factor(TARA_Killing$Timepoint, levels = rev(levels(TARA_Killing$Timepoint)))
TARA_Killing$`Specific Killing`<-NULL

TARA_Killing <- TARA_Killing %>%
  filter(HIV != "HUU")



# Update the Viral_Load column for PID CP018 at Timepoint entry
TARA_Killing <- TARA_Killing %>%
  mutate(`viral load` = ifelse(PID == "CP018" & Timepoint == "Entry", 176970, `viral load`))



# Standardize the 'viral load' variable because of high values

TARA_Killing$`viral load` <- scale(TARA_Killing$`viral load`)

###### HUT Prediction (Untreated - HUT) ############
# Round the values in the 'specific killing' column to two decimal points
TARA_Killing$HUT78_Specific_Killing <- as.numeric(TARA_Killing$HUT78_Specific_Killing)
TARA_Killing$HUT78_Specific_Killing <- round(TARA_Killing$HUT78_Specific_Killing, 2)
TARA_Killing$K562_Specific_Killing <- as.numeric(TARA_Killing$K562_Specific_Killing)
TARA_Killing$K562_Specific_Killing <- round(TARA_Killing$K562_Specific_Killing, 2)

# List of columns representing flow populations, starting from column 13 onwards
flow_population_columns <- colnames(TARA_Killing)[14:(ncol(TARA_Killing)-7)] 

TARA_Killing <- TARA_Killing[!is.na(TARA_Killing$HUT78_Specific_Killing), ]

##### FLOW HUT #######
setwd("C:/Users/ammas/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/Mixed_Model_Plots/Untreated_HUT_Prediction")

# Initialize lists to store results

fixed_part_of_formula <- "HUT78_Specific_Killing ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_Killing, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05)

# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on HUT78 Specific Killing (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal()

ggsave("Untreated_Flow_Effects_on_Specific_Killing_HUT78_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')


HUT78_filtered <- HUT78 %>%
  filter(P.Value <0.05)

ggplot(HUT78_filtered, aes(x = reorder(Subset, Estimate), y = Estimate, color = Estimate)) +
  geom_segment(aes(x = reorder(Subset, Estimate), xend = reorder(Subset, Estimate), y = 0, yend = Estimate), size = 1) +
  geom_point(size = 4) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on HUT78 Specific Killing", 
       x = "Subsets", 
       y = "Effect Size") +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(family = "Calibri"), # Set font to Calibri
        axis.text.x = element_text(size = 8, family = "Calibri"), 
        axis.text.y = element_text(size = 9, family = "Calibri"), 
        plot.title = element_text(hjust = 0.5, size = 15, family = "Calibri")) # Center th
ggsave("Untreated_Flow_Effects_on_Specific_Killing_HUT78_TARA_Freq_significant_only_plot_lolipop.png", width = 7, height = 4, dpi = 300,bg='white')


##### IFNg HUT #######
fixed_part_of_formula <- "`CD56dimCD16+/IFNy` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_Killing, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05)


# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/IFNy` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Untreated_Flow_Effects_on_IFNy_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')

##### CD56dimCD16+/MIP-1B HUT #######
# Initialize lists to store results

fixed_part_of_formula <- "`CD56dimCD16+/MIP-1B` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_Killing, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05)

# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/MIP-1B` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Untreated_Flow_Effects_on_MIP-1B_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')



##### CD56dimCD16+/CD107A HUT #######
# Initialize lists to store results
fixed_part_of_formula <- "`CD56dimCD16+/CD107a` ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
HUT78 <- fit_lmer_models(TARA_Killing, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05)

# Create the plot with only significant effects
ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on `CD56dimCD16+/CD107a` (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill = FALSE)  # Hide the legend for fill

ggsave("Untereated_Flow_Effects_on_CD107a_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')

###### Untreated vs K562 ############


##### FLOW HUT #######
setwd("C:/Users/ammas/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/Mixed_Model_Plots/Untreated_HUT_Prediction")

# Initialize lists to store results

fixed_part_of_formula <- "K562_Specific_Killing ~ `viral load` + Timepoint  + gender + HIV + (1 | PID)"
K562 <- fit_lmer_models(TARA_Killing, fixed_part_of_formula, flow_population_columns, p_value_threshold = 0.05)

# Create the plot with only significant effects
ggplot(filter_top_effects(K562), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on K562 Specific Killing (Significant Effects Only)", x = "Subsets", y = "Effect Size") +
  theme_minimal()

ggsave("Untreated_Flow_Effects_on_Specific_Killing_K562_TARA_Freq_significant_plot.png", width = 10, height = 8, dpi = 300,bg='white')

K562_filtered <- K562 %>%
  filter(Subset != 'NKG2A-')


K562_filtered <- K562_filtered %>%
  filter(P.Value <0.05)

ggplot(K562_filtered, aes(x = reorder(Subset, Estimate), y = Estimate, color = Estimate)) +
  geom_segment(aes(x = reorder(Subset, Estimate), xend = reorder(Subset, Estimate), y = 0, yend = Estimate), size = 1) +
  geom_point(size = 4) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title = "Effect of Flow Subsets on K562 Specific Killing", 
       x = "Subsets", 
       y = "Effect Size") +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(family = "Calibri"), # Set font to Calibri
        axis.text.x = element_text(size = 8, family = "Calibri"), 
        axis.text.y = element_text(size = 9, family = "Calibri"), 
        plot.title = element_text(hjust = 1, size = 15, family = "Calibri")) # Center th
ggsave("Flow_Effects_on_Specific_Killing_K562_TARA_Freq_significant_plot_lolipop_Untreated.png", width = 7, height = 4, dpi = 300,bg='white')


#################################################################################################
### Differential Abundance analysis
TARA_Freq_DA <- tara_Freq[, !(names(tara_Freq) %in% c("Batch", "Group","age group", "age (yrs)","Count", "P1", "P2", "P3", "P4", "Specific Killing","Target","Target/Dead"))]
TARA_Freq_DA$Timepoint <- as.factor(TARA_Freq_DA$Timepoint)
TARA_Freq_DA <- droplevels(TARA_Freq_DA)
covariate_columns <- c("PID", "gender", "HIV", "viral load", "Timepoint", "Treatment")
TARA_Freq_DA <- TARA_Freq_DA[TARA_Freq_DA$HIV != "HUU", ]
#### Effect of IL-15 ###################
setwd("C:/Users/ammas/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/IL15_Effect")

# Initialize a list to store results
results_list_paired <- list()
# Update the Viral_Load column for PID CP018 at Timepoint entry
TARA_Freq_DA <- TARA_Freq_DA %>%
  mutate(`viral load` = ifelse(PID == "CP018" & Timepoint == "Entry", 176970, `viral load`))
# Define pairs of treatments to compare
treatment_pairs <- list(
  "untreated" = "IL15",
  "CEM" = "CEM+IL15",
  "HUT78" = "HUT78+IL15",
  "K562" = "K562+IL15"
)
# Process each treatment pair
for (baseline in names(treatment_pairs)) {
  treated <- treatment_pairs[[baseline]]
  
  # Filter data for only the current pair of treatments
  pair_data <- TARA_Freq_DA[TARA_Freq_DA$Treatment %in% c(baseline, treated), ]
  
  # Create a factor to distinguish between baseline and treated groups
  pair_data$group <- factor(ifelse(pair_data$Treatment == baseline, "baseline", "treated"))
  names(pair_data)[names(pair_data) == "viral load"] <- "viral_load"
  # Check if there are enough data points to proceed
  if(nrow(pair_data) > 1) {
    # Design matrix including Timepoint as a factor and HIV status
    design <- model.matrix(~0+ group + Timepoint + HIV  +viral_load, data = pair_data)
    
    # Voom transformation
    v <- voom(t(pair_data[, !(names(pair_data) %in% c("PID", "HIV", "gender", "Timepoint","viral_load", "Treatment", "group"))]), 
              design, plot = FALSE)
    
    # Fit the linear model
    corfit <- duplicateCorrelation(v, design, block = pair_data$PID)
    fit <- lmFit(v, design, block = pair_data$PID, correlation = corfit$consensus)
    
    # Specify the correct contrast based on the baseline and treated groups
    contrast.matrix <- makeContrasts(grouptreated - groupbaseline, levels = design)
    
    # Apply contrasts
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # Store results, specifying treatment pair
    results_key <- paste(baseline, "vs", treated, sep="_")
    results_list_paired[[results_key]] <- topTable(fit2, adjust="BH", number=1000)
  }
}



## VOlcano Plots
results_list_paired$`untreated_vs_IL15`$feature<- rownames(results_list_paired$`untreated_vs_IL15`)
results_list_paired$`untreated_vs_IL15`$abs_logFC<- abs(results_list_paired$`untreated_vs_IL15`$logFC)

results_list_paired$`CEM_vs_CEM+IL15`$feature<- rownames(results_list_paired$`CEM_vs_CEM+IL15`)
results_list_paired$`CEM_vs_CEM+IL15`$abs_logFC<- abs(results_list_paired$`CEM_vs_CEM+IL15`$logFC)

results_list_paired$`HUT78_vs_HUT78+IL15`$feature<- rownames(results_list_paired$`HUT78_vs_HUT78+IL15`)
results_list_paired$`HUT78_vs_HUT78+IL15`$abs_logFC<- abs(results_list_paired$`HUT78_vs_HUT78+IL15`$logFC)

results_list_paired$`K562_vs_K562+IL15`$feature<- rownames(results_list_paired$`K562_vs_K562+IL15`)
results_list_paired$`K562_vs_K562+IL15`$abs_logFC<- abs(results_list_paired$`K562_vs_K562+IL15`$logFC)

### Untreated
library(ggrepel)
# Create the volcano plot and add labels for significant points
Untreated_volcano_plot <- ggplot(results_list_paired$`untreated_vs_IL15`, 
                                 aes(x = logFC, y = -log10(P.Value), 
                                     text = paste("Cell Type:", feature, 
                                                  "<br>LogFC:", logFC, 
                                                  "<br>Adj. P-Val:", adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot NK Cell Fold Change  IL15 vs Untreated", 
       x = "Log2 Fold Change (IL15 vs Untreated)", y = "-Log10 P-value") +
  
  # Add labels for significant points with lines connecting them
  geom_label_repel(aes(label = ifelse(adj.P.Val < 0.05, as.character(feature), "")), 
                   size = 3, color = "black", 
                   box.padding = 0.35, point.padding = 0.5, 
                   segment.color = "grey50", max.overlaps = 15)  # Adjust max.overlaps as needed

# Print the plot
print(Untreated_volcano_plot)
ggsave("Untreated_volcano_plot.png", plot = Untreated_volcano_plot, width = 10, height = 8, dpi = 300,bg='white')


### CEM

# Create the volcano plot and add labels for significant points
CEM_volcano_plot <- ggplot(results_list_paired$`CEM_vs_CEM+IL15`, 
                                 aes(x = logFC, y = -log10(P.Value), 
                                     text = paste("Cell Type:", feature, 
                                                  "<br>LogFC:", logFC, 
                                                  "<br>Adj. P-Val:", adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot NK Cell Fold Change for CEM+IL15 vs CEM ", 
       x = "Log2 Fold Change (CEM+IL15 vs CEM)", y = "-Log10 P-value") +
  
  # Add labels for significant points with lines connecting them
  geom_label_repel(aes(label = ifelse(adj.P.Val < 0.05, as.character(feature), "")), 
                   size = 3, color = "black", 
                   box.padding = 0.35, point.padding = 0.5, 
                   segment.color = "grey50", max.overlaps = 15)  # Adjust max.overlaps as needed


# Print the plot
print(CEM_volcano_plot)
ggsave("CEM_volcano_plot.png", plot = CEM_volcano_plot, width = 10, height = 8, dpi = 300,bg='white')

### HUT

# Create the volcano plot and add labels for significant points
HUT_volcano_plot <- ggplot(results_list_paired$`HUT78_vs_HUT78+IL15`, 
                                 aes(x = logFC, y = -log10(P.Value), 
                                     text = paste("Cell Type:", feature, 
                                                  "<br>LogFC:", logFC, 
                                                  "<br>Adj. P-Val:", adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), size=3,alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot NK Cell Fold Change  HUT78+IL15 vs HUT78", 
       x = "Log2 Fold Change (HUT78+IL15 vs HUT78)", y = "-Log10 P-value") +
  
  # Add labels for significant points with lines connecting them
  geom_label_repel(aes(label = ifelse(adj.P.Val < 0.05, as.character(feature), "")), 
                   size = 5, color = "black", 
                   box.padding = 0.35, point.padding = 0.5, 
                   segment.color = "grey50", max.overlaps = 15)  # Adjust max.overlaps as needed

# Print the plot
print(HUT_volcano_plot)
ggsave("HUT_volcano_plot_large.png", plot = HUT_volcano_plot, width = 10, height = 8, dpi = 300,bg='white')

### K562

# Create the volcano plot and add labels for significant points
K562_volcano_plot <- ggplot(results_list_paired$`K562_vs_K562+IL15`, 
                                 aes(x = logFC, y = -log10(P.Value), 
                                     text = paste("Cell Type:", feature, 
                                                  "<br>LogFC:", logFC, 
                                                  "<br>Adj. P-Val:", adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot NK Cell Fold Change  K562+IL15 vs K562", 
       x = "Log2 Fold Change (K562+IL15 vs K562)", y = "-Log10 P-value") +
  
  # Add labels for significant points
  geom_label_repel(aes(label = ifelse(adj.P.Val < 0.05, as.character(feature), "")), 
                   size = 3, color = "black", 
                   box.padding = 0.35, point.padding = 0.5, 
                   segment.color = "grey50", max.overlaps = 15)  # Adjust max.overlaps as needed

# Print the plot
print(K562_volcano_plot)
ggsave("K562_volcano_plot.png", plot = K562_volcano_plot, width = 10, height = 8, dpi = 300,bg='white')





######################
### Fix missing VL Value






#### Fit the mixed effects model for each treatment condition ####

# Define the treatment levels
treatments <- levels(TARA_Killing$Treatment)

######## Unstimulated ###############

#### Mixed effect model using HUT specific killing as outcome variable ####

# Determine the columns to keep
columns_to_keep <- c(1:(ncol(TARA_Killing) - 7), which(colnames(TARA_Killing) == "HUT78_Specific_Killing"))

# Subset the dataframe
TARA_Killing_HUT <- TARA_Killing[, columns_to_keep]
# Remove rows with NA in HUT78_Specific_Killing column using base R
TARA_Killing_HUT <- TARA_Killing_HUT[!is.na(TARA_Killing_HUT$HUT78_Specific_Killing), ]
TARA_Killing_HUT <- droplevels(TARA_Killing_HUT)
TARA_Killing_HUT$HUT78_Specific_Killing<- as.numeric(TARA_Killing_HUT$HUT78_Specific_Killing)
TARA_Killing_HUT$HUT78_Specific_Killing<- round(TARA_Killing_HUT$HUT78_Specific_Killing, 2)

# Standardize the 'viral load' variable because of high values

TARA_Killing_HUT$`viral load` <- scale(TARA_Killing_HUT$`viral load`)

# List of columns representing flow populations, starting from column 13 onwards
flow_population_columns <- colnames(TARA_Killing_HUT)[13:(ncol(TARA_Killing_HUT) - 1)] 

# Function to run the mixed-effects model and return the summary
run_model_flow <- function(data, flow_col) {
  flow_col_quoted <- paste0("`", flow_col, "`")
  formula <- as.formula(paste("HUT78_Specific_Killing ~", flow_col_quoted, "+ Timepoint + HIV + gender + (1 | PID)"))
  model <- lmer(formula, data = data)
  tidy(model)
}

run_model_flow <- function(data, flow_col) {
  flow_col_quoted <- paste0("`", flow_col, "`")
  formula <- as.formula(paste("HUT78_Specific_Killing ~", flow_col_quoted, "+ Timepoint + HIV + gender + (1 | PID)"))
  model <- lmer(formula, data = data)
  tidy(model)
}

run_model_CD107a <- function(data, flow_col) {
  flow_col_quoted <- paste0("`", flow_col, "`")
  formula <- as.formula(paste("CD56dimCD16+/CD107a ~", flow_col_quoted, "+ Timepoint + HIV + gender + (1 | PID)"))
  model <- lmer(formula, data = data)
  tidy(model)
}
# Collect results for each flow population
results <- list()
for (flow_col in flow_population_columns) {
  model_summary <- run_model(TARA_Killing_HUT, flow_col)
  results[[flow_col]] <- model_summary
}

# Combine results into a single dataframe
results_df <- do.call(rbind, lapply(names(results), function(flow_col) {
  df <- results[[flow_col]]
  df$FlowPop <- flow_col
  return(df)
}))

# Filter results for fixed effects (excluding intercept) and significant effects (p < 0.05)
results_df <- results_df %>%
  filter(term != "(Intercept)" & p.value < 0.05)

# Select top 10 positive and top 10 negative effects
top_positive_results <- results_df %>%
  arrange(desc(estimate)) %>%
  slice(1:20)

top_negative_results <- results_df %>%
  arrange(estimate) %>%
  slice(1:20)

top_results_df <- bind_rows(top_positive_results, top_negative_results)

# Add significance stars and color coding
top_results_df <- top_results_df %>%
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ),
  Color = ifelse(estimate > 0, "lightblue", "lightgreen"))

# Plot the results with significance stars and error bars
ggplot(top_results_df, aes(x = reorder(FlowPop, estimate), y = estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.4) +
  geom_text(aes(label = significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Top 10  Significant Effects of Flow Populations on HUT78 Specific Killing",
       x = "Flow Cytometry Gated Populations",
       y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

  
  # Mixed effect model using CD56dimCD16+/CD107a as outcome variable

# Mixed effect model using CD56dimCD16+/IFNy as outcome variable

# Mixed effect model using CD56dimCD16+/MIP-1B as outcome variable

# Covariate Analysis for HEI vs HEU, for viral load and for timepoint



for (treatment in treatments) {
  # Create a subset for the current treatment
  subset_data <- subset(TARA_Killing, Treatment == treatment)
  
  # Fit a mixed-effects model for the subset
  model_formula <- as.formula(`specific killing` ~ (HIV + Timepoint + gender + Batch)^2 + (1|PID))
  model <- lmer(model_formula, data = subset_data)
  
  # Print the summary of the model
  print(paste("Model for Treatment:", treatment))
  print(summary(model))
  
  # Optionally, save the model for later use
  assign(paste("TARA_Killing_", gsub("[+]", "", treatment), sep = ""), subset_data, envir = .GlobalEnv)
  assign(paste("TARA_Killing_", gsub("[+]", "", treatment), "_model", sep = ""), model, envir = .GlobalEnv)
}


TARA_Killing_CEM$predicted_killing <- predict(TARA_Killing_CEM_model, re.form = NA)
TARA_Killing_CEMIL15$predicted_killing <- predict(TARA_Killing_CEMIL15_model, re.form = NA)
TARA_Killing_HUT78$predicted_killing <- predict(TARA_Killing_HUT78_model, re.form = NA)
TARA_Killing_HUT78IL15$predicted_killing <- predict(TARA_Killing_HUT78IL15_model, re.form = NA)
TARA_Killing_K562$predicted_killing <- predict(TARA_Killing_K562_model, re.form = NA)
TARA_Killing_K562IL15$predicted_killing <- predict(TARA_Killing_K562IL15_model, re.form = NA)



# Example plot for one treatment subset, e.g., TARA_Killing_CEM
ggplot(TARA_Killing_HUT78IL15, aes(x = Timepoint, y = predicted_killing, group = PID, color = HIV)) +
  geom_point() + # Adds the points
  geom_line() + # Connects points from the same individual with a line
  geom_smooth(method = "lm", aes(group = 1), se = FALSE, color = "black") + # Adds a global regression line
  theme_minimal() +
  ylab("Predicted Specific Killing") +
  xlab("Timepoint")


ggplot(TARA_Killing_HUT78IL15, aes(x = Timepoint, y = `specific killing`, group = PID, color = HIV)) +
  geom_point(position = position_dodge(width = 0.2)) + # Adds individual points, slightly dodged to reduce overplotting
  geom_line(position = position_dodge(width = 0.2)) + # Connects points from the same individual
  facet_wrap(~gender) + # Creates separate plots for each gender
  theme_minimal() +
  labs(title = "Change in Specific Killing Across Timepoints by HIV Status",
       y = "Specific Killing", x = "Timepoint") +
  scale_color_brewer(palette = "Set1") # Optional: use a colorblind-friendly palette


ggplot(TARA_Killing_HUT78IL15, aes(x = Timepoint, y = `specific killing`, color = gender)) +
  geom_point(aes(shape = gender)) +
  geom_line(aes(group = PID)) +
  facet_grid(. ~ HIV) +
  theme_minimal() +
  labs(title = "Specific Killing Across Treatments and Timepoints by HIV Status", y = "Specific Killing")

TARA_Killing <- TARA_Killing %>%
  mutate(HIV_Timepoint = paste(HIV, Timepoint, sep = "_"))
heatmap_data <- TARA_Killing %>%
  group_by(Treatment, HIV_Timepoint) %>%
  summarise(mean_killing = mean(`specific killing`, na.rm = TRUE)) %>%
  ungroup()



ggplot(heatmap_data, aes(x = HIV_Timepoint, y = Treatment, fill = mean_killing)) +
  geom_tile() + # Creates the tiles for the heatmap
  scale_fill_viridis_c() + # Uses a color scale that's perceptually uniform
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotates x-axis labels for readability
  labs(title = "Heatmap of Mean Specific Killing by Treatment and HIV/Timepoint",
       x = "HIV Status and Timepoint", y = "Treatment", fill = "Mean Specific Killing")

#######






TARA_Killing_CEM <- subset(TARA_Killing, Treatment == "HUT78")
model_A <- lmer(`specific killing` ~ (HIV + Timepoint + gender)^2 + (1|PID), data = TARA_Killing_CEM)
summary_model_A <- summary(model_A)
coef_table <- summary_model_A$coefficients
coef_df <- as.data.frame(coef_table)
coef_df$Term <- rownames(coef_df)
fixed_effects_confint <- confint_model_A[rownames(coef_df), ]

# Ensure that the row names match between coef_df and the confidence intervals
# This is crucial for correct alignment
# Adding the confidence intervals to the dataframe
coef_df$LowerCI <- fixed_effects_confint[, "2.5 %"]
coef_df$UpperCI <- fixed_effects_confint[, "97.5 %"]


ggplot(coef_df, aes(x = Term, y = Estimate, ymin = LowerCI, ymax = UpperCI)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  coord_flip() +
  ylab("Effect Size (95% CI)") +
  xlab("Predictors")

confint_model_A <- confint(model_A, level = 0.95) # Default is 0.95 confidence level
# Joining the confidence intervals with the coefficients dataframe. Ensure names match or adjust as necessary.
coef_df <- cbind(coef_df, confint_model_A)

# Checking the dataframe
head(coef_df)

summary(model_A)
coef_df$Term <- rownames(coef_df)

ggplot(coef_df, aes(x = Term, y = Estimate, ymin = Estimate - 2. * Std..Error, ymax = Estimate + 2. * Std..Error)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  coord_flip() +
  ylab("Effect size (95% CI)") +
  xlab("Predictors")


TARA_Killing_CEM$predicted_killing <- predict(model_A)

ggplot(TARA_Killing_CEM, aes(x = as.factor(Timepoint), y = predicted_killing, group = PID, color = HIV)) +
  geom_point(position = position_dodge(0.1)) + # Dodge points for clarity
  geom_line(aes(group = interaction(PID, HIV)), position = position_dodge(0.1)) + # Connect points with lines
  geom_smooth(method = "lm", aes(group = HIV), se = FALSE, color = "black", linetype = "dashed") + # Add separate regression lines for each HIV status
  theme_minimal() +
  ylab("Predicted Specific Killing") +
  xlab("Timepoint")





# Summary of the model
# Extracting the summary as a dataframe
# Extract p-values for all coefficients from model summary
p_values <- summary(model)$coefficients[, "Pr(>|t|)"]


###### TREATMENT
# Extract treatment names from preds_treatment
preds_treatment <- ggpredict(model, terms = "Treatment")
treatment_names <- levels(preds_treatment$x)

# Initialize significance vector with NA
treatment_significance <- rep(NA, length(treatment_names))

# Find indices of treatment names in model coefficients
treatment_indices <- match(paste0("Treatment", treatment_names), rownames(summary(model)$coefficients))

# Remove NAs from treatment_indices
valid_indices <- which(!is.na(treatment_indices))

# Assign significance values based on extracted p-values
treatment_significance[valid_indices] <- ifelse(p_values[valid_indices] < 0.05, "*", "")


# Add significance to preds_treatment

preds_treatment$significance <- treatment_significance
preds_treatment$significance[is.na(preds_treatment$significance)] <- ""


str(preds_treatment)

# Plot with annotations (alternative method with significance stars on top of error bars and red color)
ggplot(preds_treatment, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(width = 0.4, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = significance, y = conf.high + 0.1), size = 5, color = "red", position = position_dodge(width = 0.5)) +  # Adjust y position and color
  labs(title = "Effect of Treatment on Specific Killing",
       x = "Treatment",
       y = "Predicted Specific Killing") +
  theme_minimal()


### GENDER

# Extract gender names from preds_gender
preds_gender <- ggpredict(model, terms = "gender")
gender_names <- levels(preds_gender$x)

# Initialize significance vector with empty strings
gender_significance <- rep("", length(gender_names))

# Find indices of gender names in model coefficients
gender_indices <- match(paste0("gender", gender_names), rownames(summary(model)$coefficients))

# Assign significance values based on extracted p-values
gender_p_values <- summary(model)$coefficients[gender_indices, "Pr(>|t|)"]
gender_significance[gender_p_values < 0.05] <- "*"

# Add significance to preds_gender
preds_gender$significance <- gender_significance

# Plot with annotations for gender
ggplot(preds_gender, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(width = 0.4, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = significance, y = conf.high + 0.1), size = 5, color = "red", position = position_dodge(width = 0.5)) +  
  labs(title = "Effect of Gender on Specific Killing",
       x = "Gender",
       y = "Predicted Specific Killing") +
  theme_minimal()



### Timepoint
# Extract Timepoint names from preds_timepoint
preds_timepoint <- ggpredict(model, terms = "Timepoint")
timepoint_names <- levels(preds_timepoint$x)

# Initialize significance vector with empty strings
timepoint_significance <- rep("", length(timepoint_names))

# Find indices of Timepoint names in model coefficients
timepoint_indices <- match(paste0("Timepoint", timepoint_names), rownames(summary(model)$coefficients))

# Assign significance values based on extracted p-values
timepoint_p_values <- summary(model)$coefficients[timepoint_indices, "Pr(>|t|)"]
timepoint_significance[timepoint_p_values < 0.05] <- "*"

# Add significance to preds_timepoint
preds_timepoint$significance <- timepoint_significance

# Plot with annotations for Timepoint
ggplot(preds_timepoint, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(width = 0.4, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = significance, y = conf.high + 0.1), size = 5, color = "red", position = position_dodge(width = 0.5)) +  
  labs(title = "Effect of Timepoint on Specific Killing",
       x = "Timepoint",
       y = "Predicted Specific Killing") +
  theme_minimal()

### HIV
preds_HIV <- ggpredict(model, terms = "HIV")
HIV_names <- levels(preds_HIV$x)

# Initialize significance vector with empty strings
HIV_significance <- rep("", length(HIV_names))

# Find indices of HIV names in model coefficients
HIV_indices <- match(paste0("HIV", HIV_names), rownames(summary(model)$coefficients))

# Assign significance values based on extracted p-values
HIV_p_values <- summary(model)$coefficients[HIV_indices, "Pr(>|t|)"]
HIV_significance[HIV_p_values < 0.05] <- "*"

# Add significance to preds_HIV
preds_HIV$significance <- HIV_significance

# Plot with annotations for HIV
ggplot(preds_HIV, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(width = 0.4, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = significance, y = conf.high + 0.1), size = 5, color = "red", position = position_dodge(width = 0.5)) +  
  labs(title = "Effect of HIV on Specific Killing",
       x = "HIV Status",
       y = "Predicted Specific Killing") +
  theme_minimal()





# Plotting the effect of 'viral load'
# Extract predictions for viral load
preds_viral_load <- ggpredict(model, terms = "viral load")

# Plot effects of viral load
ggplot(preds_viral_load, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
  geom_point() +
  geom_errorbar(width = 0.4) +
  labs(title = "Effect of Viral Load on Specific Killing",
       x = "Viral Load",
       y = "Predicted Specific Killing") +
  theme_minimal()

# Plotting the effect of 'Timepoint'
# Adjusting the plot for 'Timepoint' as a categorical variable
ggplot(preds_timepoint, aes(x = x, y = predicted, color = x)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.4, position = position_dodge(width = 0.5)) +
  labs(title = "Effect of Timepoint on Specific Killing", x = "Timepoint", y = "Predicted Specific Killing") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis label readability

# Adjusting the plot for 'Treatment' as a categorical variable
ggplot(preds_treatment, aes(x = x, y = predicted, color = x)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.4, position = position_dodge(width = 0.5)) +
  labs(title = "Effect of Treatment on Specific Killing", x = "Treatment", y = "Predicted Specific Killing") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve x-axis label readability

###### NK Freq Analysis #############
# Removing the 'Group' column from both dataframes
TARA_Freq$Group <- NULL

# Get all column names
column_names <- names(TARA_Freq)

# Number of columns
n_columns <- length(column_names)

# Move the last column to be the 10th column
# We create a sequence from the first to the (n_columns - 1) to exclude the last column initially,
# then cbind the last column name at the 10th position.
new_order <- c(column_names[-n_columns][1:9], column_names[n_columns], column_names[-n_columns][10:(n_columns-1)])

# Rearrange the dataframe columns based on the new order
TARA_Freq <- TARA_Freq[, new_order]

# Round the values in the 'specific killing' column to two decimal points
TARA_Freq$`specific killing` <- as.numeric(TARA_Freq$`specific killing`)
TARA_Freq$`specific killing` <- round(TARA_Freq$`specific killing`, 2)

TARA_Freq <- TARA_Freq %>%
  filter(HIV != "HUU")

TARA_Freq[, 11:ncol(TARA_Freq)] <- TARA_Freq[, 11:ncol(TARA_Freq)] / TARA_Freq$Count * 100
TARA_Freq[, 11:ncol(TARA_Freq)] <- round(TARA_Freq[, 11:ncol(TARA_Freq)], 3)
TARA_Freq <- TARA_Freq[, !apply(TARA_Freq == 0, 2, all)]

### Differential Abundance analysis
TARA_Freq_DA <- TARA_Freq[, !(names(TARA_Freq) %in% c("Batch", "age (yrs)", "Count", "P1", "P2", "P3", "P4", "specific killing"))]
TARA_Freq_DA$Timepoint <- as.factor(TARA_Freq_DA$Timepoint)
TARA_Freq_DA <- droplevels(TARA_Freq_DA)
covariate_columns <- c("PID", "gender", "HIV", "viral load", "Timepoint", "Treatment")

### Limma ####
# Initialize a list to store results
results_list <- list()

# Loop over each combination of Treatment and Timepoint
for (treatment in levels(TARA_Freq_DA$Treatment)) {
  for (timepoint in levels(TARA_Freq_DA$Timepoint)) {
    
    # Subset data for current treatment and timepoint
    current_data <- TARA_Freq_DA[TARA_Freq_DA$Treatment == treatment & TARA_Freq_DA$Timepoint == timepoint,]
    names(current_data)[names(current_data) == "viral load"] <- "viral_load"
    current_data <- na.omit(current_data)
    
    # Check if there are enough data points to proceed
    if(nrow(current_data) > 1){
      # Design matrix including intercept
      design <- model.matrix(~ 0+ HIV + gender + viral_load, data = current_data)

      # Voom transformation
      v <- voom(t(current_data[, !(names(current_data) %in% c("PID", "gender", "HIV", "viral_load", "Timepoint", "Treatment"))]), 
                design, plot = FALSE)
      
      # Fit the linear model
      fit <- lmFit(v, design)

      # Contrast for HIV positive vs. negative
      contrast.matrix <- makeContrasts(HIVpositive - HIVnegative, levels=design)
      
      # Apply contrasts
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      # Store results
      results_list[[paste(treatment, timepoint, sep="_")]] <- topTable(fit2, adjust="BH", number=1000)
      
    }
  }
}

## VOlcano Plots
results_list$HUT78_12$feature<- rownames(results_list$HUT78_12)
results_list$HUT78_12$abs_logFC<- abs(results_list$HUT78_12$logFC)



HUT78_12_volcano_plot <- ggplot(results_list$HUT78_12, aes(x = logFC, y = -log10(P.Value), text = paste("Cell Type:", feature, "<br>LogFC:", logFC, "<br>Adj. P-Val:", adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot Dendritic Cell Fold Change HEI Vs HEU", x = "Log2 Fold Change (HEI vs HEU)", y = "-Log10 P-value")

HUT78_12_volcano_plot

# Convert the ggplot object to a plotly object
HUT78_12_volcano_plot <- ggplotly(HUT78_12_volcano_plot, tooltip = "text")

# Print the interactive plot object to display it in the RMarkdown HTML output
HUT78_12_volcano_plot


#### Mixed Effects Models ####
library(tidyverse)

filter_top_effects <- function(data_frame) {
  positive_effects <- data_frame %>% 
    arrange(desc(Estimate)) %>% 
    head(35)
  
  negative_effects <- data_frame %>% 
    arrange(Estimate) %>% 
    head(35)
  
  bind_rows(positive_effects, negative_effects)
}
TARA_Freq_M <- TARA_Freq[, !(names(TARA_Freq) %in% c("Batch", "age (yrs)", "Count", "P1", "P2", "P3", "P4"))]

# Initialize an empty list to store intermediate dataframes
merge_treatment_killing <- function(main_df, treatment_df, treatment_name) {
  # Create a unique key by pasting PID and Timepoint
  main_df$key <- paste(main_df$PID, main_df$Timepoint)
  treatment_df$key <- paste(treatment_df$PID, treatment_df$Timepoint)
  
  # Rename the specific killing column in the treatment dataframe
  col_name <- paste0("specific_killing_", treatment_name)
  names(treatment_df)[names(treatment_df) == paste0("specific_killing_", treatment_name)] <- col_name
  
  # Merge based on the unique key
  enriched_df <- left_join(main_df, treatment_df[c("key", col_name)], by = "key")
  
  # Remove the key column to clean up
  enriched_df$key <- NULL
  
  return(enriched_df)
}
unique_treatments <- unique(TARA_Freq_M$Treatment)

# Applying the function for each treatment
TARA_Freq_M_enriched <- TARA_Freq_M

for(treatment in unique_treatments) {
  # Prepare the treatment-specific dataframe
  specific_treatment_df <- TARA_Freq_M %>%
    filter(Treatment == treatment) %>%
    select(PID, Timepoint, `specific killing`) %>%
    mutate(!!paste0("specific_killing_", treatment) := `specific killing`) %>%
    select(-`specific killing`) # We don't need the Treatment column anymore
  
  # Merge this treatment's killing values into the main dataframe
  TARA_Freq_M_enriched <- merge_treatment_killing(TARA_Freq_M_enriched, specific_treatment_df, treatment)
}

# Number of columns in the dataframe
num_cols <- ncol(TARA_Freq_M_enriched)

# Create a vector with the new column order
# First 8 columns + last 8 columns moved to 9th position + remaining columns
new_order <- c(1:8, (num_cols-7):num_cols, 9:(num_cols-8))

# Reorder the columns based on the new order
TARA_Freq_M_enriched <- TARA_Freq_M_enriched[, new_order]

### Subset untreated 
untreated_df <- TARA_Freq_M_enriched %>%
  filter(Treatment == "untreated")
filtered_df <- untreated_df %>%
  filter(!is.na(specific_killing_HUT78) & !is.na(specific_killing_K562))

fixed_part_of_formula <- "~ `viral load` + Timepoint  + gender + HIV + specific_killing_HUT78 + specific_killing_K562 + (1 | PID)"

## HIV Effect

# Initialize lists to store results
models_list <- list()
percentage_columns <- colnames(filtered_df[17:ncol(filtered_df)])

# Loop through each  subset, ensuring column names are correctly handled
for (subset_name in percentage_columns) {
  # Escape subset_name if it contains special characters
  escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
  
  # Construct the formula string with escaped column names
  formula_str <- paste(escaped_subset_name, fixed_part_of_formula)
  
  # Convert the string to a formula
  current_formula <- as.formula(formula_str)
  
  
  # Fit the model using the current formula
  model <- tryCatch({
    lmer(current_formula, data = filtered_df)
  }, error = function(e) {
    cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
    return(NULL)  # Return NULL if there was an error fitting the model
  })
  
  # Store the model if successfully fitted
  if (!is.null(model)) {
    models_list[[subset_name]] <- model
  }
}

# Initialize an empty data frame to store the results
results_df <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                                   Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)

# Loop through each model to extract information
for (subset_name in names(models_list)) {
  model <- models_list[[subset_name]]
  summary_model <- summary(model)
  df <- as.data.frame(summary_model$coefficients)
  df$Subset <- subset_name
  df$Effect <- rownames(df)
  results_df <- rbind(results_df, df)
}

results_df <- results_df %>% 
  rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
  select(Subset, Effect, Estimate, `Std.Error`, P.Value)

HIV_effects <- results_df %>% filter(Effect == "HIVpositive")
HIV_effects$Significance <- ifelse(HIV_effects$P.Value < 0.05, "*", "")
HIV_effects <- HIV_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

ggplot(filter_top_effects(HIV_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of HIV Status", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

Gender_effects <- results_df %>% filter(Effect == "gendermale")
Gender_effects$Significance <- ifelse(Gender_effects$P.Value < 0.05, "*", "")
Gender_effects <- Gender_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

Timepoint_effects <- results_df %>% filter(Effect == "TimepointEntry")
Timepoint_effects$Significance <- ifelse(Timepoint_effects$P.Value < 0.05, "*", "")
Timepoint_effects <- Timepoint_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))



HUT78 <- results_df %>% filter(Effect == "specific_killing_HUT78")
HUT78$Significance <- ifelse(HUT78$P.Value < 0.05, "*", "")
HUT78 <- HUT78 %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))


K562 <- results_df %>% filter(Effect == "specific_killing_K562")
K562$Significance <- ifelse(K562$P.Value < 0.05, "*", "")
K562 <- K562 %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

#########
ggplot(filter_top_effects(Gender_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Gender", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

ggplot(filter_top_effects(Timepoint_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Timepoint", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of HUT78", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

ggplot(filter_top_effects(K562), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of K562", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

#### Mixed Model for k562 and Hut78

HUT78_df <- TARA_Freq_M_enriched %>%
  filter(Treatment == "HUT78")
K562_df <- TARA_Freq_M_enriched %>%
  filter(Treatment == "K562")

h_fixed_part_of_formula <- "~ `viral load` + Timepoint  + gender + HIV + specific_killing_HUT78  + (1 | PID)"
k_fixed_part_of_formula <- "~ `viral load` + Timepoint  + gender + HIV + specific_killing_K562  + (1 | PID)"


# Initialize lists to store results
h_models_list <- list()
k_models_list <- list()

h_percentage_columns <- colnames(HUT78_df[17:ncol(HUT78_df)])
k_percentage_columns <- colnames(K562_df[17:ncol(K562_df)])

### Model for HUT78

# Loop through each  subset, ensuring column names are correctly handled
for (subset_name in h_percentage_columns) {
  # Escape subset_name if it contains special characters
  escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
  
  # Construct the formula string with escaped column names
  formula_str <- paste(escaped_subset_name, h_fixed_part_of_formula)
  
  # Convert the string to a formula
  current_formula <- as.formula(formula_str)
  
  
  # Fit the model using the current formula
  model <- tryCatch({
    lmer(current_formula, data = HUT78_df)
  }, error = function(e) {
    cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
    return(NULL)  # Return NULL if there was an error fitting the model
  })
  
  # Store the model if successfully fitted
  if (!is.null(model)) {
    h_models_list[[subset_name]] <- model
  }
}

# Initialize an empty data frame to store the results
h_results_df <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                         Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)

# Loop through each model to extract information
for (subset_name in names(h_models_list)) {
  model <- models_list[[subset_name]]
  summary_model <- summary(model)
  df <- as.data.frame(summary_model$coefficients)
  df$Subset <- subset_name
  df$Effect <- rownames(df)
  h_results_df <- rbind(h_results_df, df)
}

### Model for K562

### Model for HUT78

# Loop through each  subset, ensuring column names are correctly handled
for (subset_name in k_percentage_columns) {
  # Escape subset_name if it contains special characters
  escaped_subset_name <- paste0("`", gsub("/", "\\/", subset_name), "`")
  
  # Construct the formula string with escaped column names
  formula_str <- paste(escaped_subset_name, k_fixed_part_of_formula)
  
  # Convert the string to a formula
  current_formula <- as.formula(formula_str)
  
  
  # Fit the model using the current formula
  model <- tryCatch({
    lmer(current_formula, data = K562_df)
  }, error = function(e) {
    cat("Error in fitting model for subset:", subset_name, "\nError message:", e$message, "\n")
    return(NULL)  # Return NULL if there was an error fitting the model
  })
  
  # Store the model if successfully fitted
  if (!is.null(model)) {
    k_models_list[[subset_name]] <- model
  }
}

# Initialize an empty data frame to store the results
k_results_df <- data.frame(Subset = character(), Effect = character(), Estimate = numeric(), 
                           Std.Error = numeric(), P.Value = numeric(), stringsAsFactors = FALSE)

# Loop through each model to extract information
for (subset_name in names(k_models_list)) {
  model <- models_list[[subset_name]]
  summary_model <- summary(model)
  df <- as.data.frame(summary_model$coefficients)
  df$Subset <- subset_name
  df$Effect <- rownames(df)
  k_results_df <- rbind(k_results_df, df)
}

####################

k_results_df <- k_results_df %>% 
  rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
  select(Subset, Effect, Estimate, `Std.Error`, P.Value)

h_results_df <- h_results_df %>% 
  rename(Estimate = Estimate, Std.Error = `Std. Error`, P.Value = `Pr(>|t|)`) %>%
  select(Subset, Effect, Estimate, `Std.Error`, P.Value)


### K Results

HIV_effects <- k_results_df %>% filter(Effect == "HIVpositive")
HIV_effects$Significance <- ifelse(HIV_effects$P.Value < 0.05, "*", "")
HIV_effects <- HIV_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

Gender_effects <- k_results_df %>% filter(Effect == "gendermale")
Gender_effects$Significance <- ifelse(Gender_effects$P.Value < 0.05, "*", "")
Gender_effects <- Gender_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

Timepoint_effects <- k_results_df %>% filter(Effect == "TimepointEntry")
Timepoint_effects$Significance <- ifelse(Timepoint_effects$P.Value < 0.05, "*", "")
Timepoint_effects <- Timepoint_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

K562 <- k_results_df %>% filter(Effect == "specific_killing_K562")
K562$Significance <- ifelse(K562$P.Value < 0.05, "*", "")
K562 <- K562 %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

#########

ggplot(filter_top_effects(HIV_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of HIV Status", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill


ggplot(filter_top_effects(Gender_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Gender", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

ggplot(filter_top_effects(Timepoint_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Timepoint", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

ggplot(filter_top_effects(K562), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of K562", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

### H effect


HIV_effects <- h_results_df %>% filter(Effect == "HIVpositive")
HIV_effects$Significance <- ifelse(HIV_effects$P.Value < 0.05, "*", "")
HIV_effects <- HIV_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

Gender_effects <- h_results_df %>% filter(Effect == "gendermale")
Gender_effects$Significance <- ifelse(Gender_effects$P.Value < 0.05, "*", "")
Gender_effects <- Gender_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))

Timepoint_effects <- h_results_df %>% filter(Effect == "TimepointEntry")
Timepoint_effects$Significance <- ifelse(Timepoint_effects$P.Value < 0.05, "*", "")
Timepoint_effects <- Timepoint_effects %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))



HUT78 <- h_results_df %>% filter(Effect == "specific_killing_HUT78")
HUT78$Significance <- ifelse(HUT78$P.Value < 0.05, "*", "")
HUT78 <- HUT78 %>%
  mutate(Color = ifelse(Estimate > 0, "lightblue", "lightgreen"),
         Significance = ifelse(P.Value < 0.05, "*", ""))
#########

ggplot(filter_top_effects(HIV_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of HIV Status", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill


ggplot(filter_top_effects(Gender_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Gender", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill

ggplot(filter_top_effects(Timepoint_effects), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of Timepoint", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill



ggplot(filter_top_effects(HUT78), aes(x = reorder(Subset, Estimate), y = Estimate, fill = Color)) +
  geom_col() +
  geom_errorbar(aes(ymin = Estimate - `Std.Error`, ymax = Estimate + `Std.Error`), width = 0.4) +
  geom_text(aes(label = Significance), vjust = -0.5, colour = "red") +
  scale_fill_identity() +
  coord_flip() +
  labs(title = "Effect of HUT78", x = "Subsets", y = "Effect Size") +
  theme_minimal() +
  guides(fill=FALSE) # Hide the legend for fill


### Volcano Effect of IL15


### Limma


# Initialize a list to store results
results_list <- list()

# Define pairs of treatments to compare
treatment_pairs <- list(
  "CEM" = "CEM+IL15",
  "HUT78" = "HUT78+IL15",
  "K562" = "K562+IL15"
)
# Process each treatment pair
for (baseline in names(treatment_pairs)) {
  treated <- treatment_pairs[[baseline]]
  
  # Filter data for only the current pair of treatments
  pair_data <- TARA_Freq_DA[TARA_Freq_DA$Treatment %in% c(baseline, treated), ]
  
  # Create a factor to distinguish between baseline and treated groups
  pair_data$group <- factor(ifelse(pair_data$Treatment == baseline, "baseline", "treated"))
  names(pair_data)[names(pair_data) == "viral load"] <- "viral_load"
  # Check if there are enough data points to proceed
  if(nrow(pair_data) > 1) {
    # Design matrix including Timepoint as a factor and HIV status
    design <- model.matrix(~0+ group + Timepoint + HIV  +viral_load, data = pair_data)
    
    # Voom transformation
    v <- voom(t(pair_data[, !(names(pair_data) %in% c("PID", "HIV", "gender", "Timepoint","viral_load", "Treatment", "group"))]), 
              design, plot = FALSE)
    
    # Fit the linear model
    corfit <- duplicateCorrelation(v, design, block = pair_data$PID)
    fit <- lmFit(v, design, block = pair_data$PID, correlation = corfit$consensus)
    
    # Specify the correct contrast based on the baseline and treated groups
    contrast.matrix <- makeContrasts(grouptreated - groupbaseline, levels = design)
    
    # Apply contrasts
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # Store results, specifying treatment pair
    results_key <- paste(baseline, "vs", treated, sep="_")
    results_list[[results_key]] <- topTable(fit2, adjust="BH", number=1000)
  }
}


## VOlcano Plots
results_list$`CEM_vs_CEM+IL15`$feature<- rownames(results_list$`CEM_vs_CEM+IL15`)
results_list$`CEM_vs_CEM+IL15`$abs_logFC<- abs(results_list$`CEM_vs_CEM+IL15`$logFC)

results_list$`HUT78_vs_HUT78+IL15`$feature<- rownames(results_list$`HUT78_vs_HUT78+IL15`)
results_list$`HUT78_vs_HUT78+IL15`$abs_logFC<- abs(results_list$`HUT78_vs_HUT78+IL15`$logFC)

results_list$`K562_vs_K562+IL15`$feature<- rownames(results_list$`K562_vs_K562+IL15`)
results_list$`K562_vs_K562+IL15`$abs_logFC<- abs(results_list$`K562_vs_K562+IL15`$logFC)


CEM_volcano_plot <- ggplot(results_list$`CEM_vs_CEM+IL15`, aes(x = logFC, y = -log10(P.Value), text = paste("Cell Type:", feature, "<br>LogFC:", logFC, "<br>Adj. P-Val:", adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  theme_minimal() +
  labs(title = "Volcano Plot NK Cell Fold Change  CEM+IL15 vs CEM", x = "Log2 Fold Change (CEM+IL15 vs CEM)", y = "-Log10 P-value")

CEM_volcano_plot <- ggplotly(CEM_volcano_plot, tooltip = "text")

CEM_volcano_plot
