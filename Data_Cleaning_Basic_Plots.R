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
setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV")
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


data_cleaned$`Specific Killing`<- as.numeric(data_cleaned$`Specific Killing`)
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

######### Plot Specific Killing for each Group #############

### TARA
tara_Freq_plot <- drop_na(tara_Freq, `Specific Killing`)

# Filter out rows where Treatment is "untreated" or "cem" and HIV is "HUU"
tara_Freq_plot <- tara_Freq_plot %>%
  filter(!Treatment %in% c("untreated", "IL15"),HIV != "HUU")

tara_Freq_plot_filtered <- tara_Freq_plot %>%
  mutate(HIV = recode(HIV, "positive" = "HEI", "negative" = "HEU")) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("Entry", "12")))

# Plot with paired comparisons within each facet (HIV and Treatment)
ggplot(tara_Freq_plot_filtered, aes(x = Timepoint, y = `Specific Killing`, group = PID, color = HIV)) +
  geom_line(aes(linetype = HIV), size = 1) +  # Set line size
  geom_point(size = 2) +  # Set point size
  facet_grid(HIV ~ Treatment, scales = "free_y") +  # Splitting by both HIV and Treatment
  theme_minimal(base_size = 15) +  # Increase base font size for readability
  labs(
    title = "Specific Killing Across Treatments, Timepoints, and HIV Status (TARA)",
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

####### Paired Correlation for Specific Killing
hut78_data_TARA <- tara_Freq_plot_filtered %>%
  filter(Treatment == "HUT78")
k562_data_TARA <- tara_Freq_plot_filtered %>%
  filter(Treatment == "K562")

paired_hut78_data <- hut78_data_TARA %>%
  group_by(PID) %>%
  filter(all(c("Entry", "12") %in% Timepoint)) %>%
  ungroup()

# Correlation between Specific Killing at Entry and 12 months
cor_test_result <- cor.test(
  paired_hut78_data$`Specific Killing`[paired_hut78_data$Timepoint == "Entry"],
  paired_hut78_data$`Specific Killing`[paired_hut78_data$Timepoint == "12"],
  method = "pearson"  # Change to "spearman" if non-parametric
)

cor_test_result

setwd("C:/Users/axi313/Documents/NK_Killing_Differential_Abundance_Analysis_HIV/FasL_Plots")

# Linear regression model to assess the relationship
model <- lmer(`Specific Killing` ~ HIV + `viral load` + `CD56dimCD16+/FasL` + Timepoint + (1 | PID), data = hut78_data_TARA)
summary(model)

ggplot(hut78_data_TARA, aes(x = `CD56dimCD16+/FasL`, y = `Specific Killing`, group = Timepoint)) +
  geom_line(color = "grey", size = 1, linetype = "dashed") +  # Dashed grey lines to connect points for each timepoint
  geom_point(aes(color = HIV), size = 4, alpha = 0.8) +  # Color points based on HIV status
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", size = 1.2) +  # Black regression line with confidence interval
  labs(
    title = "HUT78 Specific Killing vs CD56dimCD16+/FasL",
    subtitle = "Faceted by Timepoint (Entry and 12 months)",
    x = "CD56dimCD16+/FasL",
    y = "Specific Killing (%)"
  ) +
  scale_color_manual(values = c("HEI" = "#D55E00", "HEU" = "#0072B2")) +  # Polished custom colors for HIV status
  facet_wrap(~ Timepoint) +  # Facet by Timepoint (Entry and 12 months)
  theme_minimal(base_size = 16) +  # Clean theme with larger base font size
  theme(
    legend.position = "top",  # Move legend to the top
    legend.title = element_blank(),  # Remove legend title for simplicity
    panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_rect(fill = "#f7f7f7", color = NA),  # Light grey background for better contrast
    plot.title = element_text(face = "bold", hjust = 0.5),  # Bold, centered title
    plot.subtitle = element_text(hjust = 0.5),  # Center the subtitle
    axis.title = element_text(face = "bold"),  # Bold axis titles
    strip.background = element_rect(fill = "#e0e0e0", color = NA),  # Light grey facet background
    strip.text = element_text(face = "bold")  # Bold facet labels
  )


ggsave("HUT_Specific_Killing_vs_CD56dimCD16_FasL_Timepoint.png", width = 8, height = 6,bg = 'white')

ggplot(k562_data_TARA, aes(x = `CD56dimCD16+/FasL`, y = `Specific Killing`, group = Timepoint)) +
  geom_line(color = "grey", size = 1, linetype = "dashed") +  # Dashed grey lines to connect points for each timepoint
  geom_point(aes(color = HIV), size = 4, alpha = 0.8) +  # Color points based on HIV status
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", size = 1.2) +  # Black regression line with confidence interval
  labs(
    title = "K562 Specific Killing vs CD56dimCD16+/FasL",
    subtitle = "Faceted by Timepoint (Entry and 12 months)",
    x = "CD56dimCD16+/FasL",
    y = "Specific Killing (%)"
  ) +
  scale_color_manual(values = c("HEI" = "#D55E00", "HEU" = "#0072B2")) +  # Polished custom colors for HIV status
  facet_wrap(~ Timepoint) +  # Facet by Timepoint (Entry and 12 months)
  theme_minimal(base_size = 16) +  # Clean theme with larger base font size
  theme(
    legend.position = "top",  # Move legend to the top
    legend.title = element_blank(),  # Remove legend title for simplicity
    panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_rect(fill = "#f7f7f7", color = NA),  # Light grey background for better contrast
    plot.title = element_text(face = "bold", hjust = 0.5),  # Bold, centered title
    plot.subtitle = element_text(hjust = 0.5),  # Center the subtitle
    axis.title = element_text(face = "bold"),  # Bold axis titles
    strip.background = element_rect(fill = "#e0e0e0", color = NA),  # Light grey facet background
    strip.text = element_text(face = "bold")  # Bold facet labels
  )

ggsave("K562_Specific_Killing_vs_CD56dimCD16_FasL_Timepoint.png", width = 8, height = 6,bg = 'white')

