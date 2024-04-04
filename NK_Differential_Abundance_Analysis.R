library(readxl)
library(dplyr)
library(lme4)
library(ggplot2)
library(ggeffects)
library(lmerTest)
library(tidyr)
library(reshape2)

########### Read in Data and Viral Data ############

data <- read_excel("NK Function Metadata.xlsx")
v_data <- read.csv("Viral_Titres.csv")

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


######### Split Data based on Killing, MFI and Frequency ############

col_names <- names(data_cleaned)
col_names_excluding_killing <- col_names[-((ncol(data_cleaned)-2):ncol(data_cleaned))]

# Identify indices for the shared first 9 columns and the last 3 columns
shared_cols <- 1:10
killing_additional_cols <- (ncol(data_cleaned)-2):ncol(data_cleaned) # Last 3 columns

# Identify columns for MFI and Freq excluding the last three columns
mfi_cols <- grep("Median", col_names_excluding_killing)
freq_cols <- grep("\\(\\%\\)", col_names_excluding_killing, perl = TRUE)

# Map these columns back to the original dataframe's column indices
mfi_indices <- match(col_names_excluding_killing[mfi_cols], col_names)
freq_indices <- match(col_names_excluding_killing[freq_cols], col_names)

# Creating the dataframes
Killing <- data_cleaned[, c(shared_cols, (ncol(data_cleaned)-2):ncol(data_cleaned))]
MFI <- data_cleaned[, unique(c(shared_cols, mfi_indices))]
Freq <- data_cleaned[, unique(c(shared_cols, freq_indices))]


### Clean Frequency Dataset prior to Splitting, and convert % to raw counts ########
names(Freq) <- gsub(" \\| Freq\\. of Parent \\(%\\)", "", names(Freq))
rename_map <- c(
  "Lymphocytes" = "P1",
  "Lymphocytes/Single Cells" = "P2",
  "Lymphocytes/Single Cells/Live" = "P3",
  "Lymphocytes/Single Cells/Live/CD3-Label-" = "P4"
)

# Loop through the rename map and assign new names to Freq dataframe
for(original_name in names(rename_map)) {
  # Check if the original_name exists in the column names of Freq
  if(original_name %in% names(Freq)) {
    # Get the index of the column with the original name
    col_index <- which(names(Freq) == original_name)
    # Assign the new name to the column
    names(Freq)[col_index] <- rename_map[original_name]
  }
}

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

# Impute NA Values
# Columns to impute (example; adjust as needed)
columns_to_impute <- c("Total NK", single_slash_cols, double_slash_cols)

# Impute NAs based on the median for the HIV status group
Freq <- Freq %>%
  mutate(across(all_of(columns_to_impute), ~if_else(is.na(.),
                                                    ave(., HIV, FUN = function(x) median(x, na.rm = TRUE)),
                                                    .),
                .names = "{.col}"))


# Step 2: Calculate values for direct subsets of "Total NK"
for(col in single_slash_cols) {
  Freq[[col]] <- Freq[["Total NK"]] * (Freq[[col]] / 100)
}


# Srep 3 Double / cols
for(col in double_slash_cols) {
  # Extract the parent subset name from the column name
  parts <- strsplit(col, "/")[[1]]
  parent_subset_name <- paste(parts[1:length(parts)-1], collapse="/")
  
  # Ensure parent_subset_name is in single_slash_cols or double_slash_cols
  if(parent_subset_name %in% names(Freq)) {
    parent_value <- Freq[[parent_subset_name]]
    Freq[[col]] <- parent_value * (Freq[[col]] / 100)
  }
}

# Update the dataframe column names to reflect the calculated values
names(Freq) <- gsub("Total NK/", "", names(Freq))



#### Split Datasets for Analysis #####

TARA_Killing <- subset(Killing, Group == "TARA")
TARA_MFI <- subset(MFI, Group == "TARA")
TARA_Freq <- subset(Freq, Group == "TARA")
Florah_Killing <- subset(Killing, Group == "Florah")
Florah_MFI <- subset(MFI, Group == "Florah")
Florah_Freq <- subset(Freq, Group == "Florah")

###### NK KILLING Analysis #############

### Cleaning Dataframe ###

# Removing the 'Group' column from both dataframes
TARA_Killing$Group <- NULL
TARA_Killing <- TARA_Killing %>%
  filter(`specific killing` != "n/a")

# Round the values in the 'specific killing' column to two decimal points
TARA_Killing$`specific killing` <- as.numeric(TARA_Killing$`specific killing`)
TARA_Killing$`specific killing` <- round(TARA_Killing$`specific killing`, 2)

TARA_Killing <- TARA_Killing %>%
  filter(HIV != "HUU")

# Remove Unneeded Columns

TARA_Killing <- TARA_Killing[, -c(ncol(TARA_Killing)-2, ncol(TARA_Killing)-1)]

TARA_Killing <- TARA_Killing[, -which(names(TARA_Killing) == "age (yrs)")]

TARA_Killing <- droplevels(TARA_Killing)


# Standardize the 'viral load' variable because of high values

TARA_Killing$`viral load` <- scale(TARA_Killing$`viral load`)
str(TARA_Killing)

TARA_Killing$Timepoint <- as.factor(TARA_Killing$Timepoint)
TARA_Killing$Timepoint <- factor(TARA_Killing$Timepoint, levels = rev(levels(TARA_Killing$Timepoint)))

#### Fit the mixed effects model for each treatment condition ####

# Define the treatment levels
treatments <- levels(TARA_Killing$Treatment)

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

