# Libraries
library(tidyverse)  # Data manipulation and visualization
library(lme4)       # Linear Mixed-Effects Models (LMM)
library(lmerTest)   # Provides p-values for LMM
library(emmeans)    # Post-hoc comparisons and marginal means

# --- CONFIGURATION ---
# Set the path to your file
file_path <- "C:/Program Files/DigiMove/DigiMove/DataAnalysisProject/Analysis_Results/hb_summary_final_semaxone.csv"

# Choose your variable of interest (e.g., "Mean_HbR_uM", "Mean_HbO_uM", "Mean_HbTot_uM")
target_variable <- "Mean_HbR_uM"

# ---------------------

# Load data using read_csv2 (handles ';' separator and ',' decimals)
df <- read_csv2(file_path)

# Data processing pipeline
df_clean <- df %>%
  # Ensure the first column is named "Subject" (fixes potential BOM encoding issues)
  rename(Subject = 1) %>% 
  # Filter for the intensities of interest
  filter(Intensity %in% c("30%", "50%")) %>%
  # Clean and transform
  mutate(
    # Create the generic analysis column dynamically
    NIRS_Value = .data[[target_variable]],
    # Convert to factors for statistical modeling
    Subject = as.factor(Subject),
    Contraction_Mode = as.factor(Contraction_Mode),
    # Remove '%' and convert Intensity to factor
    Intensity = as.factor(gsub("%", "", Intensity))
  ) %>%
  # Remove rows with missing values in our target variable
  drop_na(NIRS_Value)

# Check the first few rows
head(df_clean)

# --------------------
#Computing some descriptive statistics
stats_summary <- df_clean %>%
  group_by(Contraction_Mode, Intensity) %>%
  summarise(
    Mean = mean(NIRS_Value),
    SD = sd(NIRS_Value),
    N = n(),
    .groups = 'drop'
  )

print("Descriptive Statistics")
print(stats_summary)

# ---------------------
#Visualization 
ggplot(df_clean, aes(x = Intensity, y = NIRS_Value, fill = Contraction_Mode)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  # Add individual points to see the distribution
  geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = paste("NIRS Analysis:", target_variable),
    subtitle = "Comparison by Contraction Mode and Intensity",
    x = "Intensity (% MVC)",
    y = paste(target_variable, "(µM)"),
    fill = "Mode"
  ) +
  scale_fill_manual(values = c("Concentric" = "#4e79a7", "Eccentric" = "#f28e2b"))

# --------------------
#Statistical modeling with a Liner Mixed Model
#fix effect : Contraction_Mode, Intensity and their interaction
#randoom effect : Subject (to account for repeated measures within subjects)
# Fit the model
model <- lmer(NIRS_Value ~ Contraction_Mode * Intensity + (1|Subject), data = df_clean)

# Display ANOVA table (Type III Wald F-tests)
print("--- ANOVA Results ---")
print(anova(model))

# -------------------
# Post-hoc pairwise comparisons with Bonferroni adjustment
print("--- Pairwise Comparisons (Bonferroni Adjustment) ---")

# Compare Concentric vs Eccentric for each Intensity level
posthoc_results <- emmeans(model, pairwise ~ Contraction_Mode | Intensity, adjust = "bonferroni")

# Display the contrasts
print(posthoc_results$contrasts)