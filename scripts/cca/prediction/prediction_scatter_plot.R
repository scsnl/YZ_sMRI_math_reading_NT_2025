############################################################
# Plot Predicted Brain vs Behavioral Scores
# Author: Yuan Zhang
# Date: 2025-07-25
#
# Description:
# This script visualizes the correlation between predicted brain 
# CCA scores and predicted behavioral CCA scores for both 
# math- and reading-related CCA modes. 
# Scatter plots with regression lines are generated and saved 
# as EPS files.
############################################################

# ------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------

# Set working directory to the project root
setwd("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub")

# Clear the environment
rm(list = ls())

# Load required package for plotting
library(ggplot2)

# Define paths to the prediction result CSV files
prediction_math    <- "results/cca/prediction/prediction_math.csv"
prediction_reading <- "results/cca/prediction/prediction_reading.csv"

# ------------------------------------------------------------------
# 2. Plot for Math Prediction
# ------------------------------------------------------------------

# Load predicted brain and behavior scores for math
df <- read.csv(prediction_math)

# Define output file for math plot (EPS format)
outputf <- 'results/cca/prediction/prediction_math.eps'

# Create scatter plot with linear regression line
ggplot(df, aes(x = predicted_brain_score, y = predicted_beh_score)) + 
  geom_point(size = 3, alpha = 0.5) +          # Scatter points
  geom_smooth(method = lm, se = TRUE) +        # Regression line with confidence interval
  theme_classic()                              # Clean, classic theme

# Save the plot
ggsave(outputf, units = "in", width = 5, height = 5, device = cairo_ps)

# ------------------------------------------------------------------
# 3. Plot for Reading Prediction
# ------------------------------------------------------------------

# Load predicted brain and behavior scores for reading
df <- read.csv(prediction_reading)

# Define output file for reading plot (EPS format)
outputf <- 'results/cca/prediction/prediction_reading.eps'

# Create scatter plot with linear regression line
ggplot(df, aes(x = predicted_brain_score, y = predicted_beh_score)) + 
  geom_point(size = 3, alpha = 0.5) +          # Scatter points
  geom_smooth(method = lm, se = TRUE) +        # Regression line with confidence interval
  theme_classic()                              # Clean, classic theme

# Save the plot
ggsave(outputf, units = "in", width = 5, height = 5, device = cairo_ps)
