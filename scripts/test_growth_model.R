# ==============================================================================
# Test Script: Growth Model with Quadratic Term
# Purpose: Compare linear-only model with linear+quadratic growth model
# ==============================================================================

# Load required packages
library(tidyverse)
library(afex)
library(lme4)
library(lmerTest)
library(emmeans)

cat("\n===== GROWTH MODEL TEST SCRIPT =====\n\n")

# ------------------------------------------------------------------------------
# Step 1: Load and preprocess data (exactly as in article.Rmd)
# ------------------------------------------------------------------------------

cat("Loading data...\n")
final_data <- read_csv("data/anonymized_final_data.csv", show_col_types = FALSE)

cat("Number of participants:", length(unique(final_data$playerId)), "\n")

# Apply same transformations as article.Rmd
final_data <- final_data %>%
  mutate(
    ret_pct_na = ifelse(investment == 0, NA, return / (3 * investment)),
    roundPayoff = 3 * investment - return,
    opponent.f = factor(gameOpponent, levels = c("AI_HMM", "AI_HMM_vol")),
    investorState.f = factor(investorState, levels = c("unhappy", "neutral", "happy")),
    d_level = as.factor(d_level),
    roundNum = as.numeric(as.character(roundNum)),
    inv_scaled = as.vector(scale(investment))
  ) %>%
  dplyr::select(-c("gameOpponent", "investorState"))

# Filter to rounds 1-24 (excluding endgame round 25)
data24 <- final_data %>% filter(roundNum < 25)

cat("Observations (rounds 1-24):", nrow(data24), "\n")
cat("NA values in ret_pct_na:", sum(is.na(data24$ret_pct_na)), "\n\n")

# ------------------------------------------------------------------------------
# Step 2: Create centered and polynomial terms for roundNum
# ------------------------------------------------------------------------------

# Center roundNum to reduce multicollinearity
data24$roundNum_c <- data24$roundNum - mean(data24$roundNum)
data24$roundNum_c2 <- data24$roundNum_c^2

cat("roundNum mean (for centering):", mean(data24$roundNum), "\n")
cat("roundNum_c range:", min(data24$roundNum_c), "to", max(data24$roundNum_c), "\n\n")

# ------------------------------------------------------------------------------
# Step 3: Fit Current Model (Linear roundNum only)
# ------------------------------------------------------------------------------

cat("===== MODEL 1: LINEAR ROUND ONLY (Current Model) =====\n\n")

mod_linear <- afex::mixed(
  ret_pct_na ~ opponent.f * inv_scaled * d_level * volatile_first * roundNum +
    (1 + opponent.f + inv_scaled | playerId),
  data = data24,
  REML = TRUE,
  method = "KR"
)

cat("Linear Model ANOVA:\n")
print(anova(mod_linear))

# Extract key results
linear_anova <- anova(mod_linear)
cat("\n--- Key Effects (Linear Model) ---\n")
cat("D-level main effect: F =", round(linear_anova$`F`[which(rownames(linear_anova) == "d_level")], 3),
    ", p =", round(linear_anova$`Pr(>F)`[which(rownames(linear_anova) == "d_level")], 4), "\n")
cat("D-level x Round: F =", round(linear_anova$`F`[which(rownames(linear_anova) == "d_level:roundNum")], 3),
    ", p =", round(linear_anova$`Pr(>F)`[which(rownames(linear_anova) == "d_level:roundNum")], 4), "\n")
cat("D-level x Opponent x Investment: F =", round(linear_anova$`F`[which(rownames(linear_anova) == "opponent.f:inv_scaled:d_level")], 3),
    ", p =", round(linear_anova$`Pr(>F)`[which(rownames(linear_anova) == "opponent.f:inv_scaled:d_level")], 4), "\n")

# ------------------------------------------------------------------------------
# Step 4: Fit Growth Model with Quadratic Term
# ------------------------------------------------------------------------------

cat("\n\n===== MODEL 2: GROWTH MODEL (Linear + Quadratic) =====\n\n")

# Using centered roundNum to reduce collinearity
# Note: Full 5-way interactions with quadratic would be very complex
# Approach A: Add quadratic only for key interactions (d_level, opponent)

mod_growth <- afex::mixed(
  ret_pct_na ~ opponent.f * inv_scaled * d_level * volatile_first * roundNum_c +
    roundNum_c2 * d_level * opponent.f +  # Quadratic interactions with key factors
    (1 + opponent.f + inv_scaled | playerId),
  data = data24,
  REML = TRUE,
  method = "KR"
)

cat("Growth Model ANOVA:\n")
print(anova(mod_growth))

# Extract key results
growth_anova <- anova(mod_growth)
cat("\n--- Key Effects (Growth Model) ---\n")

# Check if roundNum_c2 is significant
if("roundNum_c2" %in% rownames(growth_anova)) {
  cat("Quadratic term (roundNum_c2): F =", round(growth_anova$`F`[which(rownames(growth_anova) == "roundNum_c2")], 3),
      ", p =", round(growth_anova$`Pr(>F)`[which(rownames(growth_anova) == "roundNum_c2")], 4), "\n")
}

if("d_level:roundNum_c2" %in% rownames(growth_anova)) {
  cat("D-level x Quadratic: F =", round(growth_anova$`F`[which(rownames(growth_anova) == "d_level:roundNum_c2")], 3),
      ", p =", round(growth_anova$`Pr(>F)`[which(rownames(growth_anova) == "d_level:roundNum_c2")], 4), "\n")
}

if("opponent.f:roundNum_c2" %in% rownames(growth_anova)) {
  cat("Opponent x Quadratic: F =", round(growth_anova$`F`[which(rownames(growth_anova) == "opponent.f:roundNum_c2")], 3),
      ", p =", round(growth_anova$`Pr(>F)`[which(rownames(growth_anova) == "opponent.f:roundNum_c2")], 4), "\n")
}

# Core 3-way interaction
if("opponent.f:inv_scaled:d_level" %in% rownames(growth_anova)) {
  cat("D-level x Opponent x Investment (CORE): F =", round(growth_anova$`F`[which(rownames(growth_anova) == "opponent.f:inv_scaled:d_level")], 3),
      ", p =", round(growth_anova$`Pr(>F)`[which(rownames(growth_anova) == "opponent.f:inv_scaled:d_level")], 4), "\n")
}

# ------------------------------------------------------------------------------
# Step 5: Model Comparison
# ------------------------------------------------------------------------------

cat("\n\n===== MODEL COMPARISON =====\n\n")

# Get underlying lmer models for comparison
lmer_linear <- mod_linear$full_model
lmer_growth <- mod_growth$full_model

# AIC/BIC comparison
cat("AIC Comparison:\n")
cat("  Linear model:", round(AIC(lmer_linear), 2), "\n")
cat("  Growth model:", round(AIC(lmer_growth), 2), "\n")
cat("  Difference (Linear - Growth):", round(AIC(lmer_linear) - AIC(lmer_growth), 2), "\n")
cat("  (Negative = Linear is better; Positive = Growth is better)\n\n")

cat("BIC Comparison:\n")
cat("  Linear model:", round(BIC(lmer_linear), 2), "\n")
cat("  Growth model:", round(BIC(lmer_growth), 2), "\n")
cat("  Difference (Linear - Growth):", round(BIC(lmer_linear) - BIC(lmer_growth), 2), "\n\n")

# ------------------------------------------------------------------------------
# Step 6: Slopes Analysis (if quadratic is significant)
# ------------------------------------------------------------------------------

cat("\n===== ROUND SLOPES BY D-LEVEL =====\n\n")

# Linear slopes from current model
linear_slopes <- emtrends(mod_linear, ~ d_level, var = "roundNum")
cat("Linear Model - Round slopes by D-level:\n")
print(summary(linear_slopes))

# Compare slopes
cat("\nSlope difference (High-D vs Low-D):\n")
print(pairs(linear_slopes))

# If growth model, also show quadratic trends
cat("\nGrowth Model - Linear round slopes by D-level:\n")
growth_linear_slopes <- emtrends(mod_growth, ~ d_level, var = "roundNum_c")
print(summary(growth_linear_slopes))

# ------------------------------------------------------------------------------
# Step 7: Summary
# ------------------------------------------------------------------------------

cat("\n\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")

cat("1. Does the quadratic term improve model fit?\n")
aic_diff <- AIC(lmer_linear) - AIC(lmer_growth)
if(aic_diff > 2) {
  cat("   YES - AIC favors growth model by", round(aic_diff, 1), "points\n")
} else if(aic_diff < -2) {
  cat("   NO - AIC favors linear model by", round(-aic_diff, 1), "points\n")
} else {
  cat("   NEGLIGIBLE difference (within 2 AIC points)\n")
}

cat("\n2. Is the quadratic term (roundNum^2) significant?\n")
if("roundNum_c2" %in% rownames(growth_anova)) {
  p_quad <- growth_anova$`Pr(>F)`[which(rownames(growth_anova) == "roundNum_c2")]
  if(p_quad < 0.05) {
    cat("   YES - p =", round(p_quad, 4), "\n")
  } else {
    cat("   NO - p =", round(p_quad, 4), "\n")
  }
}

cat("\n3. Does D-level interact with the quadratic term?\n")
if("d_level:roundNum_c2" %in% rownames(growth_anova)) {
  p_dquad <- growth_anova$`Pr(>F)`[which(rownames(growth_anova) == "d_level:roundNum_c2")]
  if(p_dquad < 0.05) {
    cat("   YES - p =", round(p_dquad, 4), "(differential trajectories)\n")
  } else {
    cat("   NO - p =", round(p_dquad, 4), "(similar curvature)\n")
  }
}

cat("\n4. Does the core finding (D x Opponent x Investment) remain significant?\n")
if("opponent.f:inv_scaled:d_level" %in% rownames(growth_anova)) {
  p_core <- growth_anova$`Pr(>F)`[which(rownames(growth_anova) == "opponent.f:inv_scaled:d_level")]
  if(p_core < 0.05) {
    cat("   YES - p =", round(p_core, 4), "\n")
  } else {
    cat("   NO - p =", round(p_core, 4), "(caution: core finding may change)\n")
  }
}

cat("\n========================================\n")
cat("Script complete. Review output above for analysis decisions.\n")
cat("========================================\n")
