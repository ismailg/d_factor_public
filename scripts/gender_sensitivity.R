# Gender Sensitivity Analysis
# Tests whether the key D × Opponent × Investment interaction holds when controlling for gender

library(tidyverse)
library(afex)

# Load data
final_data <- read_csv("data/anonymized_final_data.csv")
demographics <- read_csv("data/anonymized_demographics.csv")

# Merge gender into final_data
final_data <- final_data %>%
  left_join(demographics %>% select(playerId, Sex), by = "playerId")

# Check merge worked
cat("Gender distribution in merged data:\n")
print(table(final_data$Sex, useNA = "ifany"))

# Filter to rounds 1-24 (excluding endgame)
data24 <- final_data %>%
  filter(roundNum < 25) %>%
  filter(Sex %in% c("Male", "Female"))  # Exclude "CONSENT_REVOKED" / "Prefer not to say"

# Check sample size
n_participants <- length(unique(data24$playerId))
cat("\nParticipants with valid gender:", n_participants, "\n")

# Prepare variables (matching article.Rmd preprocessing)
data24 <- data24 %>%
  mutate(
    ret_pct_na = ifelse(investment == 0, NA, return / (3 * investment)),
    opponent.f = factor(ifelse(gameOpponent == "AI_HMM", "Human-like", "Responsive")),
    inv_scaled = scale(investment)[,1],
    d_level = factor(d_level),
    volatile_first = factor(volatile_first),
    Sex = factor(Sex)
  )

cat("\nGender × D-level crosstab:\n")
print(table(data24$Sex, data24$d_level) / 48)  # Divide by 48 rounds to get participant count

# Run model WITH gender as main effect (not interacted)
cat("\n\n=== MODEL WITH GENDER COVARIATE ===\n\n")

mod_with_gender <- afex::mixed(
  ret_pct_na ~ Sex + opponent.f*inv_scaled*d_level*volatile_first*roundNum +
    (1 + opponent.f + inv_scaled | playerId),
  data = data24,
  REML = TRUE,
  method = "KR"
)

cat("\nANOVA table:\n")
print(anova(mod_with_gender))

# Extract key interaction
anova_table <- anova(mod_with_gender)
d_opp_inv <- anova_table["d_level:opponent.f:inv_scaled", ]
cat("\n=== KEY RESULT ===\n")
cat("D × Opponent × Investment interaction:\n")
cat(sprintf("  F(%d, %.1f) = %.2f, p = %.4f\n",
            d_opp_inv$`num Df`, d_opp_inv$`den Df`, d_opp_inv$F, d_opp_inv$`Pr(>F)`))

if (d_opp_inv$`Pr(>F)` < 0.05) {
  cat("\n✓ KEY FINDING HOLDS: 3-way interaction remains significant when controlling for gender.\n")
} else {
  cat("\n✗ WARNING: 3-way interaction no longer significant with gender covariate.\n")
}

# Also check gender main effect
gender_effect <- anova_table["Sex", ]
cat("\nGender main effect:\n")
cat(sprintf("  F(%d, %.1f) = %.2f, p = %.4f\n",
            gender_effect$`num Df`, gender_effect$`den Df`, gender_effect$F, gender_effect$`Pr(>F)`))
