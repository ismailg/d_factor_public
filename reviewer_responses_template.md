# Point-by-Point Response to Reviewers

## D-Factor Study: "The Price of Dark Traits: Strategic Exploitation and Its Limitations in Repeated Interactions"

---

We thank both reviewers for their thoughtful and constructive feedback, and for their prompt review of our manuscript. Their comments have helped us strengthen the manuscript substantially. Below we provide point-by-point responses to each comment, with changes tracked in the revised manuscript.

---

# REVIEWER 1

We thank Reviewer 1 for their careful reading of the manuscript and their insightful suggestions regarding the growth model analysis, D-factor conceptualization, and manuscript organization. These comments have led to meaningful improvements in both our analytical approach and theoretical framing.

## Comment R1.1: Manuscript Organization

### R1.1a - Trust Game Explanation

**Reviewer Comment:**
> "The paper assumes familiarity with the trust game. This would be easily fixed by adding a few sentences explaining the structure (and the dominant response options)."

**Response:**
We have added a paragraph at the beginning of the Repeated Trust Game section that explains the game structure and the tension between self-interest and cooperation.

**Changes Made:**
- Added to Methods/Repeated Trust Game section (line 377): "The trust game is a two-player sequential game that captures the tension between self-interest and mutual benefit. An investor transfers part of their endowment to a trustee; this amount is multiplied (typically tripled), and the trustee decides how much to return. From a purely self-interested perspective, the trustee's dominant strategy is to keep everything, yet returning a fair share can sustain cooperation and yield higher long-term gains for both parties."

---

### R1.1b - Opponent Ratings Section

**Reviewer Comment:**
> "The whole section on opponent ratings could be reduced to a couple of sentences with the detailed results moved to the supplement. I am particularly unsure concerning the purpose of the play again ratings. Maybe the cooperative ratings could be regarded manipulation/validity check, so consider reporting this earlier (again in a single sentence)."

**Response:**
We have condensed the opponent ratings section in the main text to four summary sentences while moving the full ANOVA results and post-hoc comparisons to the Supplementary Materials. As the reviewer suggests, we now explicitly note that these ratings serve as a manipulation check: the significant main effect of opponent type on cooperativeness and trustworthiness ratings confirms that participants perceived the intended behavioral differences between the Human-like and Responsive HMMs. The main text now reports that High-D participants rated opponents as less cooperative, less trustworthy, and were less willing to play with them again, with significant main effects of D-level and opponent type. The figure showing the three rating dimensions has been retained in the main text as it provides an accessible visual summary. Full statistical details including all post-hoc pairwise comparisons are now available in the supplement.

**Changes Made:**
- Condensed Results/Opponent Ratings section to 4 summary sentences (lines 2130-2133)
- Added explicit framing of ratings as manipulation/validity check
- Added full statistical results to Supplementary Materials (lines 211-243), including ANOVA tables and post-hoc comparisons with verified statistics
- Retained figure in main text for visual summary

---

### R1.1c - Game Order Effect

**Reviewer Comment:**
> "The effect of game order is weird, and the authors should do more to explain this."

**Response:**
We have expanded the Discussion to explain why presentation order matters for High-D participants. The key mechanism is strategy transfer: when High-D individuals encounter the Human-like HMM first, they learn that exploitation is sustainable (the HMM does not rapidly punish non-reciprocation). This learned approach then carries over to their second game, creating a "head start" in exploitative behavior even when facing a different opponent type. This interpretation is consistent with research showing that people transfer learned opponent models between games (Guennouni & Speekenbrink, 2022). This finding is explicitly labeled as exploratory, given that we did not predict it a priori.

**Changes Made:**
- Expanded Discussion section (lines 2582-2591) with strategy transfer mechanism explanation
- Added citation: Guennouni & Speekenbrink (2022) on strategy transfer between games
- Labeled finding as exploratory

---

### R1.1d - Debriefing Variable as Covariate

**Reviewer Comment:**
> "It would be better to consider the debriefing variable as a covariate to see whether results change depending on whether respondents thought they played versus a human or an agent."

**Response:**
We addressed this concern through a sensitivity analysis that we believe provides a more stringent test than covariate adjustment. Rather than statistically controlling for belief while retaining all participants, we restricted the sample to only those participants who believed they were playing against a human or were unsure (N = 104), completely excluding participants who were certain they played a computer (N = 79, 43%). This approach directly tests whether our effects hold when removing potentially confounded observations entirely.

The key D × Opponent × Investment three-way interaction remained significant in this restricted sample (p = .007), demonstrating that our primary findings are robust to participants' beliefs about their opponents. This sensitivity analysis is reported in the Supplementary Materials, and we reference these results in the main text Limitations section (see also our response to R2.4).

**Changes Made:**
- Sensitivity analysis reported in Supplementary Materials (lines 336-428)
- Results referenced in Limitations section (line 2620)
- See also response to R2.4 for full details

---

## Comment R1.2: Growth Model Analysis

**Reviewer Comment:**
> "Instead of dividing the 25 rounds into three periods, better apply a growth model (excluding the final round) with linear and quadratic trends. Indeed, the authors already consider round number as predictor, so they simply should also include a quadratic trend (along with the interactions). The difference between the periods then becomes evident in different slopes for high/low D individuals (as reported at 257: significant interaction between D-factor and round number). Then the whole section on the three periods could be dropped, as this would be entirely redundant, and the figure should rather show individual rounds at the x axis."

**Response:**
We thank the reviewer for this methodologically important suggestion. We evaluated a growth model including both linear and quadratic trends for round number, with the quadratic term interacting with D-level and opponent type. A likelihood ratio test indicated that adding the quadratic terms significantly improved model fit (χ²(4) = 15.05, p = .005). However, critically, D-level did not interact with the quadratic term (F = 0.46, p = .50), indicating that both high-D and low-D participants showed similar curvature in their trajectories—they differ in their linear slopes (rate of change over rounds), not in the shape of change over time.

Given that the theoretically meaningful D × Quadratic interaction was non-significant—meaning the quadratic term does not reveal differential learning trajectories between groups—we retained the linear specification for our primary analyses. The improvement in fit from the quadratic term appears to capture overall curvature that is shared by both groups, rather than group-specific dynamics. As the reviewer correctly notes, the significant D-level × Round interaction already captures the key finding that high-D participants show steeper declines in return percentages over time (~0.17% per round for high-D vs. ~0.05% for low-D). We have clarified the rationale for our modeling approach in the Methods and provide full model comparison details including the likelihood ratio test in the Supplementary Materials.

We have also clarified that the period divisions (early/middle/late) are used for descriptive purposes in summarizing the data, while the formal statistical test of temporal dynamics uses round number as a continuous predictor in the linear mixed model.

**Changes Made:**

- Revised Statistical Analysis section (lines 803-805) to clarify that round number is modeled continuously, with note about growth model comparison
- Added Growth Model Comparison section to Supplementary Materials (lines 432-533) with full model fit statistics, quadratic term results, and justification for linear specification

---

## Comment R1.3: Theoretical Justification

**Reviewer Comment:**
> "More should be done in the theory explaining why we should expect diverging interaction patterns of low vs high D individuals depending on the interaction partner. 10.1017/jdm.2025.1 and 10.1037/xge0001232 might help in this regard."

**Response:**
We have added a new paragraph to the Introduction providing theoretical justification for the predicted D × Opponent interaction. Drawing on Hilbig et al. (2022), we explain that D is strongly associated with distrust-related beliefs (cynicism, expectations of selfish behavior in others) that serve to justify exploitative behavior. This framework predicts that high-D individuals will calibrate their behavior based on opponent feedback: exploiting predictable opponents who fail to punish non-reciprocation, while showing restraint when opponents respond adaptively. Low-D individuals, lacking these justificatory beliefs, should behave more consistently regardless of opponent type.

**Changes Made:**

- Added new paragraph to Introduction (line 68) explaining distrust-based justification mechanism
- Added citation: Hilbig, B. E., Moshagen, M., Thielmann, I., & Zettler, I. (2022). Making rights from wrongs: The crucial role of beliefs and justifications for the expression of aversive personality. Journal of Experimental Psychology: General, 151(11), 2730-2755.

---

## Comment R1.4: Strategic Exploitation Visibility

**Reviewer Comment:**
> "The strategic exploitation pattern observed for high D individuals with the 'human-like' agent simply cannot be observed with the responsive opponent, given that there are no sustained periods of similarly high or low offers. This should be made more explicit."

**Response:**
We agree this point needed clarification. We have added an explicit explanation in the Discussion noting that the exploitation pattern could not have emerged with the Responsive opponent by design: the Responsive HMM's rapid state-switching meant it did not generate sustained periods of consistently high investments that would allow gradual exploitation to be detected. In contrast, the Human-like HMM's characteristic 'sticky' mid-trust state maintained relatively stable investment levels, providing the necessary conditions for observing exploitation over time.

**Changes Made:**
- Added to Discussion (line 2578): "Importantly, this pattern could not have emerged with the Responsive opponent by design: because the Responsive HMM rapidly shifted between trust states based on immediate returns, it did not generate sustained periods of consistently high investments that would allow gradual exploitation to be detected. The Human-like HMM's characteristic 'sticky' mid-trust state, which maintained relatively stable investment levels despite fluctuations in returns, provided the necessary conditions for observing exploitation over time."

---

## Comment R1.5: D-factor Conceptualization

**Reviewer Comment:**
> "As limitation, the authors mention that they 'cannot determine which specific aspects of the D-factor (e.g., Machiavellianism, psychopathy, or narcissism)' drive the effects. Same in the introduction ('incorporating elements like Machiavellianism and Psychopathy'). This is misleading, as specific traits are not an aspect or element of D, but rather are an expression of D that typically comprise additional essentially non-aversive aspects (see e.g. 10.1016/j.copsyc.2025.102111)."

**Response:**
We thank the reviewer for this important clarification regarding D-factor conceptualization. The reviewer is correct that our original language was imprecise: D is the common core underlying dark traits, not a composite of them. Specific traits like Machiavellianism and Psychopathy are expressions of D that also include unique non-aversive aspects (e.g., strategic thinking in Machiavellianism, fearlessness in Psychopathy). We have revised all relevant passages to reflect this proper conceptualization.

**Changes Made:**
- Introduction (line 49): Changed "incorporating elements like Machiavellianism and Psychopathy" to "represents the shared variance among dark traits such as Machiavellianism and Psychopathy"; added citation to Moshagen et al. (2025)
- Introduction (line 53): Removed "its Machiavellian component" phrasing; now describes "strategic manipulation characteristic of high-D individuals"
- Limitations (lines 2631-2632): Reframed from "aspects of the D-factor" to clarify that we "measured only the common core (D)" and cannot determine whether patterns differ across "specific dark trait expressions"
- Added citation: Moshagen, M., Hilbig, B. E., & Zettler, I. (2025). Reconceptualizing ethically and socially aversive ("dark") personality traits. Current Opinion in Psychology, 66, 102111.

---

## Comment R1.6: Minor Points

### R1.6a - Broken References

**Reviewer Comment:**
> "fix reference at 39 and 423 '(Benjamin E. Hilbig et al. 2021b;'"

**Response:**
We thank the reviewer for catching this error. The issue was caused by inconsistent author name formatting in our bibliography file (some entries used "Benjamin E" while others used "Benjamin E." with a period). We have standardized all author name formats, which resolves the rendering issue.

**Changes Made:**
- Standardized author name formatting in bib/DFP2.bib

---

### R1.6b - Additional Citations

**Reviewer Comment:**
> "better cite 10.1017/jdm.2025.1 concerning the relation between D and single shot econ games; also consider 10.1177/09567976221116893 to make the point that D is a very strong contender for the prediction of pro-vs-antisociality and 10.1037/xge0001232 for the point that distrust is an important justification mechanism in high D individuals"

**Response:**
We have added all three recommended citations to the manuscript:

**Changes Made:**

- Added Hilbig & Thielmann (2025) to Introduction (line 51) supporting the link between D and behavior in single-shot economic games
- Added Hilbig et al. (2023) to Discussion (line 2620) demonstrating that D is the singular dispositional predictor of social preferences among 58 traits tested
- Added Hilbig et al. (2022) to Introduction (line 68) explaining the role of distrust-related beliefs in justifying exploitative behavior among high-D individuals

Citations added:

- Hilbig, B. E., & Thielmann, I. (2025). Toward a (more) parsimonious account of the link between 'dark' personality and social decision-making in economic games. Judgment and Decision Making, 20, e16.
- Hilbig, B. E., Thielmann, I., Zettler, I., & Moshagen, M. (2023). The dispositional essence of proactive social preferences: The dark core of personality vis-à-vis 58 traits. Psychological Science, 34(2), 201-220.
- Hilbig, B. E., Moshagen, M., Thielmann, I., & Zettler, I. (2022). Making rights from wrongs: The crucial role of beliefs and justifications for the expression of aversive personality. Journal of Experimental Psychology: General, 151(11), 2730-2755.

---

### R1.6c - Predictor Coding

**Reviewer Comment:**
> "please state how the binary predictors (D level etc) were coded"

**Response:**
We have added a statement in the Statistical Analysis section specifying that categorical predictors were effect-coded using sum-to-zero contrasts. This is the default behavior of the afex package we used for model estimation.

**Changes Made:**
- Added to Methods/Statistical Analysis: "Categorical predictors (D-factor group, Opponent type, and presentation Order) were effect-coded using sum-to-zero contrasts, such that main effects represent deviations from the grand mean."

---

### R1.6d - Gender Distribution

**Reviewer Comment:**
> "please state the gender distribution"

**Response:**
We have added explicit gender distribution information to the Participants section. The sample comprised 56% male and 44% female participants. We also report that gender distribution differed significantly between D-factor groups (High-D: 70% male; Low-D: 41% male; χ²(1) = 15.04, p < .001), which also addresses Reviewer 2's Comment 7 regarding gender balance.

**Changes Made:**
- Added to Participants section: "The sample comprised 56% male and 44% female participants... Gender distribution differed significantly between D-factor groups: the High-D group was 70% male, while the Low-D group was 41% male, χ²(1) = 15.04, p < .001."

---

### R1.6e - Session Timing

**Reviewer Comment:**
> "What was the lag between screening and the actual game? And clarify what is meant by 'multiple sessions'."

**Response:**
We have clarified the study timeline. The study comprised two sessions: (1) an initial screening survey where participants completed the D-factor measure, and (2) the main experimental task approximately one week later. The phrase "multiple sessions" in the original manuscript referred to our data collection schedule—we collected data over multiple recruitment sessions between October and November 2024—rather than to the participant experience. We have revised the wording to eliminate this ambiguity.

**Changes Made:**
- Clarified in Participants section (line 357): "These participants were then invited via Prolific to complete the main experiment approximately one week after the screening session (the study thus comprised two sessions: an initial screening survey and the main experimental task)."

---

### R1.6f - HMM Randomization

**Reviewer Comment:**
> "From the analyses it becomes clear that HMM types were randomly assigned to the study phases, but please make this explicit earlier"

**Response:**
We have made the randomization procedure explicit in the Design section. The text now clearly states: "Participants were randomly assigned to encounter either the Human-like or Responsive HMM first, with presentation order counterbalanced across participants."

**Changes Made:**
- Added explicit statement about random assignment and counterbalancing in Design section (line 370)

---

### R1.6g - "High D" Label Qualification

**Reviewer Comment:**
> "Consider to qualify the realization of low vs high D individuals in the present study based on the distribution observed in the large multinational sample reported on in 10.1073/pnas.2500830122. I am not that sure whether 'high D' is actually the proper label as this rather seems to be a moderate level."

**Response:**
We appreciate this important point about contextualizing our groups. We have added normative comparisons in the Participants section. Our High-D group (M = 2.95 on the 1–5 scale) scored approximately 1.2 SDs above the normative mean (M = 2.14, SD = 0.67), placing them around the 88th percentile; our Low-D group (M = 1.22) scored approximately 1.4 SDs below, around the 8th percentile. While we retain the "High-D" and "Low-D" labels for consistency with our extreme-groups design, we acknowledge these represent relative extremes within our sample rather than clinical or absolute categorizations.

**Changes Made:**
- Added normative context in Participants section (line 356): percentile rankings relative to population distribution

---

### R1.6h - Payoff Explanation Placement

**Reviewer Comment:**
> "delete the explanation of the lack of differences in the payoff at 243 and just dwell on this in the discussion"

**Response:**
We have streamlined the Results section to report only the statistical finding for total payoffs. The interpretation and explanation of why High-D participants did not achieve higher payoffs despite their exploitative strategy is now presented solely in the Discussion section (lines 2591-2592), where we explain that the adaptive nature of the HMM opponent reduced investments in response to lower returns from High-D participants.

**Changes Made:**
- Results section (line 1397): Now reports only the statistical test results without interpretive explanation
- Discussion (lines 2591-2592): Contains the full interpretation of the payoff finding

---

### R1.6i - P-value Formatting

**Reviewer Comment:**
> "266: replace 'p = 0.000' by p < 0.001"

**Response:**
We have corrected the p-value formatting throughout the manuscript. We created a helper function that formats p-values according to APA guidelines, displaying "p < .001" for values below 0.001 rather than "p = 0.000".

**Changes Made:**
- Added format_p() helper function to ensure consistent p-value formatting
- Updated all inline R code to use this function

---

# REVIEWER 2

We thank Reviewer 2 for their thoughtful engagement with our work and their helpful suggestions regarding the conceptual positioning of D-factor, the sensitivity analyses for participant beliefs, and the expansion of our discussion of learning dynamics. These comments have strengthened both our empirical presentation and theoretical interpretation.

## Comment R2.1: Strengths Acknowledged

**Reviewer Comment:**
> "The manuscript has a number of notable strengths. First, the use of computer-simulated Hidden Markov Model (HMM) agents to create responsive, human-like investors is innovative and allows for a controlled yet potentially ecologically valid experimental environment. Second, the focus on sustained interactions rather than one-shot encounters is a clear strength, as it more closely reflects real-world strategic behavior and allows for the examination of learning and adaptation over time. Third, the authors' attention to behavioral change across trials is commendable and aligns well with contemporary interest in dynamic decision-making processes rather than static trait–behavior associations."

**Response:**
We thank the reviewer for this positive assessment of our methodological approach. We appreciate the recognition that our HMM-based methodology and focus on dynamic interactions represents a meaningful contribution to the literature.

---

## Comment R2.2: D-factor vs. Agreeableness/Honesty-Humility

**Reviewer Comment:**
> "The manuscript centers on the D-factor, which remains a controversial construct. A persistent concern in the literature is whether D captures something meaningfully distinct from low agreeableness or from a combination of low agreeableness and low honesty–humility. It would strengthen the theoretical contribution of the manuscript if the authors more explicitly discussed how their findings might generalize to, or be interpreted through, these better-established personality dimensions. Addressing this issue would help situate the results within a broader personality framework and clarify what – if anything – is gained by invoking D specifically."

**Response:**
We thank the reviewer for raising this important point. We have added a new paragraph in the Discussion (Theoretical implications) that explicitly addresses the relationship between D and established personality dimensions. Drawing on Hilbig et al. (2021), we explain that while D and low Honesty-Humility share substantial variance (~66%), they remain functionally distinct: D better captures distrust-related beliefs and lack of empathy, whereas low Honesty-Humility better accounts for pretentiousness. We also explain how D more comprehensively captures utility-maximization accompanied by justificatory beliefs that directly map onto strategic exploitation in economic games. We cite Hilbig et al. (2023), who tested D against 58 personality traits and found that D emerged as the singular dispositional predictor of how individuals weigh their own versus others' utility. We conclude that our findings are best understood through the D framework, though they likely generalize to individuals low in Honesty-Humility given the substantial overlap.

**Changes Made:**
- Added new paragraph to Discussion/Theoretical implications section (line 2620)
- Added citation: Hilbig, B. E., Thielmann, I., Zettler, I., & Moshagen, M. (2023). The dispositional essence of proactive social preferences: The dark core of personality vis-à-vis 58 traits. Psychological Science, 34(2), 201-220.

---

## Comment R2.3: D-factor vs. Dark Triad Clarity

**Reviewer Comment:**
> "The manuscript occasionally describes D in a manner that resembles a composite or blend of the Dark Triad traits. Although D is often conceptualized as a unifying core underlying malevolent traits, it may be misleading to describe it as consisting of specific Dark Triad components (e.g., Machiavellianism, psychopathy, narcissism). For example, statements suggesting that the observed effects might be driven by particular Dark Triad traits risk conflating the theoretical status of D with that of its correlates. Clarifying this distinction would improve conceptual precision."

**Response:**
We agree with the reviewer and have revised all passages that could be misread as treating D as a composite of Dark Triad traits. The revised text now clearly positions D as the common core that dark traits share, with specific traits being expressions of D rather than components of it. See our response to R1.5 for detailed changes, which address both reviewers' concerns about this conceptual issue.

**Changes Made:**

- See R1.5 response for specific line changes

---

## Comment R2.4: HMM Effectiveness / Ecological Validity

**Reviewer Comment:**
> "The use of HMM agents is a particularly appealing feature of the study; however, it appears that these agents may not have been especially effective at mimicking human opponents. According to the authors' report, approximately 57% of participants either believed they were playing against a human or were unsure, implying that over 40% believed they were playing against a computer. This raises concerns about ecological validity and participants' engagement with the task. It may be informative to examine whether the primary effects replicate when analyses are restricted to participants who believed they were interacting with a human opponent, as perceptions of social interaction could meaningfully influence exploitative behavior."

**Response:**
We thank the reviewer for this important suggestion. We have addressed this concern in two ways. First, we added an explicit acknowledgment in the Limitations section that 43% of participants believed they were playing against a computer, which may affect the generalizability of findings to contexts where participants are certain of human interaction. Second, we conducted a sensitivity analysis restricting analyses to the 104 participants who believed they were playing against a human or were unsure (excluding those who indicated "Computer"). The key D × Opponent × Investment three-way interaction remains significant in this restricted sample (F(1, 6936.32) = 7.19, p = .007), demonstrating that our primary findings are robust to potential concerns about ecological validity. Full results are reported in the Supplementary Materials.

**Changes Made:**
- Added limitation in Discussion (line 2641) acknowledging the 43% who believed computer
- Added sensitivity analysis to Supplementary Materials (lines 248-339) with full model results for N=104 subsample
- Sensitivity analysis confirms key three-way interaction (p = .007)

---

## Comment R2.5: Extreme-Groups Limitation

**Reviewer Comment:**
> "The authors' use of an extreme-groups approach is interesting, but it necessarily excludes individuals in the middle of the D distribution. This limits the ability to draw conclusions about how D operates across its full range. I would encourage the authors to acknowledge this limitation more explicitly and to recommend that future research treat D as a continuous variable rather than dichotomizing it. Doing so would provide a more comprehensive understanding of how exploitative behavior varies as a function of D."

**Response:**
We agree with the reviewer and have added a new paragraph in the Limitations section that explicitly acknowledges this constraint. We note that while the extreme-groups design maximizes statistical power to detect group differences, it excludes individuals in the middle of the D distribution and limits our ability to draw conclusions about how D operates across its full range. We also note the possibility of non-linear relationships (e.g., moderate-D individuals might show qualitatively different patterns than those at the extremes) and explicitly recommend that future research treat D as a continuous variable.

**Changes Made:**
- Added new paragraph to Limitations (line 2639): "Second, our extreme-groups design, while maximizing statistical power to detect group differences, excludes individuals in the middle of the D distribution and limits our ability to draw conclusions about how D operates across its full range. The relationship between D and exploitation may not be linear; for instance, moderate-D individuals might show qualitatively different patterns than those at the extremes. Future research should treat D as a continuous variable to provide a more comprehensive understanding of how exploitative behavior varies as a function of D-factor scores."

---

## Comment R2.6: Power Analysis Concerns

**Reviewer Comment:**
> "I am also somewhat surprised by the claim that a sample size of 180 participants provides sufficient power to detect the targeted three-way interaction examined in the manuscript. This concern is heightened by the fact that the authors go on to test additional interactions, including a four-way interaction, despite basing their power analysis on a three-way interaction. Greater justification of the power assumptions, or a more cautious interpretation of higher-order interactions, would be warranted."

**Response:**
We appreciate the reviewer raising this important methodological concern. We have clarified the power analysis in two ways.

First, regarding power for the confirmatory three-way interaction: Our Monte Carlo simulations using the simr package specifically targeted the D × Opponent × Investment interaction, incorporating the repeated-measures structure of the design. Crucially, while our sample included 183 participants, the repeated-measures design (48 observations per participant across two games of 24 rounds each) yields approximately 8,800 total observations, providing substantially greater statistical power than 183 independent observations would suggest. We have now made this explicit in the Participants section. The fact that our confirmatory three-way interaction was indeed significant (p = .004) provides empirical confirmation that the study was adequately powered for its primary hypothesis.

Second, regarding higher-order interactions: We fully agree that four-way and five-way interactions require more statistical power than our study was designed to detect. We have now explicitly distinguished between confirmatory and exploratory analyses in the Statistical Analysis section. Specifically, we clarify that the three-way interaction (D × Opponent × Investment) was our confirmatory hypothesis as specified in the power analysis, while higher-order interactions involving presentation order and round number are treated as exploratory and should be interpreted with appropriate caution. This distinction is reiterated in the Limitations section, where we note that "the higher-order interactions involving presentation order and round number were exploratory and should be interpreted with appropriate caution given reduced statistical power for such complex effects."

**Changes Made:**
- Added sentence to Participants section (line 365) explaining that repeated-measures design yields ~8,800 observations, clarifying why N=183 provides adequate power
- Statistical Analysis section (lines 787-789) now explicitly distinguishes confirmatory (3-way) from exploratory (4-way, 5-way) analyses
- Limitations section (line 2622) acknowledges reduced power for higher-order interactions

---

## Comment R2.7: Gender Balance

**Reviewer Comment:**
> "Additionally, it would be helpful to know whether the high-D and low-D groups were balanced with respect to gender. Unequal gender distributions could complicate the interpretation of the results, particularly given known gender differences in both personality traits and economic decision-making."

**Response:**
We have now reported gender distribution by D-factor group. As the reviewer anticipated, we found a significant gender imbalance: the High-D group was 70% male (64/91), while the Low-D group was 41% male (37/91), χ²(1) = 15.04, p < .001. This pattern is consistent with the well-documented finding that males tend to score higher on dark personality traits. We acknowledge this imbalance in the Discussion/Limitations section.

**Changes Made:**
- Added to Participants section: Gender distribution by D-level group with chi-square test results
- See also response to R1.6d

---

## Comment R2.8: Analytic Strategy Clarity

**Reviewer Comment:**
> "At several points, the manuscript would benefit from greater clarity regarding the analytic strategy. The authors describe the study as employing a 2 × 2 design, yet elsewhere note that the sample size was chosen to detect a three-way interaction, and later report analyses involving four-way interactions. As a reader, it was somewhat difficult to determine which analyses corresponded to core hypotheses and which were exploratory. Streamlining the analytic approach, or more clearly distinguishing core analyses from exploratory analyses, would substantially improve readability and interpretability."

**Response:**
We thank the reviewer for this helpful suggestion. We have added a paragraph to the Statistical Analysis section that explicitly clarifies our analytic strategy. We now explain that:

1. **Confirmatory hypothesis:** The three-way interaction between D-factor, Opponent type, and Investment level was our confirmatory hypothesis, as specified in the power analysis. This tests whether high-D individuals show differential exploitation depending on opponent characteristics and investment magnitude.

2. **Control variable:** Presentation Order was included as a control variable to account for counterbalancing effects, not as a factor of primary theoretical interest.

3. **Exploratory analyses:** Higher-order interactions involving Order and Round (i.e., four-way and five-way interactions) are explicitly treated as exploratory and interpreted with appropriate caution.

This distinction is maintained throughout the Results and Discussion sections, with exploratory findings (such as the order effect) now labeled as such.

**Changes Made:**
- Added clarifying paragraph to Statistical Analysis section (lines 787-789)
- Discussion section now explicitly labels the four-way order interaction as "exploratory"
- See also response to R2.6 regarding power for higher-order interactions

---

## Comment R2.9: Learning Dynamics Interpretation

**Reviewer Comment:**
> "The manuscript raises a particularly intriguing idea that deserves further elaboration. The suggestion that high-D individuals may not possess superior social intelligence per se, but instead differ in their learning dynamics—such as being quicker to abandon pro-social norms when environmental predictability is detected—is very interesting. Expanding on this possibility would add depth to the discussion and help reframe the observed behavior as differential adaptation rather than simple strategic dominance. This distinction has important implications and could provide a fruitful direction for future research, particularly if individual learning rates can be modeled directly."

**Response:**
We thank the reviewer for highlighting this point. We have substantially expanded the learning dynamics interpretation in the Discussion. The revised text grounds the interpretation in our observed data: exploitation emerged gradually over rounds rather than immediately (suggesting learning rather than a pre-formed plan), and the strategy transfer effect points to experience-dependent adaptation—consistent with research showing that people transfer learned opponent models between games (Guennouni & Speekenbrink, 2022). We frame this in reinforcement learning terms—High-D individuals may update strategies more rapidly or hold weaker prosocial priors. We cite relevant literature on dual-process theories of cooperation (Bear & Rand, 2016) and research on altered prediction error processing in trait expressions through which D may manifest, such as psychopathy (Gregory et al., 2015; Atanassova et al., 2025), while being careful to note that D represents the common core underlying these expressions rather than any single trait. We discuss practical implications (interventions might focus on strengthening prosocial priors) and suggest that future research could test this by fitting computational models to individual trajectories.

**Changes Made:**
- Expanded Discussion paragraph on learning dynamics (lines 2552-2558) from 3 sentences to full paragraph (~230 words)
- Added three new citations to bibliography:
  - Bear & Rand (2016) PNAS - dual-process cooperation
  - Gregory et al. (2015) Lancet Psychiatry - reinforcement learning in psychopathy
  - Atanassova et al. (2025) Translational Psychiatry - foraging decisions and psychopathy

---

# Summary of Changes

| Section | Change Type | Description |
|---------|-------------|-------------|
| Introduction | Conceptual | Revised D-factor description |
| Methods | Addition | Trust game explanation |
| Methods | Addition | Predictor coding details |
| Methods | Addition | Gender distribution |
| Results | Clarification | Growth model comparison added; linear model retained |
| Results | Reorganization | Opponent ratings: full stats to supplement; summary + figure retained |
| Discussion | Addition | D vs. Agreeableness/HH discussion |
| Discussion | Addition | Learning dynamics elaboration |
| Discussion | Addition | Extreme-groups limitation |
| Supplement | Addition | Sensitivity analyses |

---

*Document prepared: [Date]*
