# Submission Guide - Royal Society Open Science

## Journal Requirements Summary

Based on the decision letter, RSOS requires:

### Required Files (all in editable format)

| File | Format | Description |
|------|--------|-------------|
| Manuscript (tracked) | .tex or .docx | Changes highlighted/bold/tracked |
| Manuscript (clean) | .tex or .docx | No highlighting, for typesetting |
| Figures | Individual files | Each figure as separate file |
| Tables | Editable file | Each table separately (N/A - no tables) |
| Captions | Editable file | All figure/table captions |
| Point-by-point response | .doc/.docx | "Essential" - must be editable |
| Supplementary material | Any | ESM files |

### Additional Requirements
- Equations must be editable text (not images)
- Vancouver style references preferred (with DOIs)
- Acknowledgements section before references
- 100-word summary for press release (entered in ScholarOne)
- Data access statement with DOIs

---

## Strategy Decision: .tex with latexdiff

### Why .tex over .docx?

| Factor | .tex | .docx |
|--------|------|-------|
| **Native format** | ✅ Already have article.tex | ❌ Requires conversion |
| **Math/equations** | ✅ Perfect rendering | ⚠️ May break |
| **Figure placement** | ✅ Preserved | ⚠️ May shift |
| **Original baseline** | ✅ In git (e9b2ca8) | ❌ Would need to create |
| **Track changes tool** | ✅ latexdiff installed | ✅ Word Compare |
| **Risk level** | Low | Medium-High |

### Potential latexdiff issues to watch:
- Major structural changes (R1.2 growth model) may produce messy diffs
- R-generated inline code changes may be marked oddly
- **Mitigation**: Test latexdiff after completing revisions, clean up if needed

---

## Submission Workflow

```
Phase 1: Content Revisions (CURRENT)
├── Complete all reviewer comments in article.Rmd
├── Update reviewer_responses_template.md as we go
└── Track progress in REVISION_TRACKER.md

Phase 2: Pre-submission Preparation
├── Final render: Rscript -e 'rmarkdown::render("article.Rmd")'
├── Verify article.pdf looks correct
├── Render supplement: Rscript -e 'rmarkdown::render("supplement.Rmd")'
└── Commit all changes

Phase 3: Generate Tracked Changes Version
├── Extract original: git show e9b2ca8:article.tex > article_original.tex
├── Run latexdiff: latexdiff article_original.tex article.tex > article_tracked.tex
├── Compile: xelatex article_tracked.tex (may need 2 runs)
├── Review article_tracked.pdf for readability
└── Clean up latexdiff output if needed

Phase 4: Extract Individual Files
├── Figures:
│   ├── Copy plots/figure1.pdf → Figure1.pdf
│   ├── Copy article_files/figure-latex/gamesPlot2-1.pdf → Figure2.pdf
│   ├── Copy article_files/figure-latex/oppExploit-1.pdf → Figure3.pdf
│   ├── Copy article_files/figure-latex/orderEffectPlot-1.pdf → Figure4.pdf
│   └── Copy article_files/figure-latex/plotRatings-1.pdf → Figure5.pdf
├── Captions: Extract from article.tex → captions.docx
└── Tables: N/A (no formal tables in manuscript)

Phase 5: Prepare Response Document
├── Convert reviewer_responses_template.md → reviewer_responses.docx
└── Review formatting, ensure all [To be completed] are filled

Phase 6: Upload to ScholarOne
├── Step 1: Attach point-by-point response
├── Step 2: Enter 100-word summary
├── Step 3: Upload files:
│   ├── article_tracked.tex (tracked changes version)
│   ├── article.tex (clean version)
│   ├── Figure1.pdf through Figure5.pdf
│   ├── captions.docx
│   ├── reviewer_responses.docx
│   └── supplement.pdf (and supplement.Rmd if required)
├── Step 6: Complete data access statement
└── Step 7: Review PDF proof and submit
```

---

## File Checklist (Pre-submission)

### Content Files
- [ ] article.Rmd - all revisions complete
- [ ] article.tex - clean version generated
- [ ] article_tracked.tex - latexdiff version generated
- [ ] article_tracked.pdf - reviewed for readability
- [ ] supplement.Rmd - updated if needed
- [ ] supplement.pdf - rendered
- [ ] bib/DFP2.bib - all new citations added

### Individual Figures
- [ ] Figure1.pdf - HMM investor behavior
- [ ] Figure2.pdf - Investment/return patterns (gamesPlot2)
- [ ] Figure3.pdf - Exploitation patterns (oppExploit)
- [ ] Figure4.pdf - Order effects (orderEffectPlot)
- [ ] Figure5.pdf - Opponent ratings (plotRatings)

### Supporting Documents
- [ ] captions.docx - all figure captions
- [ ] reviewer_responses.docx - point-by-point response
- [ ] 100-word summary - for ScholarOne entry

### Data/Ethics
- [ ] Data access statement ready
- [ ] Dataset DOI available
- [ ] Acknowledgements section in manuscript

---

## Commands Reference

```bash
# Render manuscript
Rscript -e 'rmarkdown::render("article.Rmd")'

# Generate tracked changes
git show e9b2ca8:article.tex > article_original.tex
latexdiff article_original.tex article.tex > article_tracked.tex
xelatex article_tracked.tex
xelatex article_tracked.tex  # Run twice for references

# Convert response to docx
pandoc reviewer_responses_template.md -o reviewer_responses.docx

# Extract captions (manual or script)
grep -A2 "\\caption{" article.tex > captions.txt
```

---

## Notes

- **Pre-revision baseline commit**: e9b2ca8
- **latexdiff version**: 1.3.3 (installed at /Library/TeX/texbin/latexdiff)
- **Figure count**: 5 figures, 0 tables
- **Original submission format**: .tex (compiled from article.Rmd)

---

*Created: 2026-01-12*
*Last updated: 2026-01-12*
