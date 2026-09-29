# Toxicity Metrics Comparison

## Introduction

In environmental risk assessment, toxicity data from laboratory-based concentration–response (C–R) experiments are routinely summarised into single-value metrics for use in decision-making frameworks such as species sensitivity distributions (SSDs). These summaries reduce complex biological responses to values that can be compared across species, and protective concentration thresholds for chemicals in the environment are derived from them. The choice of toxicity metric therefore affects estimates of ecological risk and the environmental guideline values based on them.

Historically, the most widely used metrics have been the no-observed-effect and lowest-observed-effect concentrations (NOEC and LOEC) and effect-based measures such as the ECx (e.g., EC10). These metrics are commonly used as inputs to SSDs, where a distribution of species-specific toxicity values is used to estimate concentrations protective of a specified proportion of species (e.g., HC5). However, the statistical properties and ecological interpretation of these metrics have been widely debated. NOEC and LOEC values are sensitive to experimental design and limited to the tested concentrations. ECx values depend on the selection of an arbitrary effect level and may not align conceptually with “no-effect” protection goals (van Straalen & Denneman, 1989; Van der Hoeven et al., 1997; Green et al., 2013).

These limitations prompted the development of alternative toxicity metrics, intended to represent the underlying biological response more closely and to have better-defined statistical properties than the NOEC and the ECx. In particular, there has been increasing emphasis on model-based approaches that use the full concentration–response relationship, rather than relying on discrete hypothesis testing or arbitrary effect thresholds.

## Scope

This repository contains code and analyses for comparing commonly used **ecotoxicological toxicity metrics** and their statistical estimation methods.

The repository is intended as an evolving framework for:

- Reproducing published methods
- Comparing numerical behaviour across metrics
- Exploring implications for ecological risk assessment and SSDs

The aim of this repository is to provide a transparent and reproducible framework for comparing these toxicity metrics using a common dataset and consistent implementation. The LECₓ, NEC, NSEC and NSEC-MSD are applied to the same data, to show how they differ in their assumptions, in their sensitivity to experimental design, and in their implications for ecological risk assessment. This forms part of a broader effort to compare toxicity metrics and their estimation methods, with future extensions to include approaches such as benchmark dose modelling and model-averaged thresholds.

## Currently included

Four estimators are currently considered:

- **LECₓ** (Lowest Effect Concentration with associated effect size), which combines hypothesis testing with model-based interpolation to retain both statistical confidence and effect magnitude [@debruyn2013lecx];
- **NEC** (No-Effect Concentration), which estimates a threshold parameter directly within a concentration–response model, assuming the existence of a true biological threshold [@fox2010nec];
- **NSEC** (No-Significant-Effect Concentration), which defines a no-effect threshold based on statistical indistinguishability from the control within a continuous model framework [@fisher2023nsec];
- **NSEC-MSD**, the LECₓ as calculated for replicated tests, written in the NSEC framing: the concentration at which the fitted curve falls below its modelled control response by the percent minimum significant difference (PMSD) of Dunnett's test.

De Bruyn and Elphick (2013) give two methods for the LECₓ. Method 1, for replicated tests, sets the effect size to the Dunnett PMSD; the NSEC-MSD is this method. Method 2, for tests without replication, sets it from the confidence band on the baseline response; the NSEC is the same calculation. The two differ only in the departure from the control that is treated as significant, and `compare_nsec_lecx.qmd` shows how that difference splits into three factors.

They differ in how they summarise a C–R relationship, and in what they assume a biologically meaningful absence of effect to be. Comparing them shows how methodological choices affect the derived toxicity values.

## Contents

- `compare_nsec_lecx.qmd`: the worked comparison of the LECₓ, NEC, NSEC and NSEC-MSD on the fish growth data of Fisher and Fox (2023).
- `R/dunnett.R`: the one-sided Dunnett critical value and minimum significant difference, shared by the two analyses.
- `case_studies/ce_nsec_analysis.R`: the NSEC and NSEC-MSD for ten single-species cerium chloride toxicity tests, for comparison with the LECₓ calculated in CETIS. The script reads the workbook `case_studies/Ce dataset - for statistical comparison.xlsx`, supplied by the data owner, and writes its results to `case_studies/output/`.
