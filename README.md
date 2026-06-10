# Toxicity Metrics Compariso

## Introduction

In environmental risk assessment, toxicity data from laboratory-based concentration–response (C–R) experiments are routinely summarised into single-value metrics for use in decision-making frameworks such as species sensitivity distributions (SSDs). These summaries aim to reduce complex biological responses into comparable units across species, enabling the derivation of protective concentration thresholds for chemicals in the environment. As a result, the choice of toxicity metric plays a central role in shaping estimates of ecological risk and environmental guideline values.

Historically, the most widely used metrics have been the no-observed-effect concentration (NOEC/LOEC) and effect-based measures such as the ECx (e.g., EC10). These metrics are commonly used as inputs to SSDs, where a distribution of species-specific toxicity values is used to estimate concentrations protective of a specified proportion of species (e.g., HC5). However, the statistical properties and ecological interpretation of these metrics have been widely debated. NOEC/LOEC values are sensitive to experimental design and limited to tested concentrations, while ECx values depend on the selection of an arbitrary effect level and may not align conceptually with “no-effect” protection goals (van Straalen & Denneman, 1989; Van der Hoeven et al., 1997; Green et al., 2013).

These limitations have motivated the development of alternative toxicity metrics that aim to better represent underlying biological responses while providing improved statistical properties. In particular, there has been increasing emphasis on model-based approaches that use the full concentration–response relationship, rather than relying on discrete hypothesis testing or arbitrary effect thresholds.

## Scope

This repository contains code and analyses for comparing commonly used **ecotoxicological toxicity metrics** and their statistical estimation methods.

The repository is intended as an evolving framework for:

- Reproducing published methods
- Comparing numerical behaviour across metrics
- Exploring implications for ecological risk assessment and SSDs

The aim of this repository is to provide a transparent and reproducible framework for comparing these toxicity metrics using a common dataset and consistent implementation. By examining how LECₓ, NEC, and NSEC behave under identical conditions, we aim to highlight differences in their assumptions, sensitivity to experimental design, and implications for their use in ecological risk assessment. This forms part of a broader effort to compare toxicity metrics and their estimation methods, with future extensions to include approaches such as benchmark dose modelling and model-averaged thresholds.

## Currently included

Three such approaches are currently considered:

- **LECₓ** (Lowest Effect Concentration with associated effect size), which combines hypothesis testing with model-based interpolation to retain both statistical confidence and effect magnitude [@debruyn2013lecx];
- **NEC** (No-Effect Concentration), which estimates a threshold parameter directly within a concentration–response model, assuming the existence of a true biological threshold [@fox2010nec];
- **NSEC** (No-Significant-Effect Concentration), which defines a no-effect threshold based on statistical indistinguishability from the control within a continuous model framework [@fisher2023nsec].

These approaches reflect different ways of summarising C–R relationships for use in environmental assessment, and embody different assumptions about what constitutes a biologically meaningful absence of effect. As such, they provide an opportunity to explore how methodological choices influence derived toxicity values.







