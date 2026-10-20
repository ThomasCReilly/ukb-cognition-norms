Normative Modelling of Cognition in the UK Biobank
===================================================

Analysis code for [PAPER TITLE / CITATION - add when available].

The pipeline:
  1. Cleans the UK Biobank cognitive test data.
  2. Builds GAMLSS normative models (age, sex, education) for 12 cognitive tests
     and generates out-of-sample z-scores with participant-level 10-fold
     cross-validation.
  3. Removes implausible extreme scores.
  4. Defines cognitive groups (healthy cognition, borderline, impaired cognition)
     and links them to all-cause dementia diagnoses.
  5. Relates cognition to plasma GFAP and NfL (Olink).
  6. Models how far ahead cognition and biomarkers predict dementia
     (GEE models, look-forward landmark models, and yearly divergence analyses).


DATA AVAILABILITY
-----------------
This repository contains code only. No UK Biobank data is included.

UK Biobank data are available to approved researchers through the UK Biobank
Access Management System (https://www.ukbiobank.ac.uk). The input files were
exported from the UK Biobank Research Analysis Platform (RAP).


REQUIREMENTS
------------
Python (developed on 3.13):

    pip install -r requirements.txt

R (developed on 4.5.2). Install the required packages from an R session:

    install.packages(c(
      "dplyr", "tidyr", "stringr", "readr",            # data handling
      "gamlss", "gamlss.dist", "gamlss.add",           # normative models
      "caret", "moments",                              # CV folds, skew/kurtosis
      "future.apply", "progressr",                     # parallel runs
      "ggplot2", "gridExtra", "svglite",               # static plots
      "plotly", "htmlwidgets", "webshot2",             # 3D surfaces (03e)
      "geepack", "broom", "survival"                   # 07* notebooks (via rpy2)
    ))

    # Optional: rasterises points inside SVGs in 03e
    install.packages("ggrastr")


INPUT FILES
-----------
Put these files in the working directory (or edit the path at the top of the
script that reads them; look for "#replace with data path").

  File                        Used in   Contents
  --------------------------  --------  -----------------------------------------
  ukb_cog_clean.csv           02        Wide cognitive test scores per instance,
                                        plus eid, Sex, Edu and Age per instance
  MP_puzzles_viewed.csv       02        Matrix pattern: number of puzzles viewed
  SD_0_60.csv                 02        Symbol digit: IDs scoring 0 or > 60 - can use unfiltered scores
  FI_0.csv                    02        Fluid intelligence: IDs scoring 0 - can use unfiltered scores
  TR_0.csv                    02        Tower rearranging: IDs scoring 0 - can use unfiltered scores
  PMem_data.csv               05a       Prospective memory results
  Caselists/ID_dateDem.txt    05b       eid and date of all-cause dementia
                                        diagnosis (tab-separated, no header)
  Caselists/id_DOB.csv        05b       eid, birth year and month
  olink_olink_instance_0.csv  06a       Olink proteomics (GFAP, NfL), instance 0
  olink_confounds.csv         06a       Age, sex, BMI and creatinine for the
                                        Olink sample

Cognitive test abbreviations:
  PMat = pairs matching          NM  = numeric memory
  TMA  = trail making A          TMB = trail making B
  MP   = matrix pattern          FI  = fluid intelligence
  RT   = reaction time           BL  = broken letters
  SD   = symbol digit            PAL = paired associate learning
  TR   = tower rearranging       PV  = picture vocabulary
  PMem = prospective memory (binary; not normatively modelled)


PIPELINE AND RUN ORDER
----------------------
Run the scripts in numerical order. Each step reads the files written by earlier steps from the working directory.

01  Data extraction on the UKB RAP (not included).

02  02_data_cleaning.ipynb
    Quality-control filters for each cognitive test (e.g. too few puzzles
    viewed, zero scores, implausibly long completion times).
    -> ukb_cog_filtered.csv

03  Normative modelling (R)
    03a_data_prep.R
        Reshapes the data to long format (one row per person, test and instance).
        -> long_data.RData
    03b_cogtest_descriptive_analysis.R
        Distribution summaries and plots for each test.
    03c_dist_test_final.R
        Fits GAMLSS models with 9 candidate distribution families per test and
        ranks them by AIC.
        -> gamlss_results/
    03d_ukb_cog_CV_final.R
        Fits the selected distribution per test with 10-fold cross-validation
        split by participant. Trained on baseline data and predicted for all
        instances.
        -> cv_results/cv_results_all_tests.csv  (z-scores and percentiles)
    03e_vis_batch_final.R
        3D normative surfaces and cross-sections. Paths can be set with the
        environment variables MODELS_ROOT, LONG_DATA_PATH and OUTDIR.

04  04_extremes_analysis.ipynb
    Tail-inflation analysis. Removes isolated extreme scores (|z| > 3 with no
    other test beyond |1.5| at the same instance).
    -> cv_results_all_tests_cleaned.csv

05  Group definition and descriptives
    05a_PMem_transform.r
        Prospective memory to long format.
        -> PMem_long.csv
    05b_group_definition.ipynb
        Flags impaired tests (1.5 SD), merges PMem and dementia diagnoses, and
        defines the groups:
          HC (healthy cognition)     0 impaired tests at every instance
          Borderline                 exactly 1 impaired test at an instance
          IC (impaired cognition)    2 or more impaired tests at one instance
        -> cv_results_all_tests_merged_with_groups.csv
    05c_IC_threshold_analysis.ipynb
        Compares impairment thresholds (1.0, 1.5, 2.0, 2.5 SD): prevalence of
        IC, and sensitivity and specificity for dementia. This supports the
        1.5 SD threshold used in 05b.
    05d_group_descriptives_clean.ipynb
        Descriptive tables by group.
    05e_UKB_phenotype_transitions_minimal.ipynb
        Transitions between groups over follow-up from baseline.

06  Olink biomarkers (GFAP, NfL)
    06a_olink_dataprep.ipynb
        Calculates eGFR (CKD-EPI 2021) and merges biomarkers with cognition.
        -> cog_biomarker_merged.csv
    06b_olink_stage1_confound.ipynb
        Residualises the biomarkers for age, sex, BMI and eGFR.
        -> stage1_adjusted_data.csv
    06c_olink_stage2_correlation.ipynb
        Association of adjusted biomarkers with cognitive z-scores.
    06d_combined_multimodal.ipynb
        Logistic model of dementia: age + composite cognition + GFAP + NfL,
        with ROC comparison against nested models.

07  Dementia prediction
    07a_geeglm_pertest_groups.ipynb
        GEE logistic models of dementia for each test and for composite
        cognition.
    07b_lookforward_dynamic_int.ipynb
    07c_lookforward_dynamic_pertest_int.ipynb
        Look-forward landmark models: risk of dementia within 1-20 years of
        a visit (composite, then per test).
    07d_divergence_yearly.ipynb
    07e_divergence_yearly_pertest.ipynb
        Nested case-control analysis by years before diagnosis, following
        You et al. (2024), to find when future cases first diverge
        cognitively (composite, then per test).
        07d -> divergence_yearly/divergence_results_pooled.csv
    07f_divergence_yearly_olink.ipynb
        The same divergence analysis for GFAP and NfL, compared against the
        cognitive results from 07d.

Several notebooks (06c, 06d, 07a) have a GROUP setting at the top
('all', 'hc' or 'ic'). Change it and re-run the notebook for each group.


NOTES
-----
- The R scripts 03c and 03d can run tests in parallel; see the USAGE notes at
  the top of each script.
- Random seeds are fixed where sampling is used (e.g. control matching in
  07d-07f).


CONTACT
-------
reilly.t@unimelb.edu.au
