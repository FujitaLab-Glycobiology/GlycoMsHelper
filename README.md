GlycoMsHelper manual
================
Weize Kong and Morihisa Fujita

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.jpg" alt="" width="80%" style="display: block; margin: auto auto auto 0;" />

# GlycoMsHelper

<!-- badges: start -->

<!-- badges: end -->

The goal of GlycoMsHelper is to relieve users from labor-intensive
glycomic MS data analysis by enabling systematic extraction of
glycan-derived spectra and generation of traceable candidate
compositions for downstream validation.

## Installation

You can install the development version of GlycoMsHelper from
[GitHub](https://github.com/FujitaLab-Glycobiology/GlycoMsHelper) with:

``` r
# install.packages("devtools")
devtools::install_github("FujitaLab-Glycobiology/GlycoMsHelper")
```

## Basic workflow of GlycoMsHelper

### STEP1:

Load the package and read the `.mzML` file.

#### Codes

``` r
library(Spectra)
library(GlycoMsHelper)

# Read the MS file
mass_spectrum_data = Spectra::Spectra(your_ms_file_path, source = MsBackendMzR())
```

### STEP2:

Check the MS (`.mzML`) file.

Check the MS (`.mzML`) file to confirm the presence of MS2 spectra. If
no MS2 spectrum is detected, the pipeline will terminate automatically
to prevent downstream errors.

#### Codes

``` r
# Check the MS file to confirm the presence of MS2 spectra. If no MS2 spectrum is detected, the pipeline will terminate automatically to prevent downstream errors.
GlycoMsHelper::MsFileChecker(mass_spectrum_data) 
```

### STEP3:

Quality control of the MS (`.mzML`) file

The `SpectrumQcFilter` function performs quality control on MS data by
evaluating three variables: peaks count (`peaksCount`), total ion
current (`totIonCurrent`), and retention time (`rtime`).

#### QC logic

- `peaksCount` & `totIonCurrent`: Scans with `peaksCount` or
  `totIonCurrent` values below the specified threshold are removed to
  exclude low-quality or noise-dominant spectrum.

- `rtime`: Spectra are filtered using a defined temporal window. Only
  MS1 and MS2 spectrum falling within the specified `rtime` range are
  retained.

#### QC method

- `mean_sd`: adaptive threshold detection based on the distribution of
  sorted variable differences.
- `quantile_prob`: thresholding based on quantiles calculated using
  stats::quantile().
- `start_end`: fixed-range filtering using predefined lower and upper
  boundaries.

#### Output

``` text
QC_result
|-- filtered_ms_data                  Filtered MS data for next step or could be exported
|-- pic_ms_varibles_filtered          Diagnostic plots showing the distributions of MS1 and MS2 variables after filtering.
|   |                                 Filtered and retained spectra are shown in different colors.
|   |                                 The upper plot shows variables ordered by retention time. 
|   |                                 the lower plot shows variables sorted in ascending order.
|   |-- MS1_peaksCount
|   |-- MS1_totIonCurrent
|   |-- MS1_rtime
|   |-- MS2_peaksCount
|   |-- MS2_totIonCurrent
|   |-- MS2_rtime
|-- pic_ms_varibles_unfiltered        Diagnostic plots showing the distributions of MS1 and MS2 variables before filtering.
|   |                                 These plots can be used to determine appropriate QC methods and thresholds.
|   |                                 The upper plot shows variables ordered by retention time.
|   |                                 the lower plot shows variables sorted in ascending order.
|   |-- MS1_peaksCount
|   |-- MS1_totIonCurrent
|   |-- MS1_rtime
|   |-- MS2_peaksCount
|   |-- MS2_totIonCurrent
|   |-- MS2_rtime
```

#### Codes

``` r
# filter low quality MS1 and MS2 spectrum based on Spectra::spectraVariables()
qc_results = GlycoMsHelper::SpectrumQcFilter(
  ms_data = mass_spectrum_data, 
  filter_method_ms1 = c(
    peaksCount = 'mean_sd', 
    totIonCurrent = 'quantile_prob', 
    rtime = 'start_end'
    ), 
  threshold_ms1 = list(
    peaksCount = 2, 
    totIonCurrent = c(0.05, 1), 
    rtime = c(8*60, 45*60)
    ), 
  filter_method_ms2 = c(
    peaksCount = 'quantile_prob', 
    totIonCurrent = 'quantile_prob', 
    rtime = 'start_end'
    ), 
  threshold_ms2 = list(
    peaksCount = c(0.1, 1), 
    totIonCurrent = c(0.05, 1), 
    rtime = c(8*60, 50*60)
    ), 
  plot_option = T
)
                          
mass_spectrum_data_filtered = qc_results$filtered_ms_data
```

#### Parameters

- **`ms_data`**  
  The MS data will be used for the quality control.

- **`filter_method_ms1` and `filter_method_ms2`**  
  Named character vectors specifying the QC method used for each MS1 and
  MS2 spectrum variable, respectively.

  - **`mean_sd`**: Adaptive threshold detection using the “knee point”
    method

    **Method:**

    1.  Sort the variable values in ascending order
    2.  Calculate consecutive differences between sorted values
    3.  Identifies significant jumps using the following criterion:
        `difference >= mean(all differences) + n * sd(all differences)`,
        where `n` is defined in `threshold_ms1` or `threshold_ms2`.
    4.  Select the position with the maximum significant jump as the
        cutoff threshold.

    **Parameter:** \|——\|——\| \| `n` \| Sensitivity of jump detection.
    Higher `n` → fewer candidate jumps; lower `n` → more candidate
    jumps. \| **Example:**

    ``` r
    filter_method_ms1 = c(peaksCount = "mean_sd")
    threshold_ms1 = list(peaksCount = 2)
    ```

    This example uses two standard deviations above the mean difference
    as the criterion for detecting a significant jump.

  - **`quantile_prob`**: Quantile-based filtering using R’s internal
    `stats::quantile()` function.

    **Method:**

    1.  Sort the variable values in ascending order
    2.  Calculate the specific quantiles to define lower and upper
        boundaries for peaksCount, totIonCurrent, or rtime (use
        `stats::quantile()` function)
    3.  Applies these calculated quantiles as lower and upper threshold
        boundaries.

    **Parameter:** \|——\|————-\| \| `c(lower bound, upper bound)` \|
    Probability vector passed to `stats::quantile()`, defining the lower
    and upper quantile bounds (set via `threshold_ms1` /
    `threshold_ms2`). \| **Example:**

    ``` r
    filter_method_ms1 = c(peaksCount = "quantile_prob")
    threshold_ms1 = list(peaksCount = c(0.1, 1))
    ```

    This uses the 10th percentile as the minimum threshold for
    `peaksCount`.

  - **`start_end`**: Fixed range filtering. This is a deterministic
    method used when you have predefined limits.

    **Method**

    1.  Directly compares variables against a hard-coded range.

    **Parameter:** \|——\|————-\| \| `c(min, max)` \| A length-2 numeric
    vector specifying the lower and upper bounds. \| **Example:**

    ``` r
    filter_method_ms1 = c(rtime = "start_end")
    threshold_ms1 = list(rtime = c(8 * 60, 45 * 60))
    ```

    This restricts the data to spectra acquired between 8 and 45
    minutes.

- **`threshold_ms1` and `threshold_ms2`**  
  The list vector for the QC method defined in `filter_method_ms1` and
  `filter_method_ms2`.

  - For `mean_sd`, the parameter is `n`.
  - For `quantile_prob`, the parameter is `c(lower_bound, upper_bound)`.
  - For `start_end`, the parameter is `c(min_value, max_value)`.

- **`plot_option`**  
  Set `plot_option = T` to include diagnostic plots in the output list.
  These plots visualize the distribution of MS1 and MS2 spectrum
  variables after filtering, which could be used to verify the QC
  results, filtered variables and left variables are in different color.
  The upper plot is variables sorted in time, the lower plot is
  variables sorted in ascending order.

### STEP4:

MS2 spectrum denoising

The `MS2SpectrumDenoising` function performs denoising on every MS2
spectrum for the MS data after QC.

#### Denoising logic

- `spline_segmentation_regression`, `spline_regression`,
  `segmentation_regression`

- `quantile_prob`: Provide dynamic noise threshold estimation for each
  of the MS2 spectrum.

- `fixed_value`: Use the fixed noise threshold for every MS2 spectrum.

#### Denoising method

- `spline_segmentation_regression`
- `spline_regression`
- `segmentation_regression`
- `quantile_prob`
- `fixed_value`

#### Output

``` text
Denoised_result
|-- denoised_ms_data                  Denoised MS data for next step or could be exported
|-- denoising_regression_info         Detailed denoised parameters for each MS2 spectrum
```

#### Codes

``` r
# MS2 spectrum denoising
ms2_denoising_info = list(
  spline_segmentation_regression = list(
    spar_start = -1.5, spar_end = 1.5, spar_step = 0.02, 
    RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1, 
    use_cv = T, top_n_to_remove = 5, 
    segmentated_non_linear_transform_fun = function(z) z^2+z
    ), 
  spline_regression = list(
    spar_start = -1.5, spar_end = 1.5, spar_step = 0.02, 
    RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1, 
    use_cv = T, top_n_to_remove = 5
    ), 
  segmentation_regression = list(
    segmentated_non_linear_transform_fun = function(z) z^2+z
    ), 
  quantile_prob = 0.05, 
  fixed_value = 30 
)

denoised_results = GlycoMsHelper::MS2SpectrumDenoising(
  ms_data = mass_spectrum_data_filtered, 
  ms2_spectrum_transform_method = 'log2_transform', 
  ms2_denoising_method = 'spline_segmentation_regression', 
  ms2_denoising_detail = ms2_denoising_info
  ) 

mass_spectrum_data_filtered_denoised = denoised_results$denoised_ms_data

# export(mass_spectrum_data_filtered_denoised,
#        backend = MsBackendMzR(),
#        file = "your_ms_file_path_and_name.mzML", BPPARAM = SerialParam())
```

#### Parameters

- **`ms_data`**  
  MS data will be used for the denoising.

- **`ms2_spectrum_transform_method`**  
  Define the non-linear signal transform method for MS2 spectrum. After
  signal transform, the noise threshold could be determined.

  - **`log2_transform`**: log2(x + 1) transform

  - **`asinh_transform`**: asinh(x) transform

  - **`non_transform`**: the signal is not non-linear transformed, may
    not suitable for the `spline_segmentation_regression`,
    `spline_regression`, and `segmentation_regression`

  - **`function(z)`**: users defined function, could be any thing.
    **Example**

    ``` r
    ms2_spectrum_transform_method = function(z) log10(z + 1)
    ```

    This uses log10(x + 1) as the non-linear transform function for MS2

- **`ms2_denoising_method`**

  - **`spline_segmentation_regression`**: automatically determine to use
    the spline-based or the segmentation regression-based method for
    denoising based on the properties of the MS2 spectrum

    **Method:**

    1.  Apply a non-linear transform to the MS2 spectrum intensities
        (defined in `ms2_spectrum_transform_method`).
    2.  remove the duplicated transformed intensity values and sort the
        transformed spectrum in ascending order.
    3.  check the top 4 transformed intensity values and their
        corresponding m/z value, the difference between the m/z value
        will be calculated, and then the cv value of the difference will
        be calculated. if the cv value is smaller than 0.05 and the m/z
        difference is smaller then 1.11, bigger than 0.09, spline-based
        denoising method will be used, or the segmentation
        regression-based denoising method will be used.
    4.  Perform either the segmentation regression or spline regression
        to fit the transformed spectrum
    5.  for the segmentation regression-based method, the noise
        threshold is determined at the breakpoint of the linear
        regression and the non-linear regression. For the spline
        regression-based method, the noise threshold is determined at
        the position with the minimum first derivative of the fitted
        curve.

    **Parameters:**

    - `spar_start`: the starting value of the spar parameter for the
      spline regression, default is -1.5
    - `spar_end`
    - `spar_step`
    - `RMSE_weight`
    - `CV_weight`
    - `D2_weight`
    - `D1_weight`
    - `use_cv`
    - `top_n_to_remove`
    - `segmentated_non_linear_transform_fun`

    **Example:**

  - `spline_regression`

    - Method:

    - Parameters:

      - `spar_start`
      - `spar_end`
      - `spar_step`
      - `RMSE_weight`
      - `CV_weight`
      - `D2_weight`
      - `D1_weight`
      - `use_cv`
      - `top_n_to_remove`

    - Example:

  - `segmentation_regression`

    - Method:

    - Parameters:

      - `segmentated_non_linear_transform_fun`

    - Example:

  - `quantile_prob`

    - Method:

    - Parameters:

    - Example:

  - `fixed_value`

    - Method:

    - Parameters:

    - Example:

- **`ms2_denoising_detail`**

``` r
# Find the spectrum likely to be glycan based on diagnostic ions 
diagnostic_frags = c(
  HexNAc =              204.08667, 
  HexNAc_ProA =         441.2708, 
  dHex_HexNAc_ProA =    587.3287, 
  HexNAc_HexNAc_ProA =  644.3502, 
  Hex_ProA =            400.2442, 
  Hex_Hex_ProA =        562.2970, 
  dHex_Hex_ProA =       546.3021, 
  Hex =                 163.06007, 
  Bi_HexNAc =           407.16607, 
  Bisecting =           1009.4824, 
  Bisecting_dHex =      1155.5403, 
  Hex_HexNAc_ProA =     603.3236
)

diagnostic_results = GlycoMsHelper::FindSpectrumByDiagnosticFragments(
  ms_data = mass_spectrum_data_filtered_denoised, 
  ms_data_raw = mass_spectrum_data_filtered, 
  diagnostic_frags_list = diagnostic_frags, 
  diagnostic_frags_exp = 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA', 
  ppm_val = 100
  )

# export(diagnostic_results$selected_ms_data,
#        backend = MsBackendMzR(),
#        file = "your_ms_file_path_and_name.mzML", BPPARAM = SerialParam())


likely_glycan_spectrum_info = diagnostic_results$spectrum_info
```

### STEP4:

Construct the glycan lib

``` r
molecular_formula_all = c(
  # monosaccharides
  Hex = 'C6H10O5',
  HexNAc = 'C8H13N1O5',
  dHex = 'C6H10O4',
  Neu5Ac = 'C11H17N1O8',
  HexA = 'C6H8O6',
  Neu5Gc = 'C11H17N1O9',
  Pentose = 'C5H8O4',
  KDN = 'C9H14O8', 
  EtNP = 'C2H6N1O3P1', 
  AHM = 'C6H8O4', 
  # label
  ProA = 'C13H23N3O1', 
  AB = 'C7H10N2O1', 
  PA = 'C5H8N2', 
  # adduct
  H = 'H1',
  Na = 'Na1',
  K = 'K1',
  Li = 'Li1',
  Mg = 'Mg1'
  )

N_glycan_lib = GlycoMsHelper::ConstructGlycanLibrary(
  glycan_type = 'N_glycan', 
  min_charge_state = 1, 
  max_charge_state = 3, 
  derivatization_type = 'ProA', 
  adduct_type = c('H', 'Na', 'K'), 
  min_total_monosaccharides_num = 3, 
  max_total_monosaccharides_num = 21, 
  min_Hex_num = 1,      max_Hex_num = 12, 
  min_HexNAc_num = 2,   max_HexNAc_num = 10, 
  min_dHex_num = 0,     max_dHex_num = 4, 
  min_Neu5Ac_num = 0,   max_Neu5Ac_num = 4, 
  min_HexA_num = 0,     max_HexA_num = 2, 
  min_Neu5Gc_num = 0,   max_Neu5Gc_num = 0, 
  min_Pentose_num = 0,  max_Pentose_num = 0, 
  min_KDN_num = 0,      max_KDN_num = 0
  )

N_glycan_library = N_glycan_lib$glycan_monosaccharides_library

# Optional: Identify the monoisotopic peak and the most abundant M+1 isotopologue for each glycan in the library, then compute their abundance ratio. 
N_glycan_library_iso_info = GetMonoisoAndIsotopologueRatio(glycan_lib = N_glycan_library, 
                                          molecular_names = colnames(N_glycan_lib$monosaccharides_adduct_num),
                                          molecular_formula_list = molecular_formula_all, 
                                          monosaccharides_names = colnames(N_glycan_lib$monosaccharides_combination), 
                                          threshold_iso_probalility = 0.01) 
```

### STEP5:

Match the likely glycan spectrum to glycan lib

``` r
likely_glycan_spectrum_matching_result = GlycoMsHelper::FindPossibleGlycanComposition(
  spectrum_info = likely_glycan_spectrum_info, 
  glycan_lib = N_glycan_library, 
  max_precursor_mz_ppm = 150, 
  max_possible_candidates_num = 5
  )

# write.csv(likely_glycan_spectrum_matching_result, file = 'your_csv_file_path_and_name.csv')
```

### STEP6:

Find the candidate glycan composition based on isotopics distribution

``` r
final_glycan_spectrum_matching_result = GlycoMsHelper::ValidateGlycanCompositionByIsotopePattern(
  spectrum_matching_info = likely_glycan_spectrum_matching_result, 
  molecular_names = colnames(N_glycan_lib$monosaccharides_adduct_num), 
  molecular_formula_list = molecular_formula_all, 
  ms_data = mass_spectrum_data_filtered, 
  ms1_window_left = 1, 
  ms1_window_right = 2, 
  bin_width = 0.05, 
  threshold_iso_probalility = 0.01
  )

# write.csv(final_glycan_spectrum_matching_result, file = 'your_csv_file_path_and_name.csv')
```

### STEP7:

``` r
ms2_spectrum_similarity_info = GetMS2SpectrumSimilarityScore(ms_data = mass_spectrum_data_filtered, 
                                         spectrum_matching_result = final_glycan_spectrum_matching_result, 
                                         glycan_composition_str = 'Hex3HexNAc4dHex1', 
                                         adduct_type = c(H = 2, Na = 0, K = 0), 
                                         bin_width = 0.3, 
                                         ms2_range_start = 100, 
                                         ms2_range_end = 2200) 
```
