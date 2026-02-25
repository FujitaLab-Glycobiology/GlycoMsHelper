
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GlycoMsHelper

<!-- badges: start -->
<!-- badges: end -->

The goal of GlycoMsHelper is to relieve you from labor-intensive
glycomic MS data analysis.

## Installation

You can install the development version of GlycoMsHelper from
[GitHub](https://github.com/FujitaLab-Glycobiology/GlycoMsHelper) with:

``` r
# install.packages("devtools")
devtools::install_github("FujitaLab-Glycobiology/GlycoMsHelper")
```

## Example

This is a basic example which shows you the basic workflow for
GlycoMsHelper:

### STEP1:

Load the package and read the `.mzML` file.

``` r
library(Spectra)
library(GlycoMsHelper)

# Read the MS file
mass_spectrum_data = Spectra::Spectra(your_ms_file_path, source = MsBackendMzR())
```

### STEP2:

Check the `.mzML` file.

``` r
GlycoMsHelper::MsFileChecker(mass_spectrum_data)
```

### STEP3:

Get the MS2 spectrum info likely to be glycan based on diagnostic ions.

``` r
# Optional: filter low quality MS1 and MS2 spectrums based on Spectra::spectraVariables()
qc_results = GlycoMsHelper::SpectrumQcFilter(
  ms_data = mass_spectrum_data, 
  filter_method_ms1 = c(peaksCount = 'mean_sd', totIonCurrent = 'quantile_prob', rtime = 'start_end'), 
  threshold_ms1 = list(peaksCount = 2, totIonCurrent = c(0.05, 1), rtime = c(8*60, 45*60)), 
  filter_method_ms2 = c(peaksCount = 'quantile_prob', totIonCurrent = 'quantile_prob', rtime = 'start_end'), 
  threshold_ms2 = list(peaksCount = c(0.1, 1), totIonCurrent = c(0.05, 1), rtime = c(8*60, 50*60)), 
  plot_option = T
  )
                          
mass_spectrum_data_filtered = qc_results$filtered_ms_data


# MS2 spectrum de-noising
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
  AB = 'C7H10N2', 
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
