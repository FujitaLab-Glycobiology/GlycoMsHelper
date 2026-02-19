library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)
library('Spectra')
library(segmented)
library(enviPat)

# source('SubFunctions.R')


#' Check MS2 Spectra and Precursor Charge Information
#'
#' This function validates whether the input Spectra object contains MS2 scans
#' and ensures that precursor charge information is available for downstream
#' glycan composition analysis.
#'
#' @param ms_data A \code{Spectra} object (typically loaded from an .mzML file by Spectra::Spectra()).
#'
#' @return The function returns nothing if the check passes, but issues a success
#'   message. It throws an error if MS2 spectra or charge states are missing.
#'
#' @importFrom Spectra filterMsLevel precursorCharge
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'raw_data' is a Spectra object
#' MsFileChecker(raw_data)
#' }
MsFileChecker = function(ms_data) {
  ms2_data = Spectra::filterMsLevel(ms_data, 2)

  if (length(ms2_data) == 0) {
    stop("The provided file contains no MS2 spectra.", call. = FALSE)
  }

  charges = Spectra::precursorCharge(ms2_data)

  if (all(is.na(charges)) || all(charges == 0)) {
    stop("Error: The .mzML file lacks 'precursorCharge' information.
         Downstream analysis (like find glycan composition) requires charge states.",
         call. = FALSE)
  } else {
    message("Check passed: MS2 charge state information is present.")
  }
}








# default setting
# define the glycan type
glycan_type_default = 'N_glycan'    # 'O_glycan', 'GPI', 'GSL'
# N_glycan library rules:
# 1. if there is Neu5Ac or HexA or Neu5Gc, Hex num >= 3
# 2. if there is Neu5Ac or HexA or Neu5Gc, HexNAc num >= 3
# 3. Hex, HexNAc, dHex, Neu5Ac, HexA are in the N-glycan library

# GSL library rules
# 1. if dHex > 0, Hex num >= 3
# 2. Hex num >= 1
# 3. if HexNAc > 0, Hex num >= 2


# define the derivatization type
derivatization_type_default = 'ProA'      #'AB', 'PA', 'user_defined'


# define the adduct type
adduct_type_default = c('H', 'Na', 'K')


# define the glycan library
# min <= mono <= max
min_Hex_num_default = 1;      max_Hex_num_default = 12
min_HexNAc_num_default = 2;   max_HexNAc_num_default = 10
min_dHex_num_default = 0;     max_dHex_num_default = 4
min_Neu5Ac_num_default = 0;    max_Neu5Ac_num_default = 4

min_HexA_num_default = 0;     max_HexA_num_default = 2;

min_Neu5Gc_num_default = 0;     max_Neu5Gc_num_default = 4

min_Pentose_num_default = 0;  max_Pentose_num_default = 0
min_KDN_num_default = 0;      max_KDN_num_default = 0

# define the total mono number
min_total_monosaccharides_num_default = 3
max_total_monosaccharides_num_default = 22


# define the charge state range
# currently only suppor the positive mode
# min_charge_state <= charge_state <= max_charge_state
min_charge_state_default = 1
max_charge_state_default = 3



# # define the glycan modification
# modification_list = c(
#                       Phosphate	79.9663
#                       Sulphate	79.9568
#                       Trifluoroacetic acid	112.9850391
#                       Amidation	74.04801282
#                       Methylation (CH2)	14.01565006414
# )
#' Construct a Glycan Theoretical Library
#'
#' This function generates a theoretical glycan library based on specified glycan types (N-glycan, GPI, etc.),
#' monosaccharide ranges, adduct types, and derivatization labels. It calculates monoisotopic weights
#' and m/z values for all possible combinations.
#'
#' @param glycan_type A character string specifying the glycan class.
#'   Supported: \code{'N_glycan'}, \code{'GPI'}, \code{'GSL'}, \code{'O_glycan'}.
#' @param min_charge_state Minimum charge state (positive integer).
#' @param max_charge_state Maximum charge state (positive integer).
#' @param monosaccharides_additional_customized_list A named numeric vector of
#'   customized monosaccharides and their monoisotopic weights.
#'   Example: \code{c(EtNP = 123.01, AHM = 161.05)}.
#' @param min_num_monosaccharides_additional_customized_list Named vector for minimum counts
#'   of customized monosaccharides.
#'   Example: \code{c(EtNP = 0, AHM = 1)}.
#' @param max_num_monosaccharides_additional_customized_list Named vector for maximum counts
#'   of customized monosaccharides.
#'   Example: \code{c(EtNP = 3, AHM = 1)}.
#' @param adduct_type Character vector of supported adducts (e.g., \code{c('H', 'Na', 'K')}).
#' @param adduct_monoisotopic_customized_list Named numeric vector for user-defined adduct weights.
#' @param adduct_monoisotopic_customized_charge_state_list Named numeric vector for
#'   user-defined adduct charge states.
#' @param derivatization_type Derivatization label type (e.g., \code{'ProA'}, \code{'AB'}, \code{'PA'}).
#' @param derivatization_monoisotopic_weight_customized Numeric weight for a custom label.
#' @param min_Hex_num,max_Hex_num Min/Max number of Hexose (Hex).
#' @param min_HexNAc_num,max_HexNAc_num Min/Max number of N-Acetylhexosamine (HexNAc).
#' @param min_dHex_num,max_dHex_num Min/Max number of Deoxyhexose (dHex).
#' @param min_Neu5Ac_num,max_Neu5Ac_num Min/Max number of N-Acetylneuraminic acid (Neu5Ac).
#' @param min_HexA_num,max_HexA_num Min/Max number of Hexuronic acid (HexA).
#' @param min_Neu5Gc_num,max_Neu5Gc_num Min/Max number of N-Glycolylneuraminic acid (Neu5Gc).
#' @param min_Pentose_num,max_Pentose_num Min/Max number of Pentose.
#' @param min_KDN_num,max_KDN_num Min/Max number of 2-keto-3-deoxy-D-glycero-D-galacto-nononic acid (KDN).
#' @param min_total_monosaccharides_num,max_total_monosaccharides_num Range for the sum of all monosaccharides.
#'
#' @return A list containing five elements:
#' \itemize{
#'   \item \code{monosaccharides_combination}: Data frame of valid monosaccharide compositions.
#'   \item \code{adduct_combination}: Data frame of valid adduct combinations and total charges.
#'   \item \code{monosaccharides_adduct_num}: Matrix of counts for each component.
#'   \item \code{monosaccharides_adduct_weight}: Matrix of weight contributions.
#'   \item \code{glycan_monosaccharides_library}: The full library with monoisotopic weights and m/z.
#' }
#'
#' @importFrom dplyr filter across where if_else mutate select bind_cols
#' @importFrom tidyr crossing
#' @export
#'
#' @examples
#' \dontrun{
#' library = ConstructGlycanLibrary(
#'   glycan_type = 'N_glycan',
#'   min_Hex_num = 3, max_Hex_num = 10,
#'   adduct_type = c('H', 'Na')
#' )
#' }
ConstructGlycanLibrary = function(glycan_type = glycan_type_default,
                                 min_charge_state = min_charge_state_default, max_charge_state = max_charge_state_default,

                             monosaccharides_additional_customized_list,
                             min_num_monosaccharides_additional_customized_list, max_num_monosaccharides_additional_customized_list,
                             # This option is used when the user wants to include monosaccharides that are NOT defined in the default 'monosaccharides_list'.
                             # Just pay attention, the name of the customized monosaccharides should be different from
                             # Hex, HexNAc, dHex, Neu5Ac, HexA, Neu5Gc, Pentose, KDN, if the customized monosaccharides name are same as Hex, HexNAc, dHex, Neu5Ac, HexA, Neu5Gc, Pentose, KDN
                             # the scripts will check to make sure that the molecular weight difference is smaller than 0.01, if not smaller than 0.01, warning will appear
                             # ex:
                             # monosaccharides_additional_customized_list = c(Hex_de-H2o = 144.0528)
                             # min_num_monosaccharides_additional_customized_list = c(Hex_de-H2o = 0)
                             # max_num_monosaccharides_additional_customized_list = c(Hex_de-H2o = 3)


                             # This option could also be used when 'glycan_type_default' is not one of 'N_glycan', 'O_glycan', 'GPI', 'GSL'
                             # For example, you want to only check the glycans that have Hex and HexNAc
                             # but attention here,if the name of the customized monosaccharides are
                             # Hex, HexNAc, dHex, Neu5Ac, HexA, Neu5Gc, Pentose, KDN,
                             # the scripts will check to make sure that the molecular weight difference between the user defined and the standard in this scripts is smaller than 0.01,
                             # if not smaller than 0.01, warning will appear
                             # ex:
                             # glycan_type = 'anything_you_want'
                             # monosaccharides_additional_customized_list = c(Hex = 162.0528, HexNAc = 203.0794)
                             # min_num_monosaccharides_additional_customized_list = c(Hex = 0, HexNAc = 0)
                             # max_num_monosaccharides_additional_customized_list = c(Hex = 2, HexNAc = 3)


                             # It can also be used to explicitly add monosaccharides that exist in
                             # 'monosaccharides_list' but are excluded by the selected 'glycan_type'.
                             # for this setting, 'min_Neu5Ac_num' and 'max_Neu5Ac_num' will be override,
                             # which means that you can NOT use 'min_Neu5Ac_num' and 'max_Neu5Ac_num' to define the monosaccharide,
                             # you have to use min_num_monosaccharides_additional_customized_list and max_num_monosaccharides_additional_customized_list
                             # ex:
                             # When glycan_type = "N_glycan", the default glycan library includes
                             # Hex, HexNAc, dHex, Neu5Ac, and HexA, along with other N-glycan-specific rules.

                             # If the user wants to additionally include Neu5Gc in the library,
                             # the following arguments can be used:
                             #   - monosaccharides_additional_customized_list
                             #   - min_num_monosaccharides_additional_customized_list
                             #   - max_num_monosaccharides_additional_customized_list

                             # example settings:
                             # monosaccharides_additional_customized_list = c(Neu5Gc = 307.0903)
                             # min_num_monosaccharides_additional_customized_list = c(Neu5Gc = 0)
                             # max_num_monosaccharides_additional_customized_list = c(Neu5Gc = 3)




                             adduct_type = adduct_type_default,
                             # 'adduct_type' could be H, Na, K, Li, Mg,
                             # if using H, Na, K, Li, Mg, just define them using 'adduct_type', no need to define by yourself
                             # for other adduct, please define using
                             # 'adduct_monoisotopic_customized_list' and 'adduct_monoisotopic_customized_charge_state_list'


                             adduct_monoisotopic_customized_list, adduct_monoisotopic_customized_charge_state_list,
                             # adduct_monoisotopic_customized_list = c(NH4 = 18, Ca = 40)
                             # adduct_monoisotopic_customized_charge_state_list = c(NH4 = 1, Ca = 2)




                             derivatization_type = derivatization_type_default,
                             derivatization_monoisotopic_weight_customized,
                             # 'derivatization_type' could be ProA, AB, PA and customized derivatization type,
                             # if using one of ProA, AB, PA, just use 'derivatization_type', no need to define by yourself
                             # if define the customized derivatization label, both 'derivatization_type' and 'derivatization_monoisotopic_weight_customized'
                             # should be provided
                             # RIGHT EXAMPLE:
                             # derivatization_type = 'customized'
                             # derivatization_monoisotopic_weight_customized = 100


                              min_Hex_num     = min_Hex_num_default,      max_Hex_num     = max_Hex_num_default,
                              min_HexNAc_num  = min_HexNAc_num_default,   max_HexNAc_num  = max_HexNAc_num_default,
                              min_dHex_num    = min_dHex_num_default,     max_dHex_num    = max_dHex_num_default,
                              min_Neu5Ac_num   = min_Neu5Ac_num_default,    max_Neu5Ac_num   = max_Neu5Ac_num_default,

                              min_HexA_num = min_HexA_num_default,        max_HexA_num = max_HexA_num_default,

                              min_Neu5Gc_num = min_Neu5Gc_num_default,      max_Neu5Gc_num = max_Neu5Gc_num_default,

                              min_Pentose_num = min_Pentose_num_default,  max_Pentose_num = max_Pentose_num_default,
                              min_KDN_num = min_KDN_num_default,          max_KDN_num = max_KDN_num_default,

                              min_total_monosaccharides_num = min_total_monosaccharides_num_default,
                              max_total_monosaccharides_num = max_total_monosaccharides_num_default
                              ) {

  # initialize variables
  adduct_list = c(H = 1.00727,
                  Na =22.989768,
                  K = 38.963707,
                  Li = 7.016003,
                  Mg = 23.985042)

  adduct_charge_state_list = c(H = 1,
                               Na = 1,
                               K = 1,
                               Li = 1,
                               Mg = 2)


  derivatization_list = c(ProA = 237.1841124,
                          AB = 122.0843983,
                          PA = 96.068761)


  monosaccharides_list = c(Hex = 162.0528,
                           HexNAc = 203.0794,
                           dHex = 146.0579,
                           Neu5Ac = 291.0954,
                           HexA = 176.03209,
                           Neu5Gc = 307.0903,
                           Pentose = 132.0423,
                           KDN = 250.0689
  )


  # check the monosaccharides number
  if (any(c(min_Hex_num, min_HexNAc_num, min_dHex_num, min_Neu5Ac_num,  min_HexA_num, min_Neu5Gc_num, min_Pentose_num, min_KDN_num, min_total_monosaccharides_num) < 0)) {
    stop("The following arguments must be >= 0:\n" ,
         "min_Hex_num, min_HexNAc_num, min_dHex_num, min_Neu5Ac_num,\n",
         "min_HexA_num, min_Neu5Gc_num, min_Pentose_num,\n",
         "min_KDN_num, min_total_monosaccharides_num.",
         call. = FALSE)
  }

  # check the monosaccharides_list
  has_customized_monosaccharides = !missing(monosaccharides_additional_customized_list) ||
    !missing(min_num_monosaccharides_additional_customized_list) ||
    !missing(max_num_monosaccharides_additional_customized_list)

  complete_customized_monosaccharides = !missing(monosaccharides_additional_customized_list) &&
    !missing(min_num_monosaccharides_additional_customized_list) &&
    !missing(max_num_monosaccharides_additional_customized_list)


  if (complete_customized_monosaccharides) {

    if (!setequal(names(monosaccharides_additional_customized_list), names(min_num_monosaccharides_additional_customized_list)) ||
        !setequal(names(monosaccharides_additional_customized_list), names(max_num_monosaccharides_additional_customized_list)) ||
        !(length(monosaccharides_additional_customized_list) == length(min_num_monosaccharides_additional_customized_list)) ||
        !(length(monosaccharides_additional_customized_list) == length(max_num_monosaccharides_additional_customized_list))) {

      stop("'monosaccharides_additional_customized_list', 'min_num_monosaccharides_additional_customized_list', ",
           "and 'max_num_monosaccharides_additional_customized_list' must match.",
           call. = FALSE)

    } else {

      monosaccharides_list[names(monosaccharides_additional_customized_list)] = monosaccharides_additional_customized_list

      }

  } else if (has_customized_monosaccharides && !complete_customized_monosaccharides) {

    stop("If providing customized monosaccharides, you must provide all of ",
         "'monosaccharides_additional_customized_list', ",
         "'min_num_monosaccharides_additional_customized_list', and ",
         "'max_num_monosaccharides_additional_customized_list'.",
         call. = FALSE)

  } else if (!has_customized_monosaccharides) {

    #monosaccharides_list = monosaccharides_list

  }


  # Check whether the monoisotopic weights match the reference list,
  # or whether their differences fall within a specified tolerance.
  standard_monosaccharides_list = c(Hex = 162.0528, HexNAc = 203.0794,
                                    dHex = 146.0579, Neu5Ac = 291.0954,
                                    HexA = 176.03209, Neu5Gc = 307.0903,
                                    Pentose = 132.0423, KDN = 250.0689)

  for (std_monos in intersect(names(monosaccharides_list), names(standard_monosaccharides_list))) {
    if (abs(monosaccharides_list[paste(std_monos)] - standard_monosaccharides_list[paste(std_monos)]) >= 0.01) {

      stop(paste("Monoisotopic weight for monosaccharide",
                 sprintf("'%s'", std_monos), " differs from the reference by >= 0.01 Da."),
           call. = FALSE)

    }

  }


  # define the glycan label
  if (derivatization_type %in% names(derivatization_list)) {

    if (!missing(derivatization_monoisotopic_weight_customized)) {

      stop("Argument 'derivatization_monoisotopic_weight_customized' is provided but 'derivatization_type' is not provided. ",
           call. = FALSE)
    }

    derivatization_monoisotopic_weight = derivatization_list[derivatization_type]

    } else {

      if (missing(derivatization_monoisotopic_weight_customized)) {

        stop("Argument 'derivatization_monoisotopic_weight_customized' is required but missing.",
             call. = FALSE)

      } else if(!is.numeric(derivatization_monoisotopic_weight_customized) ||
              derivatization_monoisotopic_weight_customized < 0 ||
              is.null(derivatization_monoisotopic_weight_customized)) {

        stop("Argument 'derivatization_monoisotopic_weight_customized' must be a single numeric value >= 0.",
             call. = FALSE)

      }

      derivatization_monoisotopic_weight = derivatization_monoisotopic_weight_customized
    }


  # define the adducts
  supported_adduct = names(adduct_list)

  all_valid_adduct = all(adduct_type %in% names(adduct_list))

  has_customized_adduct = !missing(adduct_monoisotopic_customized_list) ||
    !missing(adduct_monoisotopic_customized_charge_state_list)

  complete_customized_adduct = !missing(adduct_monoisotopic_customized_list) &&
    !missing(adduct_monoisotopic_customized_charge_state_list)


  if (!has_customized_adduct && all_valid_adduct) {

    adduct_list = adduct_list[adduct_type]
    adduct_charge_state_list = adduct_charge_state_list[adduct_type]

  } else if (!all_valid_adduct) {

    stop(sprintf("Invalid 'adduct_type' detected. Supported adducts are: %s. ", paste(supported_adduct, collapse = ", ")),
         "For customized adduct(s), please define it/them using ",
         "'adduct_monoisotopic_customized_list' and 'adduct_monoisotopic_customized_charge_state_list'.",
         call. = FALSE)

  } else if (has_customized_adduct && all_valid_adduct) {

    if (!complete_customized_adduct) {

      stop("If providing customized adduct(s), you must provide all of ",
           "'adduct_monoisotopic_customized_list' and 'adduct_monoisotopic_customized_charge_state_list'",
           call. = FALSE)

    } else {

      if (!setequal(names(adduct_monoisotopic_customized_list), names(adduct_monoisotopic_customized_charge_state_list)) ||
          !(length(adduct_monoisotopic_customized_list) == length(adduct_monoisotopic_customized_charge_state_list))) {

        stop("'adduct_monoisotopic_customized_list' and 'adduct_monoisotopic_customized_charge_state_list' MUST match.",
             call. = FALSE)

        } else if (!is.numeric(adduct_monoisotopic_customized_list) ||
                 !is.numeric(adduct_monoisotopic_customized_charge_state_list) ||
                 any(adduct_monoisotopic_customized_list < 0) ||
                 any(adduct_monoisotopic_customized_charge_state_list <= 0) ||
                 any(adduct_monoisotopic_customized_charge_state_list %% 1 != 0)
                 ) {

          stop("Argument 'adduct_monoisotopic_customized_list' must be numeric value(s) >= 0. ",
             "Argument 'adduct_monoisotopic_customized_charge_state_list' must be numeric integer(s) >= 1. ",
             call. = FALSE)

        } else if (any(names(adduct_monoisotopic_customized_list) %in% names(adduct_list)) ||
                   any(adduct_monoisotopic_customized_list %in% adduct_list)) {

          stop("Names or values in 'adduct_monoisotopic_customized_list' must NOT overlap with 'adduct_list'. ",
               "If want to add adduct like H, Na, K, Li, Mg, please use 'adduct_type'. ",
               call. = FALSE)

        } else {

          adduct_list = adduct_list[adduct_type]
          adduct_list[names(adduct_monoisotopic_customized_list)] = adduct_monoisotopic_customized_list

          adduct_charge_state_list = adduct_charge_state_list[adduct_type]
          adduct_charge_state_list[names(adduct_monoisotopic_customized_charge_state_list)] = adduct_monoisotopic_customized_charge_state_list

        }

    }

  }


  # create the glycan library
  default_minmax_map = list(
    Hex     = c(min_Hex_num,     max_Hex_num),
    HexNAc  = c(min_HexNAc_num,  max_HexNAc_num),
    dHex    = c(min_dHex_num,    max_dHex_num),
    Neu5Ac  = c(min_Neu5Ac_num,  max_Neu5Ac_num),
    HexA    = c(min_HexA_num,    max_HexA_num),
    Neu5Gc  = c(min_Neu5Gc_num,  max_Neu5Gc_num),
    Pentose = c(min_Pentose_num, max_Pentose_num),
    KDN     = c(min_KDN_num,     max_KDN_num)
  )

  # define the glycan type(glycan library)
  if (glycan_type == 'N_glycan') {

    N_glycan_std_monosaccharides = c('Hex', 'HexNAc', 'dHex', 'Neu5Ac', 'HexA', 'Neu5Gc')

    if (complete_customized_monosaccharides) {

      specific_glycan_monosaccharides_list = monosaccharides_additional_customized_list
      specific_glycan_monosaccharides_list[N_glycan_std_monosaccharides] = monosaccharides_list[N_glycan_std_monosaccharides]

      monosaccharides_combos = MakeAllMonosCombos(specific_glycan_monosaccharides_list,
                                               min_num_custom = min_num_monosaccharides_additional_customized_list,
                                               max_num_custom = max_num_monosaccharides_additional_customized_list,
                                               minmax_map = default_minmax_map)

    } else {
      # create list of sequences for expand.grid
      seq_list = lapply(N_glycan_std_monosaccharides, function(x) {
        seq(from = default_minmax_map[[x]][1], to = default_minmax_map[[x]][2], by = 1)
      })

      for (i in N_glycan_std_monosaccharides) {
        if (default_minmax_map[[i]][1] > default_minmax_map[[i]][2]) {
          stop(sprintf("Invalid range (min > max) for monosaccharide(s): %s.", i),
               call. = FALSE)
        }
      }
      names(seq_list) = paste0(N_glycan_std_monosaccharides)

      monosaccharides_combos = expand.grid(seq_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      specific_glycan_monosaccharides_list = monosaccharides_list[N_glycan_std_monosaccharides]
    }

    monosaccharides_combos = dplyr::filter(monosaccharides_combos,
                                           rowSums(dplyr::across(dplyr::where(is.numeric))) >= min_total_monosaccharides_num & rowSums(dplyr::across(dplyr::where(is.numeric))) <= max_total_monosaccharides_num)

    # apply N-glycan rules
    if ('Neu5Gc' %in% colnames(monosaccharides_combos)) {
      monosaccharides_combos = dplyr::filter(monosaccharides_combos,
                                             dplyr::if_else(Neu5Ac > 0 | HexA > 0 | Neu5Gc > 0,
                                                            Hex >= 3 & HexNAc >= 3,
                                                            TRUE)
      )
    } else {
      monosaccharides_combos = dplyr::filter(
        monosaccharides_combos, dplyr::if_else(Neu5Ac > 0 | HexA > 0, Hex >= 3 & HexNAc >= 3, TRUE)
      )
    }


    if (dim(monosaccharides_combos)[1] <= 0) {
      stop("Can NOT find any monosaccharides combinations that fullfill the criteria of N-glycan. ",
           "N-glycan rules are: ",
           "1. If Neu5Ac or HexA or Neu5Gc num >= 1, Hex num >= 3 and HexNAc num >= 3",
           call. = FALSE)
    }


  } else if (glycan_type == 'O_glycan') {

  } else if (glycan_type == 'GPI') {

    GPI_std_monosaccharides = c('Hex', 'HexNAc', 'Neu5Ac', 'Neu5Gc')

    if (complete_customized_monosaccharides) {

      specific_glycan_monosaccharides_list = monosaccharides_additional_customized_list
      specific_glycan_monosaccharides_list[GPI_std_monosaccharides] = monosaccharides_list[GPI_std_monosaccharides]

      monosaccharides_combos = MakeAllMonosCombos(specific_glycan_monosaccharides_list,
                                                  min_num_custom = min_num_monosaccharides_additional_customized_list,
                                                  max_num_custom = max_num_monosaccharides_additional_customized_list,
                                                  minmax_map = default_minmax_map)
    } else {

      stop(
        "For GPI glycan library construction, 'EtNP' and 'AHM' must be explicitly defined in 'monosaccharides_additional_customized_list'.\n",
        "Please provide:\n",
        "  monosaccharides_additional_customized_list = c(EtNP = __, AHM = __)\n",
        "  min_num_monosaccharides_additional_customized_list = c(EtNP = __, AHM = __)\n",
        "  max_num_monosaccharides_additional_customized_list = ...\n",
        call. = FALSE
      )

    }

    monosaccharides_combos = dplyr::filter(monosaccharides_combos,
                                           rowSums(dplyr::across(dplyr::where(is.numeric))) >= min_total_monosaccharides_num & rowSums(dplyr::across(dplyr::where(is.numeric))) <= max_total_monosaccharides_num)

    # apply GPI rules
    if ('Neu5Gc' %in% colnames(monosaccharides_combos)) {
      monosaccharides_combos = dplyr::filter(monosaccharides_combos, dplyr::if_else(Neu5Ac > 0 | Neu5Gc > 0, Hex >= 4, TRUE)) |>
        dplyr::filter(dplyr::if_else(EtNP > 0, Hex >= EtNP, TRUE)) |>
        dplyr::filter(dplyr::if_else(Hex == 5, HexNAc == 1, TRUE))

    } else {
      monosaccharides_combos = dplyr::filter(monosaccharides_combos, dplyr::if_else(Neu5Ac > 0, Hex >= 4, TRUE)) |>
        dplyr::filter(dplyr::if_else(EtNP > 0, Hex >= EtNP, TRUE)) |>
        dplyr::filter(dplyr::if_else(Hex == 5, HexNAc == 1, TRUE))
    }

    if (dim(monosaccharides_combos)[1] <= 0) {
      stop("Can NOT find any monosaccharides combinations that fullfill the criteria of GPI. ",
           "GPI rules are: ",
           "1. Hex num >= EtNP",
           "2. If Hex num == 5, HexNAc num must = 1",
           "3. If Neu5Ac/Neu5Gc num == 1, Hex num >= 4",
           call. = FALSE)
    }


  } else if (glycan_type == 'GSL') {

    GSL_std_monosaccharides = c('Hex', 'HexNAc', 'dHex', 'Neu5Ac', 'Neu5Gc')

    if (complete_customized_monosaccharides) {

      specific_glycan_monosaccharides_list = monosaccharides_additional_customized_list
      specific_glycan_monosaccharides_list[GSL_std_monosaccharides] = monosaccharides_list[GSL_std_monosaccharides]

      monosaccharides_combos = MakeAllMonosCombos(specific_glycan_monosaccharides_list,
                                                  min_num_custom = min_num_monosaccharides_additional_customized_list,
                                                  max_num_custom = max_num_monosaccharides_additional_customized_list,
                                                  minmax_map = default_minmax_map)
    } else {
      # create list of sequences for expand.grid
      seq_list = lapply(GSL_std_monosaccharides, function(x) {
        seq(from = default_minmax_map[[x]][1], to = default_minmax_map[[x]][2], by = 1)
      })

      for (i in GSL_std_monosaccharides) {
        if (default_minmax_map[[i]][1] > default_minmax_map[[i]][2]) {
          stop(sprintf("Invalid range (min > max) for monosaccharide(s): %s.", i),
               call. = FALSE)
        }
      }
      names(seq_list) = paste0(GSL_std_monosaccharides)

      monosaccharides_combos = expand.grid(seq_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      specific_glycan_monosaccharides_list = monosaccharides_list[GSL_std_monosaccharides]
    }

    monosaccharides_combos = dplyr::filter(monosaccharides_combos,
                                           rowSums(dplyr::across(dplyr::where(is.numeric))) >= min_total_monosaccharides_num & rowSums(dplyr::across(dplyr::where(is.numeric))) <= max_total_monosaccharides_num)

    # apply GSL rules
    monosaccharides_combos = dplyr::filter(monosaccharides_combos, dplyr::if_else(dHex >= 1, Hex >= 3, TRUE)) |>
      dplyr::filter(dplyr::if_else(HexNAc >= 3, Hex >= 2, TRUE))

    if (dim(monosaccharides_combos)[1] <= 0) {
      stop("Can NOT find any monosaccharides combinations that fullfill the criteria of GSL. ",
           "GSL rules are: ",
           "1. If dHex num >= 1, HexNAc num >= 3",
           "3. If HexNAc num >= 1, Hex num >= 2",
           call. = FALSE)
    }


  } else {
    # warning("Argument 'glycan_type' is missing, use complete library.",
    #         call. = FALSE)
  }


  # check the charge state
  if (min_charge_state > max_charge_state) {
    stop(sprintf("Invalid range (min > max) for charge state: min_charge_state = %s, max_charge_state = %s.", min_charge_state, max_charge_state),
         call. = FALSE)
  }

  # create the adduct library
  # get the adduct combos
  temp_adduct_combos <- lapply(adduct_charge_state_list, function(x) 0:floor(max_charge_state / x))
  adduct_combos_all_grids <- expand.grid(temp_adduct_combos)


  adduct_matrix <- as.matrix(adduct_combos_all_grids)
  adduct_charges <- as.numeric(adduct_charge_state_list[colnames(adduct_matrix)])
  total_charges <- adduct_matrix %*% adduct_charges

  adduct_combos_all_grids <- adduct_combos_all_grids[total_charges >= min_charge_state & total_charges <= max_charge_state, ]

  # # get the number of monos and adducts
  # mono_adduct_combos_full_num <- tidyr::crossing(monosaccharides_combos, adduct_combos_all_grids)
  # mono_adduct_combos_full_num[[derivatization_type]] <- 1


  adduct_combos_all_grids$total_charge <- total_charges[total_charges >= min_charge_state & total_charges <= max_charge_state]

  # get the monos and adducts combos
  mono_adduct_combos_full <- tidyr::crossing(monosaccharides_combos, adduct_combos_all_grids)

  weight_vector <- unlist(c(specific_glycan_monosaccharides_list, adduct_list))
  # check the columns order
  relevant_cols <- names(weight_vector)
  weight_matrix <- as.matrix(mono_adduct_combos_full[, relevant_cols])

  mono_adduct_combos_full$glycan_monoisotopic_weight <-
    as.numeric(weight_matrix %*% weight_vector) + derivatization_monoisotopic_weight

  mono_adduct_combos_full[[derivatization_type]] <- 1

  # glycan_monosaccharides_lib <- dplyr::mutate(mono_adduct_combos_full, glycan_monoisotopic_mz = glycan_monoisotopic_weight/total_charge)
  glycan_monosaccharides_lib = mono_adduct_combos_full

  # get the monoisotopic weight
  full_weight_vector <- c(weight_vector, "ProA" = unname(derivatization_monoisotopic_weight))

  calc_cols <- names(full_weight_vector)

  # matrix calculation
  weight_contribution_matrix <- as.matrix(glycan_monosaccharides_lib[, calc_cols]) %*% diag(full_weight_vector)

  mono_adduct_weight <- as.data.frame(weight_contribution_matrix)
  colnames(mono_adduct_weight) <- calc_cols

  mono_adduct_weight <- dplyr::bind_cols(mono_adduct_weight,
                                         dplyr::select(glycan_monosaccharides_lib, total_charge, glycan_monoisotopic_weight)
  )


  glycan_monosaccharides_lib <- dplyr::mutate(mono_adduct_combos_full, glycan_monoisotopic_mz = glycan_monoisotopic_weight/total_charge)




  mono_adduct_combos_full_num = dplyr::select(glycan_monosaccharides_lib, names(full_weight_vector))




  message(
    paste0(
      'Constructed ', glycan_type, ' library with ',
      nrow(monosaccharides_combos), ' monosaccharide combinations. \n',
      'Monosaccharides: ', paste(colnames(monosaccharides_combos), collapse = ", "), '. ',
      'Labeling: ', derivatization_type, '. \n',
      'Adducts: ', paste(names(adduct_list), collapse = ", "), '.'
    )
  )


  return(list(
    monosaccharides_combination = monosaccharides_combos,
    adduct_combination = adduct_combos_all_grids,
    monosaccharides_adduct_num = mono_adduct_combos_full_num,
    monosaccharides_adduct_weight = mono_adduct_weight,
    glycan_monosaccharides_library = glycan_monosaccharides_lib
    )
  )
}


# res_list = ConsructGlycanLibrary(min_charge_state = min_charge_state_default,
#
#                                   max_charge_state = 2, min_Hex_num     = 0,
#
#                                   monosaccharides_additional_customized_list = c(Hex_de_H2o = 144.0528, Neu5Gc = 307.0903),
#
#                                   min_num_monosaccharides_additional_customized_list=c(Hex_de_H2o = 0, Neu5Gc = 0),
#                                   max_num_monosaccharides_additional_customized_list=c(Hex_de_H2o = 3, Neu5Gc = 4),
#
#                                  adduct_type = c('H'),
#
#                                   adduct_monoisotopic_customized_list = c(NH4 = 18, Ca = 40),
#                                   adduct_monoisotopic_customized_charge_state_list = c(NH4 = 1, Ca = 2))
#
# monosaccharides_adduct_num = res_list$monosaccharides_adduct_combination
#
# glycan_library = res_list$glycan_monosaccharides_library
#
# hist(glycan_library$glycan_monoisotopic_weight,
#      breaks = 50,
#      main = "Distribution of Glycan Mass",
#      xlab = "Monoisotopic Weight",
#      col = "lightblue")
#


















# filter low quality MS1 and MS2 spectrums based on Spectra::spectraVariables()
# filter_method_ms1, filter_method_ms2: mean_sd, quantile_prob, start_end
# attention, for the 'mean_sd', it gives the ranges from the jump point to the max of that varible
# ex:
# filter_method_ms1 = c(peaksCount = 'mean_sd', totIonCurrent = 'mean_sd', rtime = 'start_end')
# threshold_ms1 = c(peaksCount = 2, totIonCurrent = 2, rtime  =c(10, 45))
# > Spectra::spectraVariables(ms_data)
# [1] "msLevel"                  "rtime"                    "acquisitionNum"           "scanIndex"                "dataStorage"
# [6] "dataOrigin"               "centroided"               "smoothed"                 "polarity"                 "precScanNum"
# [11] "precursorMz"              "precursorIntensity"       "precursorCharge"          "collisionEnergy"          "isolationWindowLowerMz"
# [16] "isolationWindowTargetMz"  "isolationWindowUpperMz"   "peaksCount"               "totIonCurrent"            "basePeakMZ"
# [21] "basePeakIntensity"        "ionisationEnergy"         "lowMZ"                    "highMZ"                   "mergedScan"
# [26] "mergedResultScanNum"      "mergedResultStartScanNum" "mergedResultEndScanNum"   "injectionTime"            "filterString"
# [31] "spectrumId"               "ionMobilityDriftTime"     "scanWindowLowerLimit"     "scanWindowUpperLimit"

filter_method_ms1_default = c(peaksCount = 'mean_sd', totIonCurrent = 'quantile_prob', rtime = 'start_end')
filter_method_ms2_default = c(peaksCount = 'quantile_prob', totIonCurrent = 'quantile_prob', rtime = 'start_end')

threshold_ms1_default = NULL
threshold_ms2_default = NULL
# threshold_ms1 = list(peaksCount = 2, totIonCurrent = c(0.1, 1), rtime = c(10*60, 45*60))
# threshold_ms2 = list(peaksCount = c(0.1, 1), totIonCurrent = c(0.1, 1), rtime = c(10*60, 45*60))

plot_option_default = T



#' Quality Control Filtering for MS1 and MS2 Spectra
#'
#' This function filters raw mass spectrometry data by analyzing spectra variables
#' (such as peak count, total ion current, and retention time). It supports three
#' different filtering logic modes and provides visual feedback of the filtering process.
#'
#' @param ms_data A \code{Spectra} object containing the mass spectrometry data.
#' @param filter_method_ms1 A named character vector specifying the filtering
#'   method for MS1. Names must match \code{Spectra::spectraVariables(ms_data)}.
#'   Methods: \code{'mean_sd'}, \code{'quantile_prob'}, or \code{'start_end'}.
#' @param filter_method_ms2 A named character vector specifying the filtering
#'   method for MS2. Similar to \code{filter_method_ms1}.
#' @param threshold_ms1 A named list of thresholds for MS1.
#'   \itemize{
#'     \item For \code{'mean_sd'}: A single numeric value (n times SD).
#'     \item For \code{'quantile_prob'}: A numeric vector of length 2 (lower and upper probabilities).
#'     \item For \code{'start_end'}: A numeric vector of length 2 (min and max values).
#'   }
#' @param threshold_ms2 A named list of thresholds for MS2. Similar to \code{threshold_ms1}.
#' @param plot_option Logical, if \code{TRUE} (default), returns diagnostic plots
#'   showing the distribution of variables before and after filtering.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{filtered_ms_data}: The filtered \code{Spectra} object.
#'   \item \code{pic_ms_varibles_unfiltered}: A list of ggplot objects showing raw distributions.
#'   \item \code{pic_ms_varibles_filtered}: (Optional) A list of ggplot objects showing distributions after filtering.
#' }
#'
#' @details
#' The function calculates thresholds dynamically if not provided.
#' For example, \code{'mean_sd'} identifies a "jump point" based on mean and standard deviation
#' to remove low-quality scans at the beginning or end of an LC-MS run.
#'
#' @importFrom Spectra spectraVariables filterMsLevel
#' @importFrom utils capture.output
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' filter_methods <- c(peaksCount = 'mean_sd', totIonCurrent = 'quantile_prob')
#' thresholds <- list(peaksCount = 2, totIonCurrent = c(0.1, 0.9))
#'
#' qc_result <- SpectrumQcFilter(
#'   ms_data = raw_spectra,
#'   filter_method_ms1 = filter_methods,
#'   threshold_ms1 = thresholds
#' )
#'
#' # Access filtered data
#' clean_data <- qc_result$filtered_ms_data
#' }
SpectrumQcFilter = function(ms_data,
                         filter_method_ms1 = filter_method_ms1_default,
                         filter_method_ms2 = filter_method_ms2_default,
                         threshold_ms1 = threshold_ms1_default,
                         threshold_ms2 = threshold_ms2_default,
                         plot_option = plot_option_default) {

  if (#!missing(filter_method_ms1) &&
      !is.null(filter_method_ms1)) {
    invalid_vars = setdiff(names(filter_method_ms1), Spectra::spectraVariables(ms_data[ms_data$msLevel == 1]))
    if (length(invalid_vars) > 0) {
      stop(sprintf(
        "Invalid filter_method_ms1 variable(s): %s. Allowed variables are: %s.",
        paste(invalid_vars, collapse = ", "),
        paste(Spectra::spectraVariables(ms_data), collapse = ", ")
      ), call. = FALSE)
    }
  } else {
    stop("Argument 'filter_method_ms1' is required and must not be NULL. ",
         call. = FALSE)
  }

  if (#!missing(filter_method_ms2) &&
      !is.null(filter_method_ms2)) {
    invalid_vars = setdiff(names(filter_method_ms2), Spectra::spectraVariables(ms_data[ms_data$msLevel == 2]))
    if (length(invalid_vars) > 0) {
      stop(sprintf(
        "Invalid filter_method_ms2 variable(s): %s. Allowed variables are: %s.",
        paste(invalid_vars, collapse = ", "),
        paste(Spectra::spectraVariables(ms_data), collapse = ", ")
      ), call. = FALSE)
    }
  } else {
    stop("Argument 'filter_method_ms2' is required and must not be NULL. ",
         call. = FALSE)
  }
















  # default plot
  pic_ms1_var_unfiltered = list()
  pic_ms2_var_unfiltered = list()

  # ms1
  pic_ms1_var_unfiltered = PlotUnfilteredMsVaribles(msdata = ms_data, ms_level = 1, filter_method = filter_method_ms1)

  # ms2
  pic_ms2_var_unfiltered = PlotUnfilteredMsVaribles(msdata = ms_data, ms_level = 2, filter_method = filter_method_ms2)


  pic_ms_var_unfiltered = c(pic_ms1_var_unfiltered, pic_ms2_var_unfiltered)
  #
  # if (length(pic_ms_var_unfiltered) > 0) {
  #   return(pic_ms_var_unfiltered)
  # }



  # get threshold and do filter
  ms_clean = ms_data


  # ms1
  threshold_value_ms1 = GetQcFilterThreshold(msdata = ms_data, ms_level = 1,
                                             filter_method = filter_method_ms1, filter_threshold = threshold_ms1)

  # ms2
  threshold_value_ms2 = GetQcFilterThreshold(msdata = ms_data, ms_level = 2,
                                             filter_method = filter_method_ms2, filter_threshold = threshold_ms2)

  # get ms data
  # test = Spectra::filterRanges(spectra_data, spectraVariables = c("msLevel", "peaksCount", "rtime"),
  #                              ranges = c(1,1, 10903,16495, 480,2700),
  #                              match = "all")
  ms_clean_temp = GetQcFilteredMsData(msdata = ms_data, ms_level = 1, threshold_value = threshold_value_ms1)

  ms_clean = GetQcFilteredMsData(msdata = ms_clean_temp, ms_level = 2, threshold_value = threshold_value_ms2)

  # ms_clean = Spectra::filterRanges(ms_data,
  #                                          spectraVariables = c("msLevel", names(filter_method_ms1)),
  #                                          ranges = c(1, 1, threshold_value_ms1))
  # ms_clean = Spectra::filterRanges(ms_data,
  #                                          spectraVariables = c("msLevel", names(filter_method_ms2)),
  #                                          ranges = c(2, 2, threshold_value_ms2))


  # plot
  # ms1
  pic_ms1_var_filtered = PlotFilteredMsVaribles(msdata = ms_data, ms_level = 1, filter_threshold_value = threshold_value_ms1)

  #ms2
  pic_ms2_var_filtered = PlotFilteredMsVaribles(msdata = ms_data, ms_level = 2, filter_threshold_value = threshold_value_ms2)


  pic_ms_var_filtered = c(pic_ms1_var_filtered, pic_ms2_var_filtered)



  message(
    paste0(
      "Threshold for MS1 spectrum: \n",
      paste(utils::capture.output(print(threshold_value_ms1)), collapse = "\n"), "\n",
      "Threshold for MS2 spectrum: \n",
      paste(utils::capture.output(print(threshold_value_ms2)), collapse = "\n"), "\n",
      "Original MS1 spectrum number: ", length(Spectra::filterMsLevel(ms_data, 1)),
      "  Filtered MS1 spectrum number: ", length(Spectra::filterMsLevel(ms_clean, 1)), "\n",
      "Original MS2 spectrum number: ", length(Spectra::filterMsLevel(ms_data, 2)),
      "  Filtered MS2 spectrum number: ", length(Spectra::filterMsLevel(ms_clean, 2)), "\n"
    )
  )


  if (plot_option == T) {
    return(list(filtered_ms_data = ms_clean,
                pic_ms_varibles_filtered = pic_ms_var_filtered,
                pic_ms_varibles_unfiltered = pic_ms_var_unfiltered)
    )
  } else {
    return(list(filtered_ms_data = ms_clean,
                pic_ms_varibles_unfiltered = pic_ms_var_unfiltered)
    )
  }

}

















# When using 'spline_segmentation_regression',
# the script automatically determines whether the glycan in a given spectrum is completely fragmented.
# This judgment is based on the following assumption: if fragmentation is incomplete,
# the top 5 highest intensity peaks in the spectrum should correspond to isotopic peaks, meaning they have a regular m/z spacing pattern.
# For the incomplete fragmentation spectrum, spline regression is used.
# For the complete fragmentation spectrum, segmentated regression is used.

# spectrum_transform_method = c('log2_transform', 'asinh_transform', 'non_transform', function(z) log(z+1))
#                                log2(z+1)          asinh(z)              z
#
# denoising_method = c('spline_segmentation_regression', 'spline_regression', 'segmentation_regression', 'quantile_prob', 'fixed_value')
#
# segmentated_non_linear_transform_fun = c('square_transform', 'exponential_transform', function(z) z^2+z)
#                                                z^2                  2^x


ms2_denoising_detail_default = list(spline_segmentation_regression = list(spar_start = -1.5, spar_end = 1.5, spar_step = 0.02,
                                                                      RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1,
                                                                      use_cv = T, top_n_to_remove = 5,
                                                                      segmentated_non_linear_transform_fun = function(z) z^2+z
                                                                      ),

                                spline_regression = list(spar_start = -1.5, spar_end = 1.5, spar_step = 0.02,
                                                         RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1,
                                                         use_cv = T, top_n_to_remove = 5
                                                         ),

                                segmentation_regression = list(segmentated_non_linear_transform_fun = function(z) z^2+z
                                                               ),

                                quantile_prob = 0.05,

                                fixed_value = 30 # define the noise threshold
)




#' MS2 Spectrum Denoising using Regression or Thresholding Methods
#'
#' This function performs denoising on MS2 spectra using various methods, including
#' spline regression, segmented regression, and adaptive thresholding. It is
#' specifically designed to distinguish between complete and incomplete glycan
#' fragmentation to apply the most appropriate denoising strategy.
#'
#' @param ms_data A \code{Spectra} object containing MS1 and MS2 data.
#' @param ms2_spectrum_transform_method Character string or function.
#'   Methods include \code{'log2_transform'} (default), \code{'asinh_transform'},
#'   or \code{'non_transform'}. A custom function can also be provided.
#' @param ms2_denoising_method Character string specifying the denoising logic.
#'   Options: \code{'spline_segmentation_regression'}, \code{'spline_regression'},
#'   \code{'segmentation_regression'}, \code{'quantile_prob'}, or \code{'fixed_value'}.
#' @param ms2_denoising_detail A nested list containing parameters for the chosen
#'   denoising method. See details for default values and structure.
#'
#' @details
#' The \code{'spline_segmentation_regression'} method automatically detects the
#' fragmentation state. If fragmentation is incomplete (indicated by regular
#' m/z spacing of top peaks), spline regression is used. Otherwise, segmented
#' regression is applied.
#'
#'
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{denoising_regression_info}: A data frame (\code{tibble}) with
#'     regression statistics and the calculated threshold for each MS2 spectrum.
#'   \item \code{denoised_ms_data}: A \code{Spectra} object with the denoising
#'     processing added via \code{Spectra::addProcessing}.
#' }
#'
#' @importFrom Spectra filterMsLevel peaksData addProcessing
#' @importFrom dplyr distinct bind_rows
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Use default settings (Spline + Segmentation)
#' result <- MS2SpectrumDenoising(ms_data = raw_data)
#'
#' # Use a fixed intensity threshold for denoising
#' custom_detail <- ms2_denoising_detail_default
#' custom_detail$fixed_value <- 100
#' result_fixed <- MS2SpectrumDenoising(
#'   ms_data = raw_data,
#'   ms2_denoising_method = 'fixed_value',
#'   ms2_denoising_detail = custom_detail
#' )
#' }
MS2SpectrumDenoising = function(ms_data,
                                ms2_spectrum_transform_method = 'log2_transform',
                                ms2_denoising_method = 'spline_segmentation_regression',
                                ms2_denoising_detail = ms2_denoising_detail_default) {

  # checking the ms2_spectrum_transform_method
  if (is.character(ms2_spectrum_transform_method)) {
    ms2_spectrum_transform_fun = switch(ms2_spectrum_transform_method,
                           log2_transform = function(z) log2(z+1),
                           asinh_transform = function(z) asinh(z),
                           non_transform = function(z) z,
                           stop("Unknown method"))
  } else if (is.function(ms2_spectrum_transform_method)) {
    ms2_spectrum_transform_fun = ms2_spectrum_transform_method
  } else {
    stop("'ms2_spectrum_transform_method' must be a function or character", call. = FALSE)
  }

  # check the denoising
  if (!(ms2_denoising_method %in% names(ms2_denoising_detail))) {
    stop(sprintf(
      "Unsupported 'ms2_denoising_method': '%s'. Supported methods are: %s",
      ms2_denoising_method,
      paste(names(ms2_denoising_detail), collapse = ", ")
    ), call. = FALSE)
  }

  spectra_data_ms1 = Spectra::filterMsLevel(ms_data, 1)
  spectra_data_ms2 = Spectra::filterMsLevel(ms_data, 2)


  all_ms2_peaks_data = Spectra::peaksData(spectra_data_ms2)
  all_ms2_spectra_id = spectra_data_ms2[['spectrumId']]
  all_ms2_spectra_rt = spectra_data_ms2[['rtime']]


  regression_info_list = list()

  ms2_spectra_num = length(spectra_data_ms2)

  pb = utils::txtProgressBar(min = 0, max = ms2_spectra_num, style = 3)
  on.exit(close(pb))

  for (i in seq_along(spectra_data_ms2)) {

    # progress bar
    utils::setTxtProgressBar(pb, i)

    current_ms2_peaks_data = as.data.frame(all_ms2_peaks_data[[i]])
    current_ms2_peaks_data_distinct = dplyr::distinct(current_ms2_peaks_data, intensity, .keep_all = TRUE)

    current_ms2_spectra_id = all_ms2_spectra_id[i]
    current_ms2_spectra_rt = all_ms2_spectra_rt[i]

    if (ms2_denoising_method == 'spline_segmentation_regression') {

      denoising_info = GetSplineSegmentationNoise(denoising_detail = ms2_denoising_detail,
                                                  denoising_method = ms2_denoising_method,
                                                  transform_fun = ms2_spectrum_transform_fun,
                                                  spectra_distinct = current_ms2_peaks_data_distinct,
                                                  ms_id = current_ms2_spectra_id,
                                                  ms_rt = current_ms2_spectra_rt
      )

      thres_val = denoising_info$threshold_value
      new_row_spectrum_info = denoising_info$row_regression_info

      } else if (ms2_denoising_method == 'spline_regression') {

        denoising_info = GetSplineNoise(denoising_detail = ms2_denoising_detail,
                                                    denoising_method = ms2_denoising_method,
                                                    transform_fun = ms2_spectrum_transform_fun,
                                                    spectra_distinct = current_ms2_peaks_data_distinct,
                                                    ms_id = current_ms2_spectra_id,
                                                    ms_rt = current_ms2_spectra_rt
        )

        thres_val = denoising_info$threshold_value
        new_row_spectrum_info = denoising_info$row_regression_info

    } else if (ms2_denoising_method == 'segmentation_regression ') {

      denoising_info = GetSegmentationNoise(denoising_detail = ms2_denoising_detail,
                                                  denoising_method = ms2_denoising_method,
                                                  transform_fun = ms2_spectrum_transform_fun,
                                                  spectra_distinct = current_ms2_peaks_data_distinct,
                                                  ms_id = current_ms2_spectra_id,
                                                  ms_rt = current_ms2_spectra_rt
      )

      thres_val = denoising_info$threshold_value
      new_row_spectrum_info = denoising_info$row_regression_info

    } else if (ms2_denoising_method == 'quantile_prob') {

      denoising_info = GetQuantileNoise(denoising_detail = ms2_denoising_detail,
                                                  denoising_method = ms2_denoising_method,
                                                  transform_fun = ms2_spectrum_transform_fun,
                                                  spectra_distinct = current_ms2_peaks_data_distinct,
                                                  ms_id = current_ms2_spectra_id,
                                                  ms_rt = current_ms2_spectra_rt
      )

      thres_val = denoising_info$threshold_value
      new_row_spectrum_info = denoising_info$row_regression_info

    } else if (ms2_denoising_method == 'fixed_value') {

      denoising_info = GetFixedNoise(denoising_detail = ms2_denoising_detail,
                                                  denoising_method = ms2_denoising_method,
                                                  transform_fun = ms2_spectrum_transform_fun,
                                                  spectra_distinct = current_ms2_peaks_data_distinct,
                                                  ms_id = current_ms2_spectra_id,
                                                  ms_rt = current_ms2_spectra_rt
      )

      thres_val = denoising_info$threshold_value
      new_row_spectrum_info = denoising_info$row_regression_info
      }

    regression_info_list[[length(regression_info_list) + 1]] = new_row_spectrum_info

  }

  # close progress bar
  close(pb)

  reg_info = dplyr::bind_rows(regression_info_list)


  # function for the addprocess to get the denoised ms file
  MS2NoiseFilter = function(x, spectrumMsLevel, spectrumId, threshold_lookup) {
    if (spectrumMsLevel != 2L) {
      return(x)
    }

    current_thres = threshold_lookup$threshold_value[threshold_lookup$ms2_spectrum_id == spectrumId]

    if (length(current_thres) == 0 || is.na(current_thres)) {
      return(x)
    }

    if (is.null(dim(x)) || nrow(x) == 0) {
      return(x)
    }

    if (!is.null(colnames(x)) &&
        "intensity" %in% colnames(x)) {

      keep = x[, "intensity"] > current_thres
      x = x[keep, , drop = FALSE]

    } else {
      return(x)
    }

    return(x)

  }

  sps_filtered = Spectra::addProcessing(ms_data,
                                        MS2NoiseFilter,
                                        threshold_lookup = reg_info, spectraVariables = c("msLevel", "spectrumId")
                                        )

  return(
    list(
      denoising_regression_info = reg_info,
      denoised_ms_data = sps_filtered
    )
  )


}




















# define the glycan diagnostics fragments
diagnostic_frags_list_default = c(HexNAc =         204.08667,   # HexNAc
                                 HexNAc_ProA =        441.2708,    # HexNAc + ProA
                                 dHex_HexNAc_ProA =   587.3287,    # HexNAc + Deoxyhexose + ProA
                                 HexNAc_HexNAc_ProA = 644.3502,    # HexNAc + HexNAc + ProA
                                 Hex_ProA = 400.2442,
                                Hex_Hex_ProA = 562.2970,
                                  dHex_Hex_ProA = 546.3021,
                                 Hex =                163.06007,   # Hex
                                 Bi_HexNAc =          407.16607,   # HexNAc + HexNAc (Bi-HexNAc)
                                 Bisecting =          1009.4824,   # Bisecting
                                 Bisecting_dHex =     1155.5403,     # Bisecting + Deoxyhexose
                                Hex_HexNAc_ProA = 603.3236
)

# diagnostic_frags_exp = 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA'
# diagnostic_frags_exp = '(HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA'




#' Identify Glycan Spectra via Diagnostic Fragment Logic
#'
#' This function screens MS2 spectra for specific glycan diagnostic ions. It allows
#' users to define complex matching logic (e.g., AND, OR, NOT) to filter spectra
#' potentially containing glycans or specific glycan motifs.
#'
#' @param ms_data A \code{Spectra} object (typically denoised/processed) used for fragment screening.
#' @param ms_data_raw A \code{Spectra} object (typically raw) from which matched spectra will be extracted.
#' @param diagnostic_frags_list A named numeric vector of target m/z values for diagnostic fragments.
#'   Defaults to \code{diagnostic_frags_list_default}.
#' @param diagnostic_frags_exp A character string representing a logical expression.
#'   Names in this expression must match the names in \code{diagnostic_frags_list}.
#'   Example: \code{"HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA"}.
#' @param ppm_val Numeric. The mass tolerance in parts-per-million (ppm) for matching fragments.
#'
#' @details
#' The function evaluates the \code{diagnostic_frags_exp} by checking the presence of each
#' fragment in the \code{diagnostic_frags_list} within the specified \code{ppm_val} tolerance.
#' If a match is found, the function also identifies the preceding MS1 scan for that MS2 spectrum.
#'
#'
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{spectrum_info}: A data frame containing metadata for matched MS1 and MS2 scans,
#'     including boolean flags (\code{_flag}) for the presence of each diagnostic fragment.
#'   \item \code{selected_ms_data}: A \code{Spectra} object containing only the identified
#'     MS1 and MS2 spectra.
#' }
#'
#' @importFrom Spectra filterMsLevel peaksData rtime
#' @importFrom dplyr bind_rows distinct
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Define custom diagnostic fragments
#' my_frags <- c(HexNAc = 204.0867, Neu5Ac = 291.0954)
#'
#' # Filter spectra that MUST have HexNAc AND Neu5Ac
#' result <- FindSpectrumByDiagnosticFragments(
#'   ms_data = processed_data,
#'   ms_data_raw = raw_data,
#'   diagnostic_frags_list = my_frags,
#'   diagnostic_frags_exp = "HexNAc & Neu5Ac",
#'   ppm_val = 20
#' )
#'
#' # View matched spectrum IDs
#' print(result$spectrum_info$ms2_spectrum_id)
#' }
FindSpectrumByDiagnosticFragments = function(ms_data, ms_data_raw, diagnostic_frags_list = diagnostic_frags_list_default,
                                          diagnostic_frags_exp, ppm_val) {

  spectra_data_ms1 = Spectra::filterMsLevel(ms_data, 1)
  spectra_data_ms2 = Spectra::filterMsLevel(ms_data, 2)

  ms2_spectra_num = length(spectra_data_ms2)


  pb = utils::txtProgressBar(min = 0, max = ms2_spectra_num, style = 3)
  on.exit(close(pb))

  res_list = list()

  ms1_rt_all = Spectra::rtime(spectra_data_ms1)
  all_ms2_peaks_data = Spectra::peaksData(spectra_data_ms2)

  for (i in seq_along(spectra_data_ms2)) {
    # progress bar
    utils::setTxtProgressBar(pb, i)

    peaks_matrix = all_ms2_peaks_data[[i]]
    mz_values = peaks_matrix[, "mz"]

    ms2_id = spectra_data_ms2$spectrumId[i]

    presence_map <- sapply(
      names(diagnostic_frags_list), function(name) {
        target_mz <- diagnostic_frags_list[[name]]
        any(abs((mz_values - target_mz) / target_mz) * 10^6 <= ppm_val)
      }
    )

    is_match = base::eval(parse(text = diagnostic_frags_exp), envir = as.list(presence_map))


    presence_list_flag <- as.list(presence_map)
    names(presence_list_flag) <- paste0(names(presence_list_flag), "_flag")



    # extract info
    if (is_match) {

      spectra_data_target_ms2 = spectra_data_ms2[i]

      ms2_rt = spectra_data_target_ms2$rtime

      candidate_index_ms1 = which(ms1_rt_all < ms2_rt)

      if (length(candidate_index_ms1) == 0) {
        stop("No MS1 scan with retention time earlier than MS2.", call. = FALSE)
        break
      }

      closest_index_ms1 = candidate_index_ms1[which.min(ms2_rt - ms1_rt_all[candidate_index_ms1])]
      spectra_data_target_ms1 = spectra_data_ms1[closest_index_ms1]

      new_row_spectrum_info = data.frame(ms1_spectrum_id = spectra_data_target_ms1$spectrumId, ms1_retention_time = spectra_data_target_ms1$rtime,

                                    ms2_spectrum_id = spectra_data_target_ms2$spectrumId, ms2_precursor_charge = spectra_data_target_ms2$precursorCharge,
                                    ms2_precursor_mz = spectra_data_target_ms2$precursorMz, ms2_retention_time = spectra_data_target_ms2$rtime,

                                    presence_list_flag,

                                    ms1_peaks_count = spectra_data_target_ms1$peaksCount, ms1_total_ion_current = spectra_data_target_ms1$totIonCurrent,
                                    ms1_base_peak_mz = spectra_data_target_ms1$basePeakMZ,
                                    ms2_peaks_count = spectra_data_target_ms2$peaksCount, ms2_total_ion_current = spectra_data_target_ms2$totIonCurrent,
                                    ms2_base_peak_mz = spectra_data_target_ms2$basePeakMZ
      )

      res_list[[length(res_list) + 1]] = new_row_spectrum_info

    }

  }

  # close progress bar
  close(pb)

  # check res_list
  if (length(res_list) == 0) {
    message("No MS2 spectrum was identified")
    return(NULL)
  }

  # get the dataframe
  glycan_spectra_info <- dplyr::bind_rows(res_list)

  ms1_spectrum_id_list = dplyr::distinct(glycan_spectra_info, ms1_spectrum_id, .keep_all = F)

  all_glycan_spectra_id = c(ms1_spectrum_id_list$ms1_spectrum_id, glycan_spectra_info$ms2_spectrum_id)

  ms_data_filtered = ms_data_raw[ms_data_raw$spectrumId %in% all_glycan_spectra_id]

  return(
    list(
      spectrum_info = glycan_spectra_info,
      selected_ms_data = ms_data_filtered
    )
  )

}


























# 'spectrum_info': the output dataframe of FindSpectrumByDiagnosticFragments() fucntion,
# could also be dataframe contains such columns:
# "ms1_spectrum_id", "ms1_retention_time", "ms2_spectrum_id",
# "ms2_precursor_charge", "ms2_precursor_mz", "ms2_retention_time"

precursor_mz_ppm_default = 150
max_possible_candidates_num_default = 5


#' Match MS2 Precursors to Potential Glycan Compositions
#'
#' This function takes the identified glycan spectra and searches a glycan library
#' for potential compositions based on the precursor m/z and charge state. It
#' calculates mass errors in ppm and returns the top candidates for each spectrum.
#'
#' @param spectrum_info A data frame containing spectrum metadata. Required columns
#'   include: \code{"ms2_spectrum_id"}, \code{"ms2_precursor_charge"}, and \code{"ms2_precursor_mz"}.
#'   Usually the output from \code{FindSpectrumByDiagnosticFragments}.
#' @param glycan_lib A data frame acting as the reference glycan library. It must contain
#'   at least: \code{total_charge}, \code{glycan_monoisotopic_mz}, and composition details.
#' @param max_precursor_mz_ppm Numeric. The maximum allowable mass tolerance in ppm
#'   between the observed precursor and the library entry. Default is 150.
#' @param max_possible_candidates_num Numeric. The maximum number of top-ranked
#'   (by lowest ppm error) candidates to return for each MS2 spectrum. Default is 5.
#'
#' @details
#' The function filters the library by charge state, calculates the ppm error for
#' each entry, and retains candidates within the specified tolerance. Results are
#' merged back with the original \code{spectrum_info} for easy downstream analysis.
#'
#'
#'
#' @return A data frame (\code{tibble}) containing the original spectrum information
#'   joined with matched glycan candidates, their calculated \code{ppm_error},
#'   and their chemical compositions.
#'
#' @importFrom dplyr filter mutate arrange slice_head bind_rows left_join
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a glycan library and spectrum info
#' candidates <- FindPossibleGlycanComposition(
#'   spectrum_info = likely_glycan_spectrum_info,
#'   glycan_lib = N_glycan_library,
#'   max_precursor_mz_ppm = 50,
#'   max_possible_candidates_num = 3
#' )
#'
#' # Inspect the top matches for the first spectrum
#' head(candidates)
#' }
FindPossibleGlycanComposition = function(spectrum_info, glycan_lib,
                                         max_precursor_mz_ppm = precursor_mz_ppm_default,
                                         max_possible_candidates_num = max_possible_candidates_num_default) {

  required_cols = c("ms1_spectrum_id", "ms1_retention_time", "ms2_spectrum_id",
                    "ms2_precursor_charge", "ms2_precursor_mz", "ms2_retention_time")
  colnames_info = colnames(spectrum_info)
  if (!all(required_cols %in% colnames_info)) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }


  length_spectrum_info = dim(spectrum_info)[1]

  all_precursor_mz = spectrum_info[['ms2_precursor_mz']]
  all_precursor_charge = spectrum_info[['ms2_precursor_charge']]
  all_ms2_id = spectrum_info[['ms2_spectrum_id']]

  # initialize the possible_glycan_list
  possible_glycan_list = list()

  for (i in seq_along(1:length_spectrum_info)) {

    precursor_mz = all_precursor_mz[i]
    precursor_charge = all_precursor_charge[i]
    precursor_ms2_id = all_ms2_id[i]

    all_possible_glycan_comp = glycan_lib |>
      dplyr::filter(total_charge == precursor_charge) |>
      dplyr::mutate(ppm_error = abs(glycan_monoisotopic_mz - precursor_mz) / glycan_monoisotopic_mz * 1e6) |>
      dplyr::filter(ppm_error <= max_precursor_mz_ppm) |>
      dplyr::arrange(ppm_error)


    # if (dim(all_possible_glycan_comp)[1] >= max_possible_candidates_num) {
    #   top_glycan_candidate = dplyr::slice_head(all_possible_glycan_comp, n = max_possible_candidates_num)
    # } else if (dim(all_possible_glycan_comp)[1] >= 1 && dim(all_possible_glycan_comp)[1] < max_possible_candidates_num) {
    #   top_glycan_candidate = all_possible_glycan_comp
    # } else {
    #   print(precursor_ms2_id)
    # }
    #
    # top_glycan_candidate = dplyr::mutate(top_glycan_candidate, ms2_spectrum_id = rep(precursor_ms2_id, dim(top_glycan_candidate)[1]))
    #
    # possible_glycan_list[[length(possible_glycan_list) + 1]] = top_glycan_candidate




    if (nrow(all_possible_glycan_comp) > 0) {

      top_glycan_candidate = all_possible_glycan_comp |>
        dplyr::slice_head(n = max_possible_candidates_num) |>
        dplyr::mutate(ms2_spectrum_id = precursor_ms2_id)

      possible_glycan_list[[length(possible_glycan_list) + 1]] = top_glycan_candidate

    } else {

      message(paste0("No candidates found for MS2 ID: ", precursor_ms2_id))
    }





  }

  candidate_glycan_list = dplyr::bind_rows(possible_glycan_list) |>
    dplyr::left_join(spectrum_info, by = c('ms2_spectrum_id' = 'ms2_spectrum_id'))


  return(candidate_glycan_composition = candidate_glycan_list)

}

#
#
#   likely_glycan_spectrum_info = diagnostic_results$spectrum_info
#
#
#
#
# N_glycan_library = N_glycan_lib$glycan_monosaccharides_library
# N_glycan_library_charge = N_glycan_lib$adduct_combination_charge
#























#
# # find the optimal ppm
# ppm_start_default = 5
# ppm_end_default = 300
# ppm_step_default = 5
# # mannual_checked_ground_truth_glycan_ms2_id
# # optional, but works best with this
# # if you provide the mannual_checked_ground_truth_glycan_ms2_id, not need to provide the high_confident_glycan_exp
# # number of the IDs in the mannual_checked_ground_truth_glycan_ms2_id must >= 100
# # ex: mannual_checked_ground_truth_glycan_ms2_id_example = data.frame(ms2_id = c('function=2 process=0 scan=1568',
# #                                                                       'function=4 process=0 scan=1693',
# #                                                                       'function=4 process=0 scan=1694',
# #                                                                       'function=6 process=0 scan=1384',
# #                                                                       'function=4 process=0 scan=1365',
# #                                                                       'function=2 process=0 scan=1242',
# #                                                                       'function=3 process=0 scan=1324',
# #                                                                       'function=2 process=0 scan=1243'))
#
#
# # define the glycan diagnostics fragments
# diagnostic_frags_list_default = c(HexNAc =         204.08667,   # HexNAc
#                                   HexNAc_ProA =        441.2708,    # HexNAc + ProA
#                                   dHex_HexNAc_ProA =   587.3287,    # HexNAc + Deoxyhexose + ProA
#                                   HexNAc_HexNAc_ProA = 644.3502,    # HexNAc + HexNAc + ProA
#                                   Hex_ProA = 400.2442,
#                                   Hex_Hex_ProA = 562.2970,
#                                   dHex_Hex_ProA = 546.3021,
#                                   Hex =                163.06007,   # Hex
#                                   Bi_HexNAc =          407.16607,   # HexNAc + HexNAc (Bi-HexNAc)
#                                   Bisecting =          1009.4824,   # Bisecting
#                                   Bisecting_dHex =     1155.5403,     # Bisecting + Deoxyhexose
#                                   Hex_HexNAc_ProA = 603.3236
# )
#
# # diagnostic_frags_exp = 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA'
#
# # for the 'high_confident_glycan_precursor_mz', try to involve as much glycan precursor mz and corresponding charge state
# # that could possibel exist in your ms file, also, it's better if the glycan abundance is high in your ms file
# # after PpmEstimater find the ms2 spectrum likely to be glycan using 'high_confident_diagnostic_frags_exp',
# # it will check whether the precursor mz and charge of this ms2 spectrum could match the mz and charge in
# # 'high_confident_glycan_precursor_mz' and 'high_confident_glycan_precursor_charge',
# # so, please involve as much glycan precursor mz value
#
#
#
# high_confident_glycan_precursor_mz = c(
#   # 2H
#   484.7315262, 557.7604762, 565.7579262, 667.2976262, 646.7843262, 740.3265762, 768.8373262, 748.3240262,
#   841.8662762, 870.3770262, 821.3529762, 727.8107262, 943.4059762, 849.8637262, 971.9167262, 1044.945676,
#   922.8926762, 829.3504262, 914.8952262, 1073.456426, 951.4034262, 902.3793762, 808.8371262, 1016.434926,
#   1146.485376, 1024.432376, 930.8901262, 995.9216262, 1117.974626, 1003.919076, 889.8635262, 1076.948026,
#   1191.003576, 970.8899262, 1051.916326, 1132.942726, 1213.969126, 1294.995526,
#   # H+K
#   503.7097447, 576.7386947, 584.7361447, 686.2758447, 665.7625447, 759.3047947, 787.8155447, 767.3022447,
#   860.8444947, 889.3552447, 840.3311947, 746.7889447, 962.3841947, 868.8419447, 990.8949447, 1063.923895,
#   941.8708947, 848.3286447, 933.8734447, 1092.434645, 970.3816447, 921.3575947, 827.8153447, 1035.413145,
#   1165.463595, 1043.410595, 949.8683447, 1014.899845, 1136.952845, 1022.897295, 908.8417447, 1095.926245,
#   1209.981795, 989.8681447, 1070.894545, 1151.920945, 1232.947345, 1313.973745,
#
#
# )
#
#
# high_confident_glycan_precursor_charge = c(rep(2, 38), rep(2, 38))
#
#
#
# PpmEstimater = function(ppm_start, ppm_end, ppm_step,
#
#                         ms_data,
#
#                         diagnostic_fragments_list,
#
#                        high_confident_diagnostic_frags_exp,
#                        high_confident_glycan_precursor_mz,
#                        high_confident_glycan_precursor_charge,
#
#                        diagnostic_frags_exp,
#
#                        # filter_method_ms1_ppm_est,
#                        # filter_method_ms2_ppm_est,
#                        # threshold_ms1_ppm_est,
#                        # threshold_ms2_ppm_est,
#
#                         mannual_checked_ground_truth_glycan_ms2_id)
# #
# #
#
#
#   if (missing(diagnostic_fragments_list)) {
#     stop(
#       paste0(
#         "Argument 'diagnostic_fragments_list' is required.\n",
#         "Example of ProA labeled N-glycan 'diagnostic_fragments_list':\n",
#         "  list()"
#       ),
#       call. = FALSE
#     )
#   }
#
#   if (!missing(mannual_checked_ground_truth_glycan_ms2_id) &&
#       !missing(high_confident_glycan_exp)) {
#
#     stop("Argument 'mannual_checked_ground_truth_glycan_ms2_id' is provided, 'high_confident_glycan_exp' is not required",
#          call. = FALSE)
#
#   } else if (missing(mannual_checked_ground_truth_glycan_ms2_id) &&
#             !missing(high_confident_glycan_exp)) {
#
#     filter_results = SpectrumQcFilter(ms_data = ms_data,
#                                       filter_method_ms1 = ,
#                                       threshold_ms1 = ,
#
#                                       filter_method_ms2 = ,
#                                       threshold_ms2 = ,
#                                       plot_option = T)
#
#
#     ms2_data_filtered = Spectra::filterMsLevel(filter_results$filtered_ms_data, 2)
#
#
#
#
#
#
#
#     for (i in seq_along(ms2_data_filtered)) {
#       # progress bar
#       utils::setTxtProgressBar(utils::txtProgressBar(min = 0, max = length(ms2_data_filtered), style = 3), i)
#       #print(i)
#
#       current_spectra_ms2 = as.data.frame(Spectra::peaksData(ms2_data_filtered)[[i]])
#
#
#
#       current_spectra_ms2_distinct = dplyr::distinct(current_spectra_ms2, intensity, .keep_all = TRUE)
#
#     }
#
#
#
#
#
#
#   } else if (!missing(mannual_checked_ground_truth_glycan_ms2_id) &&
#              missing(high_confident_glycan_exp)) {
#
#     if (dim(mannual_checked_ground_truth_glycan_ms2_id) < 100) {
#
#       stop("number of the IDs in the 'mannual_checked_ground_truth_glycan_ms2_id' must >= 100",
#            call. = FALSE)
#
#     } else {
#
#     }
#
#   } else {
#
#     stop("Please provide either 'high_confident_glycan_exp' or 'mannual_checked_ground_truth_glycan_ms2_id'.", call. = FALSE)
#
#   }
#
#
#
#   high_confident_ms2_id = 23
#
# }
#
#
#
#
# list1_and = c(a = 100, b = 300, c=400)
# list2_or = c(b = 300, c=400)
# list3_not = c(d = 900)
#
# logical_exp = (a and b) and (b or c) and not (d)
#
#
# if (logical_exp) {
#   print{i}
# }
#
#
#
# stop("Please provide ground truth: either 'high_confident_glycan_exp' or 'mannual_checked_ground_truth_glycan_ms2_id'.", call. = FALSE)
#
#
# f <- function(x) {
#   missing(x)
# }
#

























molecular_formula_list_default = c(
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

molecular_formula_name_default = c('Hex', 'HexNAc', 'dHex', 'Neu5Ac', 'HexA', 'Neu5Gc',
                           'H', 'Na', 'K',
                           'ProA')

threshold_iso_probalility_default = 0.01

# > ms1_window_left = 1
# > ms1_window_right = 2

bin_width_default = 0.3




#' Validate Glycan Compositions using MS1 Isotope Pattern Matching
#'
#' This function validates candidate glycan compositions by comparing their theoretical
#' isotopic distributions with the observed MS1 spectra. It calculates a similarity
#' score (cosine score) between the predicted and actual isotope patterns to
#' identify the most likely glycan composition for each MS2 spectrum.
#'
#' @param spectrum_matching_info A data frame containing initial matching results.
#'   Required columns: \code{"ms1_spectrum_id"}, \code{"ms2_spectrum_id"}, and
#'   \code{"ms2_precursor_mz"}. Typically the output from \code{FindPossibleGlycanComposition}.
#' @param molecular_names A character vector of names from \code{molecular_formula_list}
#'   to be used for formula calculation (e.g., sugars and adducts).
#' @param molecular_formula_list A named character vector of molecular formulas for
#'   monosaccharides, labels, and adducts.
#' @param ms_data A \code{Spectra} object containing the raw MS1 spectra.
#' @param ms1_window_left Numeric. The left (lower) m/z offset from the precursor
#'   m/z to define the MS1 isolation window for isotope pattern extraction.
#' @param ms1_window_right Numeric. The right (upper) m/z offset from the precursor
#'   m/z to define the MS1 isolation window.
#' @param bin_width Numeric. The bin width for m/z discretization during cosine
#'   score calculation.
#' @param threshold_iso_probalility Numeric. The minimum relative abundance threshold
#'   for theoretical isotopic peaks to be included in the comparison. Default is 0.01.
#'
#' @details
#' For each MS2 spectrum with multiple potential glycan candidates, the function
#' retrieves the corresponding MS1 scan. It then simulates the theoretical isotope
#' distribution based on the candidate's molecular formula and compares it to the
#' experimental peaks within the specified m/z window.
#'
#'
#'
#' @return A data frame (\code{tibble}) updated with isotope validation metrics.
#'   Candidates with higher cosine scores represent better matches to the
#'   experimental data.
#'
#' @importFrom Spectra peaksData
#' @importFrom dplyr filter bind_rows
#' @importFrom stats setNames
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Validate candidates using a 1 Da window on the left and 2 Da on the right
#' validated_results <- ValidateGlycanCompositionByIsotopePattern(
#'   spectrum_matching_info = candidate_list,
#'   ms_data = raw_spectra,
#'   ms1_window_left = 1.0,
#'   ms1_window_right = 2.5,
#'   bin_width = 0.05
#' )
#' }
ValidateGlycanCompositionByIsotopePattern = function(spectrum_matching_info,

                                     molecular_names = molecular_formula_name_default,
                                     molecular_formula_list = molecular_formula_list_default,


                                     ms_data,
                                     ms1_window_left, ms1_window_right, bin_width,

                                     threshold_iso_probalility = threshold_iso_probalility_default) {

  required_cols = c("ms1_spectrum_id", "ms2_spectrum_id")

  colnames_info = colnames(spectrum_matching_info)

  if (!all(required_cols %in% colnames_info)) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }



  ms2_id_all = unique(spectrum_matching_info[['ms2_spectrum_id']])

  ms_id_map = stats::setNames(seq_along(ms_data), ms_data$spectrumId)


  ms2_id_num = length(ms2_id_all)



  # progress bar
  pb = utils::txtProgressBar(min = 0, max = ms2_id_num, style = 3)
  on.exit(close(pb))

  res_list = list()

  for (i in seq_along(ms2_id_all)) {

    # progress bar
    utils::setTxtProgressBar(pb, i)

    current_ms2_id = ms2_id_all[i]

    current_matching_info = dplyr::filter(spectrum_matching_info, ms2_spectrum_id == current_ms2_id)

    if (dim(current_matching_info)[1] == 1) {

      res_list[[length(res_list) + 1]] = current_matching_info
      next

    } else {

      current_ms1_id = unique(current_matching_info[['ms1_spectrum_id']])

      current_precursor_mz = unique(current_matching_info[['ms2_precursor_mz']])

      ms1_window_left_range = current_precursor_mz - ms1_window_left
      ms1_window_right_range = current_precursor_mz + ms1_window_right


      this_index = ms_id_map[current_ms1_id]

      current_ms1_data = as.data.frame(Spectra::peaksData(ms_data)[[this_index]])

      ms1_isolation_window_data = dplyr::filter(current_ms1_data, mz >= ms1_window_left_range & mz <= ms1_window_right_range)

      molecular_list = molecular_formula_list[molecular_names]


      best_matching_info = GetIsotopePatternCosineScore(matching_info = current_matching_info,
                                                            mol_list = molecular_list,
                                                            threshold_iso = threshold_iso_probalility,
                                                            ms1_window_data = ms1_isolation_window_data,
                                                            bin_w = bin_width)


      res_list[[length(res_list) + 1]] = best_matching_info

    }

  }

  # close progress bar
  close(pb)

  validate_spectrum_matching_info = dplyr::bind_rows(res_list)

  return(validate_spectrum_matching_info)


}


