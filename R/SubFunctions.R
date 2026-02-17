library('Spectra')
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(enviPat)

#============
# functions 
#============

# function 1: get the RMSE, CV, d2_smoothness, d1_total_var values to determine the spar value 

# input: 
# spline_fit_result: the fitted spline 
# x_val: the x vlaue used to fit the spline 
# y_val: the y value used to fit the spline 
# y_predict_val: the predicted y value by the fitted spline 

# output: 
# RMSE, CV, d2_smoothness, d1_total_var 

GetSplineSparEvaValue = function(spline_fit_result, x_val, y_val, y_predict_val) {
  # add more point to get more precise estimation 
  d1 = stats::predict(spline_fit_result, 
                      seq(from = x_val[1], to = x_val[length(x_val)], length.out = length(x_val)*2), deriv = 1)$y
  d2 = stats::predict(spline_fit_result, 
                      seq(from = x_val[1], to = x_val[length(x_val)], length.out = length(x_val)*2), deriv = 2)$y
  
  # RMSE
  RMSE_val = sqrt(mean((y_val - y_predict_val)^2))
  
  # cv: in the smooth.spline(): leave-one-out cross-validation(LOOCV)
  if (!is.null(spline_fit_result$cv.crit)) {
    loocv = spline_fit_result$cv.crit
  } else {
    loocv = NA
  }
  
  # first derivative 
  # first derivative smoothness: total variation 
  # (smaller changes in the derivative indicate greater smoothness)
  # because of the data was sorted by the intensity, so theoretically, 
  # d1 should always be positive value 
  d1_total_var = sum(abs(diff(d1)))
  
  # second derivative 
  # rate of change of tangent slope
  dx = diff(seq(from = x_val[1], to = x_val[length(x_val)], length.out = length(x_val)*2))
  dx = c(dx, mean(dx))
  d2_smoothness = sum((d2^2) * dx)
  
  return(list(RMSE = RMSE_val, CV = loocv, D1 = d1_total_var, D2 = d2_smoothness))
}









# function2: find the best spar parameter for the smooth::spline() used in glycan spectrum

# input: 
# x: order(start from 1, and to the length(y)), y: spectrum intensity, ln() transformed and must be sort()
# spar_start: start value for the spar setting in smooth.spline()
# spar_end: end value for the spar setting in smooth.spline()
# spar_step: step size 

# weight: weight of RMSE, CV, d2_smoothness, d1_total_var, default: 0.3, 0.3, 0.3, 0.1, respectively

# use_cv: TRUE, according to the manual of smooth.spline(), 
# when there are duplicate points in x, avoid cv = TRUE (because of the LOOCV)
# plot: whether plot the best spar fitting plot
FindSplineSpar = function(x, y,
                            spar_start, spar_end, spar_step,
                            RMSE_weight, CV_weight, D1_weight, D2_weight,
                            use_cv = TRUE, plot = TRUE) {
  stopifnot(length(x) == length(y))
  
  spar_seq = seq(spar_start, spar_end, by = spar_step)
  RMSE_w = RMSE_weight 
  CV_w = CV_weight 
  D1_w = D1_weight
  D2_w = D2_weight
  
  spar_data_info = data.frame(spar = numeric(), RMSE = numeric(), CV = numeric(), D1 = numeric(), D2 = numeric())
  
  for (i in spar_seq) {
    fit = stats::smooth.spline(x, y, spar = i, cv = use_cv)
    # if (inherits(fit, "try-error")) return(NULL)
    y_predict = stats::predict(fit, x)$y
    
    eva_value = GetSplineSparEvaValue(spline_fit_result = fit, 
                                          x_val = x, y_val = y, y_predict_val = y_predict) 
    
    new_spar_info = data.frame(spar = i, RMSE = eva_value$RMSE, CV = eva_value$CV, 
                               D1 = eva_value$D1, D2 = eva_value$D2)
    
    spar_data_info = dplyr::bind_rows(spar_data_info, new_spar_info)
  }
  
  # max min normalization 
  max_min_norm = function(v) {
    (v - min(v)) / (max(v) - min(v) + 1e-12)
  } 
  spar_data_info_norm = data.frame(spar = spar_data_info$spar, 
                                   RMSE_norm = max_min_norm(spar_data_info$RMSE), 
                                   CV_norm = max_min_norm(spar_data_info$CV), 
                                   D1_norm = max_min_norm(spar_data_info$D1), 
                                   D2_norm = max_min_norm(spar_data_info$D2))
  
  # evaluation and determine the best spar 
  if (all(is.na(spar_data_info_norm$CV_norm))) {
    RMSE_w = RMSE_w + CV_w
    CV_w = 0
    spar_data_info_norm$CV_norm = 0
  }
  spar_data_info_norm = dplyr::mutate(spar_data_info_norm, 
                                      evaluation_score = RMSE_norm*RMSE_w + CV_norm*CV_w + D1_norm*D1_w + D2_norm*D2_w)
  
  best_idx = which.min(spar_data_info_norm$evaluation_score)
  best_spar = spar_data_info_norm$spar[best_idx]
  best_fit = stats::smooth.spline(x, y, spar = best_spar, cv = use_cv)
  
  # plotting part 
  final_plot = NULL
  
  if (plot) {
    base_plot <- function(data, mapping, title, y_lab) {
      ggplot(data, mapping) +
        geom_point(color = "black") +
        geom_vline(xintercept = best_spar, linetype = "dashed", color = "red") +
        theme_minimal() +
        labs(title = title, x = "spar", y = y_lab)
    }
    
    # RMSE
    RMSE_plot = base_plot(data = spar_data_info, 
                          mapping = aes(x=spar, y=RMSE), 
                          title = "RMSE", 
                          y_lab = "RMSE")
    # CV
    if (!all(is.na(spar_data_info$CV))) {
      
      cv_plot = base_plot(data = spar_data_info, 
                          mapping = aes(x=spar, y=CV), 
                          title = "Leave-one-out cross-validation (LOOCV)", 
                          y_lab = "LOOCV")
    } else {
      
      cv_plot = ggplot() + 
        annotate("text", x=0.5, y=0.5, label="CV not computed") + theme_void()
    }
    # D2
    D2_plot = base_plot(data = spar_data_info, 
                         mapping = aes(x=spar, y=D2), 
                         title = "D2 smoothness (∫(f'' )^2)", 
                         y_lab = "D2_smoothness")
    # D1
    D1_plot = base_plot(data = spar_data_info, 
                         mapping = aes(x=spar, y=D1), 
                         title = "D1 total variation", 
                         y_lab = "D1_total_var")
    
    final_plot = (RMSE_plot + cv_plot) / (D2_plot + D1_plot) + 
      patchwork::plot_annotation(title = paste("Optimization for best spar:", round(best_spar, 4)),
                      subtitle = "Red dashed line indicates the selected values")
  }

  return(list(
    best_spar = best_spar, 
    spar_data_info = spar_data_info, 
    spar_data_info_norm = spar_data_info_norm, 
    spar_data_info_plot = final_plot
    )
  )
}





































# monos_list: named numeric vector of monosaccharide masses for specific type of glycan 
# min_num_custom/max_num_custom: min_num_monosaccharides_additional_customized_list / max_num_monosaccharides_additional_customized_list:
#   named numeric vectors for customized monosaccharides (may be missing/empty)

MakeAllMonosCombos = function(monos_list,
                              min_num_custom = NULL,
                              max_num_custom = NULL,
                              minmax_map) {
  
  monos = names(monos_list)
  
  # build min/max vectors aligned to monos
  min_vec = setNames(rep(NA_real_, length(monos)), monos)
  max_vec = setNames(rep(NA_real_, length(monos)), monos)
  
  # fill defaults (if present)
  for (nm in intersect(names(minmax_map), monos)) {
    min_vec[nm] = minmax_map[[nm]][1]
    max_vec[nm] = minmax_map[[nm]][2]
  }
  
  # fill customized (override or add)
  if (!is.null(min_num_custom)) {
    min_vec[names(min_num_custom)] = min_num_custom
  }
  if (!is.null(max_num_custom)) {
    max_vec[names(max_num_custom)] = max_num_custom
  }
  
  # validate: no missing ranges
  missing_range = monos[is.na(min_vec) | is.na(max_vec)]
  
  if (length(missing_range) > 0) {
    stop(sprintf("Missing min/max range for monosaccharide(s): %s.", 
                 paste(sprintf("'%s'", missing_range), collapse = ", ")), 
         call. = FALSE)
  }
  
  # validate: min <= max
  bad <- monos[min_vec > max_vec]
  
  if (length(bad) > 0) {
    stop(sprintf("Invalid range (min > max) for monosaccharide(s): %s.", 
                 paste(sprintf("'%s'", bad), collapse = ", ")), 
         call. = FALSE)
  }
  
  # create list of sequences for expand.grid
  seq_list <- lapply(monos, function(n) {
    seq(from = min_vec[[n]], to = max_vec[[n]], by = 1)
  })
  
  names(seq_list) = paste0(monos)  
  
  expand.grid(seq_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  seq_list = expand.grid(seq_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  return(seq_list)
}

































PlotUnfilteredMsVaribles = function(msdata, ms_level, filter_method) {
  
  ms_level_vec = msdata[["msLevel"]]
  
  pic_ms_var = list()
  
  # ms plot 
  for (i in names(filter_method)) {
    
    ms_v = msdata[[i]][ms_level_vec == ms_level]
    
    df_raw = data.frame(index = seq_along(ms_v), value = ms_v)
    df_sorted = data.frame(index = seq_along(ms_v), value = sort(ms_v))
    
    p1 = ggplot(df_raw, aes(x = index, y = value)) +
      geom_col(width = 1, fill = "black") +
      labs(x = "Index", y = i,  title = paste0('MS',ms_level , '_', i)) + 
      theme_classic() +
      theme(legend.position = "none", plot.title = element_text(size = 12))
    
    p2 = ggplot(df_sorted, aes(x = index, y = value)) +
      geom_col(width = 1, fill = "black") +
      labs(x = "Sorted_index", y = i,  title = paste0('MS',ms_level , '_', i, '_sorted')) +
      theme_classic() +
      theme(legend.position = "none", plot.title = element_text(size = 12))
    
    pic_ms_var[[paste0('MS',ms_level , '_', i)]] = p1 / p2
    
  }
  
  return(pic_ms_var)
  
}
  

  
  




#                                                      filter_method

PlotFilteredMsVaribles = function(msdata, ms_level, filter_threshold_value) {
  
  ms_level_vec = msdata[["msLevel"]]
  
  pic_ms_var = list()
  
  # ms plot 
  for (i in names(filter_threshold_value)) {
    
    ms_v = msdata[[i]][ms_level_vec == ms_level]
    
    df_raw = data.frame(index = seq_along(ms_v), value = ms_v) %>% 
      dplyr::mutate(
        status = dplyr::if_else(
          value >= filter_threshold_value[[i]][1] & value <= filter_threshold_value[[i]][2], 'keep', 'filtered', missing = "filtered"
          )
        )
    
    df_sorted = data.frame(index = seq_along(ms_v), value = sort(ms_v)) %>% 
      dplyr::mutate(
        status = dplyr::if_else(
          value >= filter_threshold_value[[i]][1] & value <= filter_threshold_value[[i]][2], 'keep', 'filtered', missing = "filtered"
        )
      )

    
    p1 = ggplot(df_raw, aes(x = index, y = value, fill = status)) + 
      geom_col(width = 1) + 
      scale_fill_manual(values = c("keep" = "black", "filtered" = "grey"),
                        name = "Status"
                        ) + 
      labs(x = "Index", y = i,  title = paste0('MS',ms_level , '_', i)) + 
      theme_classic() +
      theme(legend.position = "right", plot.title = element_text(size = 12))

    p2 = ggplot(df_sorted, aes(x = index, y = value, fill = status)) +
      geom_col(width = 1) + 
      scale_fill_manual(values = c("keep" = "black", "filtered" = "grey"),
                        name = "Status"
                        ) + 
      labs(x = "Sorted_index", y = i,  title = paste0('MS',ms_level , '_', i, '_sorted')) +
      theme_classic() +
      theme(legend.position = "none", plot.title = element_text(size = 12))
    
    
    pic_ms_var[[paste0('MS',ms_level , '_', i)]] = p1 / p2
    
  }
  
  return(pic_ms_var)
  
}












GetQcFilterThreshold = function(msdata, ms_level, filter_method, filter_threshold) {
  
  ms_level_vec = msdata[["msLevel"]]
  
  threshold_value = list()
  
  # start
  for (i in names(filter_method)) {
    # i:        "peaksCount"
    mtd_ms = filter_method[[i]]      # mtd_ms:  "mean_sd" 
    thr_ms = filter_threshold[[i]]        # thr_ms:   2
    
    ms_vars = msdata[[i]][ms_level_vec == ms_level]
    
    if (mtd_ms == 'mean_sd') {
      
      if (thr_ms == 0) {
        
        stop(paste0("Invalid ", i, "threshold value: ", thr_ms, " for ", mtd_ms, 
                    "If you do NOT want to do the filter for ", i, "do not write it in 'filter_method_ms1' or 'filter_method_ms2'. 
                    If you want to do the ", mtd_ms, "filter for ", i, ", threshold value should > 0. "),
             call. = FALSE)
        
      } else if (thr_ms < 0) {
        
        stop(paste0("Invalid value: ", i, "threshold for method 'mean_sd' should > 0 "),
             call. = FALSE)
        
      } else {
        
        ms_vars_sorted = data.frame(data_varibles = sort(ms_vars)) %>%
          dplyr::mutate(data_diff = data_varibles - dplyr::lag(data_varibles)) %>%
          dplyr::mutate(data_diff = tidyr::replace_na(data_diff, 0))
        
        ms_vars_knee = dplyr::filter(ms_vars_sorted, data_diff >= mean(ms_vars_sorted$data_diff) + thr_ms * sd(ms_vars_sorted$data_diff))
        
        if (nrow(ms_vars_knee) == 0) {
          stop(paste0("No differ (sorted) of the ", i, ">= ", mean(ms_vars_sorted$data_diff) + thr_ms * sd(ms_vars_sorted$data_diff), 
                      ". Please consider decrease the ", thr_ms, " or using other method. "), 
               call. = FALSE)
        }
        #
        max_gap = ms_vars_knee$data_varibles[which.max(ms_vars_knee$data_diff)]
        # threshold_value = c(threshold_value, max_gap, max(ms_vars))
        threshold_value[[i]] = c(max_gap, max(ms_vars))
        
      }
      
    } else if (mtd_ms == 'quantile_prob') {
      
      if (any(thr_ms < 0) || 
          any(thr_ms > 1)) {
        
        stop(sprintf("Invalid threshold for '%s' using method '%s': value must be within [0, 1].", i, mtd_ms), 
             call. = FALSE)
        
      } else if (thr_ms[1] >= thr_ms[2]) {
        
        stop(sprintf("Invalid threshold for '%s' using method '%s': expected c(start, end) with start < end, got c(%s, %s).", 
                     i, mtd_ms, thr_ms[1], thr_ms[2]), 
             call. = FALSE)
        
      } else {
        # 
        quantile_val = as.numeric(stats::quantile(ms_vars, probs = c(thr_ms), na.rm = TRUE))
        # threshold_value = c(threshold_value, quantile_val)
        threshold_value[[i]] = c(quantile_val)
        
      }
      
    } else if (mtd_ms == 'start_end') {
      
      if (thr_ms[1] >= thr_ms[2]) {
        
        stop(sprintf("Invalid threshold for '%s' using method '%s': expected c(start, end) with start < end, got c(%s, %s).", 
                     i, mtd_ms, thr_ms[1], thr_ms[2]), 
             call. = FALSE)
        
      } else if (thr_ms[1] < min(ms_vars) || 
                 thr_ms[2] > max(ms_vars)) {
        
        stop(sprintf("Invalid threshold for '%s' using method '%s': expected c(start, end) with start < end, got c(%s, %s).", 
                     i, mtd_ms, thr_ms[1], thr_ms[2]), 
             call. = FALSE)
        
      } else {
        
        # threshold_value = c(threshold_value, thr_ms)
        threshold_value[[i]] = c(thr_ms)
        
      }
      
    }
    
  }
  
  return(threshold_value)
  
}






GetQcFilteredMsData = function(msdata, ms_level, threshold_value) {
  
  ms_level_vec = msdata[["msLevel"]]
  
  idx_ms = which(ms_level_vec == ms_level)
  logic_ms_keep = rep(TRUE, length(idx_ms))
  
  for (i in names(threshold_value)) {
    
    rng = threshold_value[[i]]
    ms_vars = msdata[[i]][idx_ms]
    
    logic_ms_keep = logic_ms_keep & (ms_vars >= rng[1]) & (ms_vars <= rng[2])
    
  }

  logic_ms_all = rep(TRUE, length(ms_level_vec)) # Default to keeping everything 
  
  logic_ms_all[idx_ms] = logic_ms_keep
  filtered_msdata = msdata[logic_ms_all]
  
  return(filtered_msdata)

}







GetSplineSegmentationNoise = function(denoising_detail, denoising_method, transform_fun, spectra_distinct, ms_id, ms_rt) {
#GetSplineSegmentationNoise = function(den_detail, den_mth, trans_fun, spectra_distinct, ms_id, ms_rt) {

  # background noise level detection 
  # the intensity is transformed!!!!!
  y = transform_fun(sort(spectra_distinct$intensity))
  n_total = length(y)
  x_idx = 1:n_total

  # according to the manual annotation, it is found that segmented regression can only well describe the
  # spectrum that have good fragmentation pattern (ions are completely fragmented), 
  # however, for the spectrum that do NOT have good fragmentation pattern, which means the ions are not completely fragmented, 
  # signal ions have low intensity, for these 'NOT have good fragmentation pattern' spectrm, 
  # using a cubic smoothing spline to do the signal regression
  # use this method to determine the fragmentation is good or not good
  
  top4_peaks = dplyr::arrange(spectra_distinct, desc(intensity)) %>% 
    head(4)
  
  # for some of the isotopic patterns, monoisotopic peak have similar intensity with M+1 peak, or even M+1 peak are higher than monoisotopic
  # if sort from the highest intensity to lower intensity, the difference(abs(diff(intensity_tops))) might be 0.5, 1, 0.5
  # so here sort() is used
  intensity_tops = sort(top4_peaks$mz)
  
  intensity_tops_diff = abs(diff(intensity_tops))
  
  is_complex_pattern <- all(abs(intensity_tops_diff - 1) < 0.1) ||
    all(abs(intensity_tops_diff - 0.5) < 0.05) ||
    all(abs(intensity_tops_diff - 0.33) < 0.03)
  
  # for some of the isotopic patterns, especially for the non-proton ions, 
  # fragmentation is not complete, the precursor isotopes have the highest intensity, 
  # and the intensity of these precursor isotopes sometimes are similar
  # so if do not remove the top 5 peaks, the minimum slope point might be among these top peaks
  if (is_complex_pattern) {
    y_ss <- y[1:max(1, (n_total - denoising_detail[[denoising_method]]$top_n_to_remove))]
    x_ss <- 1:length(y_ss)
    
    method_for_regression = 'Spline_regression'
    
    # determine the best spar parameter for the stats::smooth.spline()
    best_spar_info = FindSplineSpar(x_ss, y_ss, 
                                       spar_start = denoising_detail[[denoising_method]]$spar_start, 
                                       spar_end = denoising_detail[[denoising_method]]$spar_end, 
                                       spar_step = denoising_detail[[denoising_method]]$spar_step, 
                                       RMSE_weight = denoising_detail[[denoising_method]]$RMSE_weight, 
                                       CV_weight = denoising_detail[[denoising_method]]$CV_weight, 
                                       D2_weight = denoising_detail[[denoising_method]]$D2_weight, 
                                       D1_weight = denoising_detail[[denoising_method]]$D1_weight, 
                                       use_cv = denoising_detail[[denoising_method]]$use_cv, 
                                       plot = F)
    
    try(
      {
        fit_ss = stats::smooth.spline(x_idx, y, spar = best_spar_info$best_spar, cv = denoising_detail[[denoising_method]]$use_cv)
      
      # regression spline and first-order derivative (slope)
      slope = stats::predict(fit_ss, x_ss, deriv = 1)$y
      
      threshold_idx = which.min(slope)

      # check the threshold index
      if (!is.null(threshold_idx) && threshold_idx > 0 && threshold_idx < n_total) {
        # intensity threshold
        threshold_intensity = sort(spectra_distinct$intensity)[threshold_idx]
      } else {
        stop(sprintf("No threshold for denoising was detected in %s spectrum.", ms_id), 
                   call. = FALSE)
      }
      
      fit_ss_eva_value = GetSplineSparEvaValue(
        spline_fit_result = fit_ss, 
        x_val = x_idx, 
        y_val = y, 
        y_predict_val = stats::predict(fit_ss, x_idx)$y
      ) 

      new_regression = data.frame(ms2_spectrum_id = ms_id, ms2_retention_time = ms_rt, 
                                  
                                  threshold_value = threshold_intensity, 
                                      
                               denoising_method = method_for_regression, 
                               spar = best_spar_info$best_spar, 
                               RMSE = fit_ss_eva_value$RMSE, LOOCV = fit_ss_eva_value$CV, D1 = fit_ss_eva_value$D1, D2 = fit_ss_eva_value$D2, 
                               
                               r_squared_linear = 0, adj_r_squared_linear = 0, residual_standard_error_linear = 0, 
                               slope_linear = 0, intercept_linear = 0, 
                               r_squared_non_linear = 0, adj_r_squared_non_linear = 0, residual_standard_error_non_linear = 0, 
                               slope_non_linear = 0, intercept_non_linear = 0)
      }, silent = F
   )
    
    
  } else {
    method_for_regression = 'Segmentation_regresion' 
    
    non_linear_func_expression = denoising_detail[[denoising_method]]$segmentated_non_linear_transform_fun
    
    if (is.character(non_linear_func_expression)) {
      non_linear_transform_fun = switch(non_linear_func_expression,
                             square_transform = function(z) z^2, 
                             exponential_transform = function(z) 2^z, 
                             stop("Unknown method"), call. = FALSE)
    } else if (is.function(non_linear_func_expression)) {
      non_linear_transform_fun = non_linear_func_expression
    } else {
      stop("'segmentated_non_linear_transform_fun' must be a function or character", call. = FALSE)
    }

    try({
      # segmented fitting
      fit_lm_init <- stats::lm(y ~ x_idx)
      sum_fit_lm_init = summary(fit_lm_init)
      
      # find break point
      fit_lm_seg = try(segmented::segmented(fit_lm_init, seg.Z = ~x_idx, psi = list(x_idx = length(y)/2)))
      if (inherits(fit_lm_seg, "try-error")) {

        stop("'segmented::segmented' try-error", call. = FALSE)
      }
      
      break_point = fit_lm_seg$psi[,"Est."]

      # linear and non_linear fitting
      fit_lm_linear = stats::lm(y[1:break_point] ~ x_idx[1:break_point])
      sum_fit_lm_linear = summary(fit_lm_linear)
      
      x_non_linear = x_idx[(break_point+1):length(y)]
      x_idx_non_linear = non_linear_transform_fun(x_non_linear)
      y_non_linear = y[(break_point+1):length(y)]

      fit_lm_non_linear = stats::lm(y_non_linear ~ x_idx_non_linear)
      # fit_lm_non_linear = stats::lm(log(y_non_linear) ~ x2)
      sum_fit_lm_non_linear = summary(fit_lm_non_linear)

      # summary(fit_lm_seg)
      threshold_idx = round(break_point)
      
      # check the threshold index
      if (!is.null(threshold_idx) && threshold_idx > 0 && threshold_idx < n_total) {
        # intensity threshold
        threshold_intensity = sort(spectra_distinct$intensity)[threshold_idx]
      } else {
        stop(sprintf("No threshold for denoising was detected in %s spectrum.", ms_id), 
             call. = FALSE)
      }

      new_regression = data.frame(ms2_spectrum_id = ms_id, ms2_retention_time = ms_rt, 
                                  
                                  threshold_value = threshold_intensity, 
                                  
                                  denoising_method = method_for_regression, 
                                  spar = 0, RMSE = 0, LOOCV = 0, D1 = 0, D2 = 0, 
                                  
                                  r_squared_linear = sum_fit_lm_linear$r.squared, adj_r_squared_linear = sum_fit_lm_linear$adj.r.squared, residual_standard_error_linear = sum_fit_lm_linear$sigma, 
                                  slope_linear = fit_lm_linear$coefficients[[2]], intercept_linear = fit_lm_linear$coefficients[[1]], 
                                  r_squared_non_linear = sum_fit_lm_non_linear$r.squared, 
                                  adj_r_squared_non_linear = sum_fit_lm_non_linear$adj.r.squared, 
                                  residual_standard_error_non_linear = sum_fit_lm_non_linear$sigma, 
                                  slope_non_linear = fit_lm_non_linear$coefficients[[2]], 
                                  #slope_non_linear_2 = fit_lm_non_linear$coefficients["x2_non_linear"], 
                                  intercept_non_linear = fit_lm_non_linear$coefficients[[1]])
      }, silent = F)
    
  }
  return(
    list(
      threshold_value = threshold_intensity, 
      row_regression_info = new_regression
    )
  )
}









GetSplineNoise = function(denoising_detail, denoising_method, transform_fun, spectra_distinct, ms_id, ms_rt) {

  # background noise level detection 
  # the intensity is transformed!!!!!
  y = transform_fun(sort(spectra_distinct$intensity))
  n_total = length(y)
  x_idx = 1:n_total
  
  
    y_ss <- y[1:max(1, (n_total - denoising_detail[[denoising_method]]$top_n_to_remove))]
    x_ss <- 1:length(y_ss)
    
    method_for_regression = 'Spline_regression'
    
    # determine the best spar parameter for the stats::smooth.spline()
    best_spar_info = FindSplineSpar(x_ss, y_ss, 
                                    spar_start = denoising_detail[[denoising_method]]$spar_start, 
                                    spar_end = denoising_detail[[denoising_method]]$spar_end, 
                                    spar_step = denoising_detail[[denoising_method]]$spar_step, 
                                    RMSE_weight = denoising_detail[[denoising_method]]$RMSE_weight, 
                                    CV_weight = denoising_detail[[denoising_method]]$CV_weight, 
                                    D2_weight = denoising_detail[[denoising_method]]$D2_weight, 
                                    D1_weight = denoising_detail[[denoising_method]]$D1_weight, 
                                    use_cv = denoising_detail[[denoising_method]]$use_cv, 
                                    plot = F)
    
    try(
      {
        fit_ss = stats::smooth.spline(x_idx, y, spar = best_spar_info$best_spar, cv = denoising_detail[[denoising_method]]$use_cv)
        
        # regression spline and first-order derivative (slope)
        slope = stats::predict(fit_ss, x_ss, deriv = 1)$y
        
        threshold_idx = which.min(slope)
        
        # check the threshold index
        if (!is.null(threshold_idx) && threshold_idx > 0 && threshold_idx < n_total) {
          # intensity threshold
          threshold_intensity = sort(spectra_distinct$intensity)[threshold_idx]
        } else {
          stop(sprintf("No threshold for denoising was detected in %s spectrum.", ms_id), 
               call. = FALSE)
        }
        
        fit_ss_eva_value = GetSplineSparEvaValue(
          spline_fit_result = fit_ss, 
          x_val = x_idx, 
          y_val = y, 
          y_predict_val = stats::predict(fit_ss, x_idx)$y
        ) 
        
        new_regression = data.frame(ms2_spectrum_id = ms_id, ms2_retention_time = ms_rt, 
                                    
                                    threshold_value = threshold_intensity, 
                                    
                                    denoising_method = method_for_regression, 
                                    spar = best_spar_info$best_spar, 
                                    RMSE = fit_ss_eva_value$RMSE, LOOCV = fit_ss_eva_value$CV, D1 = fit_ss_eva_value$D1, D2 = fit_ss_eva_value$D2, 
                                    
                                    r_squared_linear = 0, adj_r_squared_linear = 0, residual_standard_error_linear = 0, 
                                    slope_linear = 0, intercept_linear = 0, 
                                    r_squared_non_linear = 0, adj_r_squared_non_linear = 0, residual_standard_error_non_linear = 0, 
                                    slope_non_linear = 0, intercept_non_linear = 0)
      }, silent = F
    )
    
  return(
    list(
      threshold_value = threshold_intensity, 
      row_regression_info = new_regression
    )
  )
}




GetSegmentationNoise = function(denoising_detail, denoising_method, transform_fun, spectra_distinct, ms_id, ms_rt) {

  # background noise level detection 
  # the intensity is transformed!!!!!
  y = transform_fun(sort(spectra_distinct$intensity))
  n_total = length(y)
  x_idx = 1:n_total
  
    method_for_regression = 'Segmentation_regresion' 
    
    non_linear_func_expression = denoising_detail[[denoising_method]]$segmentated_non_linear_transform_fun
    
    if (is.character(non_linear_func_expression)) {
      non_linear_transform_fun = switch(non_linear_func_expression,
                                        square_transform = function(z) z^2, 
                                        exponential_transform = function(z) 2^z, 
                                        stop("Unknown method"), call. = FALSE)
    } else if (is.function(non_linear_func_expression)) {
      non_linear_transform_fun = non_linear_func_expression
    } else {
      stop("'segmentated_non_linear_transform_fun' must be a function or character", call. = FALSE)
    }
    
    try({
      # segmented fitting
      fit_lm_init <- stats::lm(y ~ x_idx)
      sum_fit_lm_init = summary(fit_lm_init)
      
      # find break point
      fit_lm_seg = try(segmented::segmented(fit_lm_init, seg.Z = ~x_idx, psi = list(x_idx = length(y)/2)))
      if (inherits(fit_lm_seg, "try-error")) {
        
        stop("'segmented::segmented' try-error", call. = FALSE)
      }
      
      break_point = fit_lm_seg$psi[,"Est."]
      
      # linear and non_linear fitting
      fit_lm_linear = stats::lm(y[1:break_point] ~ x_idx[1:break_point])
      sum_fit_lm_linear = summary(fit_lm_linear)
      
      x_non_linear = x_idx[(break_point+1):length(y)]
      x_idx_non_linear = non_linear_transform_fun(x_non_linear)
      y_non_linear = y[(break_point+1):length(y)]
      
      fit_lm_non_linear = stats::lm(y_non_linear ~ x_idx_non_linear)
      # fit_lm_non_linear = stats::lm(log(y_non_linear) ~ x2)
      sum_fit_lm_non_linear = summary(fit_lm_non_linear)
      
      # summary(fit_lm_seg)
      threshold_idx = round(break_point)
      
      # check the threshold index
      if (!is.null(threshold_idx) && threshold_idx > 0 && threshold_idx < n_total) {
        # intensity threshold
        threshold_intensity = sort(spectra_distinct$intensity)[threshold_idx]
      } else {
        stop(sprintf("No threshold for denoising was detected in %s spectrum.", ms_id), 
             call. = FALSE)
      }
      
      new_regression = data.frame(ms2_spectrum_id = ms_id, ms2_retention_time = ms_rt, 
                                  
                                  threshold_value = threshold_intensity, 
                                  
                                  denoising_method = method_for_regression, 
                                  spar = 0, RMSE = 0, LOOCV = 0, D1 = 0, D2 = 0, 
                                  
                                  r_squared_linear = sum_fit_lm_linear$r.squared, adj_r_squared_linear = sum_fit_lm_linear$adj.r.squared, residual_standard_error_linear = sum_fit_lm_linear$sigma, 
                                  slope_linear = fit_lm_linear$coefficients[[2]], intercept_linear = fit_lm_linear$coefficients[[1]], 
                                  r_squared_non_linear = sum_fit_lm_non_linear$r.squared, 
                                  adj_r_squared_non_linear = sum_fit_lm_non_linear$adj.r.squared, 
                                  residual_standard_error_non_linear = sum_fit_lm_non_linear$sigma, 
                                  slope_non_linear = fit_lm_non_linear$coefficients[[2]], 
                                  #slope_non_linear_2 = fit_lm_non_linear$coefficients["x2_non_linear"], 
                                  intercept_non_linear = fit_lm_non_linear$coefficients[[1]])
    }, silent = F)
    
  return(
    list(
      threshold_value = threshold_intensity, 
      row_regression_info = new_regression
    )
  )
}





GetQuantileNoise = function(denoising_detail, denoising_method, transform_fun, spectra_distinct, ms_id, ms_rt) {

  y = spectra_distinct$intensity
  n_total = length(y)
  
  quantile_prob_val = denoising_detail[[denoising_method]]
  
  quantile_value = stats::quantile(y, probs = quantile_prob_val)
  
  new_regression = data.frame(ms2_spectrum_id = ms_id, ms2_retention_time = ms_rt, 
                              
                              threshold_value = quantile_value, 
                              
                              denoising_method = 'Quantile', 
                              
                              spar = 0, RMSE = 0, LOOCV = 0, D1 = 0, D2 = 0, 
                              r_squared_linear = 0, adj_r_squared_linear = 0, residual_standard_error_linear = 0, 
                              slope_linear = 0, intercept_linear = 0, 
                              r_squared_non_linear = 0, adj_r_squared_non_linear = 0, residual_standard_error_non_linear = 0, 
                              slope_non_linear = 0, intercept_non_linear = 0)

  return(
    list(
      threshold_value = quantile_value, 
      row_regression_info = new_regression
    )
  )
}





GetFixedNoise = function(denoising_detail, denoising_method, transform_fun, spectra_distinct, ms_id, ms_rt) { 
  
  new_regression = data.frame(ms2_spectrum_id = ms_id, ms2_retention_time = ms_rt, 
                              
                              threshold_value = denoising_detail[[denoising_method]], 
                              
                              denoising_method = 'Fixed value', 
                              
                              spar = 0, RMSE = 0, LOOCV = 0, D1 = 0, D2 = 0, 
                              r_squared_linear = 0, adj_r_squared_linear = 0, residual_standard_error_linear = 0, 
                              slope_linear = 0, intercept_linear = 0, 
                              r_squared_non_linear = 0, adj_r_squared_non_linear = 0, residual_standard_error_non_linear = 0, 
                              slope_non_linear = 0, intercept_non_linear = 0)
  
  return(
    list(
      threshold_value = denoising_detail[[denoising_method]], 
      row_regression_info = new_regression
    )
  )
  
}










GetIsotopePatternCosineScore = function(matching_info, mol_list, threshold_iso, ms1_window_data, bin_w) {

  data(isotopes, package = "enviPat")
  data(adducts, package = "enviPat")
  
  ms1_window_tic = sum(ms1_window_data[['intensity']])
  ms1_window_data = dplyr::mutate(ms1_window_data, relative_intensity = intensity/ms1_window_tic)
  
  cosine_values = numeric(nrow(matching_info))
  
  for (i in seq_along(dim(matching_info)[1])) {
    
    current_row = matching_info[i, ]
    current_charge = current_row[['total_charge']]

    molecular_formula = NULL
    
    for (mol_name in names(mol_list)) {

      current_formula_num = current_row[[mol_name]] 
      
      if (current_formula_num == 0) {
        next
      }
      
      current_base_formula = mol_list[[mol_name]]
      current_comp_formula = enviPat::multiform(current_base_formula, current_formula_num)

      if (is.null(molecular_formula)) {
        molecular_formula = current_comp_formula
      } else {
        molecular_formula = enviPat::mergeform(molecular_formula, current_comp_formula)
      }
 
    }
    
    # isotopic distribution calculation
    molecular_formula = enviPat::check_chemform(isotopes, molecular_formula)
    new_molecular_formula = molecular_formula$new_formula
    
    molecular_iso_pattern = enviPat::isopattern(isotopes, new_molecular_formula, 
                                                charge = current_charge, rel_to = 2, threshold = threshold_iso, verbose = F)

    theoretical_iso_pattern = as.data.frame(molecular_iso_pattern[[new_molecular_formula]])

    mz_min = min(theoretical_iso_pattern[['m/z']], ms1_window_data[['mz']])
    mz_max = max(theoretical_iso_pattern[['m/z']], ms1_window_data[['mz']])
    
    centers = seq(from = mz_min, to = mz_max, by = bin_w)
    
    
    GetCenteredVector = function(mz_value, intensity_value, centers_value) {
      # attention!!!!!: mz_value, intensity_value must be in same order
      if (length(mz_value) != length(intensity_value)) {
        stop("mz_value and intensity_value must have the same length.")
      }
      
      center_idx_map = sapply(mz_value, function(m) {
        which.min(abs(centers_value - m))
      })
      
      center_vector = data.frame(mz = centers_value, intensity = rep(0, length(centers_value)))
      unique_map_idx = unique(center_idx_map)
      
      for (idx in unique_map_idx) {
        total_int = sum(intensity_value[center_idx_map == idx])
        center_vector$intensity[idx] = total_int
      }
      
      return(center_vector)
      
    }
    
    ms1_vector = GetCenteredVector(mz_value = ms1_window_data[['mz']], intensity_value = ms1_window_data[['relative_intensity']], centers_value = centers)
    theoretical_iso_vector = GetCenteredVector(mz_value = theoretical_iso_pattern[['m/z']], intensity_value = theoretical_iso_pattern[['abundance']], centers_value = centers)

    
    cosine_sim = function(a, b) {
      sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    }
    
    
    similarity_value = cosine_sim(ms1_vector$intensity, theoretical_iso_vector$intensity)
    
    
    cosine_value[i] = similarity_value

  }
  
  #matching_info = dplyr::slice_max(matching_info, cosine_value, n = 1, with_ties = T)
  
  matching_info <- matching_info %>%
    dplyr::mutate(score = cosine_value) %>% 
    dplyr::slice_max(score, n = 1, with_ties = FALSE)
  
  return(matching_info)
  

}


































