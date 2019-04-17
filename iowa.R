library(readr)
library(readxl)
library(lmom)

stations_table <- read_csv("C:/Users/q0hecgsk/Downloads/NIST_Extreme_Winds/stations_analyzed.csv")
iowa_stations <- stations_table[stations_table$STATE == "IA",]
iowa_stations <- iowa_stations[complete.cases(iowa_stations),]
unique_usaf_id <- unique(iowa_stations$USAF)
setwd("C:/Users/q0hecgsk/Downloads/NIST_Extreme_Winds/final_qc_data")

parse_NIST_excel_file <- function(USAF_id) {
  fname <- paste0("station_matrix_",USAF_id,".xlsx")
  df <- read_excel(fname, sheet = "Sheet1", skip = 6)
  colnames(df) <- c("DateTime", "WindSpeedRaw", "WindSpeedStandard", "WindDirection",
                    "NonThunderstorm", "IndependentNonThunderstorm",
                    "Thunderstorm", "IndependentThunderstorm",
                    "AllIndependent", "Tropical",
                    "AMaxNT", "AMaxT", "AMaxAll")
  df
}

data_list <- lapply(unique_usaf_id, parse_NIST_excel_file)

extract_thunderstorm_pds <- function(NIST_df, threshold = 42) {
  x <- NIST_df[NIST_df$IndependentThunderstorm == 1,]
  x[x$WindSpeedStandard >= threshold,]
}

t_pds <- lapply(data_list, extract_thunderstorm_pds)

compute_NIST_df_l_moments <- function(NIST_df) {
  x <- NIST_df[["WindSpeedStandard"]]
  samlmu(x, nmom = 5)
}

t_pds_lmom <- lapply(t_pds, compute_NIST_df_l_moments)

lmrd_points_from_list <- function(lmom_list, new_diagram = F) {
  if(new_diagram) lmrd()
  x <- lapply(lmom_list, lmrdpoints)
}

lmrd_points_from_list(t_pds_lmom, new_diagram = T)

extract_thunderstorm_ams <- function(NIST_df) {
  NIST_df[NIST_df$AMaxT == 1,]
}

t_ams <- lapply(data_list, extract_thunderstorm_ams)
t_ams_lmom <- lapply(t_ams, compute_NIST_df_l_moments)
lmrd_points_from_list(t_ams_lmom, new_diagram = F)
  
extract_nonthunderstorm_pds <- function(NIST_df, threshold = 42) {
  x <- NIST_df[NIST_df$IndependentNonThunderstorm == 1,]
  x[x$WindSpeedStandard >= threshold,]
}
extract_nonthunderstorm_ams <- function(NIST_df) {
  NIST_df[NIST_df$AMaxNT == 1,]
}
nt_pds <- lapply(data_list, extract_nonthunderstorm_pds)
nt_pds_lmom <- lapply(nt_pds, compute_NIST_df_l_moments)
lmrd_points_from_list(nt_pds_lmom, new_diagram = T)
nt_ams <- lapply(data_list, extract_nonthunderstorm_ams)
nt_ams_lmom <- lapply(nt_ams, compute_NIST_df_l_moments)
lmrd_points_from_list(nt_ams_lmom, new_diagram = F)

compute_sample_size_from_series <- function(df_list) {
  unlist(lapply(df_list, nrow))
}
t_ams_sizes <- compute_sample_size_from_series(t_ams)
t_pds_sizes <- compute_sample_size_from_series(t_pds)
t_ev_per_year <- t_pds_sizes / t_ams_sizes

plot_NIST_time_series <- function(NIST_df) {
  with(NIST_df, plot(WindSpeedStandard ~ DateTime))
}

gev_xi <- function(q0, alpha, lambda, kappa) {
  return(q0 + (alpha/kappa) * (1 - lambda ^ (-kappa)))
}
gev_alpha <- function(alpha, lambda, kappa) {
  return(alpha * lambda ^ (-kappa))
}

regional_l_moments_from_list <- function(df_list) {
  get_t <- function(lmom_df) {
    unlist(lapply(lmom_df, function(x) unname(x["l_2"] / x["l_1"])))
  }
  get_t3 <- function(lmom_df) {
    unlist(lapply(lmom_df, function(x) unname(x["t_3"])))
  }
  get_t4 <- function(lmom_df) {
    unlist(lapply(lmom_df, function(x) unname(x["t_4"])))
  }
  n <- compute_sample_size_from_series(df_list)
  x_lmom <- lapply(df_list, compute_NIST_df_l_moments)
  t_vec <- get_t(x_lmom)
  t3_vec <- get_t3(x_lmom)
  t4_vec <- get_t4(x_lmom)
  t_wt <- weighted.mean(t_vec, n)
  t3_wt <- weighted.mean(t3_vec, n)
  t4_wt <- weighted.mean(t4_vec, n)
  c("l_1" = 1, "l_2" = t_wt, "t_3" = t3_wt, "t_4" = t4_wt)
}

t_pds_kappa <- pelkap(regional_l_moments_from_list(t_pds))
t_mean_storm_rate <- weighted.mean(t_ev_per_year, t_ams_sizes)

compute_theoretical_ams_distribution <- function(pds, rate) {
  pds_kappa <- pelkap(regional_l_moments_from_list(pds))
  theoretical_xi <- gev_xi(pds_kappa["xi"], pds_kappa["alpha"], rate, pds_kappa["k"])
  theoretical_alpha <- gev_alpha(pds_kappa["alpha"], rate, pds_kappa["k"])
  c("xi" = unname(theoretical_xi), "alpha" = unname(theoretical_alpha), "k" = unname(pds_kappa["k"]), "h" = 1/rate)
}

test_case <- t_ams[[4]]
test_case_lmom <- samlmu(test_case$WindSpeedStandard, nmom = 5)
test_case_l1 <- unname(test_case_lmom["l_1"])
test_case_unitized_series <- unname(test_case$WindSpeedStandard)/test_case_l1
