# Heading --------------------------------------------------------
# title: "Emerging water limitation for global ecosystem productivity driven by increased solar radiation"
# author: "Liyao Yu"
# date: "2026-04-16"
# Description: The script used for analysis of remote sensing (RS) data and for visualization in Fig. 2, 3, S2-S6, and S8-S9

# load dependent packages
library(terra)
library(tidyterra)
library(tidyverse)
library(data.table)
library(patchwork)
library(rnaturalearth)
library(sf)
library(tigris)
library(caret)
library(randomForest)
library(fastshap)
library(glue)
library(ggpmisc)
source("User_Defined_Functions.R")

# define land boundaries for the maps
land <- ne_download(scale = "medium", type = "land", category = "physical", returnclass = "sf")

# define map themes for Fig. 2 and 3
theme_map_Yu <- coord_sf(ylim = c(-60, 90), expand = FALSE) +  # cut to 90°N–60°S
  scale_x_continuous(
    name = NULL,
    breaks = seq(-180, 180, by = 60),
    labels = function(x) {
      sapply(x, function(xi) {
        if (xi < 0) paste0(abs(xi), "°W")
        else if (xi > 0) paste0(xi, "°E")
        else "0°"
      })
    }
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = seq(-60, 90, by = 30),
    labels = function(y) {
      sapply(y, function(yi) {
        if (yi < 0) paste0(abs(yi), "°S")
        else if (yi > 0) paste0(yi, "°N")
        else "0°"
      })
    }
  ) +
  guides(
    fill = guide_colorbar(
      ticks = TRUE,
      ticks.colour = "black",
      frame.colour = NA,
      barwidth = unit(0.3, "cm"),
      barheight = unit(3, "cm"),
      title.vjust = 1,
      label.position = "right",
      nbin = 30                     
    )
  ) +
  theme_Yu

# Read LCSPP data (Photosynthesis proxy) -------------------------------------------------------------------
files <- list.files(path = "LCSPP_resampled", pattern = "\\.nc$", full.names = TLUE)

LCSPP <- files %>%
  lapply(function(f) {
    r <- rast(f)
    r
  }) %>%
  do.call(c, .) %>% 
  `names<-`(as.character(1982:2023))

# set coordinate system
crs(LCSPP) <- "EPSG:4326"

# Read CRU Precipitation data (Pre) ----------------------------------------------------------------------
dir <- list.files("CRU/Pre", pattern = "\\.nc$", full.names = TLUE)

raster_list <- list()
year_list <- c()

for (file in dir) {
  year <- str_extract(file, "(?<=pre\\.)\\d{4}")
  
  if (is.na(year) || as.numeric(year) < 1982) next
  
  message(paste0("--- Reading CRU Pre data for ", year, " ---"))
  
  yearly_raster <- rast(file) %>% sum(na.rm = TLUE)
  names(yearly_raster) <- year
  
  raster_list[[length(raster_list) + 1]] <- yearly_raster
  year_list <- c(year_list, year)
}

CRU_Pre <- rast(raster_list)
names(CRU_Pre) <- year_list
crs(CRU_Pre) <- "EPSG:4326"

rm(raster_list)
rm(yearly_raster)

# Read CRU Solar radiation data (SW) ----------------------------------------------------------------------
dir <- list.files("CRU/SW", pattern = "\\.nc$", full.names = TLUE)

raster_list <- list()
year_list <- c()

for (file in dir) {
  year <- str_extract(file, "(?<=dswrf\\.)\\d{4}")
  
  if (is.na(year) || as.numeric(year) < 1982) next
  
  message(paste0("--- Reading CRU SW data for ", year, " ---"))
  
  yearly_raster <- rast(file) %>% sum(na.rm = TLUE)
  names(yearly_raster) <- year
  
  raster_list[[length(raster_list) + 1]] <- yearly_raster
  year_list <- c(year_list, year)
}

CRU_SW <- rast(raster_list)
names(CRU_SW) <- year_list
crs(CRU_SW) <- "EPSG:4326"

rm(raster_list)
rm(yearly_raster)

# Read Vegetation type and climate zone data -----------------------------------------------------------------
IGBP <- rast("IGBP.nc")
Koppen <- rast("koppen_geiger_0p00833333.tif")

# Resample SW, Pre, and IGBP according to LCSPP -----------------------------------------------------------------------
CRU_Pre_resampled <- CRU_Pre %>% 
  resample(LCSPP, method = "bilinear")

CRU_SW_resampled <- CRU_SW %>% 
  resample(LCSPP, method = "bilinear")

IGBP_resampled <- IGBP %>% resample(LCSPP[["1982"]], method = "near") # for `cat` values, use `near` resampling, the same for Koppen

Koppen_resampled <- Koppen %>% resample(LCSPP[["1982"]], method = "near")

# Calc LUE and PUE, their IAV, and LWII --------------------------------------------------------
LUE <- LCSPP / CRU_SW
LUE[is.infinite(LUE)] <- NA

PUE <- LCSPP / CRU_Pre
PUE[is.infinite(PUE)] <- NA

year_list_int <- year_list %>% as.integer()

LUE_CV_list <- list()
PUE_CV_list <- list()

for (i in 5:length(year_list_int)) { # calc using a 5-year moving window
  idx <- (i - 4):i
  
  LUE_window <- LUE[[idx]]
  
  LUE_mean <- app(LUE_window, mean, na.rm = TLUE)
  LUE_sd   <- app(LUE_window, sd,   na.rm = TLUE)
  
  LUE_CV <- LUE_sd / LUE_mean
  LUE_CV[is.infinite(LUE_CV)] <- NA
  
  names(LUE_CV) <- year_list[i]
  
  LUE_CV_list[[i - 4]] <- LUE_CV
}

LUE_CV_stack <- rast(LUE_CV_list)

for (i in 5:length(year_list_int)) {
  idx <- (i - 4):i
  
  PUE_window <- PUE[[idx]]
  
  PUE_mean <- app(PUE_window, mean, na.rm = TLUE)
  PUE_sd   <- app(PUE_window, sd,   na.rm = TLUE)
  
  PUE_CV <- PUE_sd / PUE_mean
  PUE_CV[is.infinite(PUE_CV)] <- NA
  
  names(PUE_CV) <- year_list[i]
  
  PUE_CV_list[[i - 4]] <- PUE_CV
}

PUE_CV_stack <- rast(PUE_CV_list)

diff_CV <- PUE_CV_stack - LUE_CV_stack 


# Calc the shift of CV between two periods ----------------------------------------------------------------------
period_1 <- 2000:2011
period_2 <- 2012:2023

shift_CV <- mean(diff_CV[[which(names(diff_CV) %in% as.character(period_2))]], na.rm = T) -
  mean(diff_CV[[which(names(diff_CV) %in% as.character(period_1))]], na.rm = T)

# filter out outliers
p99 <- quantile(values(shift_CV), probs = 0.99, na.rm = TRUE)
p1 <- quantile(values(shift_CV), probs = 0.01, na.rm = TRUE)
shift_CV <- mask(shift_CV, shift_CV > p99, maskvalue = TRUE)
shift_CV <- mask(shift_CV, shift_CV < p1, maskvalue = TRUE)


# Calc the shift and mean of SW and Pre -----------------------------------------------------------
shift_SW <- mean(CRU_SW_resampled[[which(names(CRU_SW_resampled) %in% as.character(period_2))]], na.rm = T) -
  mean(CRU_SW_resampled[[which(names(CRU_SW_resampled) %in% as.character(period_1))]], na.rm = T)

shift_Pre <- mean(CRU_Pre_resampled[[which(names(CRU_Pre_resampled) %in% as.character(period_2))]], na.rm = T) -
  mean(CRU_Pre_resampled[[which(names(CRU_Pre_resampled) %in% as.character(period_1))]], na.rm = T)

mean_SW <- mean(CRU_SW_resampled, na.rm = T)
mean_Pre <- mean(CRU_Pre_resampled, na.rm = T)

# Integrate all data -------------------------------------------- 
raster_sum <- c(shift_CV, shift_SW, shift_Pre,
                mean_SW, mean_Pre,
                IGBP_resampled, Koppen_resampled) %>% 
  `names<-` (c("shift_CV", "shift_SW", "shift_Pre", 
               "mean_SW", "mean_Pre",
               "IGBP", "Koppen"))

df_raster_sum <- raster_sum %>% 
  as.data.frame(xy = T) %>% 
  filter(if_all(everything(), Negate(is.na)))

# Convert IGBP values to category labels
IGBP_levels <- 1:16
IGBP_labels <- c("ENF", "EBF", "DNF", "DBF", "MF", "CS", "OS", "WSA", 
                 "SAV", "GRA", "WET", "CRO", "URB", "CVM", "SI", "BAR")

df_raster_sum$IGBP_cat <- factor(df_raster_sum$IGBP, 
                                 levels = IGBP_levels, 
                                 labels = IGBP_labels)

df_raster_sum$Koppen <- factor(df_raster_sum$Koppen)

df_raster_sum <- df_raster_sum %>% 
  filter(!(IGBP_cat %in% c("URB", "SI", "BAR")) & !is.na(IGBP_cat)) %>% 
  mutate(IGBP_cat = case_when(
    IGBP_cat %in% c("CS", "OS") ~ "SH",
    IGBP_cat %in% c("WSA") ~ "SAV",
    .default = IGBP_cat
  ))


# Plot Fig. 2a shift of LWII -----------------
fig_diff_CV_shift <-
  df_raster_sum %>%
  filter(!is.na(shift_CV)) %>%
  ggplot() +
  geom_sf(data = land, color = "grey", size = 0.05) +
  geom_tile(aes(x = x, y = y, fill = shift_CV)) +
  scale_fill_gradient2(
    low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0,
    breaks = seq(-0.3, 0.2, by = 0.1),
    labels = scales::label_number(accuracy = 0.1),
    limits = c(min(df_raster_sum$shift_CV, na.rm = T), 
               max(df_raster_sum$shift_CV, na.rm = T))
  ) +
  labs(fill = expression(Delta * "LWII")) +
  theme_map_Yu

ggsave("Fig_LWII_shift.pdf", width = 16, height = 8, units = "cm")

# sub panel
total_n <- df_raster_sum %>%
  filter(!is.na(shift_CV)) %>%
  nrow()

df_raster_sum %>% 
  filter(df_raster_sum >= 0.02 | df_raster_sum <= -0.02) %>%
  mutate(category = ifelse(df_raster_sum >= 0.02, "Light", "Water")) %>%
  count(category) %>%
  ggplot(aes(x = category, y = n, fill = category)) +
  geom_col(color = NA) +
  labs(x = "Limitation", y = "Count of pixels") +
  scale_fill_manual(values = c("#3C5488", "#E64B35")) +
  annotate("text",
           x = "Light",
           y = 9000,
           size = rel(1),
           label = "Delta*'LWII ='~'± 0.02'",
           parse = TRUE) +
  annotate("text",
           x = "Light",
           y = 9500,
           size = rel(1),
           label = "Threshold:") +
  theme_Yu

ggsave("Fig_LWII_count_0.02.pdf", width = 4, height = 4, units = "cm")  

# suppl. summary for hotspot regions
# define region boundary first
mean_df_raster_sum_au <- df_raster_sum %>% # australia
  filter(x > 113 & x < 154 & y > -39 & y < -10) %>%
  summarise(mean_shift_CV = mean(shift_CV, na.rm = T),
            mean_shift_SW = mean(shift_SW, na.rm = T))

mean_df_raster_sum_wna <- df_raster_sum %>% # west north america
  filter(x > -125 & x < -102 & y > 20 & y < 70) %>%
  summarise(mean_shift_CV = mean(shift_CV, na.rm = T),
            mean_shift_SW = mean(shift_SW, na.rm = T))

mean_df_raster_sum_sa <- df_raster_sum %>% # south africa
  filter(x > 16 & x < 33 & y > -35 & y < -26) %>%
  summarise(mean_shift_CV = mean(shift_CV, na.rm = T),
            mean_shift_SW = mean(shift_SW, na.rm = T))


# Plot Fig. 2b ---------------------------------------
day_to_second_coef <- 12 * 3600 # need to convert the unit of SW from W m-2 to J year-1

fig_SW_shift <-
  df_raster_sum %>%
  filter(!is.na(shift_SW)) %>%
  mutate(shift_SW = shift_SW / day_to_second_coef) %>% 
  ggplot() +
  geom_sf(data = land, color = "grey", size = 0.05) +
  geom_tile(aes(x = x, y = y, fill = shift_SW)) +
  scale_fill_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0,
    limits = c(min(df_raster_sum$shift_SW / day_to_second_coef, na.rm = T), 
               max(df_raster_sum$shift_SW / day_to_second_coef, na.rm = T))               
  ) +
  labs(fill = expression(Delta * "SW")) +
  theme_map_Yu

ggsave("Fig_SW_shift.pdf", width = 16, height = 8, units = "cm")

# subpanel
df_raster_sum %>% 
  filter(diff_CV < -0.02) %>% 
  mutate(diff_sign = ifelse(diff_SW > 0, "SW > 0", "SW < 0")) %>%
  count(diff_sign) %>%                          # count pixels
  mutate(pct = n / sum(n) * 100) %>%            # convert to %
  ggplot(aes(x = diff_sign, y = pct)) + 
  geom_col(fill = c("#0571b0", "#ca0020"), 
           color = NA, width = 0.6) +
  labs(x = " ", y = "Percent of area (%)") +
  theme_Yu

ggsave("Fig_SW_count_neg_LWII.pdf", width = 4, height = 4, units = "cm") 


# Random forest model ---------------------------------
set.seed(31)
df <- df_raster_sum[complete.cases(df_raster_sum), ] %>% 
  # filter((IGBP_cat %in% c("DBF"))) %>%  # model for different PFTs
  select(shift_CV, shift_SW, shift_Pre, mean_SW, mean_Pre, Koppen)
train <- sample(1:nrow(df), nrow(df) / 2)
df_test <- df[-train, ]

# hyperparameter tuning for random forest
rf_grid = expand.grid(
  mtry = c(2, 3, 4, 5)
)

control <- trainControl(
  method = "cv", 
  number = 5
) # k-fold, or number = 10

rf_model <- train(
  shift_CV ~ .,
  data = df[train, ],
  method = "rf",
  trControl = control,
  tuneGrid = rf_grid,
  ntree = 500
)

# the best hyper-parameters
rf_model$bestTune

rf_df <- randomForest(
  shift_CV ~ .,
  data = df[train, ],
  mtry = rf_model$bestTune$mtry,
  ntree = 500,
  importance = TRUE
)

randomForest::importance(rf_df)

varImpPlot(rf_df)

shift_CV_hat <- predict(rf_df,
                        newdata = df[-train, ])

mean((shift_CV_hat - df_test$shift_CV) ^ 2)

data.frame(x=shift_CV_hat, y=df_test$shift_CV) %>% 
  ggplot(aes(x, y)) + 
  geom_point() + 
  stat_smooth(method = "lm") +
  xlab("sim") + ylab("obs") +
  theme_Yu

# calculate r2 for train and test
(r2_rf_train <- cor(df$shift_CV, predict(rf_df, df[, -1]))^2)
(r2_rf_test <- cor(df_test$shift_CV, predict(rf_df, df_test[, -1]))^2)

# calculate variable importance
imp <- randomForest::importance(rf_df) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Variable") %>% 
  select(Variable, `%IncMSE`) %>%
  arrange(desc(`%IncMSE`)) %>% 
  mutate(imp_percentage=`%IncMSE`/sum(`%IncMSE`)*100,
         Variable=factor(Variable))


# Define a prediction function for the random forest model
predict_function <- function(model, newdata) {
  predict(model, newdata = newdata)
}

# Calculate SHAP values
shap_values_rf <- fastshap::explain(rf_df, 
                                    X = df[, -which(names(df) == "shift_CV")], 
                                    pred_wrapper = predict_function, 
                                    nsim = 50)

# Convert SHAP values to a data frame
shap_values_df <- as.data.frame(shap_values_rf) %>% 
  reshape2::melt()


# Plot Fig. 3a ---------------------------------------
ggplot(imp, aes(x = Variable, y = imp_percentage)) +
  geom_col(fill = "#E19979") +
  labs(x = "Variables", y = "Relative importance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
  theme_Yu
  
ggsave("Fig_RF_Importance.pdf", width = 6, height = 6, units = "cm")


# Plot Fig. 3b and S3 ---------------------------------------
var <- "shift_SW"

shap_partial <- data.frame(
  df %>% pull(var),
  shap_long[shap_values_df$variable == var, ] %>% pull("value")
) %>% 
  setNames(c(var, "value"))

year_to_second_coef <- 365 * 12 * 3600

shap_partial %>% 
  sample_frac(0.2) %>% # adjust for speed of plotting
  ggplot(aes(x = shift_SW, y = value)) + # shift_SW / year_to_second_coef
  geom_point(color = "darkgrey", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.2) +
  stat_smooth(color = "#E19979", method = "lm", se = T) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~~")),
    formula = y ~ x, size = 2, parse = TRUE) +
  labs(x = expression(Delta * "SW (W m-2)"), y = "SHAP Value") +
  ylim(-0.08, 0.08) +
  theme_Yu

ggsave("Fig_RF_Shap_SW.pdf", width = 6, height = 6, units = "cm")


# Read future climate map ------------------------------------------------------------
# SW and Pre
base_dir <- "CMIP6"

model_dirs <- list.files(base_dir, full.names = TRUE)
model_dirs <- model_dirs[file.info(model_dirs)$isdir]

Pre_list <- list()
SW_list <- list()

for (model_path in model_dirs) {
  model_name <- basename(model_path)
  
  ssp_dirs <- list.files(model_path, full.names = TRUE)
  ssp_dirs <- ssp_dirs[file.info(ssp_dirs)$isdir]
  
  for (ssp_path in ssp_dirs) {
    ssp_name <- basename(ssp_path)
    
    nc_files <- list.files(ssp_path, pattern = "\\.nc$", full.names = TRUE)
    
    # Separate pr and sw files
    pr_files <- nc_files[grepl("_pr_", basename(nc_files))]
    sw_files <- nc_files[grepl("_rsds_", basename(nc_files))]
    
    # Process pr files
    for (f in pr_files) {
      fname <- basename(f)
      year <- gsub(".*_([0-9]{4})_[0-9]{4}\\.nc$", "\\1", fname)
      key <- paste(model_name, toupper(ssp_name), "pr", year, sep = "_")
      print(key)
      Pre_list[[key]] <- mean(rast(f))
    }
    
    # Process sw files
    for (f in sw_files) {
      fname <- basename(f)
      year <- gsub(".*_([0-9]{4})_[0-9]{4}\\.nc$", "\\1", fname)
      key <- paste(model_name, toupper(ssp_name), "sw", year, sep = "_")
      print(key)
      SW_list[[key]] <- mean(rast(f))
    }
  }
}

# filter land area
SW_list_masked <- lapply(SW_list, function(r) {
  mask(r, land)
})

Pre_list_masked <- lapply(Pre_list, function(r) {
  mask(r, land)
})

# Koppen
folder_path <- "Koppen"
files <- list.files(folder_path, pattern = "^Koppen_.*\\.tif$", full.names = TRUE)

Koppen_list <- setNames(
  lapply(files, rast),
  tools::file_path_sans_ext(basename(files))
)

Koppen_list_resampled <- lapply(Koppen_list, function(x) {
  resample(x, LCSPP[["2000"]], method = "near")
})

names(Koppen_list_resampled) <- names(Koppen_list)

# define labels
models <- paste0("Model", sprintf("%02d", 1:5))
ssps <- paste0("SSP", sprintf("%02d", 1:3))
decades <- paste0("Decade", sprintf("%02d", 1:8))

# Create all combinations in the expected order
new_names <- c()

for (model in models) {
  for (ssp in ssps) {
    for (decade in decades) {
      new_names <- c(new_names, paste(model, ssp, decade, sep = "_"))
    }
  }
}

names(SW_list_masked) <- new_names
names(Pre_list_masked) <- new_names

# integrate decades and models
ssp <- sub(".*_(SSP\\d{2})_.*", "\\1", new_names)
decade <- as.integer(sub(".*_Decade(\\d{2})", "\\1", new_names))

merged_decade <- ((decade - 1) %/% 2) + 1
merged_decade <- sprintf("Decade%02d", merged_decade)

group_keys <- paste0(ssp, "_", merged_decade)

unique_groups <- unique(group_keys)

mean_layers <- lapply(unique_groups, function(g) {
  idx <- which(group_keys == g)
  mean(SW_list_masked[[idx]])
})

SW_list_sum <- rast(mean_layers)

mean_layers <- lapply(unique_groups, function(g) {
  idx <- which(group_keys == g)
  mean(Pre_list_masked[[idx]])
})

Pre_list_sum <- rast(mean_layers)

names(SW_list_sum) <- unique_groups
names(Pre_list_sum) <- unique_groups

# conver J s-1 to J yr-1 for SW
SW_list_sum <- SW_list_sum * 24 * 3600 * 365

# convert kg m-2 s-1 to mm yr-1 for Pre
Pre_list_sum <- Pre_list_sum * 24 * 3600 * 365


# merge future SW and Pre into main dataset 
df_raster_sum_future <- df_raster_sum %>% 
  mutate(x = round(x, 2), y = round(y, 2))

for (ssp in 1:3) {
  for (decade in 1:4) {
    name_r <- glue("SSP0{ssp}_Decade0{decade}")
    r1 <- SW_list_sum[[name_r]] %>% 
      `names<-`(glue("SW_{name_r}")) %>% 
      as.data.frame(xy = T) %>% 
      mutate(x = round(x, 2), y = round(y, 2))
    r2 <- Pre_list_sum[[name_r]] %>% 
      `names<-`(glue("Pre_{name_r}")) %>% 
      as.data.frame(xy = T) %>% 
      mutate(x = round(x, 2), y = round(y, 2))
    
    df_raster_sum_future <- df_raster_sum_future %>% 
      left_join(r1, by = c("x", "y")) %>% 
      left_join(r2, by = c("x", "y"))
  }
}

df_raster_sum_future <- df_raster_sum_future %>% 
  mutate(
    shift_SW_1_SSP1 = SW_SSP01_Decade01 - mean_SW,
    shift_SW_1_SSP2 = SW_SSP02_Decade01 - mean_SW,
    shift_SW_1_SSP3 = SW_SSP03_Decade01 - mean_SW,
    shift_Pre_1_SSP1 = Pre_SSP01_Decade01 - mean_Pre,
    shift_Pre_1_SSP2 = Pre_SSP02_Decade01 - mean_Pre,
    shift_Pre_1_SSP3 = Pre_SSP03_Decade01 - mean_Pre,
    
    shift_SW_2_SSP1 = SW_SSP01_Decade02 - SW_SSP01_Decade01,
    shift_SW_2_SSP2 = SW_SSP02_Decade02 - SW_SSP02_Decade01,
    shift_SW_2_SSP3 = SW_SSP03_Decade02 - SW_SSP03_Decade01,
    shift_Pre_2_SSP1 = Pre_SSP01_Decade02 - Pre_SSP01_Decade01,
    shift_Pre_2_SSP2 = Pre_SSP02_Decade02 - Pre_SSP02_Decade01,
    shift_Pre_2_SSP3 = Pre_SSP03_Decade02 - Pre_SSP03_Decade01,
    
    shift_SW_3_SSP1 = SW_SSP01_Decade03 - SW_SSP01_Decade02,
    shift_SW_3_SSP2 = SW_SSP02_Decade03 - SW_SSP02_Decade02,
    shift_SW_3_SSP3 = SW_SSP03_Decade03 - SW_SSP03_Decade02,
    shift_Pre_3_SSP1 = Pre_SSP01_Decade03 - Pre_SSP01_Decade02,
    shift_Pre_3_SSP2 = Pre_SSP02_Decade03 - Pre_SSP02_Decade02,
    shift_Pre_3_SSP3 = Pre_SSP03_Decade03 - Pre_SSP03_Decade02,
    
    shift_SW_4_SSP1 = SW_SSP01_Decade04 - SW_SSP01_Decade03,
    shift_SW_4_SSP2 = SW_SSP02_Decade04 - SW_SSP02_Decade03,
    shift_SW_4_SSP3 = SW_SSP03_Decade04 - SW_SSP03_Decade03,
    shift_Pre_4_SSP1 = Pre_SSP01_Decade04 - Pre_SSP01_Decade03,
    shift_Pre_4_SSP2 = Pre_SSP02_Decade04 - Pre_SSP02_Decade03,
    shift_Pre_4_SSP3 = Pre_SSP03_Decade04 - Pre_SSP03_Decade03
  )

# Plot Fig. 3c ---------------------------------------
# Future change in SW (e.g. until 2040 under SSP370)
fig_diff_SW_shift_future_period1_ssp2 <-
  df_raster_sum_future %>%
  mutate(shift_SW_1_SSP2 = shift_SW_1_SSP2 / year_to_second_coef) %>% 
  ggplot() +
  geom_sf(data = land, color = "grey", size = 0.05) +
  geom_tile(aes(x = x, y = y, fill = shift_SW_1_SSP2)) +
  scale_fill_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0,
    limits = c(min(df_raster_sum_future$shift_SW_1_SSP2 / (year_to_second_coef * 2), na.rm = T), 
               max(df_raster_sum_future$shift_SW_1_SSP2 / (year_to_second_coef * 2), na.rm = T))  
  ) +
  labs(fill = " ") +
  theme_map_Yu

fig_diff_SW_shift_future_period1_ssp2
ggsave("Fig_diff_SW_shift_future_period1_ssp2.pdf", width = 16, height = 8, units = "cm")

# Fig.3c sub panel
# area with change in SW
df_shift_SW_future <- df_raster_sum_future %>%
  select(x, y, starts_with("shift_SW_")) %>%   # keep only x, y, and shift_SW_* cols
  pivot_longer(
    cols = starts_with("shift_SW_"),
    names_to = c("decades", "ssps"),
    names_pattern = "shift_SW_(\\d+)_(SSP\\d+)",
    values_to = "shift_SW"
  ) %>%
  mutate(decades = factor(decades),
         ssps = factor(ssps))

area_shift_SW_future <- df_shift_SW_future %>%
  filter(shift_SW > 50000000 | shift_SW < -50000000) %>%  # filter out pixels without substantial change in SW
  group_by(ssps, decades) %>%
  summarise(
    pixel_pos = sum(shift_SW > 0, na.rm = TRUE),
    pixel_neg = sum(shift_SW < 0, na.rm = TRUE),
    .groups = "drop"
  )

area_shift_SW_future %>% 
  ggplot(aes(x = decades, y = pixel_pos/nrow(area_shift_SW_future), group = ssps)) +
  geom_line(aes(color = ssps)) +
  ylab("Percent of area (%)") + xlab("Year") +
  annotate("text", x = 3, y = 0.2, 
           color = "#1b9e77", label = "SSP585") +
  annotate("text", x = 3, y = 0.15, 
           color = "#d95f02", label = "SSP370") +
  annotate("text", x = 3, y = 0.1, 
           color = "#7570b3", label = "SSP126") +
  scale_x_discrete(labels = c(2040, 2060, 2080, 2100)) +
  scale_color_manual(labels = c("SSP126", "SSP370", "SSP585"),
                     values = c("#7570b3", "#d95f02", "#1b9e77")) +
  theme_Yu

ggsave("Fig_area_shift_SW_future_.pdf", width = 6, height = 6, units = "cm")

# Future projection and Plot Fig. 3d and Fig. S4-6 ---------------------------------------
area_shift_future <- data.frame(
  ssps = NA,
  decades = NA,
  pixel = NA
)

for (ssp in 1:3) {
  for (decade in 1:4) {
    temp <- df_raster_sum_future %>% 
      mutate(shift_CV_predict = predict(rf.df,
                                        newdata = df_raster_sum_future %>% 
                                          select(!!sym(glue("shift_SW_{decade}_SSP{ssp}")), 
                                                 !!sym(glue("shift_Pre_{decade}_SSP{ssp}")), 
                                                 !!sym(glue("SW_SSP0{ssp}_Decade0{decade}")), 
                                                 !!sym(glue("Pre_SSP0{ssp}_Decade0{decade}")), 
                                                 IGBP_cat, Koppen) %>% 
                                          mutate(shift_SW = !!sym(glue("shift_SW_{decade}_SSP{ssp}")), 
                                                 shift_Pre = !!sym(glue("shift_Pre_{decade}_SSP{ssp}")), ,
                                                 mean_SW = !!sym(glue("SW_SSP0{ssp}_Decade0{decade}")), 
                                                 mean_Pre = !!sym(glue("Pre_SSP0{ssp}_Decade0{decade}"))))) %>%
      select(x, y, shift_CV_predict)
    
    temp %>% 
      ggplot() +
      geom_sf(data = land, color = "grey", size = 0.05) +
      geom_tile(aes(x = x, y = y, fill = shift_CV_predict)) +
      scale_fill_gradient2(
        low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0,
        breaks = seq(-0.5, 0.5, by = 0.1),
        labels = scales::label_number(accuracy = 0.1),
        limits = c(min(temp$shift_CV_predict, na.rm = T), max(temp$shift_CV_predict, na.rm = T))
      ) +
      labs(fill = " ") +
      theme_map_Yu
    
    ggsave(glue("Fig_shift_limitation_Period{decade}_SSP{ssp}.pdf"),
           width = 16, height = 8, units = "cm")
    
    area_shift_future <- area_shift_future %>% bind_rows(data.frame(
      ssps = ssp,
      decades = decade,
      pixel = temp %>% filter(shift_CV_predict < 0) %>% nrow(),
      pixel_threshold = temp %>% filter(shift_CV_predict < -0.02) %>% nrow()
    ))
  }
}

# Fig.3d sub panel
# area with shift to water limitation
sum_area_shift_future <- area_shift_future %>% 
  mutate(ssps = factor(ssps),
         decades = factor(decades)) %>% 
  filter(!is.na(pixel)) 

sum_area_shift_future %>% 
  ggplot(aes(x = decades, y = pixel_threshold, group = ssps)) +
  geom_line(aes(color = ssps)) +
  ylab("Pixels") + xlab("Year") +
  annotate("text", x = 4, y = 1050, 
           color = "#1b9e77", label = "SSP585") +
  annotate("text", x = 4, y = 950, 
           color = "#d95f02", label = "SSP370") +
  annotate("text", x = 4, y = 600, 
           color = "#7570b3", label = "SSP126") +
  scale_x_discrete(labels = c(2040, 2060, 2080, 2100)) +
  scale_color_manual(labels = c("SSP126", "SSP370", "SSP585"),
                     values = c("#7570b3", "#d95f02", "#1b9e77")) +
  theme_Yu

ggsave("Fig_area_shift_LWII_future_.pdf", width = 6, height = 6, units = "cm")

# correlation for sub panels Fig. 3c and Fig. 3d
corr_area_future_SW_LWII <- data.frame(
  d_SW = area_shift_SW_future$pixel_pos,
  d_LWII = sum_area_shift_future$pixel_threshold
)

cor(corr_area_future_SW_LWII$d_SW, corr_area_future_SW_LWII$d_LWII)


# Plot Fig. S2a ---------------------------------------
# Trend of LWII, SW, and Pre

# extract the time series of LWII, SW, and Pre
df_diff_CV_time_series <- diff_CV %>% 
  as.data.frame(xy = T) %>% 
  pivot_longer(cols = 3: 36, names_to = "Year", values_to = "diff_CV") %>% 
  mutate(Year = Year %>%  as.numeric())

df_SW_time_series <- CRU_SW %>% 
  as.data.frame(xy = T) %>% 
  pivot_longer(cols = 3: 40, names_to = "Year", values_to = "SW") %>% 
  mutate(Year = substr(Year, 8, 11) %>% as.numeric())

df_Pre_time_series <- CRU_Pre %>% 
  as.data.frame(xy = T) %>% 
  pivot_longer(cols = 3: 40, names_to = "Year", values_to = "Pre") %>% 
  mutate(Year = substr(Year, 9, 12) %>% as.numeric())

df_sum_time_series <- df_SW_time_series %>% 
  left_join(df_Pre_time_series, by = c("x", "y", "Year")) %>% 
  left_join(df_diff_CV_time_series, by = c("x", "y", "Year")) %>% 
  rename(lon = x,
         lat = y)

df_sum_time_series_new <- df_sum_time_series %>% 
  filter(Year >= 2000) %>% 
  mutate(Year_group = (Year - 2000) %/% 3 + 1) %>% 
  group_by(lon, lat, Year_group) %>% 
  summarise(mean_diff_CV = mean(diff_CV, na.rm = T),
            mean_SW = mean(SW, na.rm = T),
            mean_Pre = mean(Pre, na.rm = T)) %>% 
  group_by(lon, lat) %>% 
  arrange(Year_group, .by_group = TRUE) %>% 
  mutate(diff_LWII = mean_diff_CV - lag(mean_diff_CV),
         diff_SW = mean_SW - lag(mean_SW),
         diff_Pre = mean_Pre - lag(mean_Pre),
         diff_LWII_ratio = diff_LWII / mean_diff_CV,
         diff_SW_ratio = diff_SW / mean_SW,
         diff_Pre_ratio = diff_Pre / mean_Pre) %>% 
  ungroup() 

area_neg_LWII_SW_Pre <- df_sum_time_series_new %>% 
  filter(if_all(everything(), ~ !is.na(.))) %>% 
  group_by(Year_group) %>% 
  summarise(Area_neg_LWII = sum(diff_LWII < 0),
            Area_pos_SW = sum(diff_SW_ratio > 0),
            Area_neg_Pre = sum(diff_Pre_ratio < 0)) %>% 
  ungroup()

nrow_df_sum_time_series <- (df_sum_time_series_new %>% 
                              filter(if_all(everything(), ~ !is.na(.))) %>% 
                              nrow() / 7) %>% round(0)

area_neg_LWII_SW_Pre %>% 
  ggplot(aes(x = Year_group)) +
  geom_line(aes(y = Area_neg_LWII / nrow_df_sum_time_series * 100), color = "black") + 
  geom_line(aes(y = Area_pos_SW / nrow_df_sum_time_series * 100), color = "#fc8d62") +
  geom_line(aes(y = Area_neg_Pre / nrow_df_sum_time_series * 100), color = "#8da0cb") +
  geom_hline(yintercept = 50, color = "grey", linetype = "dashed", linewidth = rel(0.5)) +
  scale_x_continuous(labels = c("2000-\n2003", "2004-\n2007", "2008-\n2011",
                                "2012-\n2015", "2016-\n2019", "2020-\n2023")) +
  xlab("Year") + ylab("Percentage of area (%)") +
  theme_Yu

ggsave("Fig_area_neg_LWII_SW_Pre_.pdf", width = 10, height = 6, units = "cm")

# Plot Fig. S2b --------------------------------------
# count how many pixels with no decrease in Pre but have negative LWII and increase in SW 
df_raster_sum %>% 
  filter(shift_CV < -0.02) %>% 
  mutate(diff_sign = case_when(
    shift_SW > 0 & shift_Pre > 0 ~ "ΔSW > 0\nΔPre > 0",
    shift_SW > 0 & shift_Pre <= 0 ~ "ΔSW > 0\nΔPre < 0",
    shift_SW <= 0 & shift_Pre > 0 ~ "ΔSW < 0\nΔPre > 0",
    shift_SW <= 0 & shift_Pre <= 0 ~ "ΔSW < 0\nΔPre < 0"
  )) %>%
  count(diff_sign) %>%                          # count pixels
  mutate(pct = n / sum(n) * 100) %>% 
  ggplot(aes(x = diff_sign, y = pct)) + 
  geom_col(fill = "#ca0020", 
           color = NA, width = 0.6) +
  labs(x = " ", y = "Percent of area (%)") +
  theme_Yu

ggsave("Fig_neg_LWII_count_quartile.pdf", width = 6, height = 6, units = "cm")  

# Plot Fig. S8a ---------------------------------------
# density distribution of LWII

df_raster_sum %>% 
  ggplot(aes(shift_CV)) +
  geom_density(linewidth = rel(0.2)) +
  geom_vline(linetype = "dashed", linewidth = rel(0.2),
             xintercept = mean(df_raster_sum$shift_CV, na.rm = T)) +
  xlab("LWII") +
  theme_Yu

ggsave("Fig_LWII_density.pdf", width = 4, height = 4, units = "cm")  

# Plot Fig. S8b ---------------------------------------
# count of shift to water limiatation (negative ΔLWII) with varying thresholds (0.01-0.05)

df_raster_sum %>% 
  filter(df_raster_sum >= 0.01 | df_raster_sum <= -0.01) %>%
  mutate(category = ifelse(df_raster_sum >= 0.01, "Light", "Water")) %>%
  count(category) %>%
  ggplot(aes(x = category, y = n, fill = category)) +
  geom_col(color = NA) +
  labs(x = "Limitation", y = "Count of pixels") +
  scale_fill_manual(values = c("#3C5488", "#E64B35")) +
  annotate("text",
           x = "Light",
           y = 9000,
           size = rel(1),
           label = "Delta*'LWII ='~'± 0.01'",
           parse = TRUE) +
  annotate("text",
           x = "Light",
           y = 9500,
           size = rel(1),
           label = "Threshold:") +
  theme_Yu

ggsave("Fig_LWII_count_0.01.pdf", width = 4, height = 4, units = "cm") 

# Plot Fig. S8c ---------------------------------------
df_raster_sum %>% 
  filter(df_raster_sum >= 0.05 | df_raster_sum <= -0.05) %>%
  mutate(category = ifelse(df_raster_sum >= 0.05, "Light", "Water")) %>%
  count(category) %>%
  ggplot(aes(x = category, y = n, fill = category)) +
  geom_col(color = NA) +
  labs(x = "Limitation", y = "Count of pixels") +
  scale_fill_manual(values = c("#3C5488", "#E64B35")) +
  annotate("text",
           x = "Light",
           y = 9000,
           size = rel(1),
           label = "Delta*'LWII ='~'± 0.05'",
           parse = TRUE) +
  annotate("text",
           x = "Light",
           y = 9500,
           size = rel(1),
           label = "Threshold:") +
  theme_Yu

ggsave("Fig_LWII_count_0.05.pdf", width = 4, height = 4, units = "cm") 

# Plot Fig. S9 ---------------------------------------
yrs <- as.numeric(names(RUE))
yrs <- as.numeric(names(PUE))

LUE_CV_all_years <- app(
  LUE[[yrs >= 2000]],
  fun = function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA)
    
    q <- quantile(x, probs = c(0.05, 0.95), na.rm = TRUE) # exclude extreme values
    x_filtered <- x[x >= q[1] & x <= q[2]]
    
    if (length(x_filtered) < 2) return(NA)
    
    sd(x_filtered) / mean(x_filtered)
  }
)

PUE_CV_all_years <- app(
  PUE[[yrs >= 2000]],
  fun = function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA)
    
    q <- quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)
    x_filtered <- x[x >= q[1] & x <= q[2]]
    
    if (length(x_filtered) < 2) return(NA)
    
    sd(x_filtered) / mean(x_filtered)
  }
)

diff_CV_mean_post2000 <- PUE_CV_all_years - LUE_CV_all_years

names(diff_CV_mean_post2000) <- "mean" 

df_diff_CV_mean_post2000 <- as.data.frame(diff_CV_mean_post2000,
                                                xy = T) %>% 
  filter(if_all(everything(), Negate(is.na)))

fig_mean_diff_CV <- df_diff_CV_mean_post2000 %>%
  filter(!is.na(mean)) %>%
  ggplot() +
  geom_sf(data = land, color = "grey", size = 0.05) +
  geom_tile(aes(x = x, y = y, fill = mean)) +
  scale_fill_gradient2(
    low = "#ca0020", mid = "white", high = "#0571b0", midpoint = 0,
    breaks = seq(-0.3, 0.3, by = 0.05),
    labels = scales::label_number(accuracy = 0.05),
    limits = c(min(df_diff_CV_mean_post2000_clean_filtered$mean, na.rm = T),
               max(df_diff_CV_mean_post2000_clean_filtered$mean, na.rm = T))
  ) +
  labs(fill = "Mean LWII") +
  theme_map_Yu
  

fig_mean_LWII
ggsave("Fig_mean_LWII.pdf", width = 16, height = 8, units = "cm")
