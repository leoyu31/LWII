# Heading --------------------------------------------------------
# title: "Emerging water limitation for global ecosystem productivity driven by increased solar radiation"
# author: "Liyao Yu"
# date: "2026-04-16"
# Description: The script used for analysis of eddy covariance (EC) data and for visualization in Fig. 1, S1, and S7

# load dependent packages
library(tidyverse)
library(patchwork)
library(data.table)
library(readxl)
library(rstudioapi)
library(magrittr)
library(ggsci)
library(ggpubr)
library(zoo)
library(ggpmisc)
library(tidyplots)
source("User_Defined_Functions.R")

# Read annual flux data ------------------------------------------------
file_path <- file.path(paste0(getwd(), "/Data_YY/Gross/"))

temp = list.files(file_path, pattern = '*.csv')

for (i in 1:length(temp)){
  print(paste0("--- Reading data of ", temp[i], " --- ", i))
  
  assign(temp[i],
         as.data.frame(lapply(data.table::fread(paste0(file_path, temp[i]), 
                                                         header = TRUE), 
                                       as.numeric))
  )
}

# add a col of site name
for (i in 1:length(temp)){
  assign(temp[i], 
         transform(get(temp[[i]]), 
                   Site = substr(temp[[i]], 1, 6)))
}

# Merge all sites data into one dataframe
dfs = sapply(.GlobalEnv, is.data.frame)

data_YY <- do.call(bind_rows, mget(names(dfs)[dfs]))

# add site info of PFT and location
Site_info <- data.table::fread("Site.info.csv")


Site_loc <- Site_info %>% select(Site, Long, Lat, PFT) %>% 
  rename(lat = Lat,
         long = Long) %>% 
  mutate(long = ifelse(long < 0, long + 360, long))


# label sites as in Europe or North America
library(rnaturalearth)
library(sf)

# Load world map with continent info
world <- ne_countries(scale = "medium", returnclass = "sf")

site_continent <- site.loc %>% 
  mutate(Continent = mapply(get_continent, lat, long)) %>% 
  select(Site, Continent)

site_continent$Continent[site.continent$Site == "GF-Guy"] <- "South America" # French Guiana is wrongly classified into Europe
site_continent$Continent[!site.continent$Continent %in% c("Europe", "North America")] <- "Other"

# select cols of interest and filter data
data_YY_filtered <- data_YY %>% 
  mutate(Year = TIMESTAMP) %>% 
  select(Site, Year, GPP_DT_VUT_REF, LE_F_MDS,  SW_IN_F, VPD_F, TA_F, P_F, CO2_F_MDS) %>% 
  rename(GPP = GPP_DT_VUT_REF,
         LE = LE_F_MDS,
         SW = SW_IN_F,
         TA = TA_F,
         Precip = P_F,
         VPD = VPD_F,
         CO2 = CO2_F_MDS) %>% 
  filter(GPP > 0, LE > 0, VPD > 0, SW > 0, Precip > 0) %>% 
  left_join(Site_loc, by = "Site")


Sites_long <- data_YY_filtered %>% 
  group_by(Site) %>% 
  summarise(nYear = n()) %>% 
  ungroup() %>% 
  filter(nYear > 14) %>% 
  pull(Site)

Site_loc <- Site_loc %>% 
  filter(Site %in% Sites_long)

# Calculate annual light use efficiency (LUE) and precipitation use efficiency (PUE) ------------------
data_YY_filtered <- data_YY_filtered %>% 
  mutate(LUE = GPP / SW,
         PUE = GPP / Precip)


# Calculate inter-annual variability (IAV) -----------------------
# Calculate the rolling mean and sd in 5-year moving window 
data_YY_rolling_mean <- data_YY_filtered %>%
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(
    across(c(LUE, PUE),
           ~ roll_stats(.x, 5),
           .names = "{.col}.{.fn}.5"),
    
    across(c(SW, Precip),
           ~ rollapply(.x, 5, mean, na.rm = TRUE, align = "right", fill = NA),
           .names = "{.col}_mean.rolling.5")
  )


# Calculate LWII using CV_PUE - CV_LUE ---------------------
data_cv <- data_YY_rolling_mean %>%  
  select(Site,
         Year,
         LUE.rolling.CV.5,
         WUE.rolling.CV.5,
         SW_mean.rolling.5,
         Precip_mean.rolling.5) %>% 
  pivot_longer(cols = 3:6,
               names_to = "RUE",
               values_to = "CV") %>% 
  mutate(RUE = str_extract(RUE, "^[^.]*")) %>% 
  pivot_wider(names_from = RUE,
              values_from = CV) %>% 
  left_join(Site_info, by = "Site") %>% 
  mutate(diff_CV = PUE - LUE) %>% 
  filter(diff_CV < 0.5 & diff_CV > -0.5) # exclude outliers for `US-Ich` in 2011 - 2012

# Calculate the shift in LWII -----------------------------------
data_cv_shift <- data_cv %>% 
  filter(!is.na(diff_CV)) %>% 
  group_by(Site) %>% 
  mutate(nYear = n(),
         row_ID = row_number()) %>% 
  summarise(shift_CV = mean(diff_CV[row_ID >= floor(nYear / 2)], na.rm = T) -
              mean(diff_CV[row_ID < floor(nYear / 2)], na.rm = T),
            shift_SW = mean(SW_mean[row_ID >= floor(nYear / 2)], na.rm = T) -
              mean(SW_mean[row_ID < floor(nYear / 2)], na.rm = T),
            shift_Pre = mean(Precip_mean[row_ID >= floor(nYear / 2)], na.rm = T) -
              mean(Precip_mean[row_ID < floor(nYear / 2)], na.rm = T)) %>% 
  left_join(Site_info %>% select(Site, Long, Lat, PFT), by = "Site") %>% 
  mutate(Continent = case_when(
    Long >= -170 & Long <= -35 ~ "AM",
    Long >= -25 & Long <= 60 ~ "EU",
    TRUE ~ "Other"
  )) %>% select(-Long, -Lat) %>% 
  filter(!(Site %in% c("IL-Yat", "GF-Guy")))

data_cv_shift %>% 
  group_by(Continent) %>% 
  summarise(mean_shift_CV = mean(shift_CV, na.rm = T),
            median_shift_CV = median(shift_CV, na.rm = T))

# Plot Fig. 1b ------------------------------------
data_cv_shift %>% 
  ggplot(aes(x = Continent, y = shift_CV, group = Continent, fill = Continent)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.2) +
  geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA) +
  geom_jitter(shape = 21, color = "black",
              width = 0.2, size = rel(0.5), stroke = rel(0.5)) +
  scale_fill_manual(values = c("AM" = "#F49C9B",
                               "EU" = "#ABCEE2")) +
  scale_x_discrete(labels = c("North America", "Europe")) +
  ylab(expression(Delta * "LWII")) +
  theme_Yu

ggsave("EC_shift_LWII.pdf", width = 6, height = 4, units = "cm")  

# Plot Fig. 1c ----------------------------------
data_cv_shift %>% 
  group_by(Continent) %>% 
  summarise(mean_shift_SW = mean(shift_SW, na.rm = T),
            median_shift_SW = median(shift_SW, na.rm = T))

data_cv_shift %>% 
  ggplot(aes(x = Continent, y = shift_SW, group = Continent, fill = Continent)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.2) +
  geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA) +
  geom_jitter(shape = 21, color = "black",
              width = 0.2, size = rel(0.5), stroke = rel(0.5)) +
  scale_fill_manual(values = c("AM" = "#F49C9B",
                               "EU" = "#ABCEE2")) +
  scale_x_discrete(labels = c("North America", "Europe")) +
  ylab(expression(Delta * "SW (" * mu * "mol" ~ m^-2 ~ s^-1 * ")")) +
  theme_Yu

ggsave("EC_shift_SW.pdf", width = 6, height = 4, units = "cm")


# Plot Fig. 1d (NL-Loo) --------------------------
color_1 <- "#f1a340"
color_2 <- "#998ec3"

data_cv_shift %>%  
    filter(Site == "NL-Loo") %>% 
    ggplot(aes(x = Year, )) +
    geom_point(aes(y = shift_CV), color = "grey", shape = 1) +
    stat_smooth(aes(y = shift_CV), color = "black", linewidth = rel(0.5),
                method = "lm", se = FALSE) +
    stat_poly_eq(
      aes(y = RUE_diff, label = paste(..eq.label.., ..p.value.label.., sep = "~~~~")),
      formula = y ~ x, size = 2, parse = TRUE) +
    ylab("LWII") +
    theme_Yu +
  data_cv_shift %>%  
    filter(Site == "NL-Loo") %>% 
    ggplot(aes(x = Year, )) +
    geom_point(aes(y = SW_mean.rolling.5), color = colorspace::lighten(color_1, amount = 0.5), shape = 1) +
    stat_smooth(aes(y = SW_mean.rolling.5), color = color_1, linewidth = rel(0.5),
                method = "lm", se = FALSE) +
    stat_poly_eq(
      aes(y = SW_mean.rolling.5, label = paste(..eq.label.., ..p.value.label.., sep = "~~~~")),
      formula = y ~ x, size = 2, parse = TRUE) +
    scale_y_continuous(
      name = NULL,
      sec.axis = sec_axis(~ ., name = expression("SW [W" ~ m^-2 * "]"))
    ) +
    theme_Yu +
  data_cv_shift %>%  
    ggplot(aes(x = Year, )) +
    geom_point(aes(y = Precip_mean.rolling.5), color = colorspace::lighten(color_2, amount = 0.5), shape = 1) +
    stat_smooth(aes(y = Precip_mean.rolling.5), color = color_2, linewidth = rel(0.5),
                method = "lm", se = FALSE) +
    stat_poly_eq(
      aes(y = Precip_mean.rolling.5, label = paste(..eq.label.., ..p.value.label.., sep = "~~~~")),
      formula = y ~ x, size = 2, parse = TRUE) +
    scale_y_continuous(
      name = NULL,
      sec.axis = sec_axis(~ ., name = expression("Precipitation [mm]"))
    ) +
    theme_Yu

ggsave("Fig_NL_Loo_LWII_SW_Pre.pdf", width = 21, height = 4, units = "cm")



# Plot Fig. S1 (site map) ------------------------------------------------
library(maps)
library(ggmap)
library(mapdata)

# Load site info data
Site_info <- readxl::read_excel(path = file.path("Site_info/Ameriflux.ICOS.Long.xlsx"))

Site_info$PFT[which(Site_info$PFT %in% c("OSH", "CSH"))] <- "SH"
Site_info$PFT[which(Site_info$PFT %in% c("WSA", "SAV"))] <- "SAV"
Site_info$PFT[which(Site_info$PFT %in% c("ENF", "DNF"))] <- "NF"
Site_info$PFT[which(Site_info$PFT %in% c("CVM"))] <- "CRO"
Site_info$PFT[which(Site_info$PFT %in% c("WAT", "SNO"))] <- "NO"

Site_info$PFT <- factor(Site_info$PFT, levels = c("EBF", "NF", "DBF", "MF", "SAV", "GRA", "CRO", "SH", "WET", "NO"))

# Create a world map using ggplot2
theme_Yu <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

world_map <- ggplot() +
  borders("world",
          xlim = c(-155, 50),   # Limit longitude to North America + Europe
          ylim = c(0, 80),    # Limit latitude range
          colour = "#e8e8e8",
          fill = "lightgrey") +
  geom_point(data = Site_info,
             aes(x = Long, y = Lat, color = PFT, shape = PFT),
             size = 1.5, stroke = 1) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#8c510a", "#33a02c",
                                "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                                "#6a3d9a"),
                     labels = c("EBF", "NF", "DBF", "MF", "SAV", "GRA",
                                "CRO", "SH", "WET")) +
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 20, 21, 22, 23),
                     labels = c("EBF", "NF", "DBF", "MF", "SAV", "GRA",
                                "CRO", "SH", "WET")) +
  coord_fixed(xlim = c(-155, 50), ylim = c(0, 80), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_Yu

# Print the map
print(world_map)
ggsave("Site.map.pdf", width = 10, height = 5, units = "cm")

# Plot Fig. S7 -----------------------------------
# Compare annual LUE with annual mean of daily LUE
dir_data_DD <- list.files("Data_DD", pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)

dir_data_DD <- dir_data_DD[!grepl("/FLUXNET/|/Gross/", dir_data_DD)]

site_id <- Site_info
dir_data_DD <- dir_data_DD[site_id %in% Site_info]

df_data_DD <- map_dfr(dir_data_DD, function(f) {
  site <- str_extract(f, "[A-Za-z]{2}-[A-Za-z0-9]+")
  
  df <- read_csv(f, show_col_types = FALSE) %>%
    select(TIMESTAMP, SW_IN_F, GPP_DT_VUT_REF) %>%
    mutate(Site = site)
  
  return(df)
})

df_data_DD <- df_data_DD %>% 
  filter(SW_IN_F > 0 & GPP_DT_VUT_REF > 0) %>% 
  mutate(LUE_DD = GPP_DT_VUT_REF / SW_IN_F,
         Year = substr(TIMESTAMP, 1, 4) %>% as.integer())

data_DD_summary <- df_data_DD %>% 
  group_by(Site, Year) %>% 
  summarise(LUE_DD = mean(LUE_DD, na.rm = T)) %>% 
  left_join(data_YY_filtered, by = c("Site", "Year")) %>% 
  rename(LUE_YY = LUE) %>% 
  select(Site, Year, LUE_DD, LUE_YY)

r2 <- cor(data_DD_summary$LUE_DD, data_DD_summary$LUE_YY)

data_DD_summary %>% 
  ggplot(aes(LUE_YY, LUE_DD)) +
  geom_point(color = "grey", alpha =0.5, fill = NA, shape = 21) +
  stat_smooth(color = "black", linewidth = rel(0.5),
              method = "lm", se = F) +
  annotate("text", x = 2, y = 0.07, size = 2,
           label = paste0("R² = ", round(r2, 2))) +
  xlab(expression("Annual LUE (g " ~ MJ^-1 * ")")) + 
  ylab(expression("Annual mean of daily LUE (mol" ~ mol^-1 * ")")) +
  theme_Yu

ggsave("Fig_LUE_YY_DD.pdf", width = 6, height = 6, units = "cm")