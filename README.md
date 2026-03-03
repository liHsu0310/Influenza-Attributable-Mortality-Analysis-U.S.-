---
title: "Technical Report: Influenza-Attributable Mortality Analysis (U.S.)"
author: "Li Hsu"
date: "2026-03-03"
output:
  html_document: default
---
```{r setup}
###Set up

# 1. Install packages 
packages <- c(
  "readr",
  "here",
  "remotes",
  "dplyr",
  "ggplot2",
  "stringr",
  "patchwork",
  "rmarkdown",
  "knitr"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 2. Load Libraries
library(remotes)
library(readr)
library(here)
library(dplyr)  
library(ggplot2) 
library(stringr)
library(patchwork)
library(rmarkdown)
library(knitr)


# 3. Install wonderapi from GitHub
if (!require("wonderapi")) {
  remotes::install_github("socdataR/wonderapi", force = TRUE)
}

library(wonderapi)

```

```{r data-acquisition}
###Data-acquisition

# 1. Define a helper function to download and query
get_cdc_data <- function(gist_url, filename, db_code) {
  # Ensure the directory exists and define path
  dest <- here::here("Data", filename)
  # Download the XML from GitHub Gist
  curl::curl_download(
    url = gist_url,
    destfile = dest,
    quiet = TRUE,
    handle = curl::new_handle()
  )
  message(paste("XML saved to:", dest))
  # Query CDC WONDER using the specific database code (D77 or D157)
  # Using send_query as per your workflow
  data <- send_query(db_code, dest)
  return(data)
}

# 2. Get Data
# Yearly, 1999-2020
yearly_url_1999 <- "https://gist.githubusercontent.com/liHsu0310/1850fe49102aca41fb7ad92f147a61c9/raw/9ca7774894b99f728e1a3a9d33e434cd1a41a3f7/flu_yearly_query.xml"
Flu_Yr_Raw_1999 <- get_cdc_data(yearly_url_1999, "flu_yearly_query_1999.xml", "D77")

# Respect the 15-second CDC limit
message("Waiting 15 seconds for CDC rate limit...")
Sys.sleep(15)

# Monthly, 1999-2020 
monthly_url_1999 <- "https://gist.githubusercontent.com/liHsu0310/27d2f03d629887a52c4ef9d8babfa8b7/raw/843e351081e0e90bba54beb5eeb42c0df756b9f2/flu_monthly_query_1999-2020.xml"
Flu_Mon_Raw_1999 <- get_cdc_data(monthly_url_1999, "flu_monthly_query_1999.xml", "D77")

# Respect the 15-second CDC limit
message("Waiting 15 seconds for CDC rate limit...")
Sys.sleep(15)

# Yearly, 2021-2024
yearly_url_2021 <- "https://gist.githubusercontent.com/liHsu0310/efffbdfce6ff4b57fb353c214a7bb5d6/raw/ea72ce0797e00b4aae0d725c2c8d4fbc509ab4a5/flu_yearly_query_2021-2024.xml"
Flu_Yr_Raw_2021 <- get_cdc_data(yearly_url_2021, "flu_yearly_query_2021.xml", "D157")
# Respect the 15-second CDC limit
message("Waiting 15 seconds for CDC rate limit...")
Sys.sleep(15)

# Monthly, 2021-2024
monthly_url_2021 <- "https://gist.githubusercontent.com/liHsu0310/c8832b7f88863731f7b1d3fd1215e674/raw/28c99f6c2d9a9d195449361eefe59609c87ffef2/flu_monthly_query_2021-2024.xml"
Flu_Mon_Raw_2021 <- get_cdc_data(monthly_url_2021, "flu_monthly_query_2021.xml", "D157")

# 3. Combine data
# Standardize column types and combine
Flu_Yr_Raw <- bind_rows(
  Flu_Yr_Raw_1999 %>% mutate(across(everything(), as.character)),
  Flu_Yr_Raw_2021 %>% mutate(across(everything(), as.character))
) 
# Sort chronologically
Flu_Yr_Raw <- Flu_Yr_Raw %>% arrange(Year, `Ten-Year Age Groups`)

# Standardize column types and combine
Flu_Mon_Raw <- bind_rows(
  Flu_Mon_Raw_1999 %>% mutate(across(everything(), as.character)),
  Flu_Mon_Raw_2021 %>% mutate(across(everything(), as.character))
) 

# Sort chronologically
Flu_Mon_Raw <- Flu_Mon_Raw %>% arrange(Year, `Ten-Year Age Groups`)

```


```{r data-cleaning}

###Data-cleaning 


# 1. Clean Yearly Data (to get the Population 'Master List')
pop_lookup <- Flu_Yr_Raw %>%
  filter(`Ten-Year Age Groups` != "Not Stated") %>%
  filter(`Year` != "Not Stated") %>%
  # Ensure Population is numeric
  # parse_number() automatically handles commas and non-numeric characters
  mutate(Population = readr::parse_number(as.character(Population))) %>%
  select(Year, `Ten-Year Age Groups`, Population) %>%
  distinct()

# 2. Process Monthly Data
cleaned_data <- Flu_Mon_Raw %>%
  # Drop "Not Stated"
  filter(`Ten-Year Age Groups` != "Not Stated") %>%
  # Handle Suppression: Convert "Suppressed" string to NA
  mutate(Deaths = readr::parse_number(Deaths, na = c("Suppressed", "Missing", "")))
  
# 3. Extract Month and Year from the "Jan., 1999" string
cleaned_data <- cleaned_data %>%  
    mutate(
    # Pull the first 3 letters (e.g., "Jan")
    Month_Abbr = str_extract(Month, "[A-Za-z]{3}"),
    # Pull the 4-digit year
    Year_Extracted = as.numeric(str_extract(Month, "\\d{4}")),
    # Map abbreviation to number (Jan = 1, Jul = 7, etc.)
    Month_Num = match(Month_Abbr, month.abb)
  ) 

# 4. Join with Population from the yearly file
cleaned_data <- cleaned_data %>%  
  left_join(pop_lookup, by = c("Year", "Ten-Year Age Groups"))
  
# 5. Define the Flu Season (July 1 to June 30)
cleaned_data <- cleaned_data %>%  
  mutate(
    Season = ifelse(Month_Num >= 7, 
                    paste0(Year_Extracted, "-", Year_Extracted + 1), 
                    paste0(Year_Extracted - 1, "-", Year_Extracted))
  ) 
    
# 6. Aggregate by Season and Age 
age_group_data <- cleaned_data %>%
  group_by(Season, `Ten-Year Age Groups`) %>%
  summarise(
    Total_Deaths = sum(Deaths, na.rm = TRUE),
    All_Suppressed = all(is.na(Deaths)),
    Avg_Population = mean(Population.y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Handle the NA requirement
    Total_Deaths = ifelse(All_Suppressed, NA, Total_Deaths),
    # Calculate Rate per 100k
    Mortality_Rate = (Total_Deaths / Avg_Population) * 100000
  ) %>%
  # Filter out incomplete boundary seasons
  filter(!Season %in% c("1998-1999", "2024-2025")) %>%
  select(-All_Suppressed)

# 7. Aggregate by Season ONLY (To get the "Grand Total")
seasonal_totals <- age_group_data %>%
  group_by(Season) %>%
  summarise(
    Total_Deaths = sum(Total_Deaths, na.rm = TRUE),
    Avg_Population = sum(Avg_Population, na.rm = TRUE),
    `Ten-Year Age Groups` = "All Ages", # Label for the total row
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate Rate per 100k
    Mortality_Rate = (Total_Deaths / Avg_Population) * 100000
  ) %>%
  # Filter out incomplete boundary seasons
  filter(!Season %in% c("1998-1999", "2024-2025")) 

```

```{r visualization-comparison, fig.width=12, fig.height=10}

###Visualization-comparison 


# 1. Create the plot (all age)
seasonal_totals <- seasonal_totals %>%
  mutate(Season = factor(Season)) # Ensure chronological order

total_plot <- ggplot(seasonal_totals, aes(x = Season, y = Mortality_Rate, group = 1)) +
  # Add a single bold line and points for the peaks
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(color = "black", size = 2) + 
  
  # Styling the X-axis as requested (small and tilted)
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "gray80")
  ) +
  
  # Labels
  labs(
    title = "National Influenza Mortality Trend (All Ages)",
    subtitle = "Aggregated mortality rate per 100,000 population by flu season",
    x = "Flu Season (July - June)",
    y = "Mortality Rate per 100,000"
  )

total_plot


# 2. By age --- Plot A: All groups plotted on a shared scale to demonstrate comparative magnitude ---
# Define the correct chronological age order
age_order <- c("< 1 year", "1-4 years", "5-14 years", "15-24 years", 
               "25-34 years", "35-44 years", "45-54 years", 
               "55-64 years", "65-74 years", "75-84 years", "85+ years")

# Apply this order to your data
age_group_data <- age_group_data %>%
  # Filter out 'All Ages' if you only want specific groups in the facet
  filter(`Ten-Year Age Groups` != "All Ages") %>%
  mutate(`Ten-Year Age Groups` = factor(`Ten-Year Age Groups`, levels = age_order))

combined_plot <- ggplot(age_group_data, aes(x = Season, 
                                       y = Mortality_Rate, 
                                       group = `Ten-Year Age Groups`, 
                                       color = `Ten-Year Age Groups`)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) + # Points help identify the seasonal peaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Influenza Mortality Rates by Age Group (1999-2024)",
    subtitle = "All groups plotted on a shared scale to demonstrate comparative magnitude",
    x = "Flu Season (July - June)",
    y = "Mortality Rate per 100,000",
    color = "Age Group"
  ) +
 scale_color_brewer(palette = "PuOr")

combined_plot


# 3. --- Plot B: Fixed Scales ---

# Create a vector of seasons to display as labels
# This selects every 3rd season so the axis isn't crowded
all_seasons <- levels(factor(age_group_data$Season))
keep_breaks <- all_seasons[seq(1, length(all_seasons), by = 3)]

plot_fixed <- ggplot(age_group_data, aes(x = Season, y = Mortality_Rate, group = `Ten-Year Age Groups`)) +
  geom_line(color = "firebrick", linewidth = 0.8) +
  # Declutter the X-axis here
  scale_x_discrete(breaks = keep_breaks) +
  facet_wrap(~ `Ten-Year Age Groups`, scales = "free_x") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
    axis.text.y = element_text(size = 8),
    axis.line = element_line(color = "gray80"),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "A: Fixed Y-Axis (Comparative Magnitude)",
    subtitle = "Showing every 3rd season on x-axis",
    y = "Mortality per 100k", x = NULL
  )

# 4. --- Plot C: Free Scales ---
plot_free <- ggplot(age_group_data, aes(x = Season, y = Mortality_Rate, group = `Ten-Year Age Groups`)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  # Declutter the X-axis here too
  scale_x_discrete(breaks = keep_breaks) +
  facet_wrap(~ `Ten-Year Age Groups`, scales = "free") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 8),
    axis.line = element_line(color = "gray80"),
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "B: Free Y-Axis (Internal Patterns)",
    subtitle = "Showing every 3rd season on x-axis",
    y = "Mortality per 100k", x = "Flu Season"
  )

# Combine them
comparison_plot <- plot_fixed / plot_free
comparison_plot

# 5. Save the plot
# Save the National Trend Plot
ggsave(
  filename = here::here("Graphs", "total_plot.png"),
  plot = total_plot,
  width = 10, height = 6, dpi = 300
)

# Save the Combined Faceted Plot 
ggsave(
  filename = here::here("Graphs", "combined_plot.png"),
  plot = combined_plot,
  width = 12, height = 8, dpi = 300
)

# Save the Comparison Plot
ggsave(
  filename = here::here("Graphs", "comparison_plot.png"),
  plot = comparison_plot,
  width = 10, height = 15, dpi = 300
)



```



