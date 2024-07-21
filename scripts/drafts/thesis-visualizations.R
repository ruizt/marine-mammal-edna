library(here)
library(readxl)
library(tidyverse)
library(ggplot2)
library(kableExtra)
library(VennDiagram)

## Initial Set Up --------------------------------------------------------------

load(here("rslt", "comp", "cv-eta-2024-04-28.RData"))

# for density estimates
load(here("data", "ceta-density-processed-2024-05-02.RData"))

# load(here("data", "edna-imputed-2024-05-02.RData")) # don't need

# for edna_agg, edna_samples, edna_imputed, etc
load(here("data", "ncog-18s-processed-2024-05-02.RData"))

file_path <- "rslt/tbl/fit-tables-2024-04-28.xlsx"

sheet_names <- excel_sheets(file_path)

asv_bm <- read_excel(file_path, sheet = sheet_names[1])
asv_bp <- read_excel(file_path, sheet = sheet_names[2])
asv_mn <- read_excel(file_path, sheet = sheet_names[3])
fit_summary <- read_excel(file_path, sheet = sheet_names[4])

count(group_by(asv_bm, p), name = 'bm') |>
  full_join(count(group_by(asv_bp, p), name = 'bp'), by = 'p') |>
  full_join(count(group_by(asv_mn, p), name = 'mn'), by = 'p')

# join asvs across species
asv_sel <- full_join(asv_bm, asv_bp, 
                     by = c('short.id', 'd', 'p', 'c', 'o', 'f', 'g'),
                     suffix = c('.bm', '.bp')) |>
  full_join(asv_mn, 
            by = c('short.id', 'd', 'p', 'c', 'o', 'f', 'g')) |>
  rename(coef.mn = coef) 
# taxonomic ranks: domain, phylum, class, order, family, genus

# average weight is 0.0495
asv_bm %>% summary()

asvs_total <- asv_sel |>
  mutate(coef.bm = log2(coef.bm/100 + 1),
         coef.bp = log2(coef.bp/100 + 1),
         coef.mn = log2(coef.mn/100 + 1)) %>%   
  mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
  dplyr::select(short.id, phylum, starts_with('coef')) |>
  rename_with(~str_remove(.x, 'coef.')) |>
  pivot_longer(bm:mn, names_to = 'species', values_to = 'coef') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                          labels = c('blue', 'fin','humpback'))) 


data <- asv_sel |> 
  dplyr::select(starts_with('coef'), everything()) |>
  rename_with(~str_remove(.x, 'coef.')) |> 
  pivot_longer(bm:mn, names_to = 'species', values_to = 'coef') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                          labels = c('Blue', 'Fin','Humpback'))) |>
  drop_na(coef) |>
  mutate(coef = log2(coef/100 + 1)) 

# marine mammal abundance ------------------------------------------------------
# original density estimates
density_estimates |>
  mutate(date = ym(str_remove(cruise, 'X')),
         season = str_to_title(season)) |>
  pivot_longer(bp:mn, names_to = 'species', values_to = 'density') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                          labels = c('Blue', 'Fin','Humpback'))) |>
  ggplot(aes(y = density, x = date, fill = season)) + # or fill = as.factor(qtr)
  geom_bar(stat="identity") +
  facet_wrap(~species, scales = "free_y", ncol = 1) +
  scale_y_continuous(breaks = c(1, 10, 100), 
                     #trans = "log2", 
                     label = c(1, 10, 100)) +
  labs(title = "Scaled Sightings over Time", 
       subtitle = "per 1000 km effort", x = "", y = "") +
  guides(fill = guide_legend(title = "Season")) +
  theme_bw()

# seasonal de-trended data ( log[y_{i}/g_{y}(i)] )
log_density_estimates_adj |>
  left_join(density_estimates %>% select(cruise, season), by = "cruise") %>% 
  mutate(date = ym(str_remove(cruise, 'X')),
         season = str_to_title(season)) |>
  pivot_longer(bp:mn, names_to = 'species', values_to = 'density') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                          labels = c('Blue', 'Fin','Humpback')),
       #  density = exp(density)
         ) |>
  ggplot(aes(y = density, x = date, fill = season)) +
  geom_bar(stat="identity") +
  facet_wrap(~species, scales = "free_y", ncol = 1) +
 # scale_y_continuous(breaks = c(0.5, 5, 10), label = c(0.5, 5, 10)) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4), label = c(-4, -2, 0, 2, 4)) +
  labs(title = "Seasonally De-trended Sightings over Time", 
       subtitle = "per 1000 km effort", x = "", y = "") +  
  guides(fill = guide_legend(title = "Season")) +  
  theme_bw()

## Boxplots
# original density estimates
density_estimates_imputed |>
  mutate(date = ym(str_remove(cruise, 'X')),
         season = str_to_title(season)) |>
  pivot_longer(bp:mn, names_to = 'species', values_to = 'density') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                          labels = c('Blue', 'Fin','Humpback'))) |>
  ggplot(aes(y = density, fill = species, x = season)) +
  scale_y_continuous(trans = "log2") +
  geom_boxplot() +  
  labs(title = "Estimated Density by Season and Whale Species", 
       x = "", y = "") + 
  guides(fill = guide_legend(title = "Species")) +  
  scale_fill_manual(values = c("#4663cb", "#9fa4b9", "#363945")) +
  theme_bw()


# seasonal de-trended data ( log[y_{i}/g_{y}(i)] )
log_density_estimates_adj |>
  left_join(density_estimates %>% select(cruise, season), by = "cruise") %>% 
  mutate(date = ym(str_remove(cruise, 'X')),
         season = str_to_title(season)) |>
  pivot_longer(bp:mn, names_to = 'species', values_to = 'density') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                          labels = c('Blue', 'Fin','Humpback')),
  ) |>
  ggplot(aes(y = density, fill = species, x = season)) +
  geom_boxplot() +  
  labs(title = "Seasonally-Detrended Density by Season and Whale Species", 
       x = "", y = "") + 
  guides(fill = guide_legend(title = "Species")) +  
  scale_fill_manual(values = c("#4663cb", "#9fa4b9", "#363945")) +
  theme_bw()



# Eta selection ----------------------------------------------------------------
# inspect to manually choose eta (0.6 seems reasonable across the board)
cv_out |>
  filter(metric %in% c('df', 'rsq', 'spe')) |>
  group_by(eta, species, metric) |>
  summarize(mean = mean(value),
            sd = sd(value)) |>
  mutate(sd = if_else(metric == 'spe', NA, sd)) |>
  ggplot(aes(x = eta, y = mean)) +
  facet_wrap(~metric + species, scales = 'free_y', nrow = 3) +
  geom_path() +
  geom_ribbon(aes(ymin = mean - sd,
                  ymax = mean + sd),
              alpha = 0.3) + 
  geom_vline(aes(xintercept = eta.approx), 
             data = tibble(species = c('bm', 'bp', 'mn'),
                           eta.approx = c(0.65, 0.55, 0.63)),
             linetype = 'dotdash')

eta_sel_approx <- tibble(species = c('bm', 'bp', 'mn'),
                         eta.approx = c(0.65, 0.55, 0.63))

knitr::kable(eta_sel_approx)

# export to latex (and optionally place on clipboard)
xtable::xtable(eta_sel_approx) %>% 
  clipr::write_clip()

## Venn Diagrams --------------------------------------------------------------

create_venn_diagram <- function(data, column = p) {
  
  data <- data |> 
    dplyr::select(starts_with('coef'), everything()) |>
    rename_with(~str_remove(.x, 'coef.')) |> 
    pivot_longer(bm:mn, names_to = 'species', values_to = 'coef') |>
    mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                            labels = c('Blue', 'Fin','Humpback'))) |>
    select(short.id, !!column, species, coef) |>
    drop_na() 
  
  # Calculate the number of unique phylum in each species
  species_group_count <- tapply(data[column], data$species, function(x) length(unique(x)))
  
  # Create a list of vectors of phylum for each species
  species_group <- lapply(unique(data$species), function(s) unique(data[[column]][data$species == s]))
  
  # Create the Venn diagram to return
  venn_plot <- venn.diagram(
    x = species_group,
    category.names = unique(data$species),
    filename = NULL,
    disable.logging = TRUE,
    cex = 2,
    cat.cex = 2,
    cat.dist = 0.1,
    margin = 0.12
  )
  
  grid.newpage() 
  
  grid.draw(venn_plot)
}  

create_venn_diagram(data = asv_sel, column = "p")

create_venn_diagram(data = asv_sel, column = "c")

create_venn_diagram(data = asv_sel, column = "o")

create_venn_diagram(data = asv_sel, column = "short.id")

# https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf

## ASV coef plots --------------------------------------------------------------

# individual ASVs for one species, colored by phylum
asv_bm %>% 
  arrange(p, desc(coef)) %>% 
  ggplot(aes(y = short.id, x = coef, fill = p)) +
  geom_bar(position = "dodge",
           stat = "identity") + 
  theme_bw() +
  labs(y = "asv") +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank()) 

# for one species, colored by phylum
asv_bm %>% 
  arrange(p, desc(coef)) %>% 
  ggplot(aes(y = p, x = coef)) +
  geom_bar(aes(fill = short.id),
           position = "dodge",
           stat = "identity") + guides(fill="none") +
  theme_bw() +
  labs(y = "asv") +
  theme(legend.position = "none",
      #  axis.text.y = element_blank(),
        axis.ticks = element_blank()) 

# all species, faceted by phylum
asv_sel |>
  mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
  dplyr::select(short.id, phylum, starts_with('coef')) |>
  rename_with(~str_remove(.x, 'coef.')) |>
  pivot_longer(bm:mn, names_to = 'species', values_to = 'coef') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), 
                          labels = c('blue', 'fin','humpback'))) |>
  drop_na() |>
  filter(phylum != "Unknown", phylum != "uncultured") |>
  group_by(phylum) |>
  filter(n_distinct(species) > 1) |>
  ungroup() |>
  ggplot(aes(x = coef, y = short.id)) +
  geom_col(aes(fill = species), position = "dodge") + 
  theme_bw() +
  labs(y = "asv") +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(~phylum, scales = "free_y")

## Trevor's plots --------------------------------------------------------------

# visualize
pal <- colorRampPalette(c('#649dfa', '#02204f'), bias = 1)

# for all species
asvs_total |>
  ggplot(aes(y = phylum)) +
  geom_col(aes(x = coef.trans, fill = phylum),
           width = 1,
           position = position_dodge2(preserve = 'single',
                                      padding = 0.1)) +
  facet_wrap(~species) +
  scale_fill_manual(values = pal(26)) +
  guides(fill = guide_none()) +
  scale_x_continuous(n.breaks = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1, 
                                          color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.1, 
                                          color = 'grey')) +
  labs(x = 'percent change in median density', 
       y = '',
       title = '')

# for one species
asv_mn |>
  mutate(coef = log2(coef/100 + 1)) %>% 
  mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
  mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x))),
         x = 100*(2^coef - 1)) |>
  ggplot(aes(y = phylum)) +
  geom_col(aes(x = x, fill = phylum),
           width = 1,
           position = position_dodge2(preserve = 'single',
                                      padding = 0.1)) +
  scale_fill_manual(values = pal(19)) +
  guides(fill = guide_none()) +
  scale_x_continuous(n.breaks = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1, 
                                          color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.1, 
                                          color = 'grey')) +
  labs(x = 'percent change in density', 
       y = '',
       title = 'humpback')


## Residual vs Fits Graphs -----------------------------------------------------
loopreds_sel <- cv_out |>
  left_join(eta_sel_approx, by = 'species') |>
  group_by(species) |>
  slice_min(abs(eta - eta.approx)) |>
  dplyr::select(-eta.approx) |>
  filter(metric %in% c('p', 'pe', 'spe'))

p1 <- loopreds_sel |>
  filter(metric == 'p') |>
  bind_cols(obs = c(bm, bp, mn)) |>
  ggplot(aes(x = obs, y = value)) +
  geom_point() +
  facet_wrap(~species) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  labs(x = '',
       y = 'prediction')


# p2 <- tibble(fitted = c(fitted_bm, fitted_bp, fitted_mn),
#        observed = c(bm, bp, mn),
#        species = rep(c('blue', 'fin', 'humpback'), 
#                      times = c(length(bm), length(bp), length(mn)))) |>
#   ggplot(aes(x = observed, y = fitted)) +
#   facet_wrap(~species) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_bw(base_size = 18) +
#   labs(x = 'scaled sighting logratio',
#        y = 'fitted value') +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linewidth = 0.1, 
#                                         color = 'black'),
#         panel.background = element_rect(fill = "#fffffeff"),
#         plot.background = element_rect(fill = "#fffffeff"),
#         text = element_text(color = "#49485a"))
# 
# library(patchwork)
# 
# fig <- p1 + p2 + plot_layout(nrow = 2)

