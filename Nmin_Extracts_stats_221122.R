#' ---
#' date: "`r format(Sys.time(), '%d-%m-%Y')`"
#' title: "Stats - Nmin Extract Data"
#' subtitle: "Take 2 - fitting with nested data instead. Update with correct C data."
#' author: Taleta Bailey
#' output:
#'   html_document:
#'     toc: true
#' ---


#' ## Set up
#+ setup, results = 'hide'
# Set up ------------


# Load data (from other scripts in same project)
# source('Nmin_Extracts.R')

packages_stats <- c('tidyverse', 'nlme', 'emmeans', 'multcomp', 'easystats','broom', 'patchwork', 'kableExtra')
lapply(packages_stats, library, character.only = T)

# Standard error function
se2 <- function(x) sd(x, na.rm = T)/sqrt(sum(!is.na(x)))


# Custom version of str() with levels cut off at 2
str2 <- \(x) {
  str(x, max.level = 2)
}


#' ## Data

# Data ---------------

inc_statdat <- read.csv("Nmin_statdat_240325.csv")
str2(inc_statdat)

inc_statdat <- inc_statdat %>% 
  mutate(across(c(where(is.character), Rep, day), factor))



# Plot for initial look
inc_statdat %>% 
  # filter(analyte %in% c('NO3')) %>% 
  ggplot(aes(x = Fraction, y = ugg_change, col = Rep)) + # units converted to ug/g fraction material
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.5) + 
  facet_grid(analyte~Soil, scales = 'free', space = 'free_x') + ggtitle('ug/g-fraction change')

inc_statdat %>% 
  # filter(analyte %in% c('NO3')) %>% 
  ggplot(aes(x = Fraction, y = mggC_change, col = Rep)) + # this one still in mg/gC
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.5) + 
  facet_grid(analyte~Soil, scales = 'free', space = 'free_x') + ggtitle('mg/g-C change')


# Print summary for each variable
summary(inc_statdat)


# '-----------------------'----------------------------------------------------

#' # ug/g fraction data
# ug/g fraction  --------------------

# ***************************
# UPDATES TO MAKE: 
# ***************************

#' Is the change in extract analyte concentration different between soils and fractions?

#' Response: change in NO3, NH4, TN, TOC, FAA, DON over 14 day incubation
#' Fixed effects: Soil + Fraction
#' Possible random effects: Rep, Soil ...

#' ## Steps in data analysis 
#' 1. Nest data
#' 2. Make plots of nested data
#' 3a. Use map to fit anova model 
#' 3b. Check assumptions
#' 4a. Refit with new gls model
#' 4b. Check assumptions
#' 5. Post-hoc 

#' ## 1. Nest data
## 1. Nest data ---------------

# Make nested data.frame from inc_statdat, so now each analyte is separate 'dataframe'
inc_nest_ugg <- inc_statdat %>% 
  group_by(analyte) %>% 
  nest()

# E.g. 
inc_nest_ugg$data[[1]] # Just the NH4 data

#' ## 2. Visualisation - Make plots of nested data
## 2. Visualisation --------------------------

# inc_nest_ugg <- inc_nest_ugg %>% 
#   mutate(plots = 
#           map2(data, analyte, ~ ggplot(., aes(x = Soil, y = ugg, fill = Fraction)) + 
#                      geom_boxplot() + theme(legend.position = 'top') + 
#          ggtitle(.y)
#        ))
# 
# inc_nest_ugg$plots[1:6]

#' Function to make visualisation plots - box plot, histogram and interaction plot
#' When mapping with nested data, the input will the nested dataframe for each analyte. 
# Resp is analyte/response variable, y is response unit, title is
viz <- function(resp, y, title){ 
  bp <- ggplot(resp, aes(x = Soil, y = .data[[y]], fill = Fraction)) + 
    geom_boxplot() + theme(legend.position = 'top')
  
  hg <- ggplot(resp, aes(x=.data[[y]])) + geom_histogram(bins = 20)
  
  interact <- ggplot(resp, aes(x = Soil, y = .data[[y]], col = Fraction, group = Fraction)) + 
    geom_line(stat = 'summary', fun = 'mean')  + theme(legend.position = 'top') + ggtitle('Interaction')
  
  bp + hg + interact + plot_annotation(title = paste('Viz data for', title)) # patchwork to stick plots together
}

# test run
viz(resp = inc_nest_ugg$data[[6]], y = 'ugg_change', title = inc_nest_ugg$analyte[[6]])
  # works 

# Now, apply to each analyte using pmap. pmap takes list of variables. 
inc_nest_ugg <- inc_nest_ugg %>% 
    mutate(viz =
             pmap(list(resp = data, # data 'column' in nested data
                       y = 'ugg_change', # used to id column in data[]
                       title = analyte) # for labelling graphs
                  , viz))
#+ uggviz, fig.width = 10
inc_nest_ugg$viz[1:6] # amazing


#' ## 3a. Fit two-way anovas
#' Fit two-way anova with interaction between Soil and Fraction. 
#' Shouldn't need to pre-define function as can be directly defined within map. 

## 3a. Two-way anova fit -----------------
inc_nest_ugg <- inc_nest_ugg %>% 
  mutate(anova_fit = 
           map(data, ~lm(ugg_change ~ Fraction * Soil, data = .))) # So much neater than lists.. 

# # Anova with error term
# inc_nest_ugg <- inc_nest_ugg %>% 
#   mutate(anova_fit2 = 
#            map(data, ~aov(ugg_change ~ Fraction * Soil + Error(Soil), data = .)))

# Anova without interaction
inc_nest_ugg <- inc_nest_ugg %>% 
  mutate(anova_fit3 = 
           map(data, ~lm(ugg_change ~ Fraction + Soil, data = .)))

# # Anova without interaction, with error
# inc_nest_ugg <- inc_nest_ugg %>% 
#   mutate(anova_fit4 = 
#            map(data, ~aov(ugg_change ~ Fraction + Soil + Error(Soil), data = .)))

# Compare models
# Interaction vs no interaction
map2_df(inc_nest_ugg$anova_fit3, inc_nest_ugg$anova_fit, anova, .id = 'list ref') # For all vars, interaction sig. 

# With error vs
# pmap(list(inc_nest_ugg$anova_fit, inc_nest_ugg$anova_fit2,inc_nest_ugg$anova_fit3,inc_nest_ugg$anova_fit4), anova)


inc_nest_ugg %>% 
  dplyr::select(analyte, anova_fit) %>% 
  mutate(glance_aov = 
           map(anova_fit, glance)) %>% 
  # select(analyte, glance_aov) %>% 
  unnest(glance_aov)

# Re-write this so more similar to ^^ ie. map and unnest so has analyte id
map_df(inc_nest_ugg$anova_fit, anova, .id = 'List ref') %>% tidy


#' ## 3b. Check model assumptions
#' Next function takes models, and produces check_model() plots. This will also work for gls() later. 

## 3b. Check assumptions ----------------------
check_plots <- function(mod, title){
  check_model(mod, check = c("qq", 'normality', 'linearity', 'homogeneity', 'pp_check')) %>%  # Function from easy_stats package
    plot() + plot_annotation(title = paste('Check_model for ', title)) # piping to plot allows adding title with plot_annotation from patchwork
}
# # test run
# check_plots(inc_nest_ugg$anova_fit[[1]], inc_nest_ugg$analyte[[1]]) # works

# Make plots for all anova models and store in same df
inc_nest_ugg <- inc_nest_ugg %>% 
  mutate(anova_check = 
    map2(anova_fit, analyte, ~check_plots(.x, .y)))

# Print plots - now all labeled with which analyte
inc_nest_ugg$anova_check

#' The only good way to label these plots seems to be to individually call them within Markdown 
#' NO3
# inc_nest_ugg[inc_nest_ugg$analyte == 'NO3', ]$anova_check

#' NH4
# inc_nest_ugg[inc_nest_ugg$analyte == 'NH4', ]$anova_check

#' DON2
# inc_nest_ugg[inc_nest_ugg$analyte == 'DON2', ]$anova_check

#' FAAN
# inc_nest_ugg[inc_nest_ugg$analyte == 'FAAN', ]$anova_check

#' TC
# inc_nest_ugg[inc_nest_ugg$analyte == 'TC', ]$anova_check

#' TN
# inc_nest_ugg[inc_nest_ugg$analyte == 'TN', ]$anova_check


#' ## Post-hoc - anovas
## Post-hoc ----------------
inc_nest_ugg %>% 
  mutate(comp = 
           map(anova_fit, estimate_contrasts, contrast = 'Fraction', at = 'Soil'))

# Estimate marginal means and comparisons
inc_nest_ugg <- inc_nest_ugg %>%
  mutate(comp = 
           map2(anova_fit, data, ~emmeans(.x, ~Fraction | Soil, data = .y))) %>% 
  mutate(cld =
           map(comp, multcomp::cld, Letters = letters)) 
  
# Tidy output of comparisons within soil
inc_nest_ugg  %>% 
    dplyr::select(analyte, cld) %>% 
    unnest(cld)

#' NO3
inc_nest_ugg[inc_nest_ugg$analyte == 'NO3', ]$cld

#' NH4
inc_nest_ugg[inc_nest_ugg$analyte == 'NH4', ]$cld

#' DON2
inc_nest_ugg[inc_nest_ugg$analyte == 'DON2', ]$cld

#' FAAN
inc_nest_ugg[inc_nest_ugg$analyte == 'FAAN', ]$cld

#' TC
inc_nest_ugg[inc_nest_ugg$analyte == 'TC', ]$cld

#' TN
inc_nest_ugg[inc_nest_ugg$analyte == 'TN', ]$cld



## 4a. gls() fit -----------------
inc_nest_ugg <- inc_nest_ugg %>% 
  mutate(gls_fit = 
           map(data, ~gls(ugg_change ~ Fraction * Soil, weights = varIdent(form = ~ 1 | Soil), data = .)))

inc_nest_ugg$gls_fit[[1]]

map_df(inc_nest_ugg$gls_fit, anova)


## 4b. Check assumptions....... doesn't work yet -----------------
# Function to plot homogeneity of variances

check_hs <- function(mod, modname){
  plot(mod, abline = 0, main = paste(modname)) 
}

# Plots to check heteroscedasticity
# At least for the residuals vs. fitted. 
inc_nest_ugg %>% 
  mutate(gls_hs =
           map2(gls_fit, analyte, check_hs)
  ) %>% 
  pull(gls_hs)


## Post-hoc -----------------------
inc_nest_ugg <- inc_nest_ugg %>% 
  mutate(gls_comp = 
           map2(gls_fit, data, ~emmeans(.x, ~Fraction | Soil, data = .y))) %>% 
  mutate(gls_cld =
           map(gls_comp, multcomp::cld, Letters = letters))

#' NO3
inc_nest_ugg[inc_nest_ugg$analyte == 'NO3', ]$gls_cld

#' NH4
inc_nest_ugg[inc_nest_ugg$analyte == 'NH4', ]$gls_cld

#' DON2
inc_nest_ugg[inc_nest_ugg$analyte == 'DON2', ]$gls_cld

#' FAAN
inc_nest_ugg[inc_nest_ugg$analyte == 'FAAN', ]$gls_cld

#' TC
inc_nest_ugg[inc_nest_ugg$analyte == 'TC', ]$gls_cld

#' TN
inc_nest_ugg[inc_nest_ugg$analyte == 'TN', ]$gls_cld

# '-----------------------'----------------------------------------------------

#' # mg/g-C data
# mg/g-C  --------------------
#' Is the change in extract analyte concentration per g C different between soils and fractions?

#' Response: change in NO3, NH4, sumN, TOC, FAA, DON over 14 day incubation
#' Fixed effects: Soil + Fraction

#' ## Steps in data analysis 
#' 1. Nest data - done for ug/g fraction 
#' 2. Make plots of nested data
#' 3a. Use map to fit anova model 
#' 3b. Check assumptions
#' 4a. Refit with new gls model
#' 4b. Check assumptions
#' 5. Choose models - anova or gls fit
#' 6. Post-hoc 

#' ## 1. Nest data
## 1. Nest data ---------------

#Already done for ug/g fraction analysis
# Make copy for mggC analysis with just analyte and data columns
inc_nest_mggC <- inc_nest_ugg %>% 
  dplyr::select(analyte, data)

# E.g. 
inc_nest_mggC$data[[1]] # Just the NH4 data, mggC_change now column of interest


#' ## 2. Visualisation - Make plots of nested data
## 2. Visualisation --------------------------

#' Function already made- viz()

# test run
# viz(resp = inc_nest_mggC$data[[6]], y = 'mggC_change', title = inc_nest_ugg$analyte[[6]])
# works 

# Now, apply to each analyte using pmap. pmap takes list of variables. 
inc_nest_mggC <- inc_nest_mggC %>% # Make new variable for mggC analysis
  mutate(viz =
           pmap(list(resp = data, # data 'column' in nested data
                     y = 'mggC_change', # used to id column in data[]
                     title = analyte) # for labelling graphs
                , viz))
#+ mggcviz, fig.width = 10
inc_nest_mggC$viz[1:6] # amazing


#' - NH4 data has evident negative skew, `soil*fraction` interaction also likely.
#' - NO3 appears to have positive skew, `soil*fraction` interaction highly likely
#' - TN removed, replaced with sumN. TN appears relatively close to normal, interaction likely
#' - TC slight positive skew, interaction likely 
#' - FAAN possibly negative skew, interaction likely
#' - DON2 negative skew, interaction likely.
#' - sumN (sum of NH4, NO3, DON2 and FAAN) has slight negative skew. 
#' Unequal variances likely across Fraction, as well as soil? Should maybe test both.
#' `soil*fraction` interaction likely.   


#' ## 3a. Fit two-way anovas
#' Fit two-way anova with interaction between Soil and Fraction. 
#' Shouldn't need to pre-define function as can be directly defined within map. 

## 3a. Two-way anova fit -----------------
# With and without interaction

inc_nest_mggC <- inc_nest_mggC %>% 
  mutate(anova_fit = 
           map(data, ~lm(mggC_change ~ Fraction * Soil, data = .))) %>%  # So much neater than lists.. 
  mutate(anova_fit3 = 
           map(data, ~lm(mggC_change ~ Fraction + Soil, data = .))) %>%  # Without interaction
  mutate(anova_fitpos = 
           map(data, ~lm(mggC_change + (abs(min(mggC_change)) + 1) ~ Fraction * Soil, data = .)))


#+ warnings = FALSE
# Compare models
# Interaction vs no interaction
inc_nest_mggC %>% 
  dplyr::select(analyte, anova_fit, anova_fitpos) %>% 
  mutate(compare = map(map2(anova_fit, anova_fitpos, anova), tidy)) %>% 
  unnest(compare)

#' All show significant effect of interaction term

# /*
map2(inc_nest_mggC$anova_fit3, inc_nest_mggC$anova_fit, anova) # For all vars, interaction sig. 
# */

# Look at model quality metrics
inc_nest_mggC %>% 
  dplyr::select(analyte, anova_fit) %>% 
  mutate(glance_aov = 
           map(anova_fit, glance)) %>% 
  # select(analyte, glance_aov) %>% 
  unnest(glance_aov)

#/*
# Anova of main effects and interactions
inc_nest_mggC %>% 
  dplyr::select(analyte, anova_fit) %>% 
  mutate(aov_tidy = 
           map(anova_fit, anova)) %>%
  mutate(aov_tidy = map(aov_tidy, tidy)) %>% 
  unnest(aov_tidy)
#*/

#' ## 3b. Check model assumptions
#' Next function takes models, and produces check_model() plots. This will also work for gls() later. 

## 3b. Check assumptions ----------------------

# # test run
# check_plots(inc_nest_mggC$anova_fit[[1]], inc_nest_mggC$analyte[[1]]) # works

#+ mggC-check-plots, include = FALSE
# Make plots for all anova models and store in same df
inc_nest_mggC<- inc_nest_mggC %>% 
  mutate(anova_check = 
           map2(anova_fit, analyte, check_plots))

#+
inc_nest_mggC$anova_check

# Run check heteroscedasticity 
inc_nest_mggC %>% 
  mutate(checkhs = map_dbl(anova_fit, ~check_heteroscedasticity(.x))) %>% 
  dplyr::select(analyte, checkhs) %>% 
  dplyr::filter(checkhs <.05) # find which has significant heteroscedasticity detected


#' - NH4 and TC appear to have significant heteroscedasticity. Could try transformation. 
#' - DON normality of residuals has skewed tails.. could also try transforming
#' - sumN doesn't look too bad


#' Going to run gls() fit for just NH4, TC and DON to avoid transforming data. 
#' Could also try Update NH4, TC and DON models with transform and compare - preferring not to since transforming can bias data, also have negatives so is difficult..

# /*not yet, need to finish model selection
#' ## Post-hoc - anovas
#' ## Post-hoc ----------------
#' inc_nest_mggC %>%
#'   mutate(comp =
#'            map(
#'              anova_fit,
#'              estimate_contrasts,
#'              contrast = 'Fraction',
#'              at = 'Soil'
#'            ))
#' 
#' # Estimate marginal means and comparisons
#' inc_nest_mggC <- inc_nest_mggC %>%
#'   mutate(comp =
#'            map2(anova_fit, data, ~ emmeans(.x, ~ Fraction |
#'                                              Soil, data = .y))) %>%
#'   mutate(cld =
#'            map(comp, multcomp::cld, Letters = letters))
#' 
#' # Tidy output of comparisons within soil
#' inc_nest_mggC  %>%
#'   dplyr::select(analyte, cld) %>%
#'   unnest(cld)
#' 
#' #' NO3
#' inc_nest_mggC[inc_nest_mggC$analyte == 'NO3',]$cld
#' 
#' #' NH4
#' inc_nest_mggC[inc_nest_mggC$analyte == 'NH4',]$cld
#' 
#' #' DON2
#' inc_nest_mggC[inc_nest_mggC$analyte == 'DON2',]$cld
#' 
#' #' FAAN
#' inc_nest_mggC[inc_nest_mggC$analyte == 'FAAN',]$cld
#' 
#' #' TC
#' inc_nest_mggC[inc_nest_mggC$analyte == 'TC',]$cld
#' 
#' #' TN
#' inc_nest_mggC[inc_nest_mggC$analyte == 'TN',]$cld
# */

#' ## 4a. gls()  fit
#' Attempting fitting gls() to data so can include weights argument which should account for heteroscedasticity. 
#' For proper comparison, need to also fit gls without weights to determine if improved fit. 


## 4a. gls() fit -----------------
inc_nest_mggC <- inc_nest_mggC %>% 
  mutate(gls_fit = 
           map(data, ~gls(mggC_change ~ Fraction * Soil, weights = varIdent(form = ~ 1 | Soil), data = .))) %>% 
  mutate(gls_fit_now = 
           map(data, ~gls(mggC_change ~ Fraction * Soil, data = .)))  # gls() without weights

# inc_nest_mggC[inc_nest_mggC$analyte == 'NH4',]$gls_fit


# Compare performance between model with and without weights
inc_nest_mggC %>% 
  mutate(gls_comp = map2(gls_fit, gls_fit_now, anova)) %>% 
  unnest(gls_comp) %>% 
  dplyr::select(analyte, call:`p-value`) %>% 
  subset(.$`p-value` <.05)

#' - Significantly improved fit with `weights` assigned to soil for NH4, TC, and DON
#' - Can stick with anova for NO3, FAA and sumN.   
#'  

#' ## 4b. Check assumptions
#' 
## 4b. Check assumptions.......  -----------------
# Function to plot homogeneity of variances and normality of residuals
# Maybe one day I'll change this to ggplots.. 


check_gls <- function(mod, modname){
  hs <- plot(mod, abline = 0, main = paste(modname))
  qq <- qqnorm(mod, abline = c(0,1), main = paste(modname))
  return(list(hs, qq))
  
}

# Test run
# check_gls(inc_nest_mggC$gls_fit[[1]], inc_nest_mggC$analyte[[1]])


# Plots to check heteroscedasticity and normality of residuals
# At least for the residuals vs. fitted. 
inc_nest_mggC <- inc_nest_mggC %>% 
  mutate(gls_check =
           map2(gls_fit, analyte, check_gls)
  )

# Print check plots - compare anova and gls fits
  # Ideally the two plots would be beside each other, but I can't figure out how.. 
inc_nest_mggC[inc_nest_mggC$analyte %in% c('NH4', 'TC', 'DON2'),]$gls_check
inc_nest_mggC[inc_nest_mggC$analyte %in% c('NH4', 'TC', 'DON2'),]$anova_check


#' DON plots  not better with gls
#'  

#' ## 5. Collate selected models and results
#' Choosing Anova fit for NO3, sumN and FAAN, gls fit for NH4, TC and DON2.  
#' May be easiest to work with gls and anova in separate objects, then try to merge together into tidy output.

## 5. Collect selected models -------------------------
# Anova will be easiest to generate tidy output
inc_nest_mggC.aov <- inc_nest_mggC %>% 
  dplyr::select(-contains('gls')) %>% 
  dplyr::filter(analyte %in% c('NO3', 'sumN', 'FAAN')) %>% 
  mutate(mod_aov = map(anova_fit, anova)) %>% 
  mutate(mod_tidy = map(mod_aov, tidy)) %>% 
  dplyr::select(analyte, mod_tidy) %>% 
  unnest(mod_tidy)

# anova.gls has different output that can't be as easily tidied using broom
# Thinking may be able to just row bind with anova, since is kind of in dataframe?

# Perform anova and extract output to tidy df for gls
inc_nest_mggC.gls <- inc_nest_mggC %>% 
  dplyr::select(-contains('anova')) %>% 
  dplyr::filter(!(analyte %in% c('NO3', 'sumN', 'FAAN'))) %>% 
  mutate(mod_aov = map(gls_fit, anova)) %>% 
  mutate(mod_aov = map(mod_aov, rownames_to_column, var = 'term')) %>% 
  mutate(mod_aov = map(mod_aov, ~rename(.x, c(df = numDF, statistic = `F-value`, p.value = `p-value`)))) %>% 
  dplyr::select(analyte, mod_aov) %>% 
  unnest(mod_aov)


# Bind tidy dfs together with id for which model was used for which variables
inc_nest_mggC.aovres <- bind_rows('anova' = inc_nest_mggC.aov, 'gls' = inc_nest_mggC.gls, .id = 'model')

#' ### ANOVA results - mg/gC main and interaction effects 
print(inc_nest_mggC.aovres, n = 100) 


#' ## 6. Post-hoc - comparisons of emmeans

## 6. Post-hoc -----------------------
# Select models, put into one list-col
inc_nest_mggC <- inc_nest_mggC %>% 
  mutate(mods = ifelse(analyte %in% c('NH4', 'TC', 'DON2'), gls_fit, anova_fit))

# Make emmeans object for each chosen model. Will use later to make contrasts and clds
# Also add contrasts for Fraction in Soil (Frac_comp) and Soil in Fraction (Soil_comp_)
inc_nest_mggC <- inc_nest_mggC %>% 
  mutate(emms = 
           map2(mods, data, ~emmeans(.x, spec = c('Soil', 'Fraction'), data = .y))) 
  


#' ### Comparisons - Fraction in Soil
### Fraction in Soil ------------------
# Compare emmeans for fraction in soil - make contrasts (Frac_comp) and build clds (Frac_cld)
inc_nest_mggC <- inc_nest_mggC %>% 
  mutate(Frac_comp = map(emms, ~contrast(.x, by = 'Soil', method = 'pairwise'))) %>% 
  mutate(Frac_cld =
           map(emms, multcomp::cld, by = 'Soil', Letters = letters))

# Print contrasts output - *Contrasts only within each level of soil for each analyte!!*

comps_mggC_Frac <- inc_nest_mggC %>% 
  dplyr::select(analyte, Frac_comp) %>% 
  mutate(Frac_comp = map(Frac_comp, tidy)) %>% 
  unnest(Frac_comp) %>% 
  dplyr::select(-null.value)

comps_mggC_Frac %>% 
  kbl(caption = "Fraction comparisons - mg/gC", digits = 3) %>% 
  kable_styling(c('hover', 'condensed'), full_width = F) %>% 
  scroll_box(height = "400px")

# Unnest so all emmeans and cld in one df - NOTE - letters are only comparisons for Fraction, within soils for each analyte
cld_mggC_Frac <- inc_nest_mggC %>% 
  dplyr::select(analyte, Frac_cld) %>% 
  unnest(Frac_cld)

#' *NOTE* - letters are only comparisons for `Fraction`, within `Soil`s for each analyte
# print(comps_mggC_Frac, n = 100)
cld_mggC_Frac %>% 
  kbl(caption = "Fraction cld - mg/gC", digits = 3) %>% 
  kable_styling(c('hover', 'condensed'), full_width = F)%>% 
  scroll_box(height = "500px")


#' ### Comparisons - Soil in Fraction
### Soil in Fraction-------------------
# Compare emmeans for soil in fraction - because of interaction, cannot look at soils across fractions
inc_nest_mggC <- inc_nest_mggC %>% 
  mutate(Soil_comp = map(emms, ~contrast(.x, by = 'Fraction', method = 'pairwise'))) %>% 
  mutate(Soil_cld =
           map(emms, multcomp::cld, by = 'Fraction', Letters = letters))

# Print contrasts output - * Contrasts only within each level of Fraction for each analyte!*
comps_mggC_Soil <- inc_nest_mggC %>% 
  dplyr::select(analyte, Soil_comp) %>% 
  mutate(Soil_comp = map(Soil_comp, tidy)) %>% 
  unnest(Soil_comp)

comps_mggC_Soil %>% 
  kbl(caption = "Soil Contrasts - mg/gC", digits = 3) %>% 
  kable_styling(c('hover', 'condensed'), full_width = F)%>% 
  scroll_box(height = "500px")

# Unnest so all emmeans and cld in one df - NOTE - letters are only comparisons for soils, within fractions for each analyte
cld_mggC_Soil <- inc_nest_mggC %>% 
  dplyr::select(analyte, Soil_cld) %>% 
  unnest(Soil_cld)

#' *NOTE* - letters are only comparisons for `Soil`, within `Fraction` for each analyte
cld_mggC_Soil %>% 
  kbl(caption = "Soil cld - mg/gC", digits = 3) %>% 
  kable_styling(c('hover', 'condensed'), full_width = F)%>% 
  scroll_box(height = "400px")

# print(comps_mggC_Soil, n = 100)



# '---------------------'-------------------------------------------
# /*
#' #+ fixing-gls, eval = FALSE
#' 
#' # gls() troubles ---------------------
#' 
#' # Check model worked for individual gls() model, but not when used within map()
#' tmpdat <- inc_statdat %>%
#'   filter(analyte == 'NO3')
#' 
#' tmp <-  gls(ugg_change ~ Fraction * Soil, weights = varIdent(form = ~ 1 | Soil), data = tmpdat)
#' 
#' tmp %>%
#'   check_heteroskedasticity() %>% plot()
#' 
#' check_model(tmp) %>% plot() + plot_annotation(title = 'NO3') # This works
#' check_model(inc_nest_ugg$gls_fit[[2]]) #but this doesn't, error 'check_model() not implemented for gls()'
#' 
#' inc_nest_ugg$gls_fit[[2]]$call$data <- substitute(inc_nest_ugg$data[[2]]) # replacing data in gls fit fixes above command..
#' 
#' check_heteroscedasticity(inc_nest_ugg$gls_fit[[2]]) %>% plot()
#' plot(tmp)
#' 
#' qqnorm(tmp, ~resid(., type = 'n'))
#' qqnorm(inc_nest_ugg$gls_fit[[2]], abline = c(0,1))
#' 
#' # # So, for this to work with gls(), need to replace the data in the model object, as is . due to map
#' # check_plots <- function(mod, dat){
#' #   mod <- mod$call$data <- substitute(dat$data)
#' #   return(mod)
#' #   # check_model(mod, check = c("qq", 'normality', 'linearity', 'homogeneity', 'pp_check')) # Function from easy_stats package
#' # }
#' 
#' inc_nest_ugg %>%
#'   mutate(update_gls =
#'            map2(gls_fit, data, check_plots)) %>%
#'   pull(update_gls)
#' 
#' # inc_nest_ugg$gls_fit
#' 
#' #' ## 4a. Fit gls()
#' #' Fit model with nlme::gls().
#' #' Same structure as 2-way anova but with unequal variance structure to account for heteroscedasticity.
#' 
#' ## 4a. gls() fit -----------------
#' # gls() has a hard time keeping the data, so write wrapper to do so
#' gls_ugg <- function(dat){
#'   mod <- gls(ugg_change ~ Fraction * Soil, weights = varIdent(form = ~ 1 | Soil), data = dat)
#'   mod$call$data <- substitute(dat)
#'   mod
#' }
#' 
#' # Apply to nested data
#' inc_temp <- inc_nest_ugg %>%
#'   mutate(gls_fit =
#'            map(data, gls_ugg))
#' 
#' inc_temp$gls_fit[[1]]
#' check_model(inc_temp$gls_fit[[1]], data = inc_temp$data[[1]])
#' 
#' check_collinearity(inc_temp$gls_fit[[1]], data = inc_temp$data[[1]])
#' 
#' qqnorm(inc_temp$gls_fit[[1]], ~resid(., type = 'n') | Soil)
#' 
#' 
#' inc_nest_ugg <- inc_nest_ugg %>%
#'   mutate(gls_fit =
#'            map(data, ~gls(ugg_change ~ Fraction * Soil, weights = varIdent(form = ~ 1 | Soil), data = .)))
#' 
#' # Compare gls to lm  # I think actually, should compre gls() with and without weights
#' # Using easystats::compare_performance.
#' map2_df(inc_nest_ugg$anova_fit, inc_nest_ugg$gls_fit, compare_performance, .id = 'List ref')  # For each, gls gives lower AIC indicating better fit
#' 
#' ## 4b. Check assumptions ----------------
#' ## NEXT: Make function to generate plots for gls() because check_model has stopped working...
#' 
#' # Function to plot homogeneity of variances
#' 
#' check_hs <- function(mod, modname){
#'   plot(mod, abline = 0, main = paste(modname))
#' }
#' 
#' check_hs(inc_nest_ugg$gls_fit[[2]], modname = inc_nest_ugg$analyte[[2]])
#' 
#' qqnorm(inc_nest_ugg$gls_fit[[1]], abline = c(0,1))
#' 
#' # Make plots for all gls models and store in same df
#' check_model(inc_nest_ugg$gls_fit[[1]], data = inc_nest_ugg$data[[1]])
#' 
#' plot(inc_nest_ugg$gls_fit[[1]], ugg_change ~ fitted(.))
#' 
#' inc_nest<- inc_nest_ugg %>%
#'   mutate(gls_check =
#'            map(gls_fit, check_model))
#' 
#' map(inc_nest_ugg$gls_fit, check_model)
#' lapply(inc_nest_ugg$gls_fit, check_model, check = c("qq", 'normality', 'linearity', 'pp_check'))
#' 
#' ## This works for heteroscedasticity ------------
#' # At least for the residuals vs. fitted.
#' inc_nest_ugg %>%
#'   mutate(gls_hs =
#'     map2(gls_fit, analyte, check_hs)
#'   ) %>%
#'   pull(gls_hs)

# */


