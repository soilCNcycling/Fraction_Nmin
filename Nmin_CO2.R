# 19/07/2022
# N mineralisation experiment - CO2 data

# Standard error function
se2 <- function(x) sd(x, na.rm = T)/sqrt(sum(!is.na(x)))

# Function to convert character cols to factors
chtofc <- function(df){
  df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], 
                                         as.factor)
  str(df)
  return(df)
}

# Processing CO2 - from recorded readings

# Set up ------------

# ***********************************************************************************
# * <--- **** MAKE SURE PLYR, AOMISC, MASS PACKAGES ARE NOT LOADED!!! ******* ----> *
# ***********************************************************************************

# unload all 'other' packages. They cause problems...
lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(
    paste0('package:', pkgs),
    character.only = T,
    unload = T,
    force = T
  ))

# Load tidyverse. Should be only package required for first bit
library('tidyverse')


# set default ggplot theme
old <- theme_set(theme_bw(base_size = 12) + theme(panel.grid = element_blank()))

# object for ggplot scales
graphcol <- c("grey50", rev(c(colorspace::sequential_hcl(5, palette = 'Purple-Yellow')[1:3], "#A1DAB4")))
graphcol2 <- c(soilpalettes::soil_palette('redox', 5)[c(1,2,4,5)], 'Grey50')
ggscalecols <- list(scale_color_manual(values = graphcol2),
                    scale_fill_manual(values = alpha(graphcol2, 0.5))
)


# Load data -----------
# CO2 dat
CO2 <- read_csv('CO2_rawdata.csv')
str(CO2)

# Change 'Ctrl' to 'Control' in fraction column
CO2[CO2$Fraction %in% 'Ctrl', 'Fraction'] <- 'Control'

# Change 'Fraction' and 'Soil' to factors
CO2$Soil <- factor(CO2$Soil, levels = c('Ctrl', 'CC', 'CV', 'CW', 'NV'))
# Change levels in Fraction
CO2$Fraction <- factor(CO2$Fraction, levels = paste(unique(CO2$Fraction))[c(5, 1:4)])


# Inc setup dat
# Read in incubation set up data - C per incubation tube
soildat <- read_csv('Nmin_Inc setup data_230116.csv')
str(soildat)

soildat[sapply(soildat, is.character)] <- lapply(soildat[sapply(soildat, is.character)], 
                                                     as.factor)

# Change levels in Soil column
soildat$Soil <- factor(soildat$Soil, levels = c('Preinc', 'Ctrl', 'CC', 'CV', 'CW', 'NV'))

# Change levels of Fraction column
soildat$Fraction <- factor(soildat$Fraction, levels = c("Control", "Coarse", "Fine","Whole soil", "Frac. Whole soil"))


# Calibration to standards --------------

# read in standards
stds <- read_csv('CO2_Stands.csv')
str(stds)
# Reshape so dates in columns, then calculate average for each ml CO2
avg_stds <- stds %>% 
  pivot_longer(cols = !c(1:2), names_to = 'date', values_to = 'reading') %>%
  group_by(mL_CO2) %>%
  dplyr::summarise(avgread = mean(reading))

# Plot average - standard curve
ggplot(avg_stds, aes(x = mL_CO2, y = avgread)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F)

# Remove highest point (10ml) - loses linearity, don't need for samples. 
avg_stds <- subset(avg_stds, mL_CO2 != 10)

# Re-plot
ggplot(avg_stds, aes(x = mL_CO2, y = avgread)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F) + 
  coord_cartesian(xlim = c(0, 8), ylim = c(0,5))
  # Looks like better fit

# Use average standards to calculate standard curve regression
calibration <- lm(avgread ~ mL_CO2, data = avg_stds) 
summary(calibration)

# Extract linear function values from calibration model
int <- summary(calibration)$coefficients[1,1]
slope <- summary(calibration)$coefficients[2,1]

# Calculate background CO2. y = 0, y = mx+c, c/m
bg_mlCO2 <- int/slope

# Add background ml CO2 to standards - calculate 'true' ml CO2
avg_stds$adj_mlCO2 <- avg_stds$mL_CO2 + bg_mlCO2

# Re-plot
ggplot(avg_stds, aes(x = adj_mlCO2, y = avgread)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F) + 
  coord_cartesian(xlim = c(0, 8), ylim = c(0,5))

# Repeat calibration regression
adj_calibration <- lm(avgread ~ adj_mlCO2, data = avg_stds) 
summary(adj_calibration)
adj_int <- summary(adj_calibration)$coefficients[1,1]
adj_slope <- summary(adj_calibration)$coefficients[2,1]
  # Slope is same, intercept now = 0

# Now need to write function to calculate ml CO2, given reading, and apply to start and end readings within CO2 dataframe
sampleCO2 <- function(reading,int = adj_int, slope = adj_slope){
  mlCO2 <- (reading - int)/slope
  return(mlCO2)
}
# Test it out
sampleCO2(3.4) # Function takes just reading as input, since intercept and slope assigned defaults
 
# Now, need to apply function to 'Start_reading' and 'End_reading'
CO2_calc <- CO2 %>%
  mutate(across(ends_with('reading'), .fns = list(mlCO2 = sampleCO2)))
  # Columns added with start and end readings converted to ml CO2. Probably need to change col names to something easier. 

# CO2 Measurements and calculations --------------

# ***********************************************************************************
# * <--- **** MAKE SURE PLYR, AOMISC, MASS PACKAGES ARE NOT LOADED!!! ******* ----> *
# ***********************************************************************************

# Calculate change in CO2 end-start
CO2_calc$delta_mlCO2 <- ifelse(CO2_calc$End_reading_mlCO2 - CO2_calc$Start_reading_mlCO2 > 0, CO2_calc$End_reading_mlCO2 - CO2_calc$Start_reading_mlCO2, 0) 

# Identify negative 'delta_mlCO2'
CO2_calc %>% 
  dplyr::filter(End_reading_mlCO2 < Start_reading_mlCO2) %>% 
  dplyr::select(Tube_no, Treat_ID, End_Day, delta_mlCO2, ends_with('reading'))
 # none :) above ifelse() code worked to make 0 instead

# # Calculate average background delta_mlCO2
# mean(CO2_calc$delta_mlCO2[CO2_calc$Treat_ID == 'Blank'])

# Add Day 0, CO2 = 0 points
# Make data frame with all grouping data, and End_Day = 0
day0 <- data.frame(Sample_no = 0, 
                   CO2_calc[CO2_calc$End_Day==1, 2:10], 
                   End_Day = 0)
# Add column names from CO2_total that contain 'CO2' and assign the values to 0
day0[grep('CO2', names(CO2_calc), value =T)] <- 0 # This line of code, amazing. Took way too long to find such a neat solution. 

# Merge data frames
CO2_calc <-bind_rows(CO2_calc, day0)

# Plot delta ml CO2 just to look at data
ggplot(CO2_calc, aes(x = End_Day, y = delta_mlCO2, col = as.factor(Rep))) + 
  geom_point(shape = 21, stroke = 1) + 
  facet_wrap(~Treat_ID) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 14, 2))

CO2_calc %>%
  filter(Treat_ID == 'CC_F') %>%
  dplyr::select(Tube_no, Treat_ID, End_Day, delta_mlCO2) %>% 
  arrange(End_Day) %>% 
  print(n = 100)
# Tube no. 22 is the dodgy one. Consistently lower than other reps, likely leaking jar.  

# Filter to remove tube 22
CO2_calc <- CO2_calc %>% 
  dplyr::filter(Tube_no != 22)
# Replotted and looks better
# May also need to remove CV_F rep 3 - looks a bit low

# Subtract control from treatments
CO2_calc %>% 
  dplyr::filter(Treat_ID == 'Ctrl') %>% 
  group_by(End_Day) %>% 
  dplyr::summarise(mean(delta_mlCO2))

CO2_calc <- CO2_calc %>% 
  group_by(End_Day) %>% 
  mutate(mlCO2_minusblk = delta_mlCO2-mean(delta_mlCO2[Treat_ID == 'Ctrl'])) %>% 
  mutate(mlCO2_minusblk = ifelse(mlCO2_minusblk < 0 , 0, mlCO2_minusblk))

# Check to make sure no negatives
subset(CO2_calc[-c(1,2,10:23)], mlCO2_minusblk <0 & !(Treat_ID %in% c('Blank', 'Ctrl')))
  # Just to check there are no negatives - nothing returned

# Convert ml CO2 to mg CO2-C using molar volume CO2 (L/mol)
CO2_mv <- 24.642 # Calculated using avg. temp and pressure recorded on sampling days

CO2_calc <- CO2_calc %>%
  mutate(mgCO2C = mlCO2_minusblk/CO2_mv*12.01)


# Calculate cumulative sum on CO2 for each closure, within each rep
CO2_total <- CO2_calc %>% 
  dplyr::filter(!(Treat_ID %in% c('Blank', 'Ctrl'))) %>% 
  group_by(Tube_no) %>% 
  arrange(End_Day) %>% 
  mutate(cumulativeCO2 = cumsum(mgCO2C))


# Plot cumulative CO2 to see if it looks right
ggplot(CO2_total, aes(x = End_Day, y = cumulativeCO2, col = as.factor(Rep))) + 
  geom_point(shape = 21, stroke = 1, fill = alpha('white', 0.25)) + 
  facet_grid(Soil~Fraction) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 14, 2))
  # Not sure about CV_F rep 3.. But rep 1 is also variable in other other direction


## Calculate mg CO2-C per g amendment C ------------

# Merge inc set up with CO2_total. Only want 'Applied_mgC' col. 
CO2_pergC <- left_join(CO2_total, soildat[,c(names(soildat[1:7]), 'fract_Ccont', 'fract_mgC')])

# Reset Soil and Fraction cols as factors
# Change levels in Soil column
CO2_pergC$Soil <- factor(CO2_pergC$Soil, levels = c('Ctrl', 'CC', 'CV', 'CW', 'NV'))
# Change levels in Fraction
CO2_pergC$Fraction <- factor(CO2_pergC$Fraction, levels = paste(unique(CO2_pergC$Fraction)))


# # Alternative calc, to check
# CO2_pergC1 <- CO2_pergC1 %>% 
#   mutate(mgCO2C_pergC = mgCO2C/(fract_mgC/1000))
# 
# CO2_pergC1 <- CO2_pergC1 %>% 
#   dplyr::group_by(Tube_no) %>% 
#   dplyr::arrange(End_Day) %>% 
#   mutate(cumulativeCO2_pergC2 = cumsum(mgCO2C_pergC))
# 
# CO2_pergC1 %>%
#   filter(Tube_no ==1) %>%
#   select(3,8, 26:31)


# Calculate cumulative CO2 per g added C
CO2_pergC <- CO2_pergC %>% 
  mutate(cumulativeCO2_pergC = cumulativeCO2/(fract_mgC/1000))

# Graph
ggplot(CO2_pergC, aes(x = End_Day, y = cumulativeCO2_pergC, col = Fraction, fill = Fraction)) + 
  geom_line(stat = 'smooth', se = F, linetype = 'dashed') + 
  geom_point(stat = 'summary', fun = mean, shape = 21, size = 2) + 
  facet_grid(~Soil) + 
  labs(x = 'Time (Days)', y = bquote('Cumulative' ~CO[2]~'-C emitted (mg' ~CO[2]~'-C g added-C'^'-1'*')')) +
  scale_x_continuous(breaks = seq(0, 14, 7)) + 
  ggscalecols

# ggsave('Nmin_rough cumulativeCO2_230116.png', width = 20, height = 10, unit = 'cm')


# Mean & sd cumulative CO2 over time
CO2_cumulative_avg <- CO2_pergC %>% 
  group_by(Soil, Fraction, End_Day) %>% 
  dplyr::summarise(avgCO2C = mean(cumulativeCO2_pergC), avgCO2_se = se2(cumulativeCO2_pergC)) %>%
  print(n = 200)

CO2_cumulative_avg %>% ggplot(aes(x = End_Day, y = avgCO2C, col = Fraction, shape = Soil)) + 
  stat_smooth(geom='line',linetype = 'dashed', se = F, size = 0.6, alpha = 0.4) + 
  geom_point()

CO2_cumulative_avg %>% 
  dplyr::filter(End_Day == '14') %>% 
  ggplot(aes(x = fct_reorder(Fraction, avgCO2C), y = avgCO2C, col = Fraction)) + 
  geom_point() + 
  geom_errorbar(aes(ymax = avgCO2C+avgCO2_se, ymin = avgCO2C-avgCO2_se)) + 
  facet_grid(~Soil) + 
  scale_color_manual(values = graphcol2)

CO2_cumulative_avg %>% 
  dplyr::filter(End_Day ==14)

# Reconstruct coarse+fine --------------------------

# Load fraction data - mass %
soilfrac_CN <- read.csv('Nmin_Fraction CN_230116.csv')
chtofc(soilfrac_CN)

# Calcualte CO2-C per g fraction, and multiply by fraction mass % in whole soil
CO2_pergC2 <- CO2_pergC %>% 
  dplyr::select(-contains('Start')) %>% # could remove more columns, but difficult to do tidy
  left_join(soilfrac_CN[!(names(soilfrac_CN) %in% c('C_mgg', 'N_mgg', 'Treat_ID'))], by = c('Soil', 'Fraction')) %>% 
  mutate(mgCO2_pergfract = cumulativeCO2_pergC * (fract_Ccont/1000))

# Add together coarse and fine reconstructed 
CO2_recon <- CO2_pergC2 %>% 
  mutate(mgCO2_recon = mgCO2_pergfract * frac_mass_pc) %>% 
  ungroup() %>% 
  dplyr::select(Soil, Fraction, End_Day, Rep, mgCO2_recon) %>% 
  pivot_wider(values_from = mgCO2_recon, names_from = Fraction) %>% 
  mutate(coarsefine = Coarse + Fine) %>% 
  pivot_longer(cols = Coarse:coarsefine, values_to = 'mgCO2_pergfract', names_to = 'Fraction')

# Rowbind coarsefine with CO2_pergC2
CO2_pergC3 <- bind_rows(CO2_pergC2, CO2_recon[CO2_recon$Fraction == 'coarsefine',])

chtofc(CO2_pergC3)

# Graph
ggplot(CO2_pergC3, aes(x = End_Day, y = mgCO2_pergfract, col = Fraction, fill = Fraction)) + 
  geom_line(stat = 'smooth', se = F, linetype = 'dotted') + 
  geom_point(stat = 'summary', fun = mean, shape = 21, size = 2) + 
  facet_grid(~Soil) + 
  labs(x = 'Time (Days)', y = bquote('Cumulative' ~CO[2]~'-C emitted (mg' ~CO[2]~'-C g Fraction)')) +
  scale_x_continuous(breaks = seq(0, 14, 7)) + 
  scale_colour_manual(values = c('Coarse' = '#BB281E',
                                 'Fine' = '#D56936',
                                 'Whole soil' = '#42B8D7',
                                 'Frac. Whole soil' = '#16879C',
                                 'coarsefine' = 'Grey30'))+ 
  scale_fill_manual(values = alpha(c('Coarse' = '#BB281E',
                                 'Fine' = '#D56936',
                                 'Whole soil' = '#42B8D7',
                                 'Frac. Whole soil' = '#16879C',
                                 'coarsefine' = 'Grey30'), 0.5))


# Then, average reps 
# Or alternative - pivot wider, add coarse and fine, then pivot longer again?

CO2_pergC2 %>% 
  dplyr::filter(Soil == 'CC', Fraction %in% c('Coarse', 'Fine'))
  

# Model fitting ---------------

# Some package options

# Didn't end up using this package ? - ed. - did use this package
# install.packages('devtools')
# devtools::install_github("onofriAndreaPG/aomisc") # MASS package causes problems with some dplyr functions. 

packages_CO2mods <- c('modelr', 'broom', 'aomisc')

lapply(packages_CO2mods, library, character.only = T)
# library('dplyr')

## Test fitting with single treatment -------------
# Tried a few options for self starting asymptotic regression. But can be done with base stats

# Make subset data to test out model fitting
testdat <- CO2_pergC %>% 
  dplyr::select(1:10, End_Day, cumulativeCO2_pergC, fract_mgC) %>%
  dplyr::filter(Treat_ID == "CC_C") 

# Not self start, define formula of model within nls(). Try various starting values.
  # Curious to see if similar output to SelfStart achieved, as can then use with geom_smooth(method = 'nls')
nls(cumulativeCO2_pergC ~ a*(1 - exp(-k*End_Day)), data = testdat, 
    start = list(a = max(testdat$cumulativeCO2_pergC), k = 1))


# self start using 'aomisc' package
mod <- nls(cumulativeCO2_pergC ~ NLS.asymReg(End_Day, init, m, plateau), data = testdat)
summary(mod)

fittest <- data.frame(End_Day = seq(0,14, length=50))

testgrid <- testdat %>% 
  data_grid(End_day = seq_range(End_Day, 28))

plot(resid(mod))
ggplot(testdat) + 
  geom_point(aes(End_Day, cumulativeCO2_pergC)) + 
  geom_line(aes(x = End_Day, y = predict(mod, fittest)), data = fittest)+ 
  ggtitle("mod1")

# NLS.negExp() from {aomisc} package
mod2 <- nls(cumulativeCO2_pergC ~ NLS.negExp(End_Day, a, c), data = testdat)
# I think a = C input, c = k rate constant, but not completely sure... 
# This is the expression for NLS.negExp()
# negExp.fun <- function(predictor, a, c) {
#   x <- predictor
#   a * (1 - exp (- c * x))
# 

summary(mod2)

plot(mod2) # Honestly, not great. Residuals increase with fitted values. 
hist(resid(mod2)) # Distribution ok ish
plot(Tube_no ~ resid(mod2), abline = 0, data = testdat)
# mix of mostly negative and positive residuals for Tube_no suggests should be included in model. Which is how is fitted to full data..?

mod2.3 <- nlme::gnls(cumulativeCO2_pergC ~ NLS.negExp(End_Day, a, c), data = testdat, correlation = nlme::corAR1())
mod2.4 <- nlme::gnls(cumulativeCO2_pergC ~ NLS.negExp(End_Day, a, c), data = testdat)


plot(nlme::ACF(mod2.4, maxLag = 10), alpha = 0.01) # Suggests autocorrelation evident


anova(mod2.4, mod2.3) # Adding correlation matrix doesn't improve model fit? P>0.05

ggplot(testdat) + 
  geom_point(aes(End_Day, cumulativeCO2_pergC)) + 
  geom_line(aes(x = End_Day, y = predict(mod2, fittest)), data = fittest)+ 
  ggtitle("mod2")

fittest %>% 
  add_predictions(mod2) %>% 
  ggplot() +
  geom_point(aes(End_Day, cumulativeCO2_pergC), data = testdat) + 
  geom_line(aes(End_Day, pred))

tidy(mod2)
glance(mod2)

# Maybe don't need aosmic package. Can use stats::SSasympOrigin() to self start, which sets intercept as 0. 
mod2.2 <- nls(cumulativeCO2_pergC ~ SSasympOrig(End_Day, Asym, lrc), data = testdat)
# SSasmyporig(input, asymp, lrc) - Asymp*(1 - exp(-exp(lrc)*input)) - what's the second exp()? Do I need to take exp(lrc)? for k value? # estimated log(k) to ensure positive 
summary(mod2.2)
plot(mod2.2)
exp(summary(mod2.2)$coef[2, 1]) # To get value of k 

ggplot(testdat) +
  geom_point(aes(End_Day, cumulativeCO2_pergC)) + 
  geom_line(aes(x = End_Day, y = predict(mod2.2, fittest)), data = fittest)+ 
  ggtitle("mod2.2")


# LM - just to see? 
mod3 <- lm(cumulativeCO2_pergC ~ End_Day, data = testdat)

ggplot(testdat) + 
  geom_point(aes(End_Day, cumulativeCO2_pergC)) + 
  geom_line(aes(x = End_Day, y = predict(mod3, fittest)), data = fittest)+ 
  ggtitle("linear mod")


## Apply model to all treatments CO2 data -------------

# Use nested data - nested by Tube_no - so each rep is fitted individually
CO2_nest <- CO2_pergC %>% 
  dplyr::select(c(1:10, End_Day, cumulativeCO2_pergC, fract_mgC)) %>% 
  group_by(Tube_no, Soil, Fraction, Treat_ID, Rep) %>% # Tube no 
  nest()

CO2_nest$data[[1]]

# write function to fit model to each data frame

CO2mod <- function(df){
  nls(cumulativeCO2_pergC ~ NLS.negExp(End_Day, a, k), data = df)
}

# Use purrr::map to apply CO2mod to each df. This adds another col to nested data with models for each tube. 
CO2_nest <- CO2_nest %>% 
  mutate(model = map(data, CO2mod))


# Add residuals for each model
CO2_nest <- CO2_nest %>% 
  mutate(resids = map2(data, model, add_residuals))

# Unnest residuals
CO2mod_resids <- unnest(CO2_nest, resids)
# Plot residuals
ggplot(CO2mod_resids, aes(End_Day, resid, col = as.factor(Rep))) + 
  geom_point() + 
  facet_wrap(~Treat_ID)

# Add predictions for each model
CO2_preds <- CO2_nest %>% 
    mutate(preds = map2(data, model, add_predictions))

CO2mod_preds <- unnest(CO2_preds, preds)


# Graph fitted models --------------

# Plot all reps - data and model predictions
ggplot(CO2mod_preds, aes(End_Day, pred, col = Fraction, fill = Fraction, group = Tube_no)) + 
  geom_line(stat = 'smooth', method = 'loess', se = F, alpha = 0.75, linetype = 5) +
  # geom_smooth(method = 'loess', se = F, alpha = 0.75) + # function to add model 
  geom_point(aes(y=cumulativeCO2_pergC, x = End_Day), size = 1.5, shape =21) + 
  facet_grid(Fraction~Soil) + 
  labs(x = 'Time (Days)', y = bquote('Cumulative' ~CO[2]~'-C (mg' ~CO[2]~'-C g added-C'^'-1'*')')) +
  scale_color_manual(values = graphcol2) +
  scale_fill_manual(values = alpha(graphcol2, 0.5)) + 
  scale_x_continuous(breaks = seq(0, 14, 7))


# Plot average - data and model predictions
ggplot(CO2mod_preds, aes(End_Day, pred, col = Fraction, fill = Fraction)) + 
  geom_line(stat = 'smooth', method = 'loess', se = F, alpha = 0.75) +
  # geom_smooth(method = 'loess', se = F, alpha = 0.75) + # function to add model 
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', aes(y=cumulativeCO2_pergC, x = End_Day), width =0, size = 0.25) + 
  geom_point(stat = 'summary', fun= 'mean', aes(y=cumulativeCO2_pergC, x = End_Day), size = 1.5, shape =21) + 
  facet_grid(~Soil) + 
  labs(x = 'Time (Days)', y = bquote('Cumulative' ~CO[2]~'-C (mg' ~CO[2]~'-C g added-C'^'-1'*')')) +
  scale_color_manual(values = graphcol2) +
  scale_fill_manual(values = alpha(graphcol2, 0.5)) + 
  scale_x_continuous(breaks = seq(0, 14, 7))

# Plot all reps - data and model predictions. Using 'proper' nls() fit in geom_smooth, rather than plotting predictions and using loess smooth.  
ggplot(CO2_pergC, aes(x = End_Day, y = cumulativeCO2_pergC, col = Fraction, fill = Fraction, group = Tube_no)) + 
  geom_line(stat = 'smooth', se = F, alpha = 0.75, linetype = 5, 
            method = 'nls', formula =  y ~ a*(1 - exp(-k*x)), 
            method.args = list(start = list(a = max(CO2_pergC$cumulativeCO2_pergC), k = 0.5))) +
  # geom_smooth(method = 'loess', se = F, alpha = 0.75) + # function to add model 
  geom_point(aes(y=cumulativeCO2_pergC, x = End_Day), size = 1.5, shape =21) + 
  facet_grid(Fraction~Soil) + 
  labs(x = 'Time (Days)', y = bquote('Cumulative' ~CO[2]~'-C (mg' ~CO[2]~'-C g added-C'^'-1'*')')) +
  scale_color_manual(values = graphcol2) +
  scale_fill_manual(values = alpha(graphcol2, 0.5)) + 
  scale_x_continuous(breaks = seq(0, 14, 7))

## MS fig -----------
# Plot average - data and model predictions. Using 'proper' nls() fit in geom_smooth, rather than plotting predictions and using loess smooth.  
  # This definitely looks better than the loess fit.
CO2_pergC %>% 
  # dplyr::filter(Soil != 'NV') %>% 
ggplot(aes(End_Day, cumulativeCO2_pergC, col = Fraction, fill = Fraction)) + 
  geom_line(stat = 'smooth', se = F, alpha = 0.75, fullrange = F,
            method = 'nls', formula =  y ~ a*(1 - exp(-k*x)), 
            method.args = list(start = list(a = max(CO2_pergC$cumulativeCO2_pergC), k = 0.5))) +
  geom_errorbar(stat = 'summary', fun.data = 'mean_se', aes(y=cumulativeCO2_pergC, x = End_Day), width =0, linewidth = 0.25) + 
  geom_point(stat = 'summary', fun= 'mean', aes(y=cumulativeCO2_pergC, x = End_Day), size = 1.5, shape =21) + 
  facet_grid(~Soil) +
  labs(x = 'Time (Days)', y = bquote('Cumulative' ~CO[2]~'-C(mg' ~CO[2]~'-C g-C'^'-1'*')')) +
  scale_color_manual(values = graphcol2) +
  scale_fill_manual(values = alpha(graphcol2, 0.5)) + 
  scale_x_continuous(breaks = seq(0, 14, 7))


# ggsave(paste0('Graphs/MS/Nmin_cumulativeCO2_', format(Sys.Date(), "%y%m%d"), '.png'), height = 8, width = 19, unit = 'cm')
# ggsave(paste('Graphs/MS/Fig4','.pdf', sep = ""),height = 8, width = 19, units = 'cm')


# , labeller = as_labeller(c('CC' = 'Continuous cotton', 'CV' = 'Cotton-vetch', 'CW' = 'Cotton-wheat', 'NV' = 'Native veg'))
# ggsave(paste0('Graphs/Nmin_pres_CO2_', format(Sys.Date(), "%y%m%d"), '.png'), height = 8, width = 20, unit = 'cm')


# Extract k and a values from models -------------

# Turn coefficients (a, k) into tidy data
CO2_nest <- CO2_nest %>% 
  mutate(coefs = map(model, tidy)) # tidy makes tibble where each row has information about model coefficients. 

CO2_modcoefs <- CO2_nest %>% 
  dplyr::select(-c(data, model, resids)) %>% 
  unnest(coefs) %>%  # all outputs from model in columns: term (a or k), estimate (of a or k), se, t stat, p val.
  dplyr::arrange(term)
  
CO2_modcoefs_stat <- CO2_modcoefs %>% 
pivot_wider(names_from = term, values_from = c(estimate, std.error, statistic, p.value)) # Now have a and k estiamtes in separate cols. 

# CO2_nest %>% 
#   mutate(glance = map_df(model, glance)) %>% 
#   pull(glance)

# function for formating, pasting and rounding
roundform <- \(x, r, n) format(round(x, digits = r), nsmall = n)


# Mean and se
CO2_modcoefs_avg <- CO2_modcoefs %>% 
  group_by(Soil, Fraction, term) %>% 
  dplyr::summarise(across(contains('estimate'), list(mean = mean, se = se2))) %>% 
  mutate(mean_se = ifelse(term =='a', paste0(roundform(estimate_mean, 1, 1), ' (', roundform(estimate_se, 1, 1),')'), 
                           paste0(roundform(estimate_mean, 2, 2), ' (',roundform(estimate_se, 2, 2),')'))) %>% 
  dplyr::arrange(term) %>% 
  print(n = 100)

# CO2_modcoefs_avg %>% write.csv(paste0('Nmin_CO2_coefs_', format(Sys.Date(), '%y%m%d'), '.csv'))


# Run statistical analysis on k (rate constant) and a (size of labile C pool)

# Packages 
lapply(c('easystats', 'emmeans'), library, character.only = T)

# k values ----------------

# k values
CO2_modcoefs %>%  
  dplyr::filter(term == 'k') %>% 
  ggplot(aes(x = fct_reorder(Fraction, estimate, .fun = 'mean'), y = estimate, col = Fraction)) + 
  geom_point() + 
  geom_point(stat = 'summary', fun = 'mean', size = 3, alpha = 0.5) + 
  geom_errorbar(stat = 'summary', fun.data = 'mean_se') +
  facet_grid(~Soil) +
  ggscalecols + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Anova k -------------------------------

# Test normallity
shapiro.test(CO2_modcoefs_stat$estimate_k) # No evidence that significant;y different to normal 
hist(CO2_modcoefs_stat$estimate_k)

str(CO2_modcoefs_stat, max.level = 2)

# Fit linear model for anova
mod_kest1 <- lm(estimate_k ~ Fraction + Soil, data = CO2_modcoefs_stat)
summary(mod_kest1)
check_model(mod_kest1) # Looks good
check_heteroscedasticity(mod_kest1) # No heteroscedasticity evident
check_normality(mod_kest1) # non-normality detected...

# Also check interaction
mod_kest2 <- lm(estimate_k ~ Fraction * Soil, data = CO2_modcoefs_stat)

summary(mod_kest2)

anova(mod_kest1, mod_kest2) # P > 0.05, interaction term not significant. Proceed with main effects mod

mod_kest1.2 <- update(mod_kest1, log(estimate_k) ~ Fraction + Soil)
summary(mod_kest1.2)
check_model(mod_kest1.2)
check_normality(mod_kest1.2) # log transform didn't help.. visual doesn't look overly far from normality

# Anova mod1
anova(mod_kest1)
  # Significant effect of Soil and Fraction


## Post-hoc k --------------------
kest_emm <- emmeans(mod_kest1, specs = list('Soil', 'Fraction'))

map(kest_emm, plot, comparisons = T)

# Fraction and Soil main effect
map(kest_emm, contrast, method = 'pairwise')

pairs(kest_emm) #%>% tidy %>% arrange(adj.p.value) %>% print(n =30)

pwpp(kest_emm)
pwpm(kest_emm)

# Differences between soils, within fractions
multcomp::cld(kest_emm, Letter = letters)
# Differences between fractions, within soils
multcomp::cld(emmeans(mod_kest1, ~ Fraction | Soil), Letter = letters)

map(kest_emm, multcomp::cld, Letters = letters)

estimate_contrasts(mod_kest1, contrast = 'Fraction')
estimate_contrasts(mod_kest1, contrast = 'Soil')


# a values -----------

CO2_modcoefs %>%  
  dplyr::filter(term == 'a') %>% 
  ggplot(aes(x = fct_reorder(Fraction, estimate, .fun = 'mean'), y = estimate, col = Fraction)) + 
  geom_point() + 
  geom_point(stat = 'summary', fun = 'mean',size = 3, alpha = 0.5) + 
  geom_errorbar(stat = 'summary', fun.data = 'mean_se') +
  facet_grid(~Soil) +
  ggscalecols + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


## Anova a --------------

shapiro.test(CO2_modcoefs_stat$estimate_a) # No evidence that different from normal
hist(CO2_modcoefs_stat$estimate_a)

str(CO2_modcoefs_stat, max.level = 2)

# Fit linear model for anova
mod_aest1 <- lm(estimate_a ~ Fraction + Soil, data = CO2_modcoefs_stat)
summary(mod_aest1)
check_model(mod_aest1) # Looks good
check_heteroscedasticity(mod_aest1) # No heteroscedascticity evident

# Also check interaction
mod_aest2 <- lm(estimate_a ~ Fraction * Soil, data = CO2_modcoefs_stat)

summary(mod_aest2)

anova(mod_aest1, mod_aest2) # P < 0.05, interaction term significant

# check assumptions
check_model(mod_aest2)
check_heteroscedasticity(mod_aest2)
check_normality(mod_aest2) # Non-normality of residuals detected

# Try transformation
mod_aest2.2 <- update(mod_aest2, log(estimate_a) ~ Fraction * Soil)
check_normality(mod_aest2.2) # Still non-normality
check_model(mod_aest2.2) # Not seeing much improvement - go with untransformed


# Anova mod1
anova(mod_aest2)
# Significant main effect of fraction, interaction significant. Soil main effect not significant. 


## Post-hoc a --------------------
aest_emm <- emmeans(mod_aest2, pairwise ~ Fraction | Soil)

plot(aest_emm, comparisons = T)

contrast(aest_emm, by = 'Soil', method = 'pairwise')
estimate_contrasts(mod_aest2, contrast = 'Fraction', at = 'Soil')

# Comparisons for Fraction within Soil
multcomp::cld(aest_emm, Letter = letters)

# No significnat effect of soil, only look at fraction comparisons
# contrast(aest_emm, by = 'Fraction', method = 'pairwise')
# estimate_contrasts(mod_aest2, contrast = 'Soil', at = 'Fraction')

# contrast(aest_emm, 'pairwise', simple = 'each', combine = T, adjust = 'Tukey')

pwpp(aest_emm)
pwpm(aest_emm)

#Comparisons for Soil 
multcomp::cld(emmeans(mod_aest2, ~  Soil | Fraction), Letter = letters) #%>%  tidy()# %>%
  # dplyr::select(Soil, Fraction, estimate) %>% 
  # pivot_wider(names_from = Soil, values_from = estimate, )

estimate_contrasts(mod_aest2, contrast = 'Fraction', at = 'Soil')
# estimate_contrasts(mod_aest2, contrast = 'Soil')

