library(Rcpp)
library(dplyr)
library(ggplot2)

sourceCpp("functions/main___mcmc.cpp")


# -----------------------------------------------------------------------------
# Reading Data ----------------------------------------------------------------

Y = read.table(
  file = paste0("datasets/woak.mat"),
  header = TRUE, quote = "\"") %>% 
  as.matrix()
colnames(Y) = NULL

# Initial value for the GP
Y = cbind(Y, (8 / 3) * exp( (Y[, 1]^2) / 30) +
            (4 / 3) * exp(-((Y[, 2] - 7)^2) / 12) - 2)


# -----------------------------------------------------------------------------
# Running MCMC ----------------------------------------------------------------

mcmc(
  Y = Y,
  H = 10,
  W = 10,

  path_to_save = "results/",
  niter = 5500,
  niter_to_save = 500,
  seed = 1,

  lambdastar_estimated = TRUE,
  lambdastar_init = 3,
  lambdastar_prior_a =  1,
  lambdastar_prior_b = 0.1,

  link_function = 'probit',

  gaussian_processes_estimated = TRUE,
  gaussian_processes_means = 0,
  gaussian_processes_variances = 4,
  gaussian_processes_ranges = 0.5,
  gaussian_processes_exponent = 3/2,

  additional_locations__lp_if = FALSE,
  additional_locations = matrix(ncol = 2),

  intensity_function_estimated = TRUE,
  intensity_function_grid_size = 101 * 101

)


# -----------------------------------------------------------------------------
# Treating results ------------------------------------------------------------

# Lambdastar ------------------------------------------------------------------
all_files = list.files(path = "results", pattern = "chain__lambdastar")
txt_files = all_files[grepl('.txt', all_files)] %>% sort()

chain__lambdastar = NULL
for (f in txt_files) {
  aux = read.table(file = file.path("results", f)) %>% as.matrix()
  chain__lambdastar = rbind(chain__lambdastar, aux)
}

saveRDS(object = chain__lambdastar,
        file = "results/chain__lambdastar.rds")

# Intensity Function ----------------------------------------------------------
all_files = list.files(path = "results", pattern = "intensity_function")
txt_files = all_files[
  grepl('.txt', all_files) & !grepl("squared", all_files)] %>% sort()

chain__intensity_function = array(dim = c(101 * 101, 3, length(txt_files)))
for (i in 1:length(txt_files)) {
  chain__intensity_function[ , , i] = read.table(
    file = file.path("results", txt_files[i])) %>% as.matrix()
}

saveRDS(object = chain__intensity_function,
        file = "results/chain__intensity_function.rds")


# -----------------------------------------------------------------------------
# Graphs ----------------------------------------------------------------------

lam_df <- chain__lambdastar %>% 
  data.frame() %>% 
  rename(lambdastar = V1) %>% 
  mutate(iteration = row_number()) %>% 
  filter(iteration > 500)
lam_df %>% 
  ggplot(aes(x = iteration, y = lambdastar)) +
  geom_line() +
  theme_minimal()


int_df <- chain__intensity_function[ , , 2:dim(chain__intensity_function)[3]] %>% 
  apply(c(1, 2), mean) %>% 
  data.frame() %>% 
  rename(x = X1, y = X2, intesity_function = X3) 
int_df %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = intesity_function)) +
  scale_fill_gradientn(colors = heat.colors(100, rev = TRUE),
                       limits = c(0, 13)) +
  geom_point(data = Y %>% data.frame(), aes(x = X1, y = X2), size = 0.5) +
  coord_fixed() +
  theme_minimal()
ggsave("intensity_full.pdf")

save(int_df, lam_df, file = "full.RData")
