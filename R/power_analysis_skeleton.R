# This is a generalizable skeleton for simulation based power analyses.

library(tidyverse)
library(magrittr)

#### Core simulation functions ----

# The first function here defines how we're simulating data. As an example I
# simulate the power to detect differences between two normally distributed
# groups with two varying parameters: mean shift and sd of the second group. The
# pwr package actually has functionality for calculating something this simple
# exactly, but this is just for demonstration. In practice you'd run a
# simulation when the outcome of the statistical test is less straightforward or
# when the you understand the data simulation process but not how that would
# impact the output measurements you'd run the test on. Other parameters that
# one might change: strength of a linear association, number of measurements per
# group, etc.

sim_once = function(mean_shift, sd2){
  # Return a dataframe of observations

  # Accurately describing the generative process is critical to the success of the
  # simulation.
  tibble(group_id = c(rep('A', 10),
                      rep('B', 10)),
         value = c(rnorm(10),
                   rnorm(10, mean = mean_shift, sd = sd2)))

}

sim_n_times = function(mean_shift, sd2, n_sim){
  # This just runs the above function n_sim times and returns it as a labelled
  # dataframe
  tibble(iteration = 1:n_sim,
         iter_result = map(iteration, ~sim_once(mean_shift, sd2)))
}

#### Set up parameter value grid ----

# I'm going to have mean_shift take 10 values from 0 to 2 and sd2 take 10 values
# from .2 to 2. I'm only doing 100 simulations of each grid point (10x10x100 =
# 1e5 simulations). 1000 at each point would be preferable. It's easy to see how
# additional varying parameters can rapidly increase the computational cost. You
# also need to ensure that your grid covers a range between A) where detecting a difference
# is impossible to B) the smallest parameter values where the power is essentially
# 100% (I haven't done that in this case). This command sets up the grid.

sim_df = expand.grid(list(mean_shift = seq(0, 2, length.out = 10),
                          sd2 = seq(.2, 2, length.out = 10))) %>%
  as_tibble

#### Simulate at each grid point ----

# Here we run the simulations for each grid point. For large simulations, change the
# map2 command to
# parallel::mcmapply(sim_n_times,
#                    mean_shift, sd2,
#                    MoreArgs = list(n_sim = n_per_point),
#                    SIMPLIFY = FALSE,
#                    mc.cores = 4) # or however many

n_per_point = 100 # More is better

sim_df %<>%
  mutate(point_sims = map2(mean_shift, sd2,
                           sim_n_times,
                           n_sim = n_per_point)) %>%
  unnest

# ^ You could also rewrite this to avoid unnest() by putting rep(grid_range)
# into the command creating sim_df and just calling sim_once directly as you
# map2 over the inputs.

# save(sim_df) # is probably wise here

#### Analyze iteration results ----

# Now we run whatever statistical processing / test we plan to use on each
# iteration. I'm just going to pretend it's a t-test between the two groups, and
# immediately extract the p-value, but it could be a lm() or something fancier,
# and it might be helpful to separate out the test results and summary
# statistics into separate columns.

# The t-test in this example is quite fast, but if you're testing the efficacy
# of Bayesian methods that are themselves simulation based, this step can also
# be slow and thus potentially worth parallelizing.

run_iter_test = function(iter_result){
  t.test(value ~ group_id,
         data = iter_result)$p.value
}

sim_df %<>%
  mutate(iter_p = map_dbl(iter_result,
                          run_iter_test))

# save() here...

#### Summarise grid points ----

# Here we apply some threshold that we would use in practice to say we detected
# a difference. p < .05 for this example. Then we compute the fraction of
# iterations at each grid point that fall below that threshold. That's our power
# (except when the null hypothesis is true, in which case it's the false
# positive rate).

summary_df = sim_df %>%
  mutate(below_threshold = iter_p < .05) %>%
  group_by(mean_shift, sd2) %>%
  summarise(simulated_power = mean(below_threshold)) %>% # mean of a logical vector is just the fraction of TRUEs
  ungroup

#### visualize ----

# Here we create curves showing the effect on the simulated power as a function
# of the varying parameters. I have one subplot for each sd2 value and
# mean_shift on the x-axis, but the visualization will need to be adapted to the
# structure of the simulation. More simulations per grid point and a denser grid
# will yield smoother curves, at the cost of more computation.

summary_df %>%
  ggplot(aes(mean_shift, simulated_power)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = .8, # a common
             lty = 2,
             color = 'grey70') +
  facet_wrap('sd2', labeller = label_both) +
  labs(title = 'Conclusion: Power increases with increasing mean_shift and decreasing sd2',
       subtitle = '80% power requires mean_shift>1 at sd2=0.2, but >2 at sd2>1.8') +
  theme_light() +
  theme(strip.text = element_text(color = 'black'))



