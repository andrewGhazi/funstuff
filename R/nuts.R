
#### Libraries ----
library(plot3D)
library(gganimate)
library(magrittr)
library(tidyverse)
library(scales)
library(rstanarm)
library(rstan)

#### Motivating example ----

# This is the fit I'm going to try to imitate. A linear regression of mtcars mpg
# ~ dist. I scale the two variables to make setting the mass matrix and so on
# easier.

summary(lm(mpg ~ disp, data = mtcars))
x = scale(mtcars$disp)
y = scale(mtcars$mpg)
rstanarm::stan_glm(y ~ x, 
                   prior_aux = exponential(autoscale = FALSE),
                   prior_intercept = normal(scale = 2, autoscale = FALSE),
                   prior = normal(0, 1, autoscale = FALSE)) %>% 
  .$stanfit %>% 
  summary() %>% 
  .$summary

#### Model function ---- 

# I think log_post_theta() is all that one would need to change to apply this code to other
# datasets/models, though I can't guarantee that I didn't bake in the
# dimensionality at a few other spots.
log_post_theta = function(theta, y, x){
  alpha = theta[1]
  beta = theta[2]
  sigma = theta[3]
  
  log_prior = dnorm(alpha, mean = 0, sd = 2, log = TRUE) + # N(0,2) prior on intercept
    dnorm(beta, sd = 1, log = TRUE) + # N(0,1) prior on slope
    dexp(sigma, rate = 1, log = TRUE) # Exp(1) prior on SE
  
  ll = sum(dnorm(y, mean = alpha + beta*x, sd = sigma, log = TRUE))
  
  return(log_prior + ll)
}

gradient_th_numerical = function(theta, y, x){
  d = length(theta)
  e = .0001
  diff = rep(NA, d)
  for (k in 1:d){
    th_hi = theta
    th_lo = theta
    th_hi[k] = theta[k] + e
    th_lo[k] = theta[k] - e
    diff[k] = (log_post_theta(th_hi, y, x) - log_post_theta(th_lo, y, x)) / (2*e)
  }
  
  return(diff)
}

#### NUTS functions ----
leapfrog = function(theta, y, x, mom, eps){
  mom_half = mom + eps/2 * gradient_th_numerical(theta, y, x)
  theta_tilde = theta + eps * mom_half
  mom_tilde = mom_half + eps/2 * gradient_th_numerical(theta_tilde, y, x)
  return(c(theta_tilde, mom_tilde))
}

log_marg = function(state, y, x){
  d = length(state) / 2
  
  log_post_theta(state[1:d], y, x) - 1/2 * sum(state[(d+1):(2*d)]^2)
}


build_tree = function(theta, mom, u, direction, tree_height, eps, delta_max = 1000){
  # returns a list of 4 elements:
  # 1 state_minus, a concatenated vector of c(theta_minus, mom_minus)
  # 2 state_plus, similarly
  # 3 C,  a matrix containing valid states on each row. Possibly 0 rows
  # 4 s, stopping indicator (s = 0 means stop)
  
  d = length(theta)
  th_ind = 1:d
  mom_ind = (d+1):(2*d)
  if (tree_height == 0) {
    state_prime = leapfrog(theta, y, x, mom, direction*eps)
    log_prime = log_marg(state_prime, y, x)
    C_prime = if (u <= exp(log_prime)) matrix(state_prime, ncol = 2*d) else matrix(ncol = d*2, nrow = 0)
    s_prime = ifelse(log_prime > (log(u) - delta_max),1,0)
    return(list(state_prime, state_prime, C_prime, s_prime))
  } else {
    
    lower_state = build_tree(theta, mom, u, direction, tree_height - 1, eps)
    if (direction == -1){
      minus_state = build_tree(lower_state[[1]][1:d], lower_state[[1]][mom_ind], u, direction, tree_height - 1, eps)
      lower_state[[1]] = minus_state[[1]] # theta_minus, mom_minus
      C_dprime = minus_state[[3]]
      s_dprime = minus_state[[4]]
    } else {
      plus_state = build_tree(lower_state[[2]][1:d], lower_state[[2]][mom_ind], u, direction, tree_height - 1, eps)
      lower_state[[2]] = plus_state[[2]] # theta_minus, mom_minus
      C_dprime = plus_state[[3]]
      s_dprime = plus_state[[4]]
    }
    theta_plus = lower_state[[2]][1:d]
    theta_minus = lower_state[[1]][1:d]
    mom_plus = lower_state[[2]][mom_ind]
    mom_minus = lower_state[[1]][mom_ind]
    
    cause_one = lower_state[[4]] 
    cause_two = s_dprime
    cause_three = ifelse(drop((theta_plus - theta_minus) %*% mom_plus) >= 0, 1, 0)
    cause_four = ifelse(drop((theta_plus - theta_minus) %*% mom_minus) >= 0, 1, 0)
    
    s_prime = lower_state[[4]] * s_dprime * 
      cause_three *
      cause_four
    
    C_prime = rbind(lower_state[[3]], C_dprime)
    
    # if (tree_height > 1){C_prime = C_prime[-1,]}
    
    return(list(lower_state[[1]],
                lower_state[[2]],
                C_prime,
                s_prime,
                c(cause_one, cause_two, cause_three, cause_four)))
    
  }
}

#### Chain setup ----
inv_mass = rep(1, 3)
d = 3
M = 500
theta = matrix(nrow = M, ncol = d)
mom = matrix(nrow = M, ncol = d)
mom_ind = (d+1):(2*d)
theta_nought = c(0, -.847, .539) # lm(y ~ x) %>% summary
theta[1,1:d] = theta_nought
eps = 8e-4
C_list = vector('list', M-1)
u_list = vector('list', M-1)
selected = vector('numeric', M-1)
for (m in 2:M){
  # This is the NUTS chain
  
  r_init = inv_mass * rnorm(3)
  init_state = c(theta[m-1,], r_init)
  slice_prob = exp(log_marg(init_state, y, x))
  u = runif(1, min = 0, max = slice_prob)
  
  theta_minus = theta[m-1,]
  theta_plus = theta[m-1,]
  mom_minus = r_init
  mom_plus = r_init
  j = 0
  C = matrix(ncol = d*2 + 1,
             data = c(theta[m-1,], r_init, j))
  s = 1

  while (s == 1) {
    vj = sample(c(-1,1), size = 1) # the direction
    if (vj == -1) {
      ot = theta_minus
      om = mom_minus
      new_state = build_tree(theta_minus, mom_minus, u, vj, j, eps)
      theta_minus = new_state[[1]][1:d]
      mom_minus = new_state[[1]][mom_ind]
      C_prime = new_state[[3]]
      s_prime = new_state[[4]]
    } else {
      ot = theta_plus
      om = mom_plus
      new_state = build_tree(theta_plus, mom_plus, u, vj, j, eps)
      theta_plus = new_state[[2]][1:d]
      mom_plus = new_state[[2]][mom_ind]
      C_prime = new_state[[3]]
      s_prime = new_state[[4]]
    }
    
    C_prime = cbind(C_prime, matrix(ncol = 1, data = rep(j, nrow(C_prime))))
    
    if (s_prime == 1){
      C = rbind(C, C_prime)
    }

    
    s = s_prime * 
      as.numeric(drop((theta_plus - theta_minus) %*% mom_minus) >= 0) * 
      as.numeric(drop((theta_plus - theta_minus) %*% mom_plus) >= 0)
    j = j + 1
    
  }
  
  selected[m-1] = sample(1:nrow(C), 1)
  # prev_lp = log_post_theta()
  # prop_lp = log_post_theta() # ah screw the metropolis step
  C_list[[m-1]] = C
  u_list[[m-1]] = u
  C_draw = C[selected[m-1],,drop = TRUE]
  theta[m,] = C_draw[1:d]
  mom[m,] = C_draw[mom_ind]
  
  if(is.na(C_draw[1])) break
  
}

nuts_df = tibble(C = C_list %>% map(as_tibble) %>% map(~set_colnames(.x, c('x', 'y', 'z', 'dv1', 'dv2', 'dv3', 'j'))),
                 chosen = selected)
nuts_df %>%
  mutate(s = map2(C, chosen, ~.x[.y,])) %>%
  unnest(c(s)) %>%
  select(x:z) %>% 
  summarise_all(mean)
# ^ basically the same result as rstanarm, hooray

#### Animation functions ----
rotate = function(iter_dat, theta, var_ranges){
  # adapted from https://github.com/AckerDWM/gg3D/blob/master/R/stat_3D.R
  pmat = perspbox(z = diag(2), theta = theta, phi = 0, plot = FALSE)
  xrange = range(iter_dat$x)
  yrange = range(iter_dat$y)
  zrange = range(iter_dat$z)
  iter_dat %<>%
    mutate(
      x = rescale(x, from=var_ranges$x[[1]], to=c(0,1)),
      y = rescale(y, from=var_ranges$y[[1]], to=c(0,1)),
      z = rescale(z, from=var_ranges$z[[1]], to=c(0,1)))
  
  # V First eight points are for the segments, then points for little "x", "y", and "z" labels
  mini_axes = tibble(x = c(.5, .6, .5, .5, .58, .58, .49, .51, .63,  .5,  .5), # Start and end at center of x/y, 0 for z
                     y = c(.5, .5, .6, .5, .51, .49, .58, .58,  .5, .63,  .5),
                     z = c( 0,  0,  0, .1,   0,   0,   0,   0,   0,   0, .12))
  
  xy_dat = trans3D(iter_dat$x,
                   iter_dat$y,
                   iter_dat$z,
                   pmat)
  
  mini_axes_xy = trans3D(mini_axes$x,
                         mini_axes$y,
                         mini_axes$z,
                         pmat) %>% 
    {tibble(x = .$x,
           y = .$y,
           part = c(rep('seg', 8),
                    rep('label', 3)))}
  mini_segments = mini_axes_xy[1:8,]
  mini_segments = tibble(x = c(rep(mini_segments$x[1], 3), rep(mini_segments$x[2], 2), rep(mini_segments$x[3],2)),
                         y = c(rep(mini_segments$y[1], 3), rep(mini_segments$y[2], 2), rep(mini_segments$y[3],2)),
                         xend = mini_segments$x[2:8],
                         yend = mini_segments$y[2:8])
  mini_labels = mini_axes_xy[9:11,]
  
  
  rot_input = tibble(x = xy_dat$x,
                     y = xy_dat$y) %>%
    bind_cols(iter_dat %>% select(-x, -y, -z))
  
  res_list = list("rot_input" = rot_input,
                  "mini_segments" = mini_segments,
                  "mini_labels" = mini_labels)
  return(res_list)
  
  
}

get_frame_rows = function(frame_n, iter_dat, n_frames, id_vec){
  data_subset = dplyr::filter(iter_dat, id/max(id) <= frame_n / n_frames)
  if (nrow(data_subset) == 0){
    data_subset = iter_dat[1,]
  }
  return(data_subset)
}

plot_iter_animation = function(iter_dat, n_frames = 50, theta_max = 15){
  theta_range = seq(0, theta_max, length.out = n_frames)
  
  iter_dat %<>% mutate(id = 1:n())
    # mutate(speed = c(0, sqrt(diff(x)^2 + diff(y)^2 + diff(z)^2)), 
    #        )
  
  # iter_dat$speed[iter_dat$speed > .0021] = mean(iter_dat$speed[iter_dat$speed < .001]) #TODO make better
  
  var_ranges = iter_dat %>% 
    select(x:z) %>% 
    summarise_all(.funs = ~list(range(.x)))
  
  frame_df = tibble(frame = 1:n_frames, 
                    theta = theta_range,
                    data_to_rotate = map(frame, get_frame_rows, iter_dat = iter_dat, 
                                         n_frames = n_frames, id_vec = id))
  
  rotated_dat = frame_df %>% 
    mutate(xy_dat = map2(data_to_rotate, theta, rotate, 
                         var_ranges = var_ranges)) %>%
    mutate(mini_axes = map(xy_dat, 2),
           mini_labels = map(xy_dat, 3),
           xy_dat = map(xy_dat, 1)) %>% 
    select(-data_to_rotate) 
  
  rotated_points = rotated_dat %>% 
    arrange(-frame) %>% 
    select(-mini_axes, -mini_labels) %>% 
    unnest(c(xy_dat))
  
  rotated_mini_axes = rotated_dat %>% 
    select(-xy_dat, -mini_labels) %>% 
    unnest(c(mini_axes))
  
  rotated_labels = rotated_dat %>% 
    select(-xy_dat, -mini_axes) %>% 
    unnest(c(mini_labels)) %>% 
    mutate(label = rep(c('alpha', 'beta', 'sigma'),
                       n_frames))
  
  tree_plot = rotated_points %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = j),
               size = 5,
               alpha = .5) +
    scale_color_viridis_c(end = .92) + 
    geom_segment(data = rotated_mini_axes,
                 aes(x = x, 
                     xend = xend,
                     y = y, 
                     yend = yend)) + 
    geom_text(data = rotated_labels,
              aes(x = x, y = y, label = label)) + 
    transition_manual(frame) +
    theme_bw() + 
    labs(color = 'tree depth') + 
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

  animate(tree_plot,
          width = 1600,
          height = 900, nframes = n_frames, renderer = ffmpeg_renderer())

}


# C_list %>% map_dbl(nrow) 
# c_df = C_list[[175]] %>% 
#   as_tibble %>% 
#   set_names(c('x', 'y', 'z', 'dv1', 'dv2', 'dv3', 'j'))
# 
# iter_dat = c_df
# plot_iter_animation(c_df, n_frames = 32, theta_max = 20)

#### Join multiple iterations into one animation ----

# find the maximal window width needed, then repeatedly animate the evolutions
# Highlight the selected point in Green and pause on it for a few frames between iterations

plot_one = function(iter_df, frame_range, theta_range, var_ranges, selected, green_pause, previous_points){
  
  n_frames = length(frame_range)
  evo_frames = n_frames - green_pause
  iter_df %<>% mutate(id = 1:n())

  evo_df = tibble(frame = 1:(evo_frames), 
                  theta = theta_range[1:evo_frames],
                  global_frame = frame_range[1:evo_frames],
                  data_to_rotate = map(frame, get_frame_rows, iter_dat = iter_df, 
                                       n_frames = evo_frames, id_vec = iter_df$id))
  
  green_frames = (evo_frames + 1):n_frames
  green_df = tibble(frame = green_frames,
                    theta = theta_range[green_frames],
                    global_frame = frame_range[green_frames],
                    data_to_rotate = list(iter_df))
  
  rotated_dat = evo_df %>% 
    mutate(xy_dat = map2(data_to_rotate, theta, rotate, 
                         var_ranges = var_ranges)) %>%
    mutate(mini_axes = map(xy_dat, 2),
           mini_labels = map(xy_dat, 3),
           xy_dat = map(xy_dat, 1)) %>% 
    select(-data_to_rotate) 
  green_rotated = green_df %>% 
    mutate(xy_dat = map2(data_to_rotate, theta, rotate, 
                         var_ranges = var_ranges)) %>%
    mutate(mini_axes = map(xy_dat, 2),
           mini_labels = map(xy_dat, 3),
           xy_dat = map(xy_dat, 1)) %>% 
    select(-data_to_rotate) 
  
  rotated_points = rotated_dat %>% 
    arrange(-frame) %>% 
    select(-mini_axes, -mini_labels) %>% 
    unnest(c(xy_dat)) %>% 
    mutate(green_hold = FALSE, 
           selected_point = id == selected)
  green_points = green_rotated %>% 
    arrange(-frame) %>% 
    select(-mini_axes, -mini_labels) %>% 
    unnest(c(xy_dat)) %>% 
    mutate(green_hold = TRUE, 
           selected_point = id == selected)
  
  if (nrow(previous_points) == 0){
    prev_rotated = rotated_points[0,]
  } else {
    prev_df = tibble(frame = 1:n_frames,
                     theta = theta_range,
                     global_frame = frame_range,
                     data_to_rotate = list(previous_points))
    prev_rotated = prev_df %>% 
      mutate(xy_dat = map2(data_to_rotate, theta, rotate, 
                           var_ranges = var_ranges)) %>%
      mutate(xy_dat = map(xy_dat, 1)) %>% 
      select(-data_to_rotate) %>% 
      unnest(c(xy_dat)) %>% 
      mutate(id = 1:n())
  }
  
  rotated_mini_axes = rotated_dat %>% 
    select(-xy_dat, -mini_labels) %>% 
    unnest(c(mini_axes))
  green_axes = green_rotated %>% 
    select(-xy_dat, -mini_labels) %>% 
    unnest(c(mini_axes))
  
  rotated_labels = rotated_dat %>% 
    select(-xy_dat, -mini_axes) %>% 
    unnest(c(mini_labels)) %>% 
    mutate(label = rep(c('alpha', 'beta', 'sigma'),
                       evo_frames))
  green_labels = green_rotated %>% 
    select(-xy_dat, -mini_axes) %>% 
    unnest(c(mini_labels)) %>% 
    mutate(label = rep(c('alpha', 'beta', 'sigma'),
                       green_pause))
  
  return(list(rotated_points, rotated_mini_axes, rotated_labels,
              green_points, green_axes, green_labels, prev_rotated))
  
}

plot_multi = function(nuts_df, frames_per_iter = 20, green_pause = 3, theta_max = 25,
                      height = 400, width = 600){
  
  n_green = green_pause * nrow(nuts_df)
  
  n_frames = nrow(nuts_df) * frames_per_iter # 160
  
  theta_vec = seq(0, theta_max, length.out = n_frames)
  
  n_draw = frames_per_iter - green_pause
  
  full_space = nuts_df %>% 
    mutate(C = map(C, as_tibble),
           iter = 1:n()) %>%
    unnest(c(C)) %>%  
    set_names(c('x', 'y', 'z', 'momx', 'momy', 'momz', 'j', 'chosen', 'iter'))
  
  nuts_df$frame_range = split(1:n_frames, ceiling(seq_along(1:n_frames)/frames_per_iter))
  nuts_df$theta_range = split(theta_vec, ceiling(seq_along(1:n_frames) / frames_per_iter))
  nuts_df$selected_point = map2(nuts_df$C, nuts_df$chosen,
                                ~.x[.y,])
  nuts_df$previous = map(0:(nrow(nuts_df) -1 ), # jank af
                         ~ bind_rows(head(nuts_df, n = .x)$selected_point))
  
  var_ranges = full_space %>% 
    select(x:z) %>% 
    summarise_all(.funs = ~list(range(.x)))
  
  
  all_rotations = pmap(list(nuts_df$C, nuts_df$frame_range, nuts_df$theta_range, nuts_df$chosen, previous_points = nuts_df$previous),
                       plot_one, 
                       var_ranges = var_ranges, green_pause = green_pause) %>% 
    purrr::transpose()
  
  dot_dat = all_rotations[[1]] %>%
    bind_rows %>% 
    bind_rows(all_rotations[[4]] %>% bind_rows)
  
  prev_dots = all_rotations[[7]] %>% 
    bind_rows
  
  rotated_mini_axes = all_rotations[[2]] %>% 
    bind_rows %>% 
    bind_rows(all_rotations[[5]] %>% bind_rows)
  
  rotated_labels = all_rotations[[3]] %>% 
    bind_rows %>% 
    bind_rows(all_rotations[[6]] %>% bind_rows)
  
  tree_plot = dot_dat %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(color = j),
               size = 3,
               alpha = .5) +
    geom_point(data = dot_dat %>% filter(selected_point, green_hold),
               color = 'green',
               size = 10) +
    geom_point(data = prev_dots,
               color = 'green', 
               size = 6,
               alpha = .5) + 
    scale_color_viridis_c(end = .92,
                          breaks = unique(dot_dat$j)) + 
    geom_segment(data = rotated_mini_axes,
                 aes(x = x, 
                     xend = xend,
                     y = y, 
                     yend = yend)) + 
    geom_text(data = rotated_labels,
              aes(x = x, y = y, label = label)) + 
    transition_manual(global_frame) +
    theme_bw() + 
    labs(color = 'tree depth') + 
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
  
  animate(tree_plot,
          width = width,
          height = height, nframes = n_frames, renderer = ffmpeg_renderer())
  
  
}

#### Make some animations ----

# nuts_df = nuts_df[3:10,]
plot_multi(nuts_df[3:10,], green_pause = 4, frames_per_iter = 25, theta_max = 800, width = 1920, height = 1080)

