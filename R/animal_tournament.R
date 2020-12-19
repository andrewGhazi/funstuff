library(tidyverse)
library(rvest)
library(beepr)
library(imager)
library(janitor)


page = read_html('https://gifts.worldwildlife.org/gift-center/gifts/Species-Adoptions.aspx') 

get_url = function(i){
  n_str = str_pad(as.character(i), width = 2, pad = '0')
  css = paste0('#ctl00_MainContentHolder_UncategorizedRepeaterVariant1_ctl', n_str, '_recordimageimg')
  img_loc = page %>% html_node(css) %>% html_attr('src')
  url_base = 'https://gifts.worldwildlife.org'
  return(paste0(url_base, img_loc))
}

animals = page %>% 
  html_node('#myUL') %>% 
  html_text %>% 
  str_split(pattern = '[\n]+') %>% 
  .[[1]] %>% 
  tibble(animal = .) %>% 
  filter(!(animal %in% c('', ' '))) %>% 
  mutate(i = 0:139) %>% 
  mutate(img_url = map_chr(i, get_url)) %>% 
  mutate(animal = gsub('NEW! ', '', animal))

wwf_dir = paste0(Sys.getenv('HOME'), '/Pictures/wwf_animals')
if (!dir.exists(wwf_dir)) dir.create(wwf_dir)

dl_photo = function(url, animal){
  Sys.sleep(.25)
  dest = paste0(wwf_dir, '/', janitor::make_clean_names(animal), '.jpg')
  download.file(destfile = dest,
                url = url, 
                quiet = TRUE)
  return(dest)
}

get_dest = function(animal){
  dest = paste0(wwf_dir, '/', janitor::make_clean_names(animal), '.jpg')

  return(dest)
}

if (file.exists(paste0(wwf_dir, '/animals.RData'))) {
  load(paste0(wwf_dir, '/animals.RData'))
  animals$img = map(animals$animal,
                     get_dest)
} else {
  message('Downloading photos...')
  animals$img = map2(animals$img_url, animals$animal,
                     dl_photo)
  save(animals,
       file = paste0(wwf_dir, '/animals.RData'))
}

shuf = animals %>% 
  sample_frac(1) %>% 
  select(-img_url)

while (nrow(shuf) >= 2){
  shuf$group = ceiling(seq(1, (nrow(shuf)), length.out = nrow(shuf)) / 2)
  
  shuf$winner = FALSE
  for (i in 1:max(shuf$group)){
    rows = which(shuf$group == i)
    if (length(rows) == 1) {
      message(paste0('Uneven bracket: ', shuf$animal[rows[1]], ' progresses to next round automatically...'))
      shuf$winner[rows] = TRUE
      break
    }
    dev.new()
    par(mfrow = c(1,2))
    plot(load.image(shuf$img[[rows[1]]]), main = paste0('A: ', shuf$animal[rows[1]]), axes = FALSE)
    plot(load.image(shuf$img[[rows[2]]]), main = paste0('B: ',shuf$animal[rows[2]]), axes = FALSE)
    winner = readline('A or B?')
    if (winner %in% c('A', 'a')){
      shuf$winner[rows[1]] = TRUE
    } else{
      shuf$winner[rows[2]] = TRUE
    }
    dev.off()
  }
  
  shuf = shuf %>% filter(winner) %>% sample_frac(1)
  if (nrow(shuf) == 1){
    beepr::beep(3)
    message(paste0('Winner: ', shuf$animal[1], '!'))
    dev.new()
    plot(load.image(shuf$img[[rows[1]]]), main = paste0('Winner: ',shuf$animal[rows[1]], '!'), axes = FALSE)
  }
}
# dev.off()



