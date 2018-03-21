library(tidyverse)
library(readr)
library(magrittr)
library(igraph)
library(network)
library(GGally)
library(ggraph)
library(lubridate)
library(reshape2)
library(ggridges)
library(migest)
library(scales)

file <- "./data/2014_mrt_hourly_daily_ridership.csv"

read_csv(file) %>%
  map_df(~ifelse(.=="-", NA, .)) %>%
  mutate_if(grepl("_entry$|_exit", names(.)), as.integer) -> mrt


#### DATA CLEANING ####

# add a date element
year <- file %>% str_extract("[0-9]{4}")
mrt <- mrt %>%
  mutate(time = str_extract(time, "^[0-9]{2}:[0-9]{2}"),
         date = mdy(paste(month, day, year)),
         datetime = mdy_hm(paste(month, day, year, time)),
         dayofweek = wday(date, abbr = F, label = T), 
         hour = hour(datetime))

mrt_long <- mrt %>%
  select_at(vars(date, datetime, dayofweek, hour, north_avenue_entry:taft_exit)) %>%
  melt(id = c("date", "datetime", "dayofweek", "hour")) %>%
  mutate(location = str_replace(variable, "_entry|_exit$", ""), 
         mode = str_extract(variable, "(?<=_)entry|exit$")) %>%
  select(-variable)
  
#### EDA ####

# dates with no rides at all (caused by systematic errors probably)
mrt_long %>%
  group_by(date) %>%
  summarise(sum = sum(value)) %>%
  filter(sum == 0) %$%
  paste(date, collapse = "\n") %>%
  cat("These dates have no rides at all: \n", .)

# average utilization by hour, by station
station_names <- c("north_avenue", "quezon_avenue", "gma_kamuning",
                   "cubao", "santolan", "ortigas", "shaw_blvd",
                   "boni_avenue", "guadalupe", "buendia", "ayala_avenue",
                   "magallanes", "taft")
station_labels <- c("North \nAve", "Quezon \nAve", "GMA \nKamuning",
                   "Cubao", "Santolan", "Ortigas", "Shaw \nBlvd",
                   "Boni \nAve", "Guadalupe", "Buendia", "Ayala \nAve",
                   "Magallanes", "Taft \nAve")
util <- mrt_long %>%
  group_by(hour, location, mode) %>%
  summarise(utilization = mean(value, na.rm = T)) %>%
  ungroup() %>%
  mutate(utilization = ifelse(mode == "exit", -utilization, utilization),
         location = factor(location, 
                           levels = station_names,
                           labels = station_labels, 
                           ordered = T)) 
  
  ggplot(util) + 
  geom_bar(aes(x=1, y = utilization, fill = mode), stat = "identity") +
  facet_grid(location ~ hour) + 
  coord_flip() + 
  theme_void() +
  theme(legend.position = "bottom")
  
  # mirror-image in utilization
  pdf("./output/mirror-image.pdf")
  util %>%
    filter(hour >= 4) %>%
  ggplot() + 
    geom_line(aes(x=hour, y = abs(utilization)), stat = "identity") +
    facet_grid(location~mode) +
    theme_void() +
    theme(legend.position = "bottom")
  dev.off()

# changes in weekly demand
weekly_usage <- mrt_long %>%
  mutate(week = week(datetime)) %>%
  group_by(week, location, mode) %>%
  summarise(value = mean(value, na.rm = T)) %>%
  arrange(week) %>%
  mutate(index = 100 * (value/first(value))) %>%
  ungroup()

ggplot(weekly_usage) +
  geom_line(aes(x = week, y = value, linetype = mode, color = location)) +
  facet_wrap(~location) +
  theme(legend.position = "bottom")
  
# notice the difference between entry and exit for the terminal stations: North Ave and Taft Ave
# does this suggest that there are users who use the mrt to go from NA/TA to some other station
# but not take the train back?


#### OD Matrix Analysis ####

#### HELPER FUNCTIONS ####
equalize_sum <- function(vector1, vector2) {
  is_vector1_greater <- sum(vector1) >= sum(vector2) 
  diff <- ifelse(
                   is_vector1_greater, 
                   sum(vector1) - sum(vector2),
                   sum(vector2) - sum(vector1)
                 )
  print(sprintf("THE GAP BETWEEN ENTRANCE AND EXIT VOLUME IS: %s", diff))
  print(sprintf("ERROR RATE: %s", ifelse(is_vector1_greater, diff/sum(vector2), diff/sum(vector1))))
  prop <- ifelse(
                  rep(is_vector1_greater, length(vector1)), 
                  vector1 / sum(vector1),
                  vector2 / sum(vector2)
                )
  
  if (is_vector1_greater) {
    # increase counts of 2
    df <- tibble(vector = vector2, prop = prop*diff) %>%
      mutate(int = floor(prop), 
             resid = prop - int,
             index = 1:n())
  } else {
    # increase counts of 1
    df <- tibble(vector = vector1, prop = prop*diff) %>%
      mutate(int = floor(prop),
             resid = prop - int,
             index = 1:n())
  }
  # counts to allocate
  allocate <- diff - sum(df$int)
  df %>%
    arrange(desc(prop)) %>%
    mutate(newvector = vector + int + c(rep(1, allocate), rep(0, n()-allocate))) %>%
    arrange(index) %$%
    newvector -> newvector
  
  if (is_vector1_greater) {
    return(list(vector1, newvector))
  } else {
    return(list(newvector, vector2))
  }
}

# force that entries and exit within the same station is not allowed
restriction_mat <- matrix(1, 13, 13)
for (i in 1:13) {
  restriction_mat[i,i] <- 0
}


mintime <- 12
maxtime <- 15
selectdate <- "08-22-2014"
flowplots <- function(mrt, mintime, maxtime, selectdate) {
  mrt %>%
    filter(hour >= mintime & hour <= maxtime & date == mdy(selectdate)) %>%
    select_at(vars(north_avenue_entry:taft_exit)) %>% 
    colSums(na.rm = T) %>%
    tibble(location_mode = names(.), value = .) %>%
    mutate(location = str_replace(location_mode, "_exit|_entry", ""),
           location = factor(location, levels = station_names, labels = station_labels, ordered = T),
           mode = str_extract(location_mode, "(?<=_)exit|entry")) %>%
    select(location, mode, value) %>%
    spread(key = "mode", value = "value") %$%
    equalize_sum(entry, exit) -> net_counts
  
  
  mat <- net_counts %>%
  {cm2(.[[1]], .[[2]], m = restriction_mat)}
  
  mat <- mat[[1]] 
  rownames(mat) <- station_names
  colnames(mat) <- station_names
  # # mat <- scale(mat, center = F)
  g <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = T)
  
  edgelist <- get.edgelist(g) %>%
    as.tibble() %>%
    rename(from = V1, to = V2) %>%
    mutate(from = factor(from, levels = station_names, labels = station_labels, ordered = T),
           to = factor(to, levels = station_names, labels = station_labels, ordered = T),
           fromindex = as.numeric(from),
           toindex = as.numeric(to)) %>%
    mutate(weight = E(g)$weight,
           direction = ifelse(fromindex < toindex, "SB", "NB"),
           directionindex = ifelse(direction == "NB", 1, 2),
           color = ifelse(direction == "NB", "red", "black")) %>%
    arrange(desc(direction), from, to) %>%
    group_by(from, direction) %>%
    mutate(curvature = seq(from = 0.5, by = 0.2, length.out = n()))
  
  g <- graph_from_data_frame(edgelist)
  E(g)$color = edgelist$color
  E(g)$curvature = edgelist$curvature
  E(g)$weight = edgelist$weight
  # plot(g,
  #      vertex.size = 0.1,
  #      vertex.label = NA,
  #      edge.width = E(g)$weight*2,
  #      edge.arrow.size = 0.01,
  #      edge.color = E(g)$color,
  #      layout = matrix(c(rep(0, 13), (1:13)*10), byrow = F, ncol = 2),
  #      edge.curved=E(g)$curvature)
  
  # flow plot
  p1 <- ggraph(g, layout = "linear") +
    geom_segment(aes(y = -4, yend = 0, x = 1:13, xend = 1:13), alpha = 0.1) +
    geom_edge_arc(aes(width = weight, color = direction, alpha = weight), curvature = 1) +
    geom_text(aes(y = -4, x = 1:13, label = station_labels), hjust = 0) +
    scale_edge_width(range = c(0.5, 2.5), limits = c(0, 5000), breaks = c(0, 2500, 5000)) +
    scale_edge_alpha(guide = F) +
    scale_edge_color_discrete(guide = F) + 
    # geom_node_text(aes(label = name)) +
    theme_graph() +
    scale_x_reverse() +
    coord_flip() +
    theme(legend.position = "none") +
    labs(edge_width = "Volume")

  # heatmap
  rownames(mat) <- station_labels
  colnames(mat) <- station_labels
  heatmap_dat <- mat %>%
    apply(1, function(x) {x/sum(x)}) %>%
    t() %>%
    melt() %>%
    rename(from = Var1, to = Var2) %>%
    mutate(indexto = rep(1:13, each = 13),
           indexfrom = rep(1:13, 13),
           from = factor(from, levels = rev(station_labels), ordered = T)) %>%
    mutate(isSB = ifelse(indexfrom == indexto, NA, indexfrom < indexto),
           value = ifelse(value == 0, NA, value))
  
  p2 <- ggplot(heatmap_dat) +
    geom_tile(aes(to, from, fill = value), alpha = 0.8, color = "white", size = 2) +
    theme_classic() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 0, colour = "grey50"),
          axis.title.x = element_text(hjust = 1, margin = margin(t = -15, unit = "pt"), size = 12),
          axis.title.y = element_text(hjust = 1, margin = margin(b = 50, unit = "pt"), size = 12)) +
    scale_x_discrete(position = "top") +
    scale_fill_continuous(low = "grey90", high = "steelblue", na.value = "black", breaks = 0.05*0:5, labels = percent, limits = c(0, 0.25)) +
    scale_color_discrete(na.value = "white") +
    guides(alpha = FALSE) +
    coord_fixed() +
    labs(fill = "Percent traffic \nfrom origin \nto destination",
         x = "Destination",
         y = "Origin")

  mintime_label <- paste0(str_pad(mintime, width = 2, side = "left", pad = "0"), "00")
  maxtime_label <- paste0(str_pad(maxtime, width = 2, side = "left", pad = "0"), "00")
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(sprintf("Passenger flow between %s and %s", mintime_label, maxtime_label), fontface='bold', hjust = 1, size = 16)
  cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(0.3, 0.8), scale = c(1,1)) %>%
    {cowplot::plot_grid(title, ., ncol = 1, rel_heights = c(0.1, 1), rel_widths = c(1, 1))}
}


pdf("./output/combinedplots.pdf", width = 11, height = 6)  
selectdate <- "08-22-2014"
for (i in 5:22) {
  p <- flowplots(mrt, i, i + 1, selectdate)
  print(p)
}
dev.off()
  

#### NET UTILIZATION ####

mrt_long %>%
  spread(key = "mode", value = "value") %>%
  mutate_at(vars(entry:exit), .funs = function(x) ifelse(is.na(x), 0, x)) %>%
  group_by(hour, location) %>%
  summarise(entry = mean(entry, na.rm = T) %>% floor,
            exit = mean(exit, na.rm = T) %>% floor) %>%
  group_by(hour, location) %>%
  mutate(cum_entry = cumsum(entry),
         cum_exit = cumsum(exit)) -> net_util

net_util %>%
  select(hour, location, cum_entry, cum_exit) %>%
  group_by(location) %>% 
  nest() -> net_util_nest

#### imputation of 0.1 hour rates
spline_wrapper <- function(df) {
  Espline_entry <- splinefun(df$hour, df$cum_entry, method = "monoH.FC")
  Espline_exit <- splinefun(df$hour, df$cum_exit, method = "monoH.FC")
  
  tibble(hour = seq(0, 23, 0.1)) %>%
    mutate(entry = Espline_entry(hour),
           exit = Espline_exit(hour)) %>%
  return()
}

net_util_nest %>%
  mutate(impute = map(data, spline_wrapper)) %>%
  unnest(impute) %>%
  group_by(location) %>%
  mutate(hour = seq(0, 23, 0.1)) -> impute_util

impute_util_long <- gather(impute_util, mode, volume, -(location:hour))
net_util_long <- select(net_util, location, hour, entry, exit) %>%
  gather(mode, volume, -(location:hour))


ggplot() + 
  geom_line(data = impute_util_long, aes(hour, volume, col = location)) +
  geom_line(data = net_util_long, aes(hour, volume, col = location), alpha = 0.5) +
  facet_wrap(~mode)



ts(net_util$cum_entry, frequency = 24, start = 0)
Espline <- splinefun(x = net_util$hour, y = net_util$cum_entry, method = "hyman")
dEdt_spline <- function(t) Espline( t , deriv = 1 )

input <- seq(0, 24, 0.1)
output <- input %>% Espline()

imputed <- tibble(hour = input,
                  cum_entry = output) %>%
  mutate(entry = c(cum_entry[1], diff(cum_entry)))

ggplot() +
  geom_line(aes(input, output), alpha = 0.01) +
  geom_point(aes(input, output)) +
  geom_line(data = net_util, aes(hour, cum_entry))



#### SAND BOX ####



  xyz2<-mat %>%
    apply(2, function(x) {x/sum(x)}) %>%
    melt() %>%
    rename(from = Var1, to = Var2) %>%
    mutate(indexto = rep(1:13, each = 13),
           indexfrom = rep(1:13, 13),
           from = factor(from, levels = rev(station_labels), ordered = T)) %>%
    mutate(isSB = ifelse(indexfrom == indexto, NA, indexfrom < indexto),
           value = ifelse(value == 0, NA, value))
  
  ggplot(xyz2) +
    geom_tile(aes(from, to, fill = value, alpha = value), color = "white", size = 2) +
    theme_classic() +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, colour = "grey50")) +
    scale_y_discrete(position = "top") +
    scale_fill_continuous(na.value = "grey10") +
    scale_color_discrete(na.value = "white") 
  
  
  
  
  

net <- as.network.matrix(mat, "adjacency", directed = T, bipartite = F, ignore.eval = F, names.eval = "weight")
get.edge.attribute(net, "weight")
ggnet2(net, edge.size = "weight")


net <- as.network.matrix(edgelist, "edgelist", directed = T, bipartite = F, ignore.eval = F)
get.edge.attribute(net, "direction")
set.edge.attribute(net, "color", ifelse(net %e% "direction" > "NB", "black", "grey75"))
ggnet2(net, edge.size = "weight", edge.color = "color") 
