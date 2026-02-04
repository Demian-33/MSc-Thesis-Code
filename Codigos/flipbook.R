# (1) densidad SN moviendose ----

# izquierda al centro
# centro: ariba y abajo y de regreso
# centro a izquierda
# de regreso
t0 <- seq(-3.5, 3.5, length=200)
lambda_neg <- seq(-3.0, 0, by=0.15)
lambda_pos <- seq(0, 3.0, by=0.15)
escala <- seq(1, 2, by=0.1)
k <- 1
x11()
setwd(r'(C:\Users\Lesau\Pictures\flipbook2)')

save_img <- function(){
  pdf(file = paste0('img-',k,'.pdf'))
  # png(file = paste0('img-',k,'.png'))
  par(mar = c(0, 0, 0, 0))
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  dev.off()
  k <<- k+1
}

for(izq in 1:length(lambda_neg)){
  y0 <- sn::dsn(x=t0, alpha=lambda_neg[izq], omega=1)
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}
for(cen in 1:length(escala)){
  y0 <- dnorm(x=t0, sd=escala[cen])
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}
for(cen in 1:length(escala)){
  y0 <- dnorm(x=t0, sd=rev(escala)[cen])
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}
for(der in 1:length(lambda_pos)){
  y0 <- sn::dsn(x=t0, alpha=lambda_pos[der], omega=1) # runif(1, 0.9, 1)
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}
# regreso
for(der in 1:length(lambda_pos)){
  y0 <- sn::dsn(x=t0, alpha=rev(lambda_pos)[der], omega=1) # runif(1, 0.9, 1)
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}
for(cen in 1:length(escala)){
  y0 <- dnorm(x=t0, sd=escala[cen])
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}
for(cen in 1:length(escala)){
  y0 <- dnorm(x=t0, sd=rev(escala)[cen])
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}
for(izq in 1:length(lambda_neg)){
  y0 <- sn::dsn(x=t0, alpha=rev(lambda_neg)[izq], omega=1)
  plot(t0, y0, type='l', axes=F, xlab='', ylab='')
  # Sys.sleep(0.1)
  save_img()
}

# (2) histograma ----
library(tibble)
library(ggplot2)
nombre <- file.path(r'(C:\Maestría_CP\Tesis\Documento-Tesis-SAOM\Documento-Tesis\Figuras\flipbook)')
nombre <- file.path(r'(C:\Users\Lesau\Pictures\flipbook)')
dir(nombre)
sprintf('%02d', 100)

t0 <- seq(-3.5, 3.5, length=200)
y0 <- sn::dsn(t0, alpha=1)
dftrue <- tibble(x0=t0, y0=y0)
n <- 200 # numero de paginas
muestra <- c()
x11()
set.seed(1996)
for(i in 1:n){
  # tomar siete muestras
  sam <- sn::rsn(n=7, alpha=1)
  muestra <- c(muestra, sam)
  dftmp <- tibble(muestra)
  g0 <- ggplot(dftmp, aes(x=muestra)) + 
    geom_histogram(aes(y=..density..), col='gray95', fill='gray70', linewidth=0.2) + 
    geom_line(dftrue, mapping=aes(x=x0, y=y0), inherit.aes = F, col='gray50', linewidth=0.1) + 
    theme_void()
  # print(g0)
  # nombre0 <- file.path(nombre, paste0('hist-',i ,'.pdf'), fsep = '\\')
  nombre0 <- file.path(nombre, paste0('hist-',sprintf('%03d', i) ,'.pdf'), fsep = '\\')
  ggsave(filename = nombre0, plot=g0, width=2.5, height=2.5, unit='cm')
  # Sys.sleep(0.1)
}

library(bayesplot)
library(dplyr)
library(ggplot2)

# (3) markov chain? ----
# ---- Parameters ----
total_frames <- 200          # total PNGs (pages)
frames_per_transition <- 3   # interpolation frames per state-to-state jump
nodes <- data.frame(         # small 4-node square (you can change)
  id = c("A","B","C"),
  x  = c(0, 1, 0),
  y  = c(0, 0, 1),
  stringsAsFactors = FALSE
)

# Transition matrix P (rows = from, cols = to). Rows must sum to 1.
P <- matrix(c(
  0.3, 0.5, 0.2,
  0.5, 0.2, 0.3,
  0.25, 0.25, 0.5
), nrow = nrow(nodes), byrow = TRUE)
rownames(P) <- nodes$id
colnames(P) <- nodes$id

# ---- simulate a chain of enough transitions ----
# number of transitions needed (each transition produces frames_per_transition frames)
n_transitions <- ceiling(total_frames / frames_per_transition)
set.seed(1)  # for reproducibility; remove for different random draws

# simulate states (length = n_transitions + 1)
states <- character(n_transitions + 1)
states[1] <- sample(nodes$id, 1)  # random start
for (t in 1:n_transitions) {
  from_idx <- which(nodes$id == states[t])
  states[t+1] <- sample(nodes$id, 1, prob = P[from_idx, ])
}
states
# ---- build interpolated frame positions ----
interp_list <- list()
frame_idx <- 1L
for (t in seq_len(n_transitions)) {
  from <- states[t]
  to   <- states[t+1]
  from_xy <- nodes %>% filter(id == from) %>% select(x, y) %>% as.numeric()
  to_xy   <- nodes %>% filter(id == to)   %>% select(x, y) %>% as.numeric()
  # sequence of fractions between 0 and 1; include endpoint only on last transition
  fracs <- seq(0, 1, length.out = frames_per_transition + 1)[1:frames_per_transition]
  for (f in fracs) {
    x <- (1 - f) * from_xy[1] + f * to_xy[1]
    y <- (1 - f) * from_xy[2] + f * to_xy[2]
    interp_list[[frame_idx]] <- data.frame(frame = frame_idx, x = x, y = y, state_from = from, state_to = to, frac = f)
    frame_idx <- frame_idx + 1L
    if (frame_idx > total_frames) break
  }
  if (frame_idx > total_frames) break
}
frames_df <- bind_rows(interp_list)

# Optionally compute small "trail" data: for each frame, keep a few previous positions
trail_length <- 3  # how many previous points to show (for fading tail)
trail_df <- do.call(rbind, lapply(seq_len(nrow(frames_df)), function(i) {
  idxs <- pmax(1, (i - (trail_length-1))):i
  data.frame(frame = i,
             trail_index = seq_along(idxs),
             x = frames_df$x[idxs],
             y = frames_df$y[idxs],
             stringsAsFactors = FALSE)
}))

# ---- prepare edges to draw (edges with P>0) ----
edges <- expand.grid(from = nodes$id, to = nodes$id, stringsAsFactors = FALSE) %>%
  left_join(as.data.frame(P) %>% mutate(from = nodes$id), by = "from") %>%
  tidyr::gather(key = "to", value = "prob", -from) %>%
  filter(prob > 0 & from != to) %>%
  left_join(nodes, by = c("from" = "id")) %>%
  rename(x = x, y = y) %>%
  left_join(nodes, by = c("to" = "id")) %>%
  rename(xend = x.y, yend = y.y, x = x.x, y = y.x) %>%
  select(from, to, prob, x, y, xend, yend)

# fix plotting limits (padding)
xmin <- min(nodes$x) - 0.3
xmax <- max(nodes$x) + 0.3
ymin <- min(nodes$y) - 0.3
ymax <- max(nodes$y) + 0.3

# ---- create and save frames ----
# dir.create("frames", showWarnings = FALSE)
x11()
for (i in seq_len(nrow(frames_df))) {
  pos <- frames_df[i, ]
  trails_i <- trail_df %>% filter(frame == i)
  p <- ggplot() +
    # edges
    geom_segment()
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
                 size = 0.4, alpha = 0.25, lineend = "round") +
    # nodes
    geom_point(data = nodes, aes(x = x, y = y), size = 2.5) +
    geom_text(data = nodes, aes(x = x, y = y, label = id), vjust = -1.4, size = 2.5) +
    # trailing small dots (fading)
    geom_point(data = trails_i, aes(x = x, y = y, alpha = rev(seq_len(nrow(trails_i)))), size = 1.5) +
    scale_alpha(range = c(0.15, 0.8), guide = "none") +
    # moving particle
    geom_point(aes(x = pos$x, y = pos$y), color = "black", size = 3.5) +
    coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    theme_void() +
    theme(
      plot.margin = grid::unit(c(0,0,0,0), "mm")
    )
  print(p)
  Sys.sleep(0.1)
  # filename <- sprintf("frames/frame_%02d.png", i)
  # Save small square suitable for footnote (adjust width/height and dpi for print)
  # ggsave(filename, p, width = 2.2, height = 2.2, units = "cm", dpi = 600)
}

cat("Saved", nrow(frames_df), "frames in ./frames/\n")

# ---- optional: combine into an animated GIF preview (requires magick) ----
if (requireNamespace("magick", quietly = TRUE)) {
  imgs <- list.files("frames", pattern = "frame_.*png$", full.names = TRUE)
  gif <- magick::image_read(imgs) %>%
    magick::image_scale("300x300") %>%
    magick::image_animate(fps = 12)  # preview speed; flipbooks ~10-16 fps
  magick::image_write(gif, "frames/preview.gif")
  cat("Wrote preview GIF: frames/preview.gif\n")
} else {
  cat("Install 'magick' to create a GIF preview (optional).\n")
}

