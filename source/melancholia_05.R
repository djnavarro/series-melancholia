
library(tidyverse)

Rcpp::sourceCpp(here::here("source", "stoneskip.cpp"))

fractals <- list(
  ambient::billow,
  ambient::fbm,
  ambient::ridged
)

generators <- list(
  ambient::gen_simplex,
  ambient::gen_worley
)

sample_list <- function(...) {
  (sample(...))[[1]]
}

skip_matrix <- function(seed = 1, grains = 100, shades = 1000, iterations = sample(x = 50:200, size = 1)) {

  set.seed(seed)

  # create the base image using ambient -------------------------------------

  #cat("making base image...\n")
  grid <- ambient::long_grid(
    x = seq(0, 1, length.out = grains),
    y = seq(0, 1, length.out = grains)
  )

  base <- ambient::fracture(
    noise     = sample_list(x = generators, size = 1),
    fractal   = sample_list(x = fractals, size = 1),
    octaves   = sample(x = 1:10, size = 1),
    frequency = sample(x = 1:10, size = 1),
    value     = "distance2",
    seed      = seed,
    x         = grid$x,
    y         = grid$y
  )

  # run cellular automaton over base image ----------------------------------

  #cat("running stepping stone automaton...\n")
  input <- round(ambient::normalise(base, to = c(1, shades)))
  input <- matrix(as.integer(input), grains, grains)
  output <- skip_stone(input, iterations)
  return(output)

}


construct_tibble <- function(mat) {

  #cat("converting to tibble...\n")

  # use the row and column names to represent co-ordinates
  rownames(mat) <- paste0("y", nrow(mat):1) # <- flip y
  colnames(mat) <- paste0("x", 1:ncol(mat))

  # convert to tibble
  tbl <- mat %>%
    as.data.frame() %>%
    rownames_to_column("y") %>%
    as_tibble()

  # reshape
  tbl <- tbl %>%
    pivot_longer(
      cols = starts_with("x"),
      names_to = "x",
      values_to = "shade"
    )

  # tidy
  tbl <- tbl %>%
    arrange(x, y) %>%
    mutate(
      x = x %>% str_remove_all("x") %>% as.numeric(),
      y = y %>% str_remove_all("y") %>% as.numeric(),
      id = row_number()
    )

  return(tbl)
}


distort <- function(tbl, its = 70) {

  #cat("adding distortion...\n")

  tbl <- tbl %>%
    mutate(
      id = 1,
      seed = sample(1000, 1),
      ind = row_number(),
      type = "stoneskip"
    )

  dat <- tbl %>%
    #   jasmines::unfold_warp(scale = 20, iterations = 3) %>%
    mutate(x = x / max(x), y = y / max(y)) %>%
    #  slice_sample(prop = .1) %>%
    #   jasmines::unfold_breeze(iterations = its, drift = 0, scale = .001) %>%
    mutate(size = (shade %% 10)) %>%
    #    filter(time > 20) %>%
    identity()
  return(dat)
}

sample_shades <- function() {
  pal <- sample(colorir::colores$palette_name, 1)
  shades <- colorir::colores$colour[colorir::colores$palette_name == pal]
  return(shades)
}

blend_shades <- function(x, y, p = .5) {
  x <- col2rgb(x)
  y <- col2rgb(y)
  z <- round(p*x + (1-p)*y)
  z <- rgb(red = z[1, ]/255, green = z[2, ]/255, blue = z[3, ]/255)
  return(z)
}

draw <- function(dat) {

  #cat("specifying plot...\n")
  its <- max(dat$time)
  shades <- sample_shades()
  #bg <- blend_shades(shades[1], "black", p = .1)
  bg <- "white"
  brd <- .1

  pic <- ggplot(dat) +
    geom_tile(
      mapping = aes(
        x = x,
        y = y,
        fill = shade
      ),
      show.legend = FALSE
    ) +
    scale_x_continuous(name = NULL, expand = c(0,0), breaks = NULL) +
    scale_y_continuous(name = NULL, expand = c(0,0), breaks = NULL) +
    scale_size_identity() +
    scale_colour_gradientn(colours = shades) +
    scale_fill_gradientn(colours = shades) +
    theme_void() +
    theme(plot.background = element_rect(fill = bg)) +
    NULL

  return(pic)
}


render <- function(pic, seed) {
  inches <- 40/3
  for(pixels in c(1000, 4000, 8000, 16000)) {
    cat(pixels, " ")
    output <- here::here(
      "image", "sys_05", pixels, paste0("melancholia_sys05_seed", seed, ".png")
    )
    ggsave(
      filename = output,
      plot = pic,
      width = inches,
      height = inches,
      dpi = pixels / inches
    )
  }
  cat("\n")
}



box <- function(x_shift, y_shift, shade_shift, seed, grains = 20) {
  skip_matrix(seed = seed, grains = grains) %>%
    construct_tibble() %>%
    distort() %>%
    mutate(
      x = x + x_shift,
      y = y + y_shift,
      shade = shade + shade_shift
    )
}

ms_paint <- function(seed) {

  set.seed(seed)
  grains <- 25

  global <- subdivide_matrix(seed = seed, grains = grains, iterations = 7) %>%
    construct_tibble() %>%
    transmute(
      x_shift = x - 1,
      y_shift = y - 1,
      shade_shift = shade * 3,
      seed = sample(1000, n())
    )

  dat <- pmap_dfr(global, box, grains = grains) %>%
    mutate(
      x = x / max(x),
      y = y / max(y),
      time = 10,
      size = 1
    )

  dat %>%
    draw() %>%
    render(seed)
}

subdivide_matrix <- function(seed, grains, iterations) {

  set.seed(seed)

  split_rectangles <- function(left, right, bottom, top) {
    vertically <- runif(1) < .5

    if(vertically == TRUE) {
      if(top - bottom < 2) {
        return(tibble(
          left = left,
          right = right,
          bottom = bottom,
          top = top
        ))
      }
      insert <- sample((bottom + 1):top, 1)
      return(tibble(
        left = c(left, left),
        right = c(right, right),
        bottom = c(bottom, insert),
        top = c(insert - 1, top)
      ))
    }

    if(vertically == FALSE) {
      if(right - left < 2) {
        return(tibble(
          left = left,
          right = right,
          bottom = bottom,
          top = top
        ))
      }
      insert <- sample((left + 1):right, 1)
      return(tibble(
        left = c(left, insert),
        right = c(insert - 1, right),
        bottom = c(bottom, bottom),
        top = c(top, top)
      ))
    }
  }

  tbl <- tibble(
    left = 1,
    right = grains,
    bottom = 1,
    top = grains
  )

  for(it in 1:iterations) {
    tbl <- pmap_dfr(tbl, split_rectangles)
  }
  tbl <- mutate(tbl, shade = sample(1000, n()))

  mat <- matrix(0, nrow = grains, ncol = grains)
  for(i in 1:nrow(tbl)) {
    with(tbl[i, ], mat[left:right, bottom:top] <<- shade) # OMG WHHHHHHHY this code is awful???????
  }
  return(mat)
}



# run it ------------------------------------------------------------------

for(seed in 1301:1600) {
  cat("seed:", seed, "\n")
  ms_paint(seed)
}




