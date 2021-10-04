
# install.packages("Rcpp")
# install.packages("here")
# install.packages("paletteer")
# install.packages("ambient")

skip_image <- function(seed = 1) {

  set.seed(seed)
  
  # import the C++ code
  Rcpp::sourceCpp(here::here("stoneskip.cpp"))
  
  # where to save the file
  filename <- here::here(paste0("stoneskip_", seed, ".png"))
  
  # plot size, grid size and number of shades are fixed
  grains <- 1000
  pixels <- 5000
  shades <- 1000
  
  # palette is sampled randomly
  pal_names <- paletteer::palettes_c_names
  pal_index <- sample(nrow(pal_names), 1)
  pal_which <- paste0(
    pal_names$package[pal_index], "::", pal_names$palette[pal_index]
  )
  palette <- paletteer::paletteer_c(pal_which, n = shades)
  
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
  
  # create the base image using ambient -------------------------------------
  
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
  
  cat("running stepping stone automaton...\n")
  iterations <- sample(x = 50:200, size = 1)
  input <- round(ambient::normalise(base, to = c(1, shades)))
  input <- matrix(as.integer(input), grains, grains) 
  output <- skip_stone(input, iterations)
  
  
  
  # write to an image file --------------------------------------------------
  
  cat("rendering image...\n")
  rast <- as.raster(matrix(palette[output], grains, grains))
  png(filename, pixels, pixels)
  op <- par(mar=c(0,0,0,0))
  plot(rast)
  dev.off()
  par(op)
  
  
  # invisibly return the output matrix --------------------------------------
  return(invisible(output))
  
}

