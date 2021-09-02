

# Function for binding two data frames with different columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Original code extracted from: https://github.com/mcruf/LGNB/blob/master/R/utilities.R

mybind <- function(x, y) {
  i <- intersect(names(x), names(y))
  ix <- setdiff(names(x), names(y))
  iy <- setdiff(names(y), names(x))
  z <- rbind(x[i], y[i])
  xna <- data.frame(lapply(iy, function(...)rep(NA, nrow(x))))
  names(xna) <- iy
  yna <- data.frame(lapply(ix, function(...)rep(NA, nrow(y))))
  names(yna) <- ix
  zx <- rbind(x[ix], yna)
  zy <- rbind(xna, y[iy])
  cbind(z, zx, zy)
}