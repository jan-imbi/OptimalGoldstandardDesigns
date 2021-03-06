library(hexSticker)
library(ggplot2)
imgurl <- here::here("sticker/screenshot.tiff")
sticker <- sticker(imgurl,
                   package = "OptimalGoldstandardDesigns", p_size = 8 * 6000 / 300,
                   s_x = 1, s_y = .85, s_width = .85,
                   p_color = "black",
                   h_fill = "white",
                   h_color = "#cdcdff",
                   filename = "man/figures/sticker.png",
                   dpi = 6000
)
