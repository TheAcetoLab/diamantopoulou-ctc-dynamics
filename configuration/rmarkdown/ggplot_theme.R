require(ggplot2)
require(cowplot)
use_font_size <- 8

theme_set(
  theme_cowplot(font_family = "Helvetica", font_size = use_font_size, rel_small = use_font_size/use_font_size, rel_tiny = use_font_size/use_font_size, rel_large = use_font_size/use_font_size) +
    theme (
      axis.line = element_line(size = rel(0.5)),
      axis.ticks = element_line(size = rel(0.5))
      )
  )

theme_mimimal_ygrid <- theme(
      strip.text.x = element_text(size=use_font_size, hjust = 0),
      strip.background = element_blank(),
      panel.grid.major.y = element_line(size = (0.2), colour="grey")
      ) 

geom_text_size <- use_font_size * 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
