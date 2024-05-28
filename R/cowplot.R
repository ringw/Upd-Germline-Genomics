replace_legend <- function(myggplot, mygetlegend) {
  myggplot <- myggplot %>% as_gtable
  myggplot$grobs[[match("guide-box", plot_component_names(myggplot))]] <- mygetlegend
  myggplot
}