# function to plot hist in ggplot
hist_ggplot <- function(x, xlabel, breaks=20) {
  # histogram data from base r
  dat_hist <- hist(x, breaks = breaks, plot = F)
  # get histogram data in dataframe
  dat_plot <- data.frame(
    mid = dat_hist$mids,
    count = dat_hist$counts)
  # smallest unit to use for width of bars
  smallest_unit <- dat_hist$mids[2] - dat_hist$mids[1]
  # ggplot
  ggplot(dat_plot, aes(x=mid, y=count)) + 
    geom_bar(width=smallest_unit, stat="identity") + 
    theme_bw() + 
    xlab(xlabel) + 
    ylab("Count")
}