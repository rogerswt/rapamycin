#
# Analysis of Derek's rapamycin data set
#
# 2017_05_09  WTR
#


source("~/git/R/tools/sourceTools.R")
source("~/git/R/cytomics_workshop/Sessions_workshop_2/rapamycin/rapa_utils.R")

# where did I put the data?
data_base = "~/git/R/cytomics_workshop/Sessions_workshop_2/workshop_data/jonesder_20150617_NZB_GC/"
pic_base = paste(data_base, "results/", sep = "")

# get the file list and class assignments
tab = read.csv(file = paste(data_base, "data_files.csv", sep = ""))

idx_control    = which(tab$class == "control")
idx_rapamycin  = which(tab$class == "rapamycin")
idx_experiment = c(idx_control, idx_rapamycin)

for (i in 1:length(idx_experiment)) {
  e_num = idx_experiment[i]
  fname = tab$filename[e_num]
  # read in the data
  ff = Get_file(fqp = paste(data_base, fname, sep = ""))

  # ditch debris
  ff_no_debris = Gate_debris(ff = ff, show = TRUE)
  pic_file = paste(pic_base,
                   tab$type[e_num], "_",
                   tab$class[e_num], "_", e_num, ".png",
                   sep = "")
  dev.print(png, pic_file, width = 600, height = 600)
}

