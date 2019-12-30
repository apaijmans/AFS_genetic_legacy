#------------------------------------------------------------
# Anneke Paijmans
# last edited: Mar 2019
#
# Headland 2009 data - Make figure for sealing effort
#                      using data from Headland 2009
#
#------------------------------------------------------------

# Load libraries
library(readxl)
library(ggplot2)
library(extrafont)


#~~ Load data

# To load all sheets in a workbook, use lapply
path <- paste("Docs/Historical info/Scans Headland 2009/Headland histrograms.xlsx")
files <- lapply(excel_sheets(path), read_excel, path = path, range=cell_cols("A:C"))

files2 <- lapply(files, function(x) x[-1,]) # in excel I added a first and last entry to each sheet to make nicer graphs. here I want to remove them again, starting with removing the first rows in each dataframe in the last

files3 <- lapply(files2, function(x) x[-nrow(x),]) # removed last row of each dataframe in the list, based on the nrows of the specific df

# Make list into 1 dataframe
t <- data.table::rbindlist(files3, idcol=T)


#~~ Make figure

# replace id with sealing ground
colnames(t)[1] <- "sealing.ground"

t$sealing.ground <- gsub(1, "South Georgia & South Sandwich", t$sealing.ground)
t$sealing.ground <- gsub(2, "South Shetlands & South Orkneys", t$sealing.ground)
t$sealing.ground <- gsub(3, "Macquarie Island, Auckland & Campbell", t$sealing.ground)
t$sealing.ground <- gsub(4, "Prince Edward, Marion & Crozet Island", t$sealing.ground)
t$sealing.ground <- gsub(5, "Kerguelen Island & Heard Island", t$sealing.ground)

t$sealing.ground <- ordered(t$sealing.ground, levels=c("South Shetlands & South Orkneys",
                                                       "South Georgia & South Sandwich",
                                                       "Prince Edward, Marion & Crozet Island",
                                                       "Kerguelen Island & Heard Island", 
                                                       "Macquarie Island, Auckland & Campbell"))


# Define colours
col <- c("#E7298A", "#6A3D9A", "#CFC855", "#B15928", "#FFA500") #"#33A02C"

# Use Martin Stoffel's GGplot theme as a base
source("Rcode/martin.R")

# Dataframe for titles for each graph with geom_text
dat_text <- data.frame(
  labs = c("(a) South Shetlands & South Orkneys", 
           "(b) South Georgia & South Sandwich", 
           "(c) Prince Edward, Marion & Crozet Island",
           "(d) Kerguelen Island & Heard Island", 
           "(e) Macquarie Island, Auckland & Campbell"),
  sealing.ground   = c("South Shetlands & South Orkneys", 
                       "South Georgia & South Sandwich",
                       "Prince Edward, Marion & Crozet Island",
                       "Kerguelen Island & Heard Island", 
                       "Macquarie Island, Auckland & Campbell"))

# Plot figure
p_sealing <- ggplot(t, aes(year, vessels), fill = sealing.ground) +
  geom_col() +
  geom_area(aes(y = Total/5, fill = sealing.ground), alpha=0.5) + # data for 2nd axis, transformed to match roughly the range of the 1st axis, aplha makes transparent
  scale_fill_manual(values=col) +
  geom_text(aes(label = labs), data = dat_text, x = 1781, y = Inf, hjust=0, vjust = 1.5) +
  #geom_ribbon(aes(ymin=0,ymax=Total, fill=sealing.ground), alpha=0.5)
  facet_wrap(~ sealing.ground, ncol=1) + 
  scale_y_continuous(name = 'Number of vessels', sec.axis = sec_axis(~.*5, name = "Cumulative number of vessels")) + #adding 2nd axis and very important, reverting the above transformation
  guides(fill=FALSE)+
  #annotate("text", xmin=1830, y= 80, label = unique(t$sealing.ground))+
  theme_martin(base_family = "Arial", highlight_family = "Arial")+
  theme(panel.spacing = unit(0.5, "lines"),
        strip.text = element_blank(),
        strip.placement = "inside",
        axis.title.y.right = element_text(vjust=3),
        axis.title.y.left = element_text(vjust=2.5))

# Check figure
p_sealing


#~~ Save figure

ggsave(filename = "Figs-Tables/sealing_effort.jpg", plot = p_sealing, width = 7.5, height = 8.7)
