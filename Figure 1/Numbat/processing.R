library(data.table)
library(dplyr)
library(glue)
library(stringr)
library(numbat)
library(ggplot2)
library(patchwork)

# read in the numbat results and spot coordinates
nb = list()
spots = list()
samples = c('sample_1','sample_2','sample_3','sample_4','sample_5','sample_6','sample_7','sample_8','sample_9','sample_10','sample_11','sample_12')
for (sample in samples) {
  nb[[sample]] = Numbat$new(glue('outs/{sample}'))
  spots[[sample]] = fread(glue('tp/{sample}/tissue_positions_list.csv'))
  colnames(spots[[sample]]) <- c("barcode","in_tissue","array_row","array_col","tissue_scaled_x","tissue_scaled_y")
}

options(repr.plot.width = 8.5, repr.plot.height = 8, repr.plot.res = 300)

#plot CNV probability on samples
lapply(
  samples,
  function(sample) {
    
    spots[[sample]] %>%
      left_join(
        nb[[sample]]$clone_post, 
        by = c('barcode' = 'cell')
      ) %>%
      filter(in_tissue == 1) %>%
      ggplot(
        aes(x = array_col, y = array_row)
      ) +
      geom_point(aes(color = p_cnv), size = 1, alpha = 0.8, pch = 16) +
      scale_color_gradient2(
        low = 'darkgreen', high = 'red3', mid = 'yellow', 
        midpoint = 0.5, limits = c(0,1), oob = scales::oob_squish
      ) +
      theme_bw() +
      ggtitle(sample) + 
      scale_y_reverse()
  }
) %>% 
  wrap_plots(guides = 'collect')


chrom_labeller <- function(chr) {
  chr[chr %in% c(19, 21, 22)] = ""
  return(chr)
}

pal = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
getPalette = colorRampPalette(pal)

cnv_colors = c("neu" = "gray",
               "neu_up" = "darkgray", "neu_down" = "gray",
               "del_up" = "royalblue", "del_down" = "darkblue",
               "loh_up" = "darkgreen", "loh_down" = "olivedrab4",
               "amp_up" = "red", "amp_down" = "tomato3",
               "del_1_up" = "royalblue", "del_1_down" = "darkblue",
               "loh_1_up" = "darkgreen", "loh_1_down" = "olivedrab4",
               "amp_1_up" = "red", "amp_1_down" = "tomato3",
               "del_2_up" = "royalblue", "del_2_down" = "darkblue",
               "loh_2_up" = "darkgreen", "loh_2_down" = "olivedrab4",
               "amp_2_up" = "red", "amp_2_down" = "tomato3",
               "del_up_1" = "royalblue", "del_down_1" = "darkblue",
               "loh_up_1" = "darkgreen", "loh_down_1" = "olivedrab4",
               "amp_up_1" = "red", "amp_down_1" = "tomato3",
               "del_up_2" = "royalblue", "del_down_2" = "darkblue",
               "loh_up_2" = "darkgreen", "loh_down_2" = "olivedrab4",
               "amp_up_2" = "red", "amp_down_2" = "tomato3",
               "bamp" = "salmon", "bdel" = "skyblue",
               "amp" = "tomato3", "loh" = "olivedrab4", "del" = "royalblue",
               "theta_up" = "darkgreen", "theta_down" = "olivedrab4",
               "theta_1_up" = "darkgreen", "theta_1_down" = "olivedrab4",
               "theta_2_up" = "darkgreen", "theta_2_down" = "olivedrab4",
               "theta_up_1" = "darkgreen", "theta_down_1" = "olivedrab4",
               "theta_up_2" = "darkgreen", "theta_down_2" = "olivedrab4",
               '0|1' = 'red', '1|0' = 'blue','major' = '#66C2A5', 'minor' = '#FC8D62')

cnv_labels = names(cnv_colors) %>%
  stringr::str_remove_all('_') %>%
  stringr::str_to_upper() %>%
  stringr::str_replace('UP', '(major)') %>%
  stringr::str_replace('DOWN', '(minor)') %>%
  stringr::str_replace('LOH', 'CNLoH') %>%
  setNames(names(cnv_colors))

numbat_plot <- function (segs) {
  chrom_labeller <- function(chr) {
    chr[chr %in% c(19, 21, 22)] = ""
    return(chr)
  }
  ggplot(segs) + geom_rect(aes(xmin = seg_start, xmax = seg_end, 
                               ymin = -0.5, ymax = 0.5, fill = cnv_state_post)) + theme_void() + 
    theme(panel.spacing = unit(1, "mm"), strip.background = element_blank(), 
          strip.text.y = element_text(angle = 0), plot.margin = margin(0, 
           0, 0, 0), legend.position = "none") + facet_grid(~CHROM, 
    space = "free_x", scales = "free") + 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )+ 
    scale_fill_manual(values = cnv_colors, labels = cnv_labels, 
                      name = "CN states") +
    ggrepel::geom_text_repel(aes(x = (seg_start +
    seg_end)/2, y = -0.5, label = str_remove(seg_cons, "\\d+")),
    min.segment.length = 0, vjust = 1, hjust = 0, direction = "x",
    segment.curvature = -0.2, segment.ncp = 3, segment.angle = 30,
    segment.inflect = TRUE, max.overlaps = 3) + scale_y_continuous(expand = expansion(add = c(0.5,
    0))) + scale_x_continuous(expand = expansion(mult = 0.05)) +
    guides(fill = "none")
}

png("numbat_plot.png", units = "in", res = 300, height = 5, width = 10, type = "cairo")
options(repr.plot.width = 14, repr.plot.height = 12, repr.plot.res = 300)
lapply(
  samples,
  function(sample) {
    numbat_plot(nb[[sample]]$segs_consensus) 
  }
) %>% 
  wrap_plots(guides = 'collect', ncol = 1)
dev.off()
