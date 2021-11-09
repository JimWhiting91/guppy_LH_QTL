#####################################################
# Plot Trinidad geo and rivers...
lib<-c("cowplot","ggplot2","data.table","viridis","patchwork","sf","ggpubr","raster","ggsn")
lapply(lib,library,character.only=T)
source("~/Exeter/five_aside/five_aside_functions.R")

# Read in the file again
trinidad <- st_read(dsn = "data/map_figure/trinidad_tobago_coastline_osm_20141020.shp")
all_rivers <- st_read(dsn = "data/map_figure/all_rivers_trinidad.shp")
oropouche_river <- st_read(dsn = "data/map_figure/oropouche_river_trinidad.shp")
yarra_river <-  st_read(dsn = "data/map_figure/yarra_river_out.shp")

# Plot full Trini
full <- ggplot() +
  geom_sf(data = trinidad)+
  theme_void()+
  xlim(c(62,60.8))+
  ylim(c(10,10.9))+
  theme(panel.border = element_rect(colour="black",fill=NA,size=3))

# Plot 
river_plot2<-ggplot() +
  geom_sf(data = trinidad)+
  xlim(c(-61.7,-60.9))+
  ylim(c(10.55,10.9))+
  theme_void()+
  geom_sf(data = all_rivers,colour="royalblue3",alpha=0.5)+
  geom_sf(data = yarra_river,colour="#9639B2",size=2)+
  geom_sf(data = oropouche_river,colour="#F5B128",size=2)+
  ggsn::scalebar(dist = 10, dist_unit = "km",
                 transform = TRUE, model = "WGS84",location="bottomleft",
                 x.min=-61.7,x.max=-60.9,
                 y.min=10.6,y.max=10.875,height = 0.03, st.size = 5,st.dist=0.05)

# Make random plot to steal legend
random_plot <- data.frame(river=c("Yarra (HP)","Quare (LP)"),
                          colour=c("#9639B2","#F5B128"))
random_plot$river_F <- factor(random_plot$river,levels=random_plot$river)
random_plot$meh <- 2:1
random_plot$meh2 <- 1:2

legend_plot <- ggplot(random_plot,aes(meh,meh2,colour=river_F))+
  geom_line(size=2)+
  scale_colour_manual(breaks = random_plot$river_F,
                      values = random_plot$colour)+
  # guides(color=guide_legend(override.aes=list(fill=NA)))+
  labs(colour="River")+
  theme(legend.key=element_blank(),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18),
        legend.position = c(0.9,0.8),
        legend.justification = "right")

river_legend <- get_legend(legend_plot)


plot.with.inset_no_mountains <-
  ggdraw() +
  draw_plot(river_plot2) +
  draw_plot(full, x = 0.1, y = .6, width = .3, height = .4)+
  draw_plot(river_legend,x=0.5,y=0.5,width=.3,height=.4)

# Save as an R object
saveRDS(plot.with.inset_no_mountains,
        "outputs/qtl_source_rivers.rds")

# Save mountains figure...
pdf("figs/FigureX_part_rivers.pdf",width=6,height=4)
plot.with.inset_no_mountains
dev.off()

