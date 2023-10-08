######------------------------------------------------
###### world map showing the sample collection
######------------------------------------------------
#https://www.datanovia.com/en/blog/how-to-create-a-map-using-ggplot2/
library(ggplot2)
library(dplyr)
require(maps)
require(viridis)
theme_set(
  theme_void()
)

world_map <- map_data("world")
unique(world_map$region)
unique(world_map$subregion)
ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")+theme_void()



some.eu.countries <- c("USA" ,"France","UK",'Israel')

some.eu.maps <- map_data("world", region = some.eu.countries)

region.lab.data <- some.eu.maps %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

ggplot(world_map, aes(x = long, y = lat)) +
  #geom_polygon(aes( group = group, fill = region))+
  geom_polygon(fill="lightgray", colour = "white")+
  geom_text(aes(label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  #scale_fill_viridis_d()+
  theme_void()+
  theme(legend.position = "none")


ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")+
  geom_text(aes(label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  #scale_fill_viridis_d()+
  theme_void()+
  theme(legend.position = "none")
##
#######
head(world_map)

world_map2= world_map
cl <- c("USA" ,"France","UK",'Israel')
world_map2$label= ifelse(world_map2$region %in% cl, world_map2$region, '')

ggplot(world_map2, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white")+theme_void()

ggplot(world_map2, aes(long, lat, label = label)) +
  ggrepel::geom_text_repel() +
  geom_point(color = ifelse(world_map2$label == "", "grey50", "red"))



library(ggplot2)
library(dplyr)

thismap = map_data("world")

# Set colors
thismap <- mutate(thismap, fill = ifelse(region %in% cl, "steelblue", "white"))

# Use scale_fiil_identity to set correct colors
g1= ggplot(thismap, aes(long, lat, fill = fill, group=group)) + 
  geom_polygon(colour="gray",linewidth = 0.05) + #ggtitle("Map of World") + 
  scale_fill_identity()+theme_void()

g2= ggplot(thismap, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="white", colour = "gray")+theme_void()

pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/worlmap_metag_collect3.pdf',
    width = 14,height = 2)
cowplot::plot_grid(g1,g2,nrow = 1)
dev.off()
##
## mydata
load('Documents/project/mmm_aging/tmpall_20230821.RData')
head(meta.crc)
table(meta.crc$country)
dim(meta.crc)
table(meta.crc$country, meta.crc$Study)
table(meta.crc$Study,meta.crc$disease)
library(dplyr)
enrichedpath=  table(meta.crc$Study) %>% as.data.frame() 
##
fmtdata= read.csv('/Users/xiaoqiangzhu/Documents/data_from_macintel/shanghai/chen_lungcancer/FMT/analysis_20230508/Combined_25_donors_from_Davar&Baruch_etal.csv',header = 
                    T)
head(fmtdata)
table(fmtdata$study)
##
enrichedpath= rbind(table(meta.crc$Study) %>% as.data.frame() , 
                    table(fmtdata$study) %>% as.data.frame()
)

dotplotenrich <- enrichedpath %>% 
  ggplot(aes(x = Var1, y = Var1, size = Freq, fill = Freq, color = Freq)) +
  geom_point(shape = 21) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_fill_distiller(palette ="RdBu", direction = -1) 
  paletteer::scale_fill_paletteer_c("viridis::plasma")
#scale_fill_distiller(palette = "RdPu", direction= 1) +
#scale_color_distiller(palette = "RdPu",direction= 1)


pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/worlmap_metag_collect2.pdf',
    width = 4,height = 4)
cowplot::plot_grid(dotplotenrich,nrow = 1)
dev.off()
