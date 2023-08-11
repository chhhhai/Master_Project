library(ggtree)
gcta<- read.table("gcta.kinship.txt")
hc<-hclust(dist(gcta))
p <- ggtree(hc,layout="circular")+
  xlim(0,5)+
  #geom_text(aes(label=node))+
  geom_highlight(node = 268,fill="red")+
  geom_highlight(node=270,fill="steelblue")+
  geom_highlight(node=271,fill="green")+
  geom_cladelabel(node=268,label="",
                  offset=1.2,barsize = 2,
                  vjust=-2,color="red")+
  geom_cladelabel(node=270,label="",
                  offset=1.2,barsize = 2,
                  hjust=1.2,color="steelblue")+
  geom_cladelabel(node=271,label="",
                  offset=1.2,barsize = 2,
                  hjust=-2.0,color="green")


ggsave('C:/Users/Admin/Desktop/kinship/kinship_poly.jpg', p, device = "jpg", dpi = 300, 
       width=6, height=6, unit = "in")