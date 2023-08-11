library('ggplot2')

dat <- read.delim("C:/Users/Admin/Desktop/admixture-plot/CV.txt",
                  header=T, 
                  row.names=1, 
                  sep="\t",
                  stringsAsFactors = FALSE,
                  check.names = FALSE)
p=ggplot(data=dat,
         aes(x=Compartment,y=CV))+
  geom_point(size=3,colour="#85BA8F")+
  labs(x="K value", y="Cross-validation error")+
  geom_line(size=1.2,colour="#85BA8F")+scale_x_continuous(breaks=seq(0, 22, 1))+scale_y_continuous(expand = c(0,0),limits = c(0.37,0.52))+theme_bw()+theme(panel.grid=element_blank())

p1 <- p+theme(axis.title.x=element_text(face="bold", color="black",size=20,family="serif"),
              axis.title.y = element_text(colour="black",size=20,face="bold",family="serif"),
              axis.text.x = element_text(colour="black",size=16,face="bold",family="serif"),
              axis.text.y = element_text(colour="black",size=16,face="bold",family="serif"))+theme(axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(size = 16,face="bold",family="serif"),legend.title = element_text(size = 16,,face="bold"))

p1
ggsave('C:/Users/Admin/Desktop/admixture-plot/test.jpg', p1, device = "jpg", dpi = 300, 
       width=10, height=8, unit = "in")