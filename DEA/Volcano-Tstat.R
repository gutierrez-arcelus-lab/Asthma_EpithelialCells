#Regression Model Plots
#Purpose: Generates plots and stats for a given a mixedmodel regression analysis
#   1. P-value barplot.
#   2. Volcano plot showing significant instances. Using a p-value threshold.
#   3. Volcano plot showing significant instances that also pass a fold-change threshold.
#Author: Daniela Fernandez Salinas
#Date: February 17th 2022
#Last Modification:
#Developer Notes: Incorporate option to modify the p_value threshold or FDR percent.

#Packages
require(ggplot2)
require(hrbrthemes)
require(gt)
require(optparse)
require(tidyverse)


parser=list(
  make_option(c("--mixedmodel","-m"),help="Table with the betas and p-values"),
  make_option(c("--outputP","-o"),help="Path where to write the plots."),
  make_option(c("--plot","-t"),help="Type of plot to generate, choices: \nh - For p-value histogram.\nt - For both volcano plot with top 10% ranked by t-statistic.\nt2 - For opposite side of the volcano.\nall - For all of the above.\nDEFAULT: t"),
  make_option(c("--number"),help="Use to provide the number of instances considered in the original dataset."),
  make_option(c("--filename","-n"),help="Indicate prefix for image filenames."),
  make_option(c("--color","-c"),help="Color to use for selected instances in the volcano plot.")
)

args=parse_args(OptionParser(option_list=parser))
if(is.character(args$plot)==FALSE){
  plot_type<-"t"
}else{
  plot_type<-args$plot
}

if(is.character(args$color)==FALSE){
  highlight_color<-"turquoise"
}else{
  highlight_color<-args$color
}

#Palettes
pastel<-c("pink","black","lightblue")
halloween<-c("orange","black","purple")
#halloween<-c("#f793a9","black","#8cc34d")

# File reading
peaks<-read.table(args$mixedmodel, header = T, row.names = 1, sep = '\t') %>% mutate(tstat=x1/Std)
if(is.character(args$filename)){
  filename<-args$filename
  }else{
  fullname<-unlist(strsplit(args$mixedmodel,"/"))
  filename<-fullname[length(fullname)]
  filename<-sub("_mixedmodel.txt","",filename)
}
print(filename)

#P-value Histogram
if(plot_type=="all" || plot_type=="h" ){
  ggplot(peaks, aes(x=p_value))+
    geom_histogram(color="#e9ecef", alpha=0.6, position="identity", fill= "#69b3a2")+
    theme_minimal()
  ggsave(paste(args$o,filename,"_histogram.png",sep=""))
}

#Volcano 
if(plot_type=="all" || plot_type=="t" ){
  p10<-round(nrow(peaks)*0.1)
  color<-c(rep("Selected",p10),rep("Non",nrow(peaks)-p10))
  print(length(color))

  peaks<-peaks %>% arrange(desc(tstat)) %>% mutate(col=factor(color,levels=c("Non","Selected"))) %>%arrange(col)
  alphas <- c("Up" = 1, "Down" = 1, "Non significant" = 0.5)
  limit<-max(abs(peaks$x1))


  ggplot(peaks,aes(x=x1, y=-log10(p_value)))+
    geom_point(aes(color= col))+
    scale_color_manual(values = c("darkgrey",highlight_color) )+
    theme_classic()+
    xlab(expression(beta))+
    ylab(bquote(.(title)~-log[10]*'(P)'))+
    #ylab("-log10(P)")+
    #geom_hline(yintercept = -log10(0.05),linetype="dashed")+
    #geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed")+
    scale_alpha_manual(values = alphas)+
    scale_x_continuous(limits = c(-limit, limit))+
    guides(color="none")+
    #ggtitle(paste("Differential Selection for",filename,sep=" "))+
    theme(
      #plot.title=element_text(hjust = 0.5,size=15),
      #plot.subtitle =element_text(hjust = 0.5,size=12),
      legend.title = element_blank(),
      axis.text=element_text(size=25), 
      legend.text=element_text(size=15),
      axis.title = element_text(size=26))

ggsave(paste(args$o,filename,"_vt.png",sep=""),width=6,height=6,dpi=600)
}

#Upside
if(plot_type=="all" || plot_type=="t2" ){

  color<-c(rep("Selected",round(nrow(peaks)*0.1)),rep("Non",round(nrow(peaks)*0.9)))
  peaks<-peaks %>% arrange(tstat) %>% mutate(col=color)
  alphas <- c("Up" = 1, "Down" = 1, "Non significant" = 0.5)
  limit<-max(abs(peaks$x1))*1.1

  ggplot(peaks, aes(x=x1, y=-log10(p_value)))+
    geom_point(aes(color= col))+
    scale_color_manual(values = c("darkgrey",highlight_color) )+
    theme_classic()+
    xlab(expression(beta))+
    ylab("-log10(P)")+
    geom_hline(yintercept = -log10(0.05),linetype="dashed")+
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed")+
    scale_alpha_manual(values = alphas)+
    scale_x_continuous(limits = c(-(limit), limit))+
    ggtitle(paste("Differential Selection for",filename,sep=" "))+
    theme(
      plot.title=element_text(hjust = 0.5,size=15),
      plot.subtitle =element_text(hjust = 0.5,size=12),
      legend.title = element_blank(),
      axis.text=element_text(size=18), 
      legend.text=element_text(size=11),
      axis.title = element_text(size=13))

ggsave(paste(args$o,filename,"_vt2.png",sep=""))
}
