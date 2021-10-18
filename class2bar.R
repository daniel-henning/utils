class2bar <- function(df, x, y, palette = "Set1", pval = T, binWidth = 1){
  require(plyr)
  require(scales)
  require(ggpubr)
  #全局变量问题解决方案一：将局部变量引入全局变量
  # list2env( list(x=x,y=y), envir=globalenv() )
  
  #解决方案二：重命名引入变量，避免全局变量出现
  countDat <- as.data.frame(group_by(count(df, c(x,y)),get(x)))
  colnames(countDat) <- c(x,"class","freq")
  #
  dt <- ddply(countDat, x, summarise,
              class = class, count = freq, Percent = count/sum(count),
              Percentage = paste0(formatC(count*100/sum(count), digits = 2), "%"))

  dt$class <- as.factor(dt$class)
  if(pval){
    if(length(unique(df[,x]))>1){
      if(nrow(df)>40){
        chisq.res <- chisq.test(table(df[,c(x, y)]))
        ggplot(data=dt, aes(get(x), Percent, fill = class)) +
          geom_col(position=position_fill(),  width = binWidth) +
          scale_y_continuous(labels=percent) +
          ggtitle(paste0("p = ", round(chisq.res$p.value, digits = 4))) +
          labs(x = x, y = "Percent", fill = y) +
          scale_fill_brewer(palette=palette) +
          geom_text(aes(label = paste0(Percentage,"(", count,")")),position = position_stack(vjust = 0.5)) +
          theme_pubr()
        
      }else{
        fisher.res <- fisher.test(table(df[,c(x, y)]))
        ggplot(data=dt, aes(get(x), Percent, fill = class)) +
          geom_col(position=position_fill(), width = binWidth) +
          scale_y_continuous(labels=percent) +
          ggtitle(paste0("p = ", round(fisher.res$p.value, digits = 4))) +
          labs(x = x, y = "Percent", fill = y) +
          scale_fill_brewer(palette=palette) +
          geom_text(aes(label = paste0(Percentage,"(", count,")")),position = position_stack(vjust = 0.5)) +
          theme_pubr()
      }
    }else{
      if(nrow(df)>40){
        ggplot(data=dt, aes(get(x), Percent, fill = class)) +
          geom_col(position=position_fill(),  width = binWidth) +
          scale_y_continuous(labels=percent) + ggtitle(" ") +
          labs(x = x, y = "Percent", fill = y) +
          scale_fill_brewer(palette=palette) +
          geom_text(aes(label = paste0(Percentage,"(", count,")")),position = position_stack(vjust = 0.5)) +
          theme_pubr()
        
      }else{
        ggplot(data=dt, aes(get(x), Percent, fill = class)) +
          geom_col(position=position_fill(),  width = binWidth) +
          scale_y_continuous(labels=percent) + ggtitle(" ") +
          labs(x = x, y = "Percent", fill = y) +
          scale_fill_brewer(palette=palette) +
          geom_text(aes(label = paste0(Percentage,"(", count,")")),position = position_stack(vjust = 0.5)) +
          theme_pubr()
      }
    }
  }else{
    ggplot(data=dt, aes(get(x), Percent, fill = class)) +
      geom_col(position=position_fill(),  width = 0.5) +
      scale_y_continuous(labels=percent) +
      labs(x = x, y = "Percent", fill = y) +
      scale_fill_brewer(palette=palette) +
      geom_text(aes(label = paste0(Percentage,"(", count,")")),position = position_stack(vjust = 0.5)) +
      theme_pubr()
  }
  
  #rm(list = c("x","y"))
}