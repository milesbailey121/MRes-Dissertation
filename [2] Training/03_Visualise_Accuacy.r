library("ggplot2","tidyverse", "ggpubr")

value <- c(0.19064893188230553,
 0.3175322001999647,
 0.2855738779076117,
 0.43418173804914134,
 0.604411642064882,
 0.4557584412934816,
 0.28118072570295183,
 0.18479911638823693,
 0.010954872465663833,
 0.07788205387146763,
 0.2694435037814652,
 0.28149994204242496,
 0.40544608648056923,
 0.2894412897589644,
 0.3965218596272431,
 0.366015963754137,
 0.2072101122146266,
 0.28168367932102306,
 0.24043768761498982,
 0.36209383218800373) 

percentage <- value * 100

mem <- c(0.19064893188230553,
          0.3175322001999647,
          0.2855738779076117,
          0.43418173804914134,
          0.604411642064882,
          0.4557584412934816,
          0.28118072570295183,
          0.18479911638823693,
          0.010954872465663833,
          0.07788205387146763)
nuc <- c(0.2694435037814652,
         0.28149994204242496,
         0.40544608648056923,
         0.2894412897589644,
         0.3965218596272431,
         0.366015963754137,
         0.2072101122146266,
         0.28168367932102306,
         0.24043768761498982,
         0.36209383218800373)


name <- c('Mx_SA TMA 5_G02',
           'Mx_SA TMA 5_G16',
           'Mx_SA TMA 5_I12',
           'Mx_SA TMA 5_I17',
           'Mx_SA TMA 5_J01',
           'Mx_SA TMA 5_J10',
           'Mx_Southampton TMA Research Block 5_G10',
           'Mx_Southampton TMA Research Block 5_G11',
           'Mx_Southampton TMA Research Block 5_H02',
           'Mx_Southampton TMA Research Block 5_H03',
           'Mx_SA TMA 5_G02',
           'Mx_SA TMA 5_G16',
           'Mx_SA TMA 5_I12',
           'Mx_SA TMA 5_I17',
           'Mx_SA TMA 5_J01',
           'Mx_SA TMA 5_J10',
           'Mx_Southampton TMA Research Block 5_G10',
           'Mx_Southampton TMA Research Block 5_G11',
           'Mx_Southampton TMA Research Block 5_H02',
           'Mx_Southampton TMA Research Block 5_H03') 


df <- data.frame(name,percentage,type)

ggplot(df, mapping = aes(x,percentage, fill = type)) + 
  geom_bar(position = "dodge", stat = "identity")+
  geom_text(aes(label = sprintf("%.1f%%", percentage)), position = position_dodge(width = 0.9), vjust = -0.5) +
  ggpubr::theme_pubr()+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(labels = c("1", "2", "3", "4", "5","6", "7", "8", "9", "10"))+
  ylab("Accuracy(IoU)")+
  xlab("Image")

ggplot(df, aes(x, percentage, fill = type)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7, size = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3) +
  ggpubr::theme_pubr() +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) +
  ylab("Accuracy (IoU)") +
  xlab("Image")+
  labs("Type")
