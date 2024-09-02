library(ggplot2)
library(dplyr)
library(stringr)

# Function to generate plot for a given file
generate_plot <- function(file_path, output_path, title, bar_color, text_size) {
  
  data <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  data$GO_Term <- sapply(strsplit(data$Term, " \\("), `[`, 1)
  data <- data %>%
    arrange(P.value) %>%
    head(10)
  
  data$GO_Term <- str_wrap(data$GO_Term, width = 30)
  data$bar_color <- ifelse(data$P.value <= 0.05, bar_color, "grey")
  
  plot <- ggplot(data, aes(x = reorder(GO_Term, Combined.Score), y = Combined.Score)) +
    geom_bar(stat = "identity", aes(fill = bar_color), color = "black") +
    geom_text(aes(label = ifelse(round(P.value, 3)>0, paste0("p=", round(P.value, 3)), paste0("p<0.001"))), 
              hjust = 0.5, vjust = -1, color = "white", size = 2, angle=90, fontface = "bold") +  # Add p-value labels with "p="
    coord_flip() +
    theme_minimal() +
    theme(aspect.ratio = 1) +
    scale_fill_identity() +
    labs(title = title,
         x = "Term",
         y = "Combined Score") +
    theme(axis.text.y = element_text(size = text_size, face = "bold"),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(output_path, plot, width = 6, height = 9, dpi = 300)
}


setwd("C:\\Users\\fopa0001\\Downloads\\OneDrive_1_04-12-2023\\enriched_terms/")

generate_plot(file_path = "./Coexpression_Predicted_GO_Cellular_Component_2018_table_mot_up.txt", 
              output_path = "./Coexpression_Predicted_GO_Cellular_Component_2018_table_mot_up.png", 
              title="Coexpression Cellular Components", bar_color = "blue3", text_size = 11)

generate_plot(file_path = "./Coexpression_Predicted_GO_Cellular_Component_2018_table_mot_down.txt", 
              output_path = "./Coexpression_Predicted_GO_Cellular_Component_2018_table_mot_down.png", 
              title="Coexpression Cellular Components", bar_color = "blue3", text_size = 11)

generate_plot(file_path = "./Coexpression_Predicted_GO_Cellular_Component_2018_table_den.txt", 
              output_path = "./Coexpression_Predicted_GO_Cellular_Component_2018_table_den.png", 
              title="Coexpression Cellular Components", bar_color = "blue3", text_size = 11)

generate_plot(file_path = "./Coexpression_Predicted_GO_Biological_Process_2018_table_mot_up.txt", 
              output_path = "./Coexpression_Predicted_GO_Biological_Process_2018_table_mot_up.png", 
              title="Coexpression Biological Process", bar_color = "blue3", text_size = 11)

generate_plot(file_path = "./Coexpression_Predicted_GO_Biological_Process_2018_table_mot_down.txt", 
              output_path = "./Coexpression_Predicted_GO_Biological_Process_2018_table_mot_down.png", 
              title="Coexpression Biological Process", bar_color = "blue3", text_size = 11)

generate_plot(file_path = "./Coexpression_Predicted_GO_Biological_Process_2018_table_den.txt", 
              output_path = "./Coexpression_Predicted_GO_Biological_Process_2018_table_den.png", 
              title="Coexpression Biological Process", bar_color="blue3", text_size=11)

generate_plot(file_path = "./BioCarta_2016_table_mot_down.txt", 
              output_path = "./BioCarta_2016_table_mot_down.png", 
              title="BioCarta Pathways", bar_color = "brown3", text_size = 8)

generate_plot(file_path = "./BioCarta_2016_table_mot_up.txt", 
              output_path = "./BioCarta_2016_table_mot_up.png", 
              title="BioCarta Pathways", bar_color="brown3", text_size=8)

generate_plot(file_path = "./BioCarta_2016_table_den.txt", 
              output_path = "./BioCarta_2016_table_den.png", 
              title="BioCarta Pathways", bar_color="brown3", text_size=8)

generate_plot(file_path = "./Panther_2016_table_mot_down.txt", 
              output_path = "./Panther_2016_table_mot_down.png", 
              title="Panther Pathways", bar_color="brown3", text_size=8)

generate_plot(file_path = "./Panther_2016_table_mot_up.txt", 
              output_path = "./Panther_2016_table_mot_up.png", 
              title="Panther Pathways", bar_color = "brown3", text_size = 8)

generate_plot(file_path = "./Panther_2016_table_den.txt", 
              output_path = "./Panther_2016_table_den.png", 
              title="Panther Pathways", bar_color="brown3", text_size=8)
