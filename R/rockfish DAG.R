library(dagitty)
library(ggdag)
library(ggplot2)
# DAG
dag <- readLines("~/My Drive/SFU REM/Research/CCIRA Rockfish/DAG August 19 2025.txt")
dag_string <- paste(dag, collapse = "\n")
ccira_dag <- dagitty::dagitty(dag_string)

print(ccira_dag)

gg_dag <- ggdag::tidy_dagitty(ccira_dag)
ggplot(data=gg_dag,aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point() +
  geom_dag_edges_diagonal() +
  geom_dag_text() +
  theme_dag()

## ggplot
ggdag(gg_dag,layout="auto") + theme_dag()


# Get node roles from dagitty
nodes <- tibble(
  name = names(dagitty::coordinates(ccira_dag)$x),
  exposure = name %in% exposures(ccira_dag),
  outcome  = name %in% outcomes(ccira_dag),
  latent   = name %in% latents(ccira_dag)
)

# Tidy DAG and join role info
tidy <- tidy_dagitty(ccira_dag)
tidy$data <- left_join(tidy$data, nodes, by = "name")

# Create a role variable
tidy$data <- tidy$data %>%
  mutate(role = case_when(
    exposure ~ "Exposure",
    outcome  ~ "Outcome",
    latent   ~ "Latent",
    TRUE     ~ "Other"
  ))

tidy$data$name <- stringr::str_wrap(tidy$data$name, width = 25)

p <- ggdag(tidy, layout = "auto",text=TRUE,node=FALSE) +
  geom_dag_point(size = 10, aes(color = role)) +
  geom_dag_text(aes(label = name), 
                 fill="white",  # keep node color visible
                 color = "black", 
                 size = 2,vjust=1.5) +
  # geom_dag_text(color = "white", size = 2) +
  geom_dag_edges(color = "#34495E", size = 0.8) +
  theme_dag() +
  scale_color_manual(
    name = "Node Type",
    values = c(
      "Exposure" = "#1f78b4",
      "Outcome"  = "#e31a1c",
      "Latent"   = "#6a3d9a",
      "Other"    = "#b2b2b2"
    )
  ) +
  theme(text = element_text(size = 11),legend.position = "top")

# Print plot
print(p)

# Save the plot with text size 11 scaling
ggsave("Figures/ccira_dag.jpeg", p, width = 7, height = 5, dpi = 600)
