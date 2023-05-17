color_scheme <- c(
    "UA" = "grey90",
    "G1" = "#566EE7",
    "G2" = "#EA8686",
    "Unassigned" = "grey90",
    "B6" = "#566EE7",
    "CAST" = "#EA8686",

    "Tumor (Liver; Reciprocal B6/CAS)" = "#7382D3",
    "Tumor (Liver; Reciprocal CAS/B6)" = "#D75430"
)

theme_job <- ggplot2::theme(
        text = ggplot2::element_text(family = "Roboto"),
        axis.text = ggplot2::element_text(size = 9),
        axis.title.x = ggtext::element_markdown(size = 12, face = "bold"),
        axis.title.y = ggtext::element_markdown(size = 12, face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        legend.position = "bottom",
        legend.key.size = ggplot2::unit(0.3, "cm"),
        legend.key.width = ggplot2::unit(0.3, "cm"),
        legend.key.height = ggplot2::unit(0.3, "cm"),
        legend.key = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.33),
        strip.background = ggplot2::element_rect(fill = "grey50"),
        strip.text = ggplot2::element_text(color = "white"),
        panel.grid.major.y = ggplot2::element_line(color = "grey75", linetype = "dotted", linewidth = ggplot2::rel(.75)),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(colour = "white"),
        plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
        axis.line = ggplot2::element_line(colour = "black", linewidth = ggplot2::rel(1))
    )

theme_anno_job <- theme_job + ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
)
