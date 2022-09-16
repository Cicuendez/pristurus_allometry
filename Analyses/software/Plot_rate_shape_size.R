
# Packages ----
libs <- c('tidyverse', 'treeio', 'phytools', 'geiger', 'ggtree', 'OUwie', 
          'patchwork', 'mvMORPH', 'ggridges', 'ggpubr')
lapply(libs, require, character.only = TRUE)

# Colors
habitat_broad_colors0 <- c(ground = '#e0710b', rock = '#6C586E', tree = '#119616')
habitat_broad_colors1 <- c(ground = '#B95F89', rock = '#6C586E', tree = '#119616')
habitat_broad_colors2 <- c(ground = '#fada5e', rock = '#6C586E', tree = '#119616')
habitat_broad_colors3 <- c(ground = '#ffd300', rock = '#6C586E', tree = '#119616')
habitat_broad_colors4 <- c(ground = '#ffd300', rock = '#50409a', tree = '#119616')
habitat_broad_colors5 <- c(ground = '#ffd300', rock = '#603e95', tree = '#119616')
habitat_broad_colors6 <- c(ground = '#ffd300', rock = '#313866', tree = '#119616')
habitat_broad_colors7 <- c(ground = '#FDB632', rock = '#C22326', tree = '#027878')
habitat_broad_colors <- c(ground = "#F1B670", rock = "#683B5E", tree = "#E93F7B")
plot(1:3, 1:3, col = habitat_broad_colors, cex = 10, pch = 16)

land_colors0 <- c(mainland = '#C56542', island = '#42978A')
land_colors <- c(mainland = '#F37338', island = '#801638')
plot(1:2, 1:2, col = land_colors0, cex = 6, pch = 16)

state_colors0 <- c(land_colors, habitat_broad_colors7)
state_colors <- c(land_colors0, habitat_broad_colors)
plot(1:5, 1:5, col = state_colors, cex = 6, pch = 16)

viridis::viridis(4)
plot(1:4, 1:4, col = viridis::viridis(4), pch = 16, cex = 8)
viridis::plasma(4)
plot(1:4, 1:4, col = viridis::plasma(4), pch = 16, cex = 8)

model_colors_shape <- c(BM1 = "#440154FF", 
                        OU1 = "#FDE725FF", 
                        BMM = "#35B779FF", 
                        OUM = "#31688EFF")
plot(1:4, 1:4, col = model_colors_shape, pch = 16, cex = 8)
model_colors_size <- c(BM1 = "#440154FF", 
                       OU1 = "#FDE725FF", 
                       BMS = "#35B779FF")
plot(1:3, 1:3, col = model_colors_size, pch = 16, cex = 8)


prismatic::check_color_blindness(state_colors)
prismatic::check_color_blindness(model_colors_shape)
prismatic::check_color_blindness(model_colors_size)

# Set ggplot2 theme
theme.clean <- function(){
  theme_minimal() +
    theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 0.5),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 11, face = "plain"),             
          axis.title.y = element_text(size = 11, face = "plain"),  
          panel.grid = element_blank(), 
          #          panel.grid.major.x = element_blank(),                                          
          #          panel.grid.minor.x = element_blank(),
          #          panel.grid.minor.y = element_blank(),
          #          panel.grid.major.y = element_blank(),  
          axis.line.y = element_line(colour = "black", size = 0.3), 
          axis.line.x = element_line(colour = "black", size = 0.3), 
          axis.ticks.length = unit(.15, "cm"),
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size = 12, vjust = 1, hjust = 0.5),
          plot.subtitle = element_text(size = 11, vjust = 1, hjust = 0.5),
          legend.text = element_text(size = 10, face = "plain"),          
          legend.position = "right")
}
theme_set(theme.clean())


# SHAPE: Import tables of mean AICc ----
meanAICc_shape_land_cons <- readRDS('objects/mvMorph/meanAICc_shape_land_cons.rds')
meanAICc_shape_land_post <- readRDS('objects/mvMorph/meanAICc_shape_land_post.rds')
meanAICc_shape_habitat_cons <- readRDS('objects/mvMorph/meanAICc_shape_habitat_broad_cons.rds')
meanAICc_shape_habitat_post <- readRDS('objects/mvMorph/meanAICc_shape_habitat_broad_post.rds')


# SHAPE: Import tables of AICc ----
aic_df_shape_land_cons <- readRDS('objects/mvMorph/aic_df_shape_land_cons.rds')
aic_df_shape_land_post <- readRDS('objects/mvMorph/aic_df_shape_land_post.rds')
aic_df_shape_habitat_cons <- readRDS('objects/mvMorph/aic_df_habitat_cons.rds')
aic_df_shape_habitat_post <- readRDS('objects/mvMorph/aic_df_habitat_post.rds')


# SHAPE: Import tables of rates (BMM models) ----
shape_rate_land_cons_BMM_df <- readRDS('objects/mvMorph/shape_land_cons_BMM_rate_df.rds')
shape_rate_land_post_BMM_df <- readRDS('objects/mvMorph/shape_land_post_BMM_rate_df.rds')
shape_rate_habitat_cons_BMM_df <- readRDS('objects/mvMorph/shape_habitat_broad_cons_BMM_rate_df.rds')
shape_rate_habitat_post_BMM_df <- readRDS('objects/mvMorph/shape_habitat_broad_post_BMM_rate_df.rds')

head(shape_rate_habitat_cons_BMM_df)
colnames(shape_rate_land_cons_BMM_df)[4] <- 'state'
colnames(shape_rate_land_post_BMM_df)[4] <- 'state'
colnames(shape_rate_habitat_cons_BMM_df)[4] <- 'state'
colnames(shape_rate_habitat_post_BMM_df)[4] <- 'state'


# SIZE: Import tables of mean AICc ----
meanAICc_size_land_cons <- readRDS('objects/OUwie/meanAICc_svl_land_cons.rds')
meanAICc_size_land_post <- readRDS('objects/OUwie/meanAICc_svl_land_post.rds')
meanAICc_size_habitat_cons <- readRDS('objects/OUwie/meanAICc_svl_habitat_cons.rds')
meanAICc_size_habitat_post <- readRDS('objects/OUwie/meanAICc_svl_habitat_post.rds')


# SIZE: Import tables of AICc ----
aic_df_size_land_cons <- readRDS('objects/OUwie/aic_df_svl_land_cons.rds')
aic_df_size_land_post <- readRDS('objects/OUwie/aic_df_svl_land_post.rds')
aic_df_size_habitat_cons <- readRDS('objects/OUwie/aic_df_svl_habitat_cons.rds')
aic_df_size_habitat_post <- readRDS('objects/OUwie/aic_df_svl_habitat_post.rds')

# SIZE: Import tables of rates (BMM models) ----
size_rate_land_cons_BMS_df <- readRDS('objects/OUwie/svl_land_cons_BMS_rate_df.rds')
size_rate_land_post_BMS_df <- readRDS('objects/OUwie/svl_land_post_BMS_rate_df.rds')
size_rate_habitat_cons_BMS_df <- readRDS('objects/OUwie/svl_habitat_cons_BMS_rate_df.rds')
size_rate_habitat_post_BMS_df <- readRDS('objects/OUwie/svl_habitat_post_BMS_rate_df.rds')



# PLOTS - PLOTS - PLOTS - PLOTS - PLOTS ----

# .====== SIZE: AICc plots ========= ----
# AICc consensus ----
meanAICc_size_land_cons
meanAICc_size_habitat_cons


# Merge land and habitat AICc for consensus
aic_df_size_land_cons
aic_df_size_habitat_cons
aic_df_size_cons <- rbind(aic_df_size_land_cons, aic_df_size_habitat_cons)
aic_plot_size_cons <- ggplot(aic_df_size_cons, aes(x = AICc, fill = model)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_histogram(color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = model_colors_size) +
  theme(legend.position = 'bottom') + 
  labs(title = 'AICc consensus tree') +
  ggsave('plots/rates/size/aic_plot_size_cons.pdf')


# AICc posterior ----
meanAICc_size_land_post
meanAICc_size_habitat_post

# Merge land and habitat AICc for posterior
aic_df_size_land_post
aic_df_size_habitat_post
aic_df_size_post <- rbind(aic_df_size_land_post, aic_df_size_habitat_post)
aic_plot_size_post <- ggplot(aic_df_size_post, aes(x = AICc, fill = model)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_density(color = 'transparent', alpha = 0.9) +
  scale_fill_manual(values = model_colors_size) +
  theme(legend.position = 'bottom') + 
  labs(title = 'AICc posterior trees', 
       subtitle = 'BODY SIZE') +
  ggsave('plots/rates/size/aic_plot_size_post.pdf')

# AICc size consensus and posterior ----
aic_plot_size_ggpubr <- ggpubr::ggarrange(aic_plot_size_cons, aic_plot_size_post, 
                                           ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ggsave('plots/rates/size/aic_plot_size_ggpubr.pdf', aic_plot_size_ggpubr)

aic_plot_size_patchwork <- aic_plot_size_cons + aic_plot_size_post + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave('plots/rates/size/aic_plot_size_patchwork.pdf', aic_plot_size_patchwork, 
       width = 8, height = 5)



# .====== SIZE: RATE plots (BMM model) ========= ----
# Rates consensus ----
size_rate_df_cons <- rbind(size_rate_land_cons_BMS_df, size_rate_habitat_cons_BMS_df)
head(size_rate_df_cons)


rate_plot_size_cons <- ggplot(data = size_rate_df_cons, aes(x = sigmasq)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_density(aes(fill = state), color = 'transparent') +
  xlim(0, 0.04) +
  scale_fill_manual(values = state_colors) +
  labs(title = 'Evolutionary rates - body size', 
       subtitle = 'Consensus tree (1,000 character maps)') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5), 
        legend.title = element_blank()) +
  ggsave('plots/rates/size/rate_plot_size_cons.pdf')







# Rates posterior ----
size_rate_df_post <- rbind(size_rate_land_post_BMS_df, size_rate_habitat_post_BMS_df)
nrow(size_rate_df_post)


rate_plot_size_post <- ggplot(data = size_rate_df_post, aes(x = sigmasq)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_density(aes(fill = state), color = 'transparent') +
  xlim(0, 0.05) +
  scale_fill_manual(values = state_colors) +
  labs(title = 'Evolutionary rates - body size', 
       subtitle = 'Posterior trees (100 trees, 100 character maps)') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5), 
        legend.title = element_blank()) +
  ggsave('plots/rates/size/rate_plot_size_post.pdf')



# .====== SHAPE: AICc plots ========= ----
# AICc consensus ----
aic_df_shape_land_cons
aic_df_shapehabitat_cons
meanAICc_shape_habitat_cons

# Merge land and habitat AICc for consensus
aic_df_shape_cons <- rbind(aic_df_shape_land_cons, aic_df_shape_habitat_cons)
aic_plot_shape_cons <- ggplot(aic_df_shape_cons, aes(x = AICc, fill = model)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_histogram(color = 'transparent', alpha = 0.8) +
  scale_fill_manual(values = model_colors_shape) +
  theme(legend.position = 'bottom') + 
  labs(title = 'AICc consensus tree',
       subtitle = 'BODY SHAPE') +
  ggsave('plots/rates/shape/aic_plot_shape_cons.pdf')

# AICc posterior ----
aic_df_shape_land_post
aic_df_shape_habitat_post

# Merge land and habitat AICc for posterior
aic_df_shape_post <- rbind(aic_df_shape_land_post, aic_df_shape_habitat_post)
aic_plot_shape_post <- ggplot(aic_df_shape_post, aes(x = AICc, fill = model)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_density(color = 'transparent', alpha = 0.8) +
  scale_fill_manual(values = model_colors_shape) +
  theme(legend.position = 'bottom') + 
  labs(title = 'AICc posterior trees', 
       subtitle = 'BODY SHAPE') +
  ggsave('plots/rates/shape/aic_plot_shape_post.pdf')

meanAICc_shape_land_post
meanAICc_shape_habitat_post

# AICc shape consensus and posterior ----
aic_plot_shape_ggpubr <- ggpubr::ggarrange(aic_plot_shape_cons, aic_plot_shape_post, 
                                           ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ggsave('plots/rates/shape/aic_plot_shape_ggpubr.pdf', aic_plot_shape_ggpubr)

aic_plot_shape_patchwork <- aic_plot_shape_cons + aic_plot_shape_post + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave('plots/rates/shape/aic_plot_shape_patchwork.pdf', aic_plot_shape_patchwork, 
       width = 8, height = 5)


# .====== SHAPE: RATE plots (BMM model) ========= ----

head(shape_rate_habitat_cons_BMM_df)


colnames(shape_rate_land_cons_BMM_df)[4] <- 'state'
colnames(shape_rate_land_post_BMM_df)[4] <- 'state'
colnames(shape_rate_habitat_cons_BMM_df)[4] <- 'state'
colnames(shape_rate_habitat_post_BMM_df)[4] <- 'state'

# Rates consensus ----
shape_rate_df_cons <- rbind(shape_rate_land_cons_BMM_df, shape_rate_habitat_cons_BMM_df)
nrow(shape_rate_land_cons_BMM_df)
nrow(shape_rate_habitat_cons_BMM_df)

ggplot(data = shape_rate_df_cons, aes(x = sigmasq)) +
  facet_wrap(trait~PC, scales = 'free') +
  geom_density(aes(fill = state), color = 'transparent') +
  scale_fill_manual(values = state_colors) +
  theme(legend.position = 'bottom')

rate_plot_cons_pc1 <- shape_rate_df_cons %>%
  filter(PC == 'PC1') %>%
  ggplot(data = ., aes(x = sigmasq)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_density(aes(fill = state), color = 'transparent') +
  scale_fill_manual(values = state_colors) +
  labs(title = 'PC1') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5))

rate_plot_cons_pc2 <- shape_rate_df_cons %>%
  filter(PC == 'PC2') %>%
  ggplot(data = ., aes(x = sigmasq)) +
  facet_grid(rows = vars(trait), scales = 'free_y', 
             labeller = labeller(.rows = label_value)) +
  geom_density(aes(fill = state), color = 'transparent') +
  scale_fill_manual(values = state_colors) +
  labs(title = 'PC2') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5))



rate_plot_shape_cons_patchwork <- rate_plot_cons_pc1 + rate_plot_cons_pc2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5), 
        legend.title = element_blank())
rate_plot_shape_cons_patchwork + 
  plot_annotation(
    title = 'Evolutionary rates body shape', 
    subtitle = 'Consensus tree (1,000 character maps)', 
    caption = ''
  ) +
  ggsave('plots/rates/shape/rate_plot_shape_cons_patchwork.pdf', 
         width = 8, height = 5)


# Rates posterior ----
shape_rate_df_post <- rbind(shape_rate_land_post_BMM_df, shape_rate_habitat_post_BMM_df)
nrow(shape_rate_land_post_BMM_df)
nrow(shape_rate_habitat_post_BMM_df)

ggplot(data = shape_rate_df_post, aes(x = sigmasq)) +
  facet_wrap(trait~PC, scales = 'free') +
  geom_density(aes(fill = state), color = 'transparent') +
  scale_fill_manual(values = state_colors) +
  theme(legend.position = 'bottom')

rate_plot_post_pc1 <- shape_rate_df_post %>%
  filter(PC == 'PC1') %>%
  ggplot(data = ., aes(x = sigmasq)) +
  facet_grid(rows = vars(trait), scales = 'free_y') +
  geom_density(aes(fill = state), color = 'transparent') +
  scale_fill_manual(values = state_colors) +
  labs(title = 'PC1') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5))

rate_plot_post_pc2 <- shape_rate_df_post %>%
  filter(PC == 'PC2') %>%
  ggplot(data = ., aes(x = sigmasq)) +
  facet_grid(rows = vars(trait), scales = 'free_y', 
             labeller = labeller(.rows = label_value)) +
  geom_density(aes(fill = state), color = 'transparent') +
  scale_fill_manual(values = state_colors) +
  labs(title = 'PC2') +
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5))



rate_plot_shape_post_patchwork <- rate_plot_post_pc1 + rate_plot_post_pc2 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5), 
        legend.title = element_blank())
rate_plot_shape_post_patchwork <- rate_plot_shape_post_patchwork + 
  plot_annotation(
    title = 'Evolutionary rates body shape', 
    subtitle = 'Posterior trees (100 trees, 100 character maps)', 
    caption = ''
  ) +
  ggsave('plots/rates/shape/rate_plot_shape_post_patchwork.pdf', 
         width = 8, height = 5)



# PLOT AICc POSTERIOR SIZE AND SHAPE ----

aic_plot_size_shape_post <- aic_plot_size_post + aic_plot_shape_post +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', 
        plot.subtitle = element_text(size = 10, vjust = 1, hjust = 0.5), 
        legend.title = element_blank())

ggsave('plots/rates/aic_plot_size_shape_post.pdf', 
       plot = aic_plot_size_shape_post,
       width = 8, height = 5)

# PLOT RATES POSTERIOR SIZE AND SHAPE ----


rate_plot_size_post + rate_plot_shape_post_patchwork

# Match plots to areas by name
design <- "A
           B"
wrap_plots(A = rate_plot_size_post, 
           B = rate_plot_shape_post_patchwork, 
           design = design)
ggsave('plots/rates/rate_plot_size_shape_post.pdf', last_plot(), 
       height = 10, width = 8)


# .____#######################__ ----
# .____#######################__ ----
# .____#######################__ ----
# Density ridges probando ----
ggplot(aic_df_land_cons, aes(x = AICc, y = model, fill = model)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2) +
  scale_fill_viridis_d() +
  theme(legend.position = 'bottom') +
  ggsave('plots/rates/shape/density_ridges_prueba.pdf')





