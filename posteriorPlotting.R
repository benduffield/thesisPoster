
QP_plot <- ggplot() + 
  geom_line(aes(x = time_pred, y = Pred_df[,1], color = "Predicted Vspt"), lwd = 1) + 
  geom_point(aes(x = tsteps_obs, y = Vspt_obs, color = "Training data")) +
  geom_line(aes(x = time_pred, y = Vspt, color = "True Vspt"), alpha = 0.5) +
  geom_line(aes(x = tsteps, y = Pred_df[,3], color = "3 standard deviations"), lty = "dotted", lwd = 1) + 
  geom_line(aes(x = tsteps, y = Pred_df[,4], color = "3 standard deviations"), lty = "dotted", lwd = 1) + 
  labs(title = "GP predictions with QP kernel",
       x = "Time (s)",
       y = "Vspt (ml)") +
  scale_color_manual(name = NULL,  # Removes legend title for color
                     values = c("Predicted Vspt" = "#007d69", 
                                "Training data" = "#003c3c",
                                "True Vspt" = "black",
                                "3 standard deviations" = "red")) +
  theme(legend.position = "none")  # Remove individual legend

SE_plot <- ggplot() + 
  geom_line(aes(x = time_pred, y = Pred_data[,1], color = "Predicted Vspt"), lwd = 1) + 
  geom_point(aes(x = tsteps_obs, y = Vspt_obs, color = "Training data")) + 
  geom_line(aes(x = time_pred, y = Vspt, color = "True Vspt"), alpha = 0.5) +
  geom_line(aes(x = tsteps, y = Pred_data[,3], colour = "3 standard deviations"), lty = "dotted", lwd = 1) +
  geom_line(aes(x = tsteps, y = Pred_data[,4], colour = "3 standard deviations"), lty = "dotted", lwd = 1) +
  labs(title = "GP predictions with SE kernel",
       x = "",
       y = "Vspt (ml)") +
  scale_color_manual(name = NULL,  # Removes legend title for color
                     values = c("Predicted Vspt" = "#007d69", 
                                "Training data" = "#003c3c", 
                                "True Vspt" = "black",
                                "3 standard deviations" = "red")) +
  theme(legend.position = "none")  # Remove individual legend

# Combine with shared legend
final_plot <- (SE_plot / QP_plot) + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(final_plot)

