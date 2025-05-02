original_data_QP$Method = "grey box mean"
DataQPp3$Method = "plus 3 s.d."
original_data_QPm3$Method = "minus 3 s.d."
noisy_data$Method = "NR"


GPStateVars = rbind.fill(original_data_QP, DataQPp3, original_data_QPm3, noisy_data)
GPStateVars[c(3001:4000),c(2:11)] = GPStateVars[c(3001:4000),c(13:22)]
View(GPStateVars)
GPStateVars = GPStateVars[,-c(13:22)]
View(GPStateVars)

noisy_data_tb_QP <- as_tibble(GPStateVars) %>%
  pivot_longer(cols = QmtQP:VpuQP, 
               names_to = "variable", 
               values_to = "vals")

p_noisy_data_tb <- ggplot(noisy_data_tb_QP) +
  geom_line(aes(x = time, y = vals, color = Method, lty = Method), lwd = 1) +
  facet_grid(vars(variable), scales = "free") +
  scale_color_manual(values = c("grey box mean" = "black", "minus 3 s.d." = "red", "plus 3 s.d." = "red", "NR" = "blue")) +
  scale_linetype_manual(values = c("grey box mean" = "solid", "minus 3 s.d." = "dotted", "plus 3 s.d." = "dotted", "NR" = "dashed")) +
  theme_minimal() + 
  labs(Title = "",
       x = "Time (s)",
       y = "")

p_noisy_data_tb

ggplot() + 
  geom_line(data = noisy_data_tb, aes(x = time, y = vals))+
  geom_line(data = noisy_data_tb_QP, aes(x = time, y = vals)) +
  facet_grid(vars(variable), scales = "free")
