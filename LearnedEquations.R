learned1 = sin(tsteps * 8.296543) + 4.008249
learned2 = sin(tsteps * 8.29644) + 4.008488
learned3 = sin(tsteps * 8.296529) + 4.0140643
learned4 = Vrv * 0.07002716+sin(sin((0.5178254-1.1208304*(tsteps+tsteps))*
                                   (-Vrv + Vrv - 1 * 3.7371519)) *1.7992761) * 1.9446352 - 1.3000957
learned5 = 4.002479 - sin(tsteps * (-8.296508))

learned6 = Vlv * 0.30795234 - (-9.534985) * sin(0.048446856 * tsteps - 0.48549715) - 19.58908

learned7 = sin(tsteps * 8.296612) + 4.0126266
learned8 = (-(-0.02852863) * Vlv + cos((tsteps - 1.9079664) * (-8.377406))) * 1.8544639
learned9 = cos(cos(tsteps * 4.2659655)) * 5.252529

learned10 = (Vlv * 0.0324214+cos(tsteps * 0.0485412-1*2.0678177)-2.0675468)*9.469835

df <- data.frame(
  t = tsteps,
  learned1, learned2, learned3, learned4, learned5,
  learned6, learned7, learned8, learned9, learned10,
  Vspt
)

# Convert to long format
df_long <- df %>%
  pivot_longer(cols = -t, names_to = "Function", values_to = "Value")

# Define colors: black for Vspt, others from a color palette
custom_colors <- c("Vspt" = "black",
                   setNames(rainbow(10), paste0("learned", 1:10)))

# Plot
ggplot(df_long, aes(x = t, y = Value, color = Function)) +
  geom_line() +
  scale_color_manual(values = custom_colors) +
  labs(title = "Learned expressions",
       x = "Time (tsteps)",
       y = "Vspt (ml)") +
  theme_minimal()


learned1MSE = mean((Vspt - learned1)^2)
learned2MSE = mean((Vspt - learned2)^2)
learned3MSE = mean((Vspt - learned3)^2)
learned4MSE = mean((Vspt - learned4)^2)
learned5MSE = mean((Vspt - learned5)^2)
learned6MSE = mean((Vspt - learned6)^2)
learned7MSE = mean((Vspt - learned7)^2)
learned8MSE = mean((Vspt - learned8)^2)
learned9MSE = mean((Vspt - learned9)^2)
learned10MSE = mean((Vspt - learned10)^2)

learned10MSE

