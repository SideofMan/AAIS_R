x=matrix(linspace(-10,10,10001), ncol=1)
y=t_mixture_pdf(x,Proposal$M,Proposal$W,Proposal$Mu,Proposal$Sigma,df)

my_df <- data.frame(x = x, y = y)
ggplot(data = my_df, aes(x,y)) +
  geom_line() +
  scale_x_continuous(limits = c(-10,10), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.3), expand = c(0,0)) +
  ggtitle("2D Plot of Gaussian 1D Mixture function") +
  theme(plot.title = element_text(hjust = 0.5))