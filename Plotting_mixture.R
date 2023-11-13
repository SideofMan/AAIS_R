x=seq(-10, 10, length = 201)
y=seq(-10, 10, length = 201)
grid=expand.grid(x=x,y=y)
Z=with(grid, t_mixture_pdf(cbind(grid$x, grid$y), Proposal$M, Proposal$W, Proposal$Mu, Proposal$Sigma, df))

library(plotly)

plot_ly()%>%
  add_surface(x = x, y = y, z = matrix(Z, nrow = length(x), ncol = length(y), byrow = TRUE), alpha = 1) %>%
  layout(title = 'My Gaussian mixture 2D function') %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1,y=1,z=0.5),
                      xaxis = list(range = c(-5.12,5.12)), 
                      yaxis = list(range = c(-5.12,5.12)),
                      zaxis = list(range = c(0,0.03))))