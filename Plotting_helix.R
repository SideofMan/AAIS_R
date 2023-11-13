samps=as.data.frame(X$Values); colnames(samps) <- c("x","y","z")

library(plotly)

plot_ly(samps, x=~x, y=~y, z=~z, mode="markers", marker=list(size=2,color="blue",symbol="circle"))%>%
  add_markers() %>%
  layout(title = 'Helix 3D toy function') %>%
  layout(scene = list(aspectmode = "manual", aspectratio = list(x=1,y=1,z=1),
                      xaxis = list(range = c(-20,20)), 
                      yaxis = list(range = c(-20,20)),
                      zaxis = list(range = c(0,50))))