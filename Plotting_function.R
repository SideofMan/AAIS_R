# myfunction can be 'r' for Rastrigin, 'g' for Gaussian mixture.
# 'q' for quadratic, 'g1' for Gaussian 1D mixture, 'o' for outerproduct

library(ggplot2)
library("plotly")

myfunction='g1'
myfunction='g'
myfunction='o'

if(myfunction=='g1'){
  x=matrix(linspace(-10,10,10001), ncol=1)
  y=exp(my_gaussian_mixture_1d(x,1)[[1]])
  
  # Plot the curve
  my_df <- data.frame(x = x, y = y)
  ggplot(data = my_df, aes(x,y)) +
    geom_line() +
    scale_x_continuous(limits = c(-10,10), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1.5), expand = c(0,0)) +
    ggtitle("2D Plot of Gaussian 1D Mixture function") +
    theme(plot.title = element_text(hjust = 0.5))
} else if(myfunction=='g'){
  # x=matrix(linspace(-10,10,201), ncol=1)
  # y=matrix(linspace(-10,10,201), ncol=1)
  # XY=meshgrid(x,y); X=XY$X; Y=XY$Y
  # Z=zeros(dim(X)[1],dim(Y)[2])
  
  # for(i in 1:dim(X)[1]){
  #   for(j in 1:dim(Y)[2]){
  #     Z[i,j]=exp(my_gaussian_mixture(matrix(c(X[i,j],Y[i,j]), nrow = 1), 1)[[1]])
  #   }
  # }
  
  x=seq(-10, 10, length = 201)
  y=seq(-10, 10, length = 201)
  grid=expand.grid(x=x,y=y)
  Z=with(grid, exp(my_gaussian_mixture(cbind(grid$x, grid$y), 1)[[1]]))
  
  plot_ly()%>%
    add_surface(x = x, y = y, z = matrix(Z, nrow = length(x), ncol = length(y), byrow = TRUE), alpha = 1) %>%
    layout(title = 'My Gaussian mixture 2D function') %>%
    layout(scene = list(aspectmode = "manual", aspectratio = list(x=1,y=1,z=0.5),
                        xaxis = list(range = c(-10,10)), 
                        yaxis = list(range = c(-10,10)),
                        zaxis = list(range = c(0,0.2))))
} else if(myfunction=='o'){
  domains=list(linspace(-50,50,1001),
               linspace(-20,20,1001),
               linspace(-20,20,1001),
               linspace(-20,20,1001),
               linspace(-20,20,1001),
               linspace(-50,50,1001),
               linspace(-15,15,1001))
  plots=list()
  functions=list(function(t) (3/5*dgamma(t+10,2,scale = 3)+2/5*dgamma(10-t,2,scale = 5)),
                 function(t) (3/4*dsn(t,3,1,5)+1/4*dsn(t,-3,3,-6)),
                 function(t) dt(t,4),
                 function(t) (1/2*dbeta(t+3,3,3)+1/2*dnorm(t,0,sqrt(1))),
                 function(t) (1/2*dexp(t,1)+1/2*dexp(-t,1)),
                 function(t) dsn(t,0,8,-3),
                 function(t) (1/8*dnorm(t,-10,sqrt(.1))+1/4*dnorm(t,0,sqrt(.15))+5/8*dnorm(t,7,sqrt(.2))))
  
  par(mfrow = c(3,3), col.axis = "white", col.lab = "white", tck = 0)
  for(i in 1:7){
    my_f=functions[[i]]
    x=domains[[i]]
    plot(x,my_f(x), typ='l')
    # plots[[i]]=plot(x,my_f(x))
  }
}