myX=X$Values
par(mfrow = c(3,3), col.axis = "white", col.lab = "white", tck = 0)
for(i in 1:7){
  x=myX[,i]
  plot(density(x, bw = 0.05), lwd = 2,
       col = "red", main = "Le density le yes")
}