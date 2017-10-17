 
# lets work with iris dataset -------------------------------------

# head(iris)

library(scales)

palette(c(alpha("gold", .6), alpha("skyblue", .6), alpha("firebrick1", .6)))

# png(filename = "slideshow/images/irisplot.png", width = 11, height = 8.5, units = "in", res= 600)
        
par(mar=c(5,5,2,1),cex.axis=1, cex.lab=1.25,las=0,mgp=c(3,1,0))   

plot(Sepal.Length ~ Petal.Length, data=iris, bg=Species, pch=21, xlim=c(0,9),
     ylim=c(0,9), cex=2, ylab="Sepal Length (mm)", xlab= "Petal Length (mm)")

legend("bottomright", legend = unique(iris$Species), pt.bg=palette(), pch=21, inset = .02, cex=1)

# dev.off()