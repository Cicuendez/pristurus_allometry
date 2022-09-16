# Example of coefficients in ANCOVA
n <- 15
x <- 1:n
y <- x * 1 + rnorm(n, 0, .001)
y2 <- x * .5+ rnorm(n, 0, .001)
y3 <- x * 1.5+ rnorm(n, 0, .001)
plot(x,y, pch=22,bg="green")
 points(x,y2, pch=22, bg="black")
 points(x,y3, pch=22, bg = 'red')

 x=rep(x,3)
 y=c(y,y2,y3)
 gp=gl(3,n) 
 
 mydat <- data.frame(x = x, y = y, gp = gp)

 fit <- lm(y~x*gp, data = mydat)  
 coef(fit)
 c(coef(fit)[1]+coef(fit)[3], coef(fit)[2]+coef(fit)[5])
 c(coef(fit)[1]+coef(fit)[4], coef(fit)[2]+coef(fit)[6])  
 
 by(mydat,gp, function(x) lm(y~x, data = x))