
x = seq(-5,5, by=.01)
a=.6
b=1
c=1.7

hilly = function(x) a/(1+(b/x)^c)
hilly_reflect_xy =  function(x) b/((a/x)-1)^(1/c)

hilly_rotate = function(theta, cx, cy){
  
}

hilly_inverse = function(y){
  if(y<a) return(hilly_reflect_xy(y))
  if(y==a) return(0)
  return(-hilly_reflect_xy(2*a-y)-2*b)#-2*b
}
  

hilly_reflect_xy(.5999)
-hilly_reflect_xy(2*a-.62)



Rval = c(seq(0,1,by=.02))#, seq(.65,1,by=.02))

plot(Rval,sapply(Rval,hilly_inverse))
points(Rval,sapply(Rval,hilly_reflect_xy),col=2)


cbind(sapply(Rval,hilly_inverse), sapply(Rval,hilly_reflect_xy))
#png("Hill_inverse_asym.png", width  = 6, height=4, units = "in", res = 400)
x = seq(-5,5, by=.05)
plot(x, hilly(x), xlim = c(-7,5), 
     ylim =c(min(hilly(x), na.rm = T), 
             max(-hilly(x)+2*a, na.rm = T)),
     main = "Invertible Hill Curve, Asymetric")
points(-x-2*b, -hilly(x)+2*a, col=4)
plot(-2*a-hilly(-x-2*b),x+2*a)
plot(-2*a+hilly(-x-2*b),-x)
points(-x, rev(2*a-rev(hilly(x))))
abline(b=1, a=a+b)

text(2.5,1,"Reflection: y = x+a+theta")
text(-4, .2, "a=0.6, theta=1, beta=1.7")
#dev.off()
x_reflect = seq(0,1, by=.01)
plot(x, hilly_reflect_xy(x),col=2)


#inverting


