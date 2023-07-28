

# Check synergy vs antagonism vs partial antagonism

library(plotly, ggplot2)


isobole_df = data.frame()
for(a_val in c(1,2)){
  for(slope_val in c(.5,1,1.5)){
    # create a parameter input matrix 
    param_matrix = as.matrix(cbind("a" = c(a_val,2), 
                                   "b" = c(.8,1.5), 
                                   "c" = c(slope_val,1),
                                   "max_R" = 5,
                                   "d" = 0))
    # specify both chems in cluster 1 of 1
    cluster_assign_vec = c(1,1)
    # create a calculator to predict response given concentration
    mix_pred = mix_function_generator(param_matrix, cluster_assign_vec)
    # create matrix of concentration values
    c_val = seq(.1, 5, by=.1)
    n_pix = length(c_val)
    c_mat = expand.grid(c_val,c_val)
    # apply calculator to concentrations
    pred_effect = apply(c_mat, MARGIN=1, FUN = function(x) mix_pred(x) )
    mix_df = data.frame(cbind(c_mat, pred_effect))
    names(mix_df) = c("C1", "C2", "effect")
    isobole_df = rbind(isobole_df, cbind(mix_df, "a"=a_val, "slope" = slope_val))
  }
}
# 2d contour plot
#pdf(file='Isoboles_RGCA.pdf',width=8, height = 6)
ggplot(isobole_df, aes(x=C1, y=C2, z=effect)) + 
  geom_contour(bins = 20) + 
  facet_grid(cols = vars(a), rows = vars(slope))
#dev.off()


# 3d plot with 2d contour plot on floor
#mix_df$effect[mix_df$effect< -1000]=NA
plot_ly(x = (c_val),y=(c_val),z=matrix(mix_df$effect, nrow = n_pix)) %>% 
  add_surface(contours= list( z = list(
    show=TRUE,
    usecolormap=TRUE,
    highlightcolor="#ff0000",
    project=list(z=TRUE)
  ), zlim = c(0,1)
  )) %>% layout(scene = list(
    camera=list(eye = list(x=2, y=.5, z=.7))))
#%>% add_surface(z = matrix(pred_regular, nrow = n_pix), opacity = .4) 
#%>% add_surface(z = matrix(pred_agonist, nrow = n_pix), opacity = .4)
beepr::beep()

plot_ly(z=matrix(mix_df$effect, nrow = n_pix), type = "contour")



# Check uniqueness of solutions 

param_matrix = as.matrix(cbind("a" = c(2,-2), 
                               "b" = c(1,1.0), 
                               "c" = c(2.5, 1.9),
                               "max_R" = 10,
                               "d" = 0))
hill_inverse_list = apply(param_matrix, 
                          MARGIN=1, 
                          function(x) do.call(hill_invs_factry,as.list(x)))

myfun = function(r) 1/hill_inverse_list[[1]](r) +1/hill_inverse_list[[2]](r)
myseq = seq(-4, 4, by=.01)
plot(myseq, 1/sapply(myseq, hill_inverse_list[[1]]), ylim = c(-10,10))
points(myseq, 1/sapply(myseq, hill_inverse_list[[2]]), col=2)
points(myseq, sapply(myseq, myfun), col=3); abline(h=1)

#computed effect
myseq[-401]
plot(myseq[-401], (-1*myseq[-401]/.6+1*1/.7)/(1+myseq[-401]/.6 + 1/.7), ylim = c(-10,20))

mix_pred = mix_function_generator(param_matrix, c(1,1))
xseq = seq(0,6, by=.01)
plot(xseq,  sapply(xseq, FUN = function(xs) mix_pred(true_conc_value =c(1, xs))))
abline(v = 1.4)

myseq0 =  seq(.01, 10, by=.01)
plot(myseq0, -1/(1+(.7/myseq0)^(1.2)),  ylim =c(-1, 2))
points(myseq0, 2/(1+(2/myseq0)^(.9)))



active_conc= c(3,3)#c(1.4, 1.4)
GCA_function = eff_response_opt(hill_inverse_list, active_conc)
GCA_function_neg = eff_response_opt(hill_inverse_list, active_conc, interval_sign =-1)
rseq = seq(-2.2, -1.5, by = 0.001)
plot(exp(rseq), sapply(rseq, FUN = function(rs) GCA_function(rs)), type = "l", xlim = c(-.3, .3), log="y", ylim = c(1e-3, 1))
lines(-exp(rseq), sapply(rseq, FUN = function(rs) GCA_function_neg(rs)), lwd=2)
iter = 1
for(xs in c(1.45, 1.5, 1.7, 1.8,2)){#1.42, 1.45, 1.484,
  iter = iter+1
  active_conc= c(1.4, xs)
  GCA_function = eff_response_opt(hill_inverse_list, active_conc)
  GCA_function_neg = eff_response_opt(hill_inverse_list, active_conc, interval_sign =-1)
  lines(exp(rseq), sapply(rseq, FUN = function(rs) GCA_function(rs)), col=iter)
  
  lines(-exp(rseq), sapply(rseq, FUN = function(rs) GCA_function_neg(rs)), col=iter, lwd=2)
}
optim_res = optimize(GCA_function,interval= c(-100,10), tol = .Machine$double.eps^0.5)
optim_res_neg = optimize(GCA_function_neg,interval= c(-100,10), tol = .Machine$double.eps^0.5)








r_guess = seq(.1, 4, by=.01)
active_conc = c(0,26)
GCA_function = eff_response_opt(hill_inverse_list, active_conc)
plot(r_guess,sapply(r_guess,GCA_function), type = "l", ylim = c(0, 10))
chem2_val = seq(1, 25, by=1)
for(acc in 1:length(chem2_val)){#, 5, 2
  active_conc = c(chem2_val[acc],10)
  GCA_function = eff_response_opt(hill_inverse_list, active_conc)
  points(r_guess,sapply(r_guess,GCA_function), type = "l", col=(1+acc))
}
legend("top", col=c(1:5), 
       legend = c("(0,15)", "(0.01, 15)", "(5,15)"), 
       lty = c(1,1,1), horiz = T)

# Plot some inverse hill functions with partial antagonists
hill_function = function(a,b,c,conc) a/(1+(b/conc)^c)
invs_111 = hill_invs_factry(.7, 1, 1, max_R = 1)

invs_514 = hill_invs_factry(.7, 1, 1.4, max_R = 5)
test_y = seq(-1,5,by=.01)
test_x = seq(-10,10,by=.1)
pred_x=sapply(test_y, invs_514)
plot(test_x,hill_function(.7,1,1.4, test_x), xlim = c(-10,10), ylim = c(-2 ,3), type = "l",
     xlab = "Concentration", ylab = "Effect", main = "General Partial Agonist")
points(pred_x, test_y, cex=.8, col=4)
lines(test_x,hill_function(.7,1,1.4, test_x), lwd=2)
abline(h = .9, lty = 2)
points(invs_514(.9), .9, col=2, lwd =2)
abline(h = .5, lty = 2)
points(invs_514(.5), .5, col=3, lwd =2)
#points(pred_x, test_y, col=3)





pos_idx = which(pred_x>0)
neg_idx = which(pred_x<=0)
plot(c(-log(-pred_x[neg_idx]),log(pred_x[pos_idx])), 
     c(test_y[neg_idx], test_y[pos_idx]))


curve_GCA_data_quantiles = t(apply(curve_data_GCA, MARGIN = 2,
                                   FUN = function(x) quantile(x, probs = c(0.05, .5, .95))))
GCA_df = data.frame(Cx_axis_values, resp_y_values, curve_GCA_data_quantiles)
names(GCA_df) = c("Cx", "resp_y", "GCA_low", "GCA_med", "GCA_up")
hill_pred_plot = hill_pred_plot+
  geom_line(data=GCA_df, aes(Cx, GCA_med,colour = "GCA"))+
  geom_ribbon(data=GCA_df, aes(ymin = GCA_low, ymax = GCA_up, fill= "GCA"),
              alpha=0.1, #fill = "red",
              color = "black", linetype = "dotted")+
  scale_colour_manual(name = "Legend",values = c("RE+DP" = "black", "GCA" ="red"))+
  scale_fill_manual(name = "Legend", values = c("RE+DP" = "green", "GCA" ="red"))







### Other stuff? ###

x = c(seq(0, 3, by=.01))

cov_mat = Matern(rdist(x), range = 1, smoothness = 1.5)
Q = solve(cov_mat)
rx = rnorm(length(x))
y_sam = t(chol(cov_mat))%*%rx
plot(x,y_sam)
abline(h = mean(y_sam), col=2)
mean_alt =sum(Q%*%y_sam)/sum(Q)
abline(h =mean_alt, col=3)
sum(rx)

plot(x, Q%*%y_sam)

plot(hist(y_sam))
abline(v = mean_alt, col=3)
abline(v = mean(y_sam), col=2)



