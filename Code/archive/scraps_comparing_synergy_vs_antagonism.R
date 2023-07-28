
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




# Check synergy vs antagonism vs partial antagonism

library(plotly)
param_matrix = as.matrix(cbind("a" = c(4,4), 
                               "b" = c(5,5), 
                               "c" = c(1.4,1),
                               "max_R" = 4,
                               "d" = 0))
cluster_assign_vec = c(1,1)
mix_pred = mix_function_generator(param_matrix, cluster_assign_vec)
mix_pred(true_conc_value =c(15, 10))

c_val = seq(0, 5, by=.05)
n_pix = length(c_val)
c_mat = expand.grid(c_val,c_val)
pred_effect = apply(c_mat, MARGIN=1, FUN = function(x) mix_pred(x) )
#pred_regular = pred_effect
# pred_synergy = pred_effect
# pred_agonist = pred_effect
mix_df = data.frame(cbind(c_mat, pred_effect))
names(mix_df) = c("c1", "c2", "effect")
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

