# calc performance criteria

tau_same <- subset(result, tau == "same")
tau_same$rb_int = (tau_same$int_mean - 100)/100
tau_same$rb_slope = (tau_same$slope_mean - 10)/10
tau_same$rb_sigma_y = (tau_same$sigma_y_mean - 10)/10
tau_same$rb_sigma_u2 = (tau_same$sigma_u2_mean - 10)/10
tau_same$rb_sigma_u2_1 = (tau_same$sigma_u2_1_mean - 10)/10
tau_same$rb_sigma_u2_2 = (tau_same$sigma_u2_2_mean - 10)/10
tau_same$rb_sigma_u2_3 = (tau_same$sigma_u2_3_mean - 10)/10
tau_same$rb_sigma_u2_4 = (tau_same$sigma_u2_4_mean - 10)/10

tau_diff <- subset(result, tau == "different")
tau_diff$rb_int = (tau_diff$int_mean - 100)/100
tau_diff$rb_slope = (tau_diff$slope_mean - 10)/10
tau_diff$rb_sigma_y = (tau_diff$sigma_y_mean - 10)/10
tau_diff$rb_sigma_u2 = (tau_diff$sigma_u2_mean - 10)/10
tau_diff$rb_sigma_u2_1 = (tau_diff$sigma_u2_1_mean - 10)/10
tau_diff$rb_sigma_u2_2 = (tau_diff$sigma_u2_2_mean - 12)/12
tau_diff$rb_sigma_u2_3 = (tau_diff$sigma_u2_3_mean - 14)/14
tau_diff$rb_sigma_u2_4 = (tau_diff$sigma_u2_4_mean - 16)/16

rbind(tau_same, tau_diff)