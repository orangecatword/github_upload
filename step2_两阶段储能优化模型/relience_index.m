% 创建一个,根据最终的储能配置,得到故障恢复率,网损率,电压偏差

f = relience_index();
[f2,Psum_loss2,Psum_load2,U2] = function2(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
[f3,Psum_loss3,Psum_load3,U3] = function3(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
[f4,Psum_loss4,Psum_load4,U4] = function4(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
[f5,Psum_loss5,Psum_load5,U5] = function5(Lc_pes1,Lc_pes2,Lc_pes3, Ees_max1,Ees_max2,Ees_max3);
[U,Pload_total] = function_standard;
% 故障恢复率
rload = 1- (Psum_load2+Psum_load3+Psum_load4+Psum_load5)/(4*Pload_total);
% 网损率
rloss = 1- (Psum_loss2+Psum_loss3+Psum_loss4+Psum_loss5)/(4*Pload_total);
% 电压偏差
deltaU = sum(sum(abs(U-U2)+abs(U-U3)+abs(U-U4)+abs(U-U5)))/4;