
meanmu_p = mean(reshape([incdeg_conv.muhat_PCE],size(degrees_conv)));
meanvar_p = mean(reshape([incdeg_conv.varhat_PCE],size(degrees_conv)));
meancv_p = mean(reshape([incdeg_conv.cvhat_PCE],size(degrees_conv)));
meant_p = mean(reshape([incdeg_conv.t],size(degrees_conv)));


meanmu_n = mean(reshape([incn_conv(:,2:end).muhat_PCE],size(oversampling_conv(:,2:end))));
meanvar_n = mean(reshape([incn_conv(:,2:end).varhat_PCE],size(oversampling_conv(:,2:end))));
meancv_n = mean(reshape([incn_conv(:,2:end).cvhat_PCE],size(oversampling_conv(:,2:end))));
meant_n = mean(reshape([incn_conv(:,2:end).t],size(oversampling_conv(:,2:end))));