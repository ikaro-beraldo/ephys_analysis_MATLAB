ripple = reshape(nanmean(results.RDS.ripple,2),1,[]);
delta_power_mpfc = reshape(nanmean(results.bands.LFP1.Delta.NREM,2),1,[]);
delta_power_ca1 = reshape(nanmean(results.bands.LFP3.Delta.NREM,2),1,[]);

nan_values = isnan(ripple) | isnan(delta_power_ca1);
ripple(nan_values) = [];
delta_power_ca1(nan_values) = [];
delta_power_mpfc(nan_values) = [];
 

subplot(2,1,1)
scatter(ripple,delta_power_ca1)
subplot(2,1,2)
scatter(ripple,delta_power_mpfc)


[R, p] = corrcoef([ripple' delta_power_ca1'])
[R, p] = corrcoef([ripple' delta_power_mpfc'])

%% Check the ripple and delta throught time for each animal

rp_final = [];
del_final = [];

for i = 1:5
     %rp = reshape(permute(results.RDS.ripple(i,:,1:5),[3 2 1]),1,[]);
     %del = reshape(permute(results.bands.LFP3.Delta.NREM(i,:,1:5),[3 2 1]),1,[]);

    rp = permute(nanmean(results.RDS.ripple(i,:,1:7),2),[1 3 2]);
    del = permute(nanmean(results.bands.LFP1.Delta.NREM(i,:,1:7),2),[1 3 2]);

    nan_values = isnan(rp) | isnan(del);
    rp(nan_values) = [];
    del(nan_values) = [];

    [R, p] = corrcoef([rp' del']);

    subplot(5,1,i)
    yyaxis left
    plot(rp,'-.')
    yyaxis right
    plot(del,'-.')
    title(sprintf('R: %0.3f; p: %0.3f',R(1,2),p(1,2)))

    rp_final = [rp_final rp];
    del_final = [del_final del];
    
end

[R, p] = corrcoef([rp_final' del_final']);
figure
scatter(rp_final,del_final)
title(sprintf('R: %0.3f; p: %0.3f',R(1,2),p(1,2)))

%% 

rp = reshape(curve_fit_results.RDS.ripple,1,[]);
del = reshape(curve_fit_results.bands.LFP3.Delta.NREM,1,[]);

nan_values = isnan(rp) | isnan(del);
rp(nan_values) = [];
del(nan_values) = [];

[R, p] = corrcoef([rp' del'])

scatter(rp,del)

subplot(5,1,i)
yyaxis left
plot(rp,'-.')
yyaxis right
plot(del,'-.')
title(sprintf('R: %0.3f; p: %0.3f',R(1,2),p(1,2)))
