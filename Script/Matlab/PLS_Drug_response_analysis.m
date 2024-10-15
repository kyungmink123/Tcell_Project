% Drug response PLS-DA prediction model

% Preparation of input matrix
Proportion  = [];
y_variable = [];
feature = {};

% transpose input data (rows as features (ratio information), columns as patient)
After_t  = Proportion';

% retrieve data w/o naive T proportion
After_t_no_naive = After_t(1:8,:);

% median centering and autoscaling before PLS analysis

After_md = After_t(:,:) - median(After_t(:,:) , 2, 'omitnan');
After_t_no_naive_md = After_t_no_naive(:,:) - median(After_t_no_naive(:,:) , 2, 'omitnan');

% quantile normalization
After_md_qnorm = quantilenorm(After_md);
After_t_no_naive_md_qnorm = quantilenorm(After_t_no_naive_md);

% autoscaling
After_md_qnorm_auto = zscore(After_md_qnorm(:,:)')';
After_t_no_naive_md_qnorm_auto = zscore(After_t_no_naive_md_qnorm(:,:)')';
y_auto      = zscore(y_variable);

% perform PLS-DA
[m,ssq,p,q,w,t,u,b] = pls(After_md_qnorm_auto' , y_auto , 10 , 1);

% Non-responder   t(y_variable(:, 1)==1,2)
figure;plot3(t(y_variable(:, 1)==1,1),t(y_variable(:, 1)==1,2),t(y_variable(:, 1)==1,3) , '.','MarkerSize',20 , 'color',  [.7 .7 .7]);
hold on

% Responder
plot3(t(y_variable(:, 1)==0,1),t(y_variable(:, 1)==0,2),t(y_variable(:, 1)==0,3), '.','MarkerSize',20  , 'color', "#EDB120");
hold on

grid on
xlabel('LV1 (60.45%)')
ylabel('LV2 (16.38%)')
zlabel('LV3 (16.58%)')
legend('Non-responder' , 'Responder')
title('PLS-DA regression model') 

% VIP calculation and null distribution generation
[err,tr, rmse,cumpressa, yp] = sub_plsda_crossval(After_md_qnorm_auto', y_auto, 9, 1, 10);

%VIP calculation
After_vip=plsVIP(w(:,1:3),ssq(1:3,4));
[After_vip_sort,ind]=sort(After_vip,'descend');
figure;bar(After_vip_sort)
	axis([0  length(After_vip_sort)+1 0 max(After_vip_sort)+0.5])
	line([0  length(After_vip_sort)+1],[1 1],'color','k','linewidth',1)
	ylabel('VIP')
	sum(After_vip_sort>1)

figure;After_vip_null = sub_multiblock_vip_null(After_md_qnorm_auto',y_auto, 1000, 2);
	axis square;xlabel('VIP');ylabel('Frequency');title('empirical distribution of VIP');

figure ; bar(sort(After_vip_null.tmp1 , 'descend') );ylim([0,3.5])

% generation of accuracy plot
err_tot = [];

figure;plot([1:10], 1-err_tot);axis square; xlabel('No. of LVs');ylabel('Accuracy');
legend({'w/o naive','with naive'})

% Check VIP null
prctile(After_vip_null.tmp1, 99)
% ans = 1.8742

prctile(After_vip_null.tmp1, 95)
%ans = 1.5815

prctile(After_vip_null.tmp1, 90)
%ans = 1.4301

% Bulid PLS model after feature selection
[m,ssq,p,q,w,t,u,b] = pls(After_md_qnorm_auto(After_vip>1 ,:)', y_auto, 4);

% Non-responder   t(y_variable(:, 1)==1,2)
figure;plot3(t(y_variable(:, 1)==1,1),t(y_variable(:, 1)==1,2),t(y_variable(:, 1)==1,3) , '.','MarkerSize',20 , 'color', [.7 .7 .7]);
hold on

% Responder
plot3(t(y_variable(:, 1)==0,1),t(y_variable(:, 1)==0,2),t(y_variable(:, 1)==0,3), '.','MarkerSize',20  , 'color', "#EDB120");
hold on

grid on
xlabel('LV1 (15.22%)')
ylabel('LV2 (0.85%)')
zlabel('LV3 (1.80%)')
legend('Non-responder' , 'Responder')
title('Prediction model') 

% Boxplot generation of scaled proportion

tmp2 = After_md_qnorm_auto(7,y_variable(:,1) == 1);
tmp3 = After_md_qnorm_auto(7,y_variable(:,1) == 0);
groups = [ones(size(tmp2)) , 2*ones(size(tmp3))];

figure
boxplot([tmp2 , tmp3] , groups)
set(gca, 'XTickLabel', {'Non-res' , 'Res'})