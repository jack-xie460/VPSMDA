function [auc] = roc_1(deci,label_y,colour)
[threshold,ind] = sort(deci,'descend');
roc_y = label_y(ind);
%cumsum函数通常用于计算矩阵中按行或列进行累加求和,cumsum累加函数默认是按照列进行计算的
stack_x = cumsum(roc_y == 0)/sum(roc_y == 0);
stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);
auc=sum((stack_x(2:length(roc_y))-stack_x(1:length(roc_y)-1)).*stack_y(2:length(roc_y)));
plot(stack_x,stack_y,colour);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
% hold on;
title(['ROC curve of (AUC = ' num2str(auc) ' )']);
end
