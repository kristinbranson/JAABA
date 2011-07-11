% confusion_matrix(1,1) = fraction of true 1 classified as 1
% confusion_matrix(1,2) = fraction of true 1 classified as 2
% confusion_matrix(1,3) = fraction of true 1 classified as unknown
% confusion_matrix(2,1) = fraction of true 2 classified as 1
% confusion_matrix(2,2) = fraction of true 2 classified as 2
% confusion_matrix(2,3) = fraction of true 2 classified as unknown
function confusion_matrix = ComputeConfusionMatrix(hs,diff_logprob,maxdiff_logprob1,mindiff_logprob2)

confusion_matrix = nan(2,3);
hsPr1 = diff_logprob < maxdiff_logprob1;
hsPr2 = diff_logprob > mindiff_logprob2;
n1 = nnz(hs==1);
n2 = nnz(hs==2);

confusion_matrix(1,1) = nnz(hsPr1(hs==1)) / n1;
confusion_matrix(1,2) = nnz(hsPr2(hs==1)) / n1;
confusion_matrix(1,3) = 1 - confusion_matrix(1,1) - confusion_matrix(1,2);

confusion_matrix(2,1) = nnz(hsPr1(hs==2)) / n2;
confusion_matrix(2,2) = nnz(hsPr2(hs==2)) / n2;
confusion_matrix(2,3) = 1 - confusion_matrix(2,1) - confusion_matrix(2,2);
