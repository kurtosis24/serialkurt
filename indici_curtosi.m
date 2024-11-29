function [K2_sn_m1,K2_dx_m1,K2_m1,K2_sn_me,K2_dx_me,K2_me,K1_sn_m1,K1_dx_m1,K1_m1,K1_sn_me,K1_dx_me,K1_me,geary,beta2,asimm,beta1] = indici_curtosi(dati)

% calcola tutti gli indici di Zenga (anche unilaterali), geary, beta_2.
% calcola anche l'indice di asimmetria e beta_1.


N = size(dati, 1);
m1 = mean(dati);
me = median(dati);
var = mean((dati-m1).^2);
scost_m1 = mean(abs(dati-m1));
scost_me = mean(abs(dati-me));
mu3 = mean((dati-m1).^3);
mu4 = mean((dati-m1).^4);

split_m1 = dati <= m1;
split_me = dati <= me;
dati_sn_m1 = dati(split_m1);
dati_dx_m1 = dati(not(split_m1));
dati_sn_me = dati(split_me);
dati_dx_me = dati(not(split_me));

perc_sn_m1 = size(dati_sn_m1, 1)/N;
perc_sn_me = size(dati_sn_me, 1)/N; %qui non uso 1/2, per verifica

scost_sn_m1 = mean(m1-dati_sn_m1);
scost_dx_m1 = mean(dati_dx_m1-m1);
scost_sn_me = mean(me-dati_sn_me);
scost_dx_me = mean(dati_dx_me-me);

delta_sn_m1 = deltaf(dati_sn_m1);
delta_dx_m1 = deltaf(dati_dx_m1);
delta_sn_me = deltaf(dati_sn_me);
delta_dx_me = deltaf(dati_dx_me);

mom2_sn_m1 = mean((m1-dati_sn_m1).^2);
mom2_dx_m1 = mean((dati_dx_m1-m1).^2);
mom2_sn_me = mean((me-dati_sn_me).^2);
mom2_dx_me = mean((dati_dx_me-me).^2);

K2_sn_m1 = perc_sn_m1 * delta_sn_m1/(2 * scost_sn_m1);
K2_dx_m1 = (1-perc_sn_m1) * delta_dx_m1/(2 * scost_dx_m1);
K2_sn_me = perc_sn_me * delta_sn_me/(2 * scost_sn_me);
K2_dx_me = (1-perc_sn_me) * delta_dx_me/(2 * scost_dx_me);

K1_sn_m1 = perc_sn_m1 * (1 - scost_sn_m1^2/mom2_sn_m1);
K1_dx_m1 = (1-perc_sn_m1) * (1 - scost_dx_m1^2/mom2_dx_m1);
K1_sn_me = perc_sn_me * (1 - scost_sn_me^2/mom2_sn_me);
K1_dx_me = (1-perc_sn_me) * (1 - scost_dx_me^2/mom2_dx_me);

K2_m1 = K2_sn_m1 + K2_dx_m1;
K2_me = K2_sn_me + K2_dx_me;
K1_m1 = K1_sn_m1 + K1_dx_m1;
K1_me = K1_sn_me + K1_dx_me;

geary = scost_m1/sqrt(var);
beta2 = mu4/(var^2);

asimm = (m1-me)/scost_me;  %attenzione alla definizione
beta1 = mu3/(var^(3/2));


end









