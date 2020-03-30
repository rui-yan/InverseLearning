clc;close all;clear all;

a_t = 2.4947e3; b_t = 765; g_t = 2.5897e3; o_t = 858;

% generate data from forward model
[d_true, K_true, G_true] = forward_model(a_t, b_t, g_t, o_t);


pt=10;

dis_a=linspace(-2300, 2300, pt);
dis_b=linspace(-700, 700, pt);
dis_g=linspace(-2400, 2400, pt);
dis_o=linspace(-800, 800, pt);

dis_a_all=repelem(dis_a,pt*pt*pt);
dis_b_all=repmat(repelem(dis_b,pt*pt), 1, pt);
dis_g_all=repmat(repelem(dis_g,pt), 1, pt*pt);
dis_o_all=repmat(dis_o,1,pt*pt*pt);
dis_all=[dis_a_all;dis_b_all;dis_g_all;dis_o_all];

F=zeros(1,size(dis_all,2));
for i=1:size(dis_all,2)
    fprintf('current is %d',i)
    a=a_t+dis_all(1,i);
    b=b_t+dis_all(2,i);
    g=g_t+dis_all(3,i);
    o=o_t+dis_all(4,i);
    [this_d1, ~, ~]  = forward_model(a, b, g, o);
    F(i) = norm(d_true - this_d1)^2;
    fprintf('\n')
end
dat=[dis_all;F];
dat=dat';  % pt samples, 5 features