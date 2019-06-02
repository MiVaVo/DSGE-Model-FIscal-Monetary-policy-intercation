    var c,k,h,b,p,r_d,r,w,lambda,g,m,gov;

varexo eps_lambda, eps_g ;

parameters
theta
beta
delta
g_bar
tau_c
gamma_lambda
gamma_g
gamma_gov
r_bar
r_bar_d
w_bar
c_bar
k_bar
h_bar
b_bar
p_bar
%m_ad
gov_bar
m_bar
w_bar_f
n_ad
al
ksi
v_ad;

theta=0.36;
beta=0.99;
B=-2.5805;
delta=0.025;
g_bar=1.2; %can vary
tau_c=0.18;
%gov_bar=0;
gamma_lambda=0.95;
gamma_g=0.95;
gamma_gov=0.95;
w_bar_f=0.03475;
al=0.03;
m_bar=1;
ksi=0.1;

r_bar=1/beta-(1-delta);
r_bar_d=g_bar/beta-1;
w_bar=(1-theta)*(theta/r_bar)^(theta/(1-theta));
p_bar=-(1+tau_c)*g_bar*B/(beta*w_bar);
c_bar=1/p_bar;
%m_ad=1/p_bar-1/(p_bar*g_bar)-r_bar_d/(p_bar*g_bar);
%b_bar=-tau_c*c_bar/m_ad;
%k_bar=(m_ad+1/p_bar)*theta/(r_bar-theta*delta);

n_ad=((1-theta)*r_bar/theta+(r_bar-delta)-al);
k_bar=(1/p_bar-tau_c)/n_ad;
h_bar=(1-theta)*r_bar*k_bar/(theta*w_bar);
%y_bar=k_bar^(theta)*h_bar^(1-theta);
v_ad=1/p_bar-(1+r_bar_d)/(p_bar*g_bar);
b_bar=(al*k_bar-tau_c*c_bar)/v_ad;
gov_bar=k_bar*al;



model(linear);
(1-ksi)*tau_c*c_bar*c+(b_bar*(b-p)/p_bar)=(b_bar*(b(-1)-p-g)/(p_bar*g_bar))+(b_bar*r_bar_d*(b(-1)-p-g+r_d(-1))/(p_bar*g_bar))+gov_bar*gov; 
%tau_c*c_bar*c+(b_bar*(b-p)/p_bar)=(b_bar*(b(-1)-p-g)/(p_bar*g_bar))+(b_bar*r_bar_d*(b(-1)-p-g+r_d(-1))/(p_bar*g_bar))+gov_bar*gov;%+
w(+1)-w+w_bar_f=beta*(r_bar+r_bar*r(+1)); %+0.03475
w+p-c(+1)-p(+1)-g(+1)=0; %+
g_bar*(w(+1)+p(+1)+g(+1)-p-w)=beta*r_bar_d*r_d;
0=w-lambda-theta*k+theta*h; %+
0=r-lambda+(1-theta)*k-(1-theta)*h;%+
p+c=0;
m=g+m(-1);
gov_bar*gov=al*(k_bar*k(+1));
k(+1)*k_bar+m_bar*(m-p)/(p_bar)+b_bar*(b-p)/p_bar=w_bar*h_bar*(w+h)+r_bar*k_bar*(r+k)+(1-delta)*k_bar*k+b_bar*(b(-1)-p-g)/(p_bar*g_bar)+b_bar*r_bar_d*(b(-1)+r_d(-1)-p-g)/(p_bar*g_bar);


lambda=gamma_lambda*lambda(-1)+eps_lambda;
g=gamma_g*g(-1)+eps_g;
%b=gamma_b*b(-1)+eps_b;
%gov=gamma_gov*gov(-1)+eps_gov;

end; 
resid(1);
check;
shocks;
var eps_g; stderr 0.01;
end;
stoch_simul(periods=300);
/* 
ksi_par=0.1:0.05:0.7;

for i=1:length(ksi_par)
ksi=ksi_par(i);
stoch_simul(periods=300,noprint,nograph);
c_eps_mat(:,i)=c_eps_g;
stoch_simul(periods=500);
end
surf(c_eps_mat(1:10,:));xlabel('ksi_par');ylabel('time');zlabel('c_eps_g');
*/
