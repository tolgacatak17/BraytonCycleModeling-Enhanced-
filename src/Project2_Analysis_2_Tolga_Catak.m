clc;
clear all;
close all;

R_uni = 8.314; %kJ/kmol*K
P_ref = 101.325; %kPa
T_ref = 298; % K
T_0 = 300; %K
P_0 = 100; %kPa

n_c = 0.8;
n_t = 0.9;

M_CO2 = 44.01; %kg/kmol
M_H2O = 18.02; %kg/kmol
M_N2 = 28.01; %kg/kmol
M_O2 = 32.00; %kg/kmol
M_CH4 = 16.04; %kg/kmol

R_CO2 = R_uni/M_CO2; %kJ/kg*K
R_H2O = R_uni/M_H2O; %kJ/kg*K
R_N2 = R_uni/M_N2; %kJ/kg*K
R_O2 = R_uni/M_O2; %kJ/kg*K
R_CH4 = R_uni/M_CH4; %kJ/kg*K

num_state=11;

Temp_max=700:50:1650;

num_iteration=length(Temp_max);

T=zeros(num_iteration,num_state);
P=zeros(num_iteration,num_state);
rel_hum=zeros(1,num_state);



h_total=zeros(num_iteration,num_state);
h_CO2=zeros(num_iteration,num_state);
h_H2O=zeros(num_iteration,num_state);
h_N2=zeros(num_iteration,num_state);
h_O2=zeros(num_iteration,num_state);
h_CH4=zeros(num_iteration,num_state);

P_sat=zeros(num_iteration,num_state);
P_CO2=zeros(num_iteration,num_state);
P_H2O=zeros(num_iteration,num_state);
P_N2=zeros(num_iteration,num_state);
P_O2=zeros(num_iteration,num_state);
P_CH4=zeros(num_iteration,num_state);

y_CO2=zeros(num_iteration,num_state);
y_H2O=zeros(num_iteration,num_state);
y_N2=zeros(num_iteration,num_state);
y_O2=zeros(num_iteration,num_state);
y_CH4=zeros(num_iteration,num_state);

s_total=zeros(num_iteration,num_state);
s_CO2=zeros(num_iteration,num_state);
s_H2O=zeros(num_iteration,num_state);
s_N2=zeros(num_iteration,num_state);
s_O2=zeros(num_iteration,num_state);
s_CH4=zeros(num_iteration,num_state);

ef_CO2=zeros(num_iteration,num_state);
ef_H2O=zeros(num_iteration,num_state);
ef_N2=zeros(num_iteration,num_state);
ef_O2=zeros(num_iteration,num_state);
ef_CH4=zeros(num_iteration,num_state);
ef_total=zeros(num_iteration,num_state);

ndot=zeros(num_iteration,num_state);

for j=1:length(Temp_max)
    for i=1:num_state
        %% STATE 1
        if i==1
    
            T(j,i) = T_0;
            P(j,i) = P_0;
            
            rel_hum(j,i) = 0.6;
            
            P_sat(j,i)=P_sat_func(T(j,i));
            P_CO2(j,i) = 0;
            P_H2O(j,i) = P_sat(j,i)*rel_hum(j,i);
            P_N2(j,i) = (P(j,i)-P_H2O(j,i))*3.76/4.76;
            P_O2(j,i) = (P(j,i)-P_H2O(j,i))*1/4.76;
            P_CH4(j,i) = 0;
            
            y_CO2(j,i)=0;
            y_H2O(j,i)=P_H2O(j,i)/P(j,i);
            y_N2(j,i)=P_N2(j,i)/P(j,i);
            y_O2(j,i)=P_O2(j,i)/P(j,i);
            y_CH4(j,i)=0;
            
            h_CO2(j,i) = 0;
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = h_O2_func(R_uni,T(j,i));
            h_CH4(j,i) = 0;
            
            s_CO2(j,i) = 0;
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = s_O2_func(R_uni,P(j,i),T(j,i),y_O2(j,i));
            s_CH4(j,i) = 0;
            
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
    
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = ef_O2_func(h_O2(j,i),s_O2(j,i),y_O2(j,i));
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
        %% STATE 2s
        if i==2
            
            comp_ratio = 4;
            P(j,i) = comp_ratio*P(j,i-1);
            
            P_CO2(j,i) = 0;
            P_H2O(j,i) =comp_ratio*P_H2O(j,i-1);
            P_N2(j,i) = comp_ratio*P_N2(j,i-1);
            P_O2(j,i) = comp_ratio*P_O2(j,i-1);
            P_CH4(j,i) = 0;
            
            y_CO2(j,i)=0;
            y_H2O(j,i)=P_H2O(j,i)/P(j,i);
            y_N2(j,i)=P_N2(j,i)/P(j,i);
            y_O2(j,i)=P_O2(j,i)/P(j,i);
            y_CH4(j,i)=0;
            
            s_total(j,i)=s_total(j,i-1);
            
            
            T(j,i) = T_find_s(s_total(j,i),P(j,i),y_CO2(j,i),y_H2O(j,i),y_N2(j,i),y_O2(j,i),y_CH4(j,i));
            
            h_CO2(j,i) = 0;
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = h_O2_func(R_uni,T(j,i));
            h_CH4(j,i) = 0;
            
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
            
            s_CO2(j,i) = 0;
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = s_O2_func(R_uni,P(j,i),T(j,i),y_O2(j,i));
            s_CH4(j,i) = 0;
            
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = ef_O2_func(h_O2(j,i),s_O2(j,i),y_O2(j,i));
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
        %% STATE 2
        if i==3
            
            h_total(j,i) = (h_total(j,i-1)-h_total(j,i-2)+h_total(j,i-2)*n_c)/n_c;
            
            P(j,i)=P(j,i-1);
            P_CO2(j,i) = P_CO2(j,i-1);
            P_H2O(j,i) =P_H2O(j,i-1);
            P_N2(j,i) = P_N2(j,i-1);
            P_O2(j,i) = P_O2(j,i-1);
            P_CH4(j,i) = P_CH4(j,i-1);
            
            y_CO2(j,i)=y_CO2(j,i-1);
            y_H2O(j,i)=y_H2O(j,i-1);
            y_N2(j,i)=y_N2(j,i-1);
            y_O2(j,i)=y_O2(j,i-1);
            y_CH4(j,i)=y_CH4(j,i-1);
            
            T(j,i) = T_find_h(h_total(j,i),y_CO2(j,i),y_H2O(j,i),y_N2(j,i),y_O2(j,i),y_CH4(j,i));
    
            h_CO2(j,i) = 0;
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = h_O2_func(R_uni,T(j,i));
            h_CH4(j,i) = 0;
            
            s_CO2(j,i) = 0;
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = s_O2_func(R_uni,P(j,i),T(j,i),y_O2(j,i));
            s_CH4(j,i) = 0;
    
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
    
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = ef_O2_func(h_O2(j,i),s_O2(j,i),y_O2(j,i));
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
        %% STATE 2L
        if i==4
            
            T(j,i) = 280; % K
            
            cp_H2O = 4.179*M_H2O; % kj/kmol*K
            
            P(j,i)=P(j,i-1);
            P_CO2(j,i) = 0;
            P_H2O(j,i) = P(j,i);
            P_N2(j,i) = 0;
            P_O2(j,i) = 0;
            P_CH4(j,i) = 0;
            
            y_CO2(j,i)=0;
            y_H2O(j,i)=1;
            y_N2(j,i)=0;
            y_O2(j,i)=0;
            y_CH4(j,i)=0;
     
            
            h_formation_liquid = -285830; % kJ/kmol
            
            h_CO2(j,i) = 0;
            h_H2O(j,i) = cp_H2O*(T(j,i)-T_ref) + h_formation_liquid; 
            h_N2(j,i) = 0;
            h_O2(j,i) = 0;
            h_CH4(j,i) = 0;
            
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
            
            s_ref_liquid = 69.95; % kJ/kmol*K
    
            s_CO2(j,i) = 0;
            s_H2O(j,i) = cp_H2O*log(T(j,i)/T_ref)-R_uni*log(y_H2O(j,i)*P_H2O(j,i)/P_ref)+s_ref_liquid;
            s_N2(j,i) = 0;
            s_O2(j,i) = 0;
            s_CH4(j,i) = 0;
    
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
            
            h0_liquid = cp_H2O*(T_0-T_ref) + h_formation_liquid;
            s0_liquid = cp_H2O*log(T_0/T_ref)-R_uni*log(y_H2O(j,i)*P_0/P_ref)+s_ref_liquid;
            
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = h_H2O(j,i)-h0_liquid-T_0*(s_H2O(j,i)-s0_liquid)+900;
            ef_N2(j,i) = 0;
            ef_O2(j,i) = 0;
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
        %% STATE 3
        if i==5
            
            P(j,i) = P(j,i-1);
            rel_hum(j,i) = 0.9;
    
            T_max = T(j,i-2);
            T_min = T(j,i-1);
            T_avg = (T_min+T_max)/2;
    
            M_air = (3.76/4.76)*M_N2+(1/4.76)*M_O2;
            w1 = (M_H2O/M_air)*(P_H2O(j,i-2)/(P(j,i-2)-P_H2O(j,i-2)));
    
            E_enter = (3.76/4.76)*h_N2(j,i-2)+(1/4.76)*h_O2(j,i-2)+w1*(M_air/M_H2O)*h_H2O(j,i-2)-w1*(M_air/M_H2O)*h_total(j,i-1);
            
            y_CO2(j,i)=0;
            y_H2O(j,i)=1;
            y_N2(j,i)=1;
            y_O2(j,i)=1;
            y_CH4(j,i)=0;
    
            error = 10^-9;
    
            while true
                h_H2O(j,i) = h_H2O_func(R_uni,T_avg);
                h_N2(j,i) = h_N2_func(R_uni,T_avg);
                h_O2(j,i) = h_O2_func(R_uni,T_avg);
                
                h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
    
                P_sat(j,i) = P_sat_func(T_avg);
                E_exit = (3.76/4.76)*h_N2(j,i)+(1/4.76)*h_O2(j,i)-(rel_hum(j,i)*P_sat(j,i))/(P(j,i)...
                    -rel_hum(j,i)*P_sat(j,i))*h_total(j,i-1)+(rel_hum(j,i)*P_sat(j,i))/(P(j,i)-rel_hum(j,i)*P_sat(j,i))*h_H2O(j,i);
    
    
                if abs(E_exit-E_enter) < error
                    T(j,i) = T_avg;
                    break
                elseif E_exit > E_enter
                    T_max = T_avg;
                    T_avg = (T_min+T_max)/2;
                elseif E_exit < E_enter
                    T_min = T_avg;
                    T_avg = (T_min+T_max)/2;
                end
            end
    
            P_H2O(j,i) = rel_hum(j,i)*P_sat_func(T(j,i));
            P_N2(j,i) = (3.76/4.76)*(P(j,i)-P_H2O(j,i));
            P_O2(j,i) = (1/4.76)*(P(j,i)-P_H2O(j,i));
    
            y_CO2(j,i)=0;
            y_H2O(j,i)=P_H2O(j,i)/P(j,i);
            y_N2(j,i)=P_N2(j,i)/P(j,i);
            y_O2(j,i)=P_O2(j,i)/P(j,i);
            y_CH4(j,i)=0;
    
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
    
            s_CO2(j,i) = 0;
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = s_O2_func(R_uni,P(j,i),T(j,i),y_O2(j,i));
            s_CH4(j,i) = 0;
    
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
            
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = ef_O2_func(h_O2(j,i),s_O2(j,i),y_O2(j,i));
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
        %% STATE 4s
        if i==6
            
            comp_ratio = 4;
            P(j,i) = comp_ratio*P(j,i-1);
            
            P_CO2(j,i) = 0;
            P_H2O(j,i) =comp_ratio*P_H2O(j,i-1);
            P_N2(j,i) = comp_ratio*P_N2(j,i-1);
            P_O2(j,i) = comp_ratio*P_O2(j,i-1);
            P_CH4(j,i) = 0;
            
            y_CO2(j,i)=0;
            y_H2O(j,i)=P_H2O(j,i)/P(j,i);
            y_N2(j,i)=P_N2(j,i)/P(j,i);
            y_O2(j,i)=P_O2(j,i)/P(j,i);
            y_CH4(j,i)=0;
            
            s_total(j,i)=s_total(j,i-1);
            
            T(j,i) = T_find_s(s_total(j,i),P(j,i),y_CO2(j,i),y_H2O(j,i),y_N2(j,i),y_O2(j,i),y_CH4(j,i));
            
            h_CO2(j,i) = 0;
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = h_O2_func(R_uni,T(j,i));
            h_CH4(j,i) = 0;
            
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
    
            s_CO2(j,i) = 0;
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = s_O2_func(R_uni,P(j,i),T(j,i),y_O2(j,i));
            s_CH4(j,i) = 0;
    
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = ef_O2_func(h_O2(j,i),s_O2(j,i),y_O2(j,i));
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
        %% STATE 4
        if i==7
            
            h_total(j,i) = (h_total(j,i-1)-h_total(j,i-2)+h_total(j,i-2)*n_c)/n_c;
            
            P(j,i)=P(j,i-1);
            P_CO2(j,i) = 0;
            P_H2O(j,i) =P_H2O(j,i-1);
            P_N2(j,i) = P_N2(j,i-1);
            P_O2(j,i) = P_O2(j,i-1);
            P_CH4(j,i) = 0;
            
            y_CO2(j,i)=0;
            y_H2O(j,i)=y_H2O(j,i-1);
            y_N2(j,i)=y_N2(j,i-1);
            y_O2(j,i)=y_O2(j,i-1);
            y_CH4(j,i)=0;
            
            T(j,i) = T_find_h(h_total(j,i),y_CO2(j,i),y_H2O(j,i),y_N2(j,i),y_O2(j,i),y_CH4(j,i));
            
            s_CO2(j,i) = 0;
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = s_O2_func(R_uni,P(j,i),T(j,i),y_O2(j,i));
            s_CH4(j,i) = 0;
    
            h_CO2(j,i) = 0;
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = h_O2_func(R_uni,T(j,i));
            h_CH4(j,i) = 0;
    
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
    
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = ef_O2_func(h_O2(j,i),s_O2(j,i),y_O2(j,i));
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
    
        %% STATE 4F
        if i==8
            
            T(j,i) = T_0; % K;
            P(j,i)=P(j,i-1);
    
            P_CO2(j,i) = 0;
            P_H2O(j,i) = 0;
            P_N2(j,i) = 0;
            P_O2(j,i) = 0;
            P_CH4(j,i) = P(j,i);
            
            y_CO2(j,i)=0;
            y_H2O(j,i)=0;
            y_N2(j,i)=0;
            y_O2(j,i)=0;
            y_CH4(j,i)=1;
            
            h_CO2(j,i) = 0;
            h_H2O(j,i) = 0;
            h_N2(j,i) = 0;
            h_O2(j,i) = 0;
            h_CH4(j,i) = h_CH4_func(R_uni,T(j,i));
            
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
            
            s_CO2(j,i) = 0;
            s_H2O(j,i) = 0;
            s_N2(j,i) = 0;
            s_O2(j,i) = 0;
            s_CH4(j,i) = s_CH4_func(R_uni,P(j,i),T(j,i),y_CH4(j,i));
    
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
    
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = 0;
            ef_H2O(j,i) = 0;
            ef_N2(j,i) = 0;
            ef_O2(j,i) = 0;
            ef_CH4(j,i) = ef_CH4_func(h_CH4(j,i),s_CH4(j,i),y_CH4(j,i))+831650;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
        %% STATE 5
        if i==9
            
            T(j,i) = Temp_max(j); % K;
            P(j,i) = P(j,i-2)*0.95; % kPa
    
            x = P_H2O(j,i-2)/(P(j,i-2)-P_H2O(j,i-2));
    
            y_CO2(j,i)=1/(7.52+1+2+2*4.76*x);
            y_H2O(j,i)=(2+2*4.76*x)/(7.52+1+2+2*4.76*x);
            y_N2(j,i)=7.52/(7.52+1+2+2*4.76*x);
            y_O2(j,i)=0;
            y_CH4(j,i)=0;
            
            P_CO2(j,i) = P(j,i)*y_CO2(j,i);
            P_H2O(j,i) = P(j,i)*y_H2O(j,i);
            P_N2(j,i) = P(j,i)*y_N2(j,i);
            P_O2(j,i) = 0;
            P_CH4(j,i) = 0;
    
            h_CO2(j,i) = h_CO2_func(R_uni,T(j,i));
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = 0;
            h_CH4(j,i) = 0;
    
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
    
            s_CO2(j,i) = s_CO2_func(R_uni,P(j,i),T(j,i),y_CO2(j,i));
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = 0;
            s_CH4(j,i) = 0;
    
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
    
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = ef_CO2_func(h_CO2(j,i),s_CO2(j,i),y_CO2(j,i));
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = 0;
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    
        %% STATE 6s
        if i==10
           
            P(j,i) = P(j,i-1)/(comp_ratio^2*0.95);
            
            P_CO2(j,i) = P_CO2(j,i-1)/(comp_ratio^2*0.95);
            P_H2O(j,i) =P_H2O(j,i-1)/(comp_ratio^2*0.95);
            P_N2(j,i) = P_N2(j,i-1)/(comp_ratio^2*0.95);
            P_O2(j,i) = 0;
            P_CH4(j,i) = 0;
            
            y_CO2(j,i)=P_CO2(j,i)/P(j,i);
            y_H2O(j,i)=P_H2O(j,i)/P(j,i);
            y_N2(j,i)=P_N2(j,i)/P(j,i);
            y_O2(j,i)=0;
            y_CH4(j,i)=0;
            
            s_total(j,i)=s_total(j,i-1);
            
            T(j,i) = T_find_s(s_total(j,i),P(j,i),y_CO2(j,i),y_H2O(j,i),y_N2(j,i),y_O2(j,i),y_CH4(j,i));
            
            h_CO2(j,i) = h_CO2_func(R_uni,T(j,i));
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = 0;
            h_CH4(j,i) = 0;
            
            h_total(j,i) = y_CO2(j,i)*h_CO2(j,i)+y_H2O(j,i)*h_H2O(j,i)+y_N2(j,i)*h_N2(j,i)+y_O2(j,i)*h_O2(j,i)+y_CH4(j,i)*h_CH4(j,i);
    
            s_CO2(j,i) = s_CO2_func(R_uni,P(j,i),T(j,i),y_CO2(j,i));
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = 0;
            s_CH4(j,i) = 0;
    
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = ef_CO2_func(h_CO2(j,i),s_CO2(j,i),y_CO2(j,i));
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = 0;
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
        %% STATE 6
        if i==11
            
            h_total(j,i) = (h_total(j,i-2)-h_total(j,i-2)*n_t+h_total(j,i-1)*n_t);
            
            P(j,i)=P(j,i-1);
    
            P_CO2(j,i) = P_CO2(j,i-1);
            P_H2O(j,i) =P_H2O(j,i-1);
            P_N2(j,i) = P_N2(j,i-1);
            P_O2(j,i) = 0;
            P_CH4(j,i) = 0;
            
            y_CO2(j,i)=y_CO2(j,i-1);
            y_H2O(j,i)=y_H2O(j,i-1);
            y_N2(j,i)=y_N2(j,i-1);
            y_O2(j,i)=0;
            y_CH4(j,i)=0;
            
            T(j,i) = T_find_h(h_total(j,i),y_CO2(j,i),y_H2O(j,i),y_N2(j,i),y_O2(j,i),y_CH4(j,i));
            
            s_CO2(j,i) = s_CO2_func(R_uni,P(j,i),T(j,i),y_CO2(j,i));
            s_H2O(j,i) = s_H2O_func(R_uni,P(j,i),T(j,i),y_H2O(j,i));
            s_N2(j,i) = s_N2_func(R_uni,P(j,i),T(j,i),y_N2(j,i));
            s_O2(j,i) = 0;
            s_CH4(j,i) = 0;
    
            h_CO2(j,i) = h_CO2_func(R_uni,T(j,i));
            h_H2O(j,i) = h_H2O_func(R_uni,T(j,i));
            h_N2(j,i) = h_N2_func(R_uni,T(j,i));
            h_O2(j,i) = 0;
            h_CH4(j,i) = 0;
    
            s_total(j,i) = y_CO2(j,i)*s_CO2(j,i)+y_H2O(j,i)*s_H2O(j,i)+y_N2(j,i)*s_N2(j,i)+y_O2(j,i)*s_O2(j,i)+y_CH4(j,i)*s_CH4(j,i);
    
            P_sat(j,i)=P_sat_func(T(j,i));
            rel_hum(j,i)=P_H2O(j,i)/P_sat(j,i);
    
            ef_CO2(j,i) = ef_CO2_func(h_CO2(j,i),s_CO2(j,i),y_CO2(j,i));
            ef_H2O(j,i) = ef_H2O_func(h_H2O(j,i),s_H2O(j,i),y_H2O(j,i));
            ef_N2(j,i) = ef_N2_func(h_N2(j,i),s_N2(j,i),y_N2(j,i));
            ef_O2(j,i) = 0;
            ef_CH4(j,i) = 0;
    
            ef_total(j,i) = ef_CO2(j,i)*y_CO2(j,i)+ef_H2O(j,i)*y_H2O(j,i)+ef_N2(j,i)*y_N2(j,i)+ef_O2(j,i)*y_O2(j,i)+ef_CH4(j,i)*y_CH4(j,i);
        end
    end
    
    
    W_turbine = 175*10^3; % kW   

    ndot(j,11) = W_turbine/(h_total(j,9)-h_total(j,11));
    ndot(j,10) = ndot(j,11);
    ndot(j,9) = ndot(j,10);
    ndot(j,8) = ndot(j,9)*y_CO2(j,9);
    ndot(j,7) = ndot(j,11)-ndot(j,8);
    ndot(j,6) = ndot(j,7);
    ndot(j,5) = ndot(j,6);
    ndot(j,4) = ndot(j,5)*(y_H2O(j,5)-y_H2O(j,3));
    ndot(j,3) = ndot(j,5)-ndot(j,4);
    ndot(j,2) = ndot(j,3);
    ndot(j,1) = ndot(j,2);

end
     

%% PROCESS
    
num_process=5;

Qdot=zeros(num_iteration,num_process);
Wdot=zeros(num_iteration,num_process);
Sdot_gen=zeros(num_iteration,num_process);
ex_flow_diff=zeros(num_iteration,num_process);
ex_dest=zeros(num_iteration,num_process);
eta_array=zeros(num_iteration,num_process);

for j=1:length(Temp_max)
    for i=1:num_process
    
    %% Qdot
    Qcycle(j) = 0;
    for i=1:num_process
        if i==1
            Qdot(j,i)=0;
        elseif i==2
            Qdot(j,i)=0;
        elseif i==3
            Qdot(j,i)=0;
        elseif i==4
            Qdot(j,i)=ndot(j,i+5)*h_total(j,i+5)-ndot(j,i+3)*h_total(j,i+3)-ndot(j,i+4)*h_total(j,i+4);
        elseif i==5
            Qdot(j,i)=0;
        end
        Qcycle(j) = Qcycle(j)+Qdot(j,i);
    end
    
    %% Wdot
    Wcycle(j) = 0;
    for i=1:num_process
        if i==1
            Wdot(j,i)=ndot(j,i)*(h_total(j,i)-h_total(j,i+2));
        elseif i==2
            Wdot(j,i)=0;
        elseif i==3
            Wdot(j,i)=ndot(j,i+2)*(h_total(j,i+2)-h_total(j,i+4));
        elseif i==4
            Wdot(j,i)=0;
        elseif i==5
            Wdot(j,i)=175*10^3;
        end
        Wcycle(j) = Wcycle(j)+Wdot(j,i);
    end
    
    %% Sdot_gen
    Sdot_cycle(j) = 0;
    for i=1:num_process
        if i==1
            Sdot_gen(j,i)=ndot(j,i)*s_total(j,i+2)-ndot(j,i)*s_total(j,i);
        elseif i==2
            Sdot_gen(j,i)=ndot(j,i+3)*s_total(j,i+3)-ndot(j,i+1)*s_total(j,i+1)-ndot(j,i+2)*s_total(j,i+2);
        elseif i==3
            Sdot_gen(j,i)=ndot(j,i+4)*s_total(j,i+4)-ndot(j,i+2)*s_total(j,i+2);
        elseif i==4
            Sdot_gen(j,i)=ndot(j,i+5)*s_total(j,i+5)-ndot(j,i+4)*s_total(j,i+4)-ndot(j,i+3)*s_total(j,i+3)-Qdot(j,i)/T(j,i+5);
        elseif i==5
            Sdot_gen(j,i)=ndot(j,i+6)*s_total(j,i+6)-ndot(j,i+4)*s_total(j,i+4);
        end
        Sdot_cycle(j) = Sdot_cycle(j) + Sdot_gen(j,i);
    end
    
    %% ex_flow_diff
    ex_flow_cycle(j) = 0;
    for i=1:num_process
        if i==1
            ex_flow_diff(j,i) = ndot(j,i+2)*ef_total(j,i+2)-ndot(j,i)*ef_total(j,i);
        elseif i==2
            ex_flow_diff(j,i) = ndot(j,i+3)*ef_total(j,i+3)-ndot(j,i+1)*ef_total(j,i+1)-ndot(j,i+2)*ef_total(j,i+2);
        elseif i==3
            ex_flow_diff(j,i) = ndot(j,i+4)*ef_total(j,i+4)-ndot(j,i+2)*ef_total(j,i+2);
        elseif i==4
            ex_flow_diff(j,i) = ndot(j,i+5)*ef_total(j,i+5)-ndot(j,i+4)*ef_total(j,i+4)-ndot(j,i+3)*ef_total(j,i+3);
        elseif i==5
            ex_flow_diff(j,i) = ndot(j,i+6)*ef_total(j,i+6)-ndot(j,i+4)*ef_total(j,i+4);
        end
        ex_flow_cycle(j) = ex_flow_cycle(j) + ex_flow_diff(j,i);
    end
    
    %% ex_dest
    ex_dest_cycle(j) = 0;
    for i=1:num_process
        if i==1
            ex_dest(j,i) = -ex_flow_diff(j,i)-Wdot(j,i)+(1-T_0/T(i))*Qdot(j,i);
        elseif i==2
            ex_dest(j,i) = -ex_flow_diff(j,i)-Wdot(j,i)+(1-T_0/T(i))*Qdot(j,i);
        elseif i==3
            ex_dest(j,i) = -ex_flow_diff(j,i)-Wdot(j,i)+(1-T_0/T(i))*Qdot(j,i);
        elseif i==4
            ex_dest(j,i) = -ex_flow_diff(j,i)-Wdot(j,i)+(1-T_0/T(i))*Qdot(j,i);
        elseif i==5
            ex_dest(j,i) = -ex_flow_diff(j,i)-Wdot(j,i)+(1-T_0/T(j,i+4))*Qdot(j,i);
        end
        ex_dest_cycle(j) = ex_dest_cycle(j) + ex_dest(j,i);
    end
    
    
    %% eta
    for i=1:num_process
        if i==1
            eta_array(j,i) = ndot(j,i)*(ef_total(j,i+1)-ef_total(j,i)) / -Wdot(j,i);
        elseif i==2
            eta_array(j,i) = 0;
        elseif i==3
            eta_array(j,i) = ndot(j,i+2)*(ef_total(j,i+4)-ef_total(j,i+2)) / -Wdot(j,i);
        elseif i==4
            eta_array(j,i) = ndot(j,i+5)*ef_total(j,i+5)/(ndot(j,i+3)*ef_total(j,i+3)+ndot(j,i+4)*ef_total(j,i+4));
        elseif i==5
            eta_array(j,i) = Wdot(j,i) / (ndot(j,i+5)*(ef_total(j,i+4)-ef_total(j,i+6))) ;
        end
    end
    eta_cycle(j) = Wcycle(j)/(ndot(j,i+3)*ef_total(j,i+3));

    end
    W_compressor(j) = -Wdot(j,1)-Wdot(j,3);
end

rel_hum(:,:) = rel_hum(:,:)*100;

% for j=1:length(Temp_max)
% 
%     %% STATE TABLE
%     
%     State = {'1';'2';'2L';'3';'4';'4F';'5';'6'};
%     Description = {'Low Pressure Compressor Inlet';'Intercooler Air Inlet';'Intercooler Water Inlet';'High Pressure Compressor Inlet';'Combustor Air Inlet';'Combustor Fuel Inlet';'Turbine Inlet';'Turbine Outlet'};
%     T_table = [T(j,1);T(j,3);T(j,4);T(j,5);T(j,7);T(j,8);T(j,9);T(j,11)];
%     P_table = [P(j,1);P(j,3);P(j,4);P(j,5);P(j,7);P(j,8);P(j,9);P(j,11)];
%     n_dot_table = [ndot(j,1);ndot(j,3);ndot(j,4);ndot(j,5);ndot(j,7);ndot(j,8);ndot(j,9);ndot(j,11)];
%     rel_hum_table = [rel_hum(j,1);rel_hum(j,3);rel_hum(j,4);rel_hum(j,5);rel_hum(j,7);rel_hum(j,8);rel_hum(j,9);rel_hum(j,11)];
%     h_table = [h_total(j,1);h_total(j,3);h_total(j,4);h_total(j,5);h_total(j,7);h_total(j,8);h_total(j,9);h_total(j,11)];
%     s_table = [s_total(j,1);s_total(j,3);s_total(j,4);s_total(j,5);s_total(j,7);s_total(j,8);s_total(j,9);s_total(j,11)];
%     ef_table = [ef_total(j,1);ef_total(j,3);ef_total(j,4);ef_total(j,5);ef_total(j,7);ef_total(j,8);ef_total(j,9);ef_total(j,11)];
%     TA1 = table(State,Description,T_table,P_table,n_dot_table,rel_hum_table,h_table,s_table,ef_table)
%     
%     
%     %% PROCESS TABLE
%     eta = 1;
%     
%     Process = {'1->2';'2 + 2L ->3';'3->4';'4 + 4F->5';'5->6';'Cycle'};
%     Description = {'Low Pressure Compressor';'Intercooler';'High Pressure Compressor';'Combustor';'Turbine';'Entire Cycle'};
%     Q_dot = [Qdot(j,1);Qdot(2);Qdot(j,3);Qdot(j,4);Qdot(j,5);Qcycle];
%     W_dot = [Wdot(j,1);Wdot(2);Wdot(j,3);Wdot(j,4);Wdot(j,5);Wcycle];
%     sigma_gen_dot = [Sdot_gen(j,1);Sdot_gen(2);Sdot_gen(j,3);Sdot_gen(j,4);Sdot_gen(j,5);Sdot_cycle];
%     del_ef = [ex_flow_diff(j,1);ex_flow_diff(2);ex_flow_diff(j,3);ex_flow_diff(j,4);ex_flow_diff(j,5);ex_flow_cycle];
%     ex_dest_table = [ex_dest(j,1);ex_dest(2);ex_dest(j,3);ex_dest(j,4);ex_dest(j,5);ex_dest_cycle];
%     eta_table = {eta_array(j,1);'-';eta_array(j,3);eta_array(j,4);eta_array(j,5);eta_cycle};
%     TA2 = table(Process,Description,Q_dot,W_dot,sigma_gen_dot,del_ef,ex_dest_table,eta_table)
% 
% end


%% Graphs


mass_fuel = ndot(:,8)*M_CH4;
eta_cycle(:) = eta_cycle(:)*100;
bwr = 100.*W_compressor/W_turbine;

figure(1)
plot(mass_fuel,bwr,'LineWidth',2)
grid on
xlabel("Fuel Mass Flow Rate (kg/s)")
ylabel("Ratio (%)")
title('Parametric Analysis 2')
hold on
plot(mass_fuel,eta_cycle,'LineWidth',2)
legend("Backwork Ratio",'2nd Law Efficiency')



%% Functions

function h=h_CO2_func(R,T_state)
   
    T_ref = 298; % K
    h_formation = -393520; %kJ/kmol

    alpha = 2.401;
    beta = 8.735*10^-3;
    gamma = -6.607*10^-6;
    delta = 2.002*10^-9;
    epsilon = 0;
    
    h_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R; 

    h = integral(h_func,T_ref,T_state) + h_formation;
end

function h=h_H2O_func(R,T_state)
    
    T_ref = 298; % K
    h_formation = -241820; %kJ/kmol

    alpha = 4.070;
    beta = -1.108*10^-3;
    gamma = 4.152*10^-6;
    delta = -2.964*10^-9;
    epsilon = 0.807*10^-12;
    
    h_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R; 

    h = integral(h_func,T_ref,T_state) + h_formation;
end

function h=h_N2_func(R,T_state)
    
    T_ref = 298; % K
    h_formation = 0; %kJ/kmol

    alpha = 3.675;
    beta = -1.208*10^-3;
    gamma = 2.324*10^-6;
    delta = -0.632*10^-9;
    epsilon = -0.226*10^-12;
    
    h_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R; 

    h = integral(h_func,T_ref,T_state) + h_formation;
end

function h=h_O2_func(R,T_state)
    
    T_ref = 298; % K
    h_formation = 0; %kJ/kmol

    alpha = 3.626;
    beta = -1.878*10^-3;
    gamma = 7.055*10^-6;
    delta = -6.764*10^-9;
    epsilon = 2.156*10^-12;
    
    h_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R; 

    h = integral(h_func,T_ref,T_state) + h_formation;
end

function h=h_CH4_func(R,T_state)
    
    T_ref = 298; % K
    h_formation = -74850; %kJ/kmol

    alpha = 3.826;
    beta = -3.979*10^-3;
    gamma = 24.558*10^-6;
    delta = -22.733*10^-9;
    epsilon = 6.963*10^-12;
    
    h_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R; 

    h = integral(h_func,T_ref,T_state) + h_formation;
end


function s = s_CO2_func(R,P_state,T_state,y)

    T_ref = 298; % K
    P_ref = 101.325; % kPa
    s_ref = 213.69; % kJ/kmol*K

    alpha = 2.401;
    beta = 8.735*10^-3;
    gamma = -6.607*10^-6;
    delta = 2.002*10^-9;
    epsilon = 0;

    s_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R./T; 

    s = integral(s_func,T_ref,T_state)-R*log(y*P_state/P_ref)+s_ref;
end

function s = s_H2O_func(R,P_state,T_state,y)

    T_ref = 298; % K
    P_ref = 101.325; % kPa
    s_ref = 188.72; % kJ/kmol*K

    alpha = 4.070;
    beta = -1.108*10^-3;
    gamma = 4.152*10^-6;
    delta = -2.964*10^-9;
    epsilon = 0.807*10^-12;

    s_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R./T; 

    s = integral(s_func,T_ref,T_state)-R*log(y*P_state/P_ref)+s_ref;
end

function s = s_N2_func(R,P_state,T_state,y)

    T_ref = 298; % K
    P_ref = 101.325; % kPa
    s_ref = 191.50; % kJ/kmol*K

    alpha = 3.675;
    beta = -1.208*10^-3;
    gamma = 2.324*10^-6;
    delta = -0.632*10^-9;
    epsilon = -0.226*10^-12;

    s_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R./T; 

    s = integral(s_func,T_ref,T_state)-R*log(y*P_state/P_ref)+s_ref;
end

function s = s_O2_func(R,P_state,T_state,y)

    T_ref = 298; % K
    P_ref = 101.325; % kPa
    s_ref = 205.03; % kJ/kmol*K

    alpha = 3.626;
    beta = -1.878*10^-3;
    gamma = 7.055*10^-6;
    delta = -6.764*10^-9;
    epsilon = 2.156*10^-12;

    s_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R./T; 

    s = integral(s_func,T_ref,T_state)-R*log(y*P_state/P_ref)+s_ref;
end

function s = s_CH4_func(R,P_state,T_state,y)

    T_ref = 298; % K
    P_ref = 101.325; % kPa
    s_ref = 186.16; % kJ/kmol*K

    alpha = 3.826;
    beta = -3.979*10^-3;
    gamma = 24.558*10^-6;
    delta = -22.733*10^-9;
    epsilon = 6.963*10^-12;

    s_func = @(T) (alpha+beta.*T+gamma.*T.^2+delta.*T.^3+epsilon.*T.^4)*R./T; 

    s = integral(s_func,T_ref,T_state)-R*log(y*P_state/P_ref)+s_ref;
end

function ef = ef_CO2_func(h_state,s_state,y_CO2)
    
    R_uni = 8.314; %kJ/kmol*K
    T_ref = 300; % K
    P_ref = 100; %kPa
    h0 = h_CO2_func(R_uni,T_ref); %kJ/kmol
    s0 = s_CO2_func(R_uni,P_ref,T_ref,y_CO2); %kJ/kmol*K

    ef = h_state-h0-T_ref*(s_state-s0);

end

function ef = ef_H2O_func(h_state,s_state,y_H2O)
    
    R_uni = 8.314; %kJ/kmol*K
    T_ref = 300; % K
    P_ref = 100; %kPa
    h0 = h_H2O_func(R_uni,T_ref); %kJ/kmol
    s0 = s_H2O_func(R_uni,P_ref,T_ref,y_H2O); %kJ/kmol*K


    ef = h_state-h0-T_ref*(s_state-s0);

end

function ef = ef_N2_func(h_state,s_state,y_N2)
    
    R_uni = 8.314; %kJ/kmol*K
    T_ref = 300; % K
    P_ref = 100; %kPa
    h0 = h_N2_func(R_uni,T_ref); %kJ/kmol
    s0 = s_N2_func(R_uni,P_ref,T_ref,y_N2); %kJ/kmol*K

    ef = h_state-h0-T_ref*(s_state-s0);

end

function ef = ef_O2_func(h_state,s_state,y_O2)
    
    R_uni = 8.314; %kJ/kmol*K
    T_ref = 300; % K
    P_ref = 100; %kPa
    h0 = h_O2_func(R_uni,T_ref); %kJ/kmol
    s0 = s_O2_func(R_uni,P_ref,T_ref,y_O2); %kJ/kmol*K


    ef = h_state-h0-T_ref*(s_state-s0);

end

function ef = ef_CH4_func(h_state,s_state,y_CH4)
    
    R_uni = 8.314; %kJ/kmol*K
    T_ref = 300; % K
    P_ref = 100; %kPa
    h0 = h_CH4_func(R_uni,T_ref); %kJ/kmol
    s0 = s_CH4_func(R_uni,P_ref,T_ref,y_CH4); %kJ/kmol*K

    ef = h_state-h0-T_ref*(s_state-s0);

end

function P_sat = P_sat_func(T)
    T = T-273.15;
    P_sat =  ((exp(34.494-(4924.99/(T+237.1))))/(T+105)^1.57)*10^-3; % kPa
end

function T = T_find_h(h,y_CO2,y_H2O,y_N2,y_O2,y_CH4)
    
    T_min = 280; %K
    T_max = 1000; %K
    T_avg = (T_max+T_min)/2; %K

    R_uni = 8.314; %kJ/kmol*K
    
    alpha1 = 2.401;
    beta1 = 8.735*10^-3;
    gamma1 = -6.607*10^-6;
    delta1 = 2.002*10^-9;
    epsilon1 = 0;

    alpha2 = 4.070;
    beta2 = -1.108*10^-3;
    gamma2 = 4.152*10^-6;
    delta2 = -2.964*10^-9;
    epsilon2 = 0.807*10^-12;

    alpha3 = 3.675;
    beta3 = -1.208*10^-3;
    gamma3 = 2.324*10^-6;
    delta3 = -0.632*10^-9;
    epsilon3 = -0.226*10^-12;

    alpha4 = 3.626;
    beta4 = -1.878*10^-3;
    gamma4 = 7.055*10^-6;
    delta4 = -6.764*10^-9;
    epsilon4 = 2.156*10^-12;

    alpha5 = 3.826;
    beta5 = -3.979*10^-3;
    gamma5 = 24.558*10^-6;
    delta5 = -22.733*10^-9;
    epsilon5 = 6.963*10^-12;

    T_ref = 298; % K
    
    h_formation1 = -393520; % kJ/kg
    h_formation2 = -242820; % kJ/kg
    h_formation3 = 0; % kJ/kg
    h_formation4 = 0; % kJ/kg
    h_formation5 = 0; % kJ/kg

    h_func1 = @(T) (alpha1+beta1.*T+gamma1.*T.^2+delta1.*T.^3+epsilon1.*T.^4)*R_uni;
    h_func2 = @(T) (alpha2+beta2.*T+gamma2.*T.^2+delta2.*T.^3+epsilon2.*T.^4)*R_uni;
    h_func3 = @(T) (alpha3+beta3.*T+gamma3.*T.^2+delta3.*T.^3+epsilon3.*T.^4)*R_uni;
    h_func4 = @(T) (alpha4+beta4.*T+gamma4.*T.^2+delta4.*T.^3+epsilon4.*T.^4)*R_uni;
    h_func5 = @(T) (alpha5+beta5.*T+gamma5.*T.^2+delta5.*T.^3+epsilon5.*T.^4)*R_uni;
    
    error = 10^-9;
    
    while true

        h1 = (integral(h_func1,T_ref,T_avg)+h_formation1)*y_CO2;
        h2 = (integral(h_func2,T_ref,T_avg)+h_formation2)*y_H2O;
        h3 = (integral(h_func3,T_ref,T_avg)+h_formation3)*y_N2;
        h4 = (integral(h_func4,T_ref,T_avg)+h_formation4)*y_O2;
        h5 = (integral(h_func5,T_ref,T_avg)+h_formation5)*y_CH4;

        h_guess = h1+h2+h3+h4+h5;

        if abs(h - h_guess) < error
            T = T_avg;
            break
        elseif h_guess < h
            T_min = T_avg;
            T_avg = (T_min+T_max)/2;
        elseif h_guess > h
            T_max = T_avg;
            T_avg = (T_min+T_max)/2;
        end
    end 
end


function T = T_find_s(s,P,y_CO2,y_H2O,y_N2,y_O2,y_CH4)

    T_min = 280; %K
    T_max = 1000; %K
    T_avg = (T_max+T_min)/2; %K

    R_uni = 8.314; %kJ/kmol*K
    
    alpha1 = 2.401;
    beta1 = 8.735*10^-3;
    gamma1 = -6.607*10^-6;
    delta1 = 2.002*10^-9;
    epsilon1 = 0;

    alpha2 = 4.070;
    beta2 = -1.108*10^-3;
    gamma2 = 4.152*10^-6;
    delta2 = -2.964*10^-9;
    epsilon2 = 0.807*10^-12;

    alpha3 = 3.675;
    beta3 = -1.208*10^-3;
    gamma3 = 2.324*10^-6;
    delta3 = -0.632*10^-9;
    epsilon3 = -0.226*10^-12;

    alpha4 = 3.626;
    beta4 = -1.878*10^-3;
    gamma4 = 7.055*10^-6;
    delta4 = -6.764*10^-9;
    epsilon4 = 2.156*10^-12;

    alpha5 = 3.826;
    beta5 = -3.979*10^-3;
    gamma5 = 24.558*10^-6;
    delta5 = -22.733*10^-9;
    epsilon5 = 6.963*10^-12;

    T_ref = 298; % K
    P_ref = 101.325; % kPa

    s_ref1 = 213.69; % kJ/kg
    s_ref2 = 188.72; % kJ/kg
    s_ref3 = 191.50; % kJ/kg
    s_ref4 = 205.03; % kJ/kg
    s_ref5 = 186.16; % kJ/kg

    s_func1 = @(T) (alpha1+beta1.*T+gamma1.*T.^2+delta1.*T.^3+epsilon1.*T.^4).*(R_uni./T);
    s_func2 = @(T) (alpha2+beta2.*T+gamma2.*T.^2+delta2.*T.^3+epsilon2.*T.^4).*(R_uni./T);
    s_func3 = @(T) (alpha3+beta3.*T+gamma3.*T.^2+delta3.*T.^3+epsilon3.*T.^4).*(R_uni./T);
    s_func4 = @(T) (alpha4+beta4.*T+gamma4.*T.^2+delta4.*T.^3+epsilon4.*T.^4).*(R_uni./T);
    s_func5 = @(T) (alpha5+beta5.*T+gamma5.*T.^2+delta5.*T.^3+epsilon5.*T.^4).*(R_uni./T);
    
    error = 10^-9;
    
    while true

        if y_CO2==0
            s1 = 0;
        else
            s1 = (integral(s_func1,T_ref,T_avg)-R_uni*log(y_CO2*P/P_ref)+s_ref1)*y_CO2;
        end   
        
        if y_O2==0
            s4 = 0;
        else
            s4 = (integral(s_func4,T_ref,T_avg)-R_uni*log(y_O2*P/P_ref)+s_ref4)*y_O2;
        end 
        
        if y_CH4==0
            s5 = 0;
        else
            s5 = (integral(s_func5,T_ref,T_avg)-R_uni*log(y_CH4*P/P_ref)+s_ref5)*y_CH4;
        end 
        
        s2 = (integral(s_func2,T_ref,T_avg)-R_uni*log(y_H2O*P/P_ref)+s_ref2)*y_H2O;
        s3 = (integral(s_func3,T_ref,T_avg)-R_uni*log(y_N2*P/P_ref)+s_ref3)*y_N2;
        
        s_guess = s1+s2+s3+s4+s5;
        
        if abs(s - s_guess) < error
            T = T_avg;
            break
        elseif s_guess < s
            T_min = T_avg;
            T_avg = (T_min+T_max)/2;
        elseif s_guess > s
            T_max = T_avg;
            T_avg = (T_min+T_max)/2;
        end
    end
end

