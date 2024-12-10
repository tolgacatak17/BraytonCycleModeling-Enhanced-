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

T=zeros(1,num_state);
P=zeros(1,num_state);
rel_hum=zeros(1,num_state);


h_total=zeros(1,num_state);
h_CO2=zeros(1,num_state);
h_H2O=zeros(1,num_state);
h_N2=zeros(1,num_state);
h_O2=zeros(1,num_state);
h_CH4=zeros(1,num_state);

P_sat=zeros(1,num_state);
P_CO2=zeros(1,num_state);
P_H2O=zeros(1,num_state);
P_N2=zeros(1,num_state);
P_O2=zeros(1,num_state);
P_CH4=zeros(1,num_state);

y_CO2=zeros(1,num_state);
y_H2O=zeros(1,num_state);
y_N2=zeros(1,num_state);
y_O2=zeros(1,num_state);
y_CH4=zeros(1,num_state);

s_total=zeros(1,num_state);
s_CO2=zeros(1,num_state);
s_H2O=zeros(1,num_state);
s_N2=zeros(1,num_state);
s_O2=zeros(1,num_state);
s_CH4=zeros(1,num_state);

ef_CO2=zeros(1,num_state);
ef_H2O=zeros(1,num_state);
ef_N2=zeros(1,num_state);
ef_O2=zeros(1,num_state);
ef_CH4=zeros(1,num_state);
ef_total=zeros(1,num_state);

for i=1:num_state
    %% STATE 1
    if i==1

        T(i) = T_0;
        P(i) = P_0;
        
        rel_hum(i) = 0.6;
        
        P_sat(i)=P_sat_func(T(i));
        P_CO2(i) = 0;
        P_H2O(i) = P_sat(i)*rel_hum(i);
        P_N2(i) = (P(i)-P_H2O(i))*3.76/4.76;
        P_O2(i) = (P(i)-P_H2O(i))*1/4.76;
        P_CH4(i) = 0;
        
        y_CO2(i)=0;
        y_H2O(i)=P_H2O(i)/P(i);
        y_N2(i)=P_N2(i)/P(i);
        y_O2(i)=P_O2(i)/P(i);
        y_CH4(i)=0;
        
        h_CO2(i) = 0;
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = h_O2_func(R_uni,T(i));
        h_CH4(i) = 0;
        
        s_CO2(i) = 0;
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = s_O2_func(R_uni,P(i),T(i),y_O2(i));
        s_CH4(i) = 0;
        
        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);
        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);

        ef_CO2(i) = 0;
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = ef_O2_func(h_O2(i),s_O2(i),y_O2(i));
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end
    %% STATE 2s
    if i==2
        
        comp_ratio = 4;
        P(i) = comp_ratio*P(i-1);
        
        P_CO2(i) = 0;
        P_H2O(i) =comp_ratio*P_H2O(i-1);
        P_N2(i) = comp_ratio*P_N2(i-1);
        P_O2(i) = comp_ratio*P_O2(i-1);
        P_CH4(i) = 0;
        
        y_CO2(i)=0;
        y_H2O(i)=P_H2O(i)/P(i);
        y_N2(i)=P_N2(i)/P(i);
        y_O2(i)=P_O2(i)/P(i);
        y_CH4(i)=0;
        
        s_total(i)=s_total(i-1);
        
        
        T(i) = T_find_s(s_total(i),P(i),y_CO2(i),y_H2O(i),y_N2(i),y_O2(i),y_CH4(i));
        
        h_CO2(i) = 0;
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = h_O2_func(R_uni,T(i));
        h_CH4(i) = 0;
        
        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);
        
        s_CO2(i) = 0;
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = s_O2_func(R_uni,P(i),T(i),y_O2(i));
        s_CH4(i) = 0;
        
        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = 0;
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = ef_O2_func(h_O2(i),s_O2(i),y_O2(i));
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end

    %% STATE 2
    if i==3
        
        h_total(i) = (h_total(i-1)-h_total(i-2)+h_total(i-2)*n_c)/n_c;
        
        P(i)=P(i-1);
        P_CO2(i) = P_CO2(i-1);
        P_H2O(i) =P_H2O(i-1);
        P_N2(i) = P_N2(i-1);
        P_O2(i) = P_O2(i-1);
        P_CH4(i) = P_CH4(i-1);
        
        y_CO2(i)=y_CO2(i-1);
        y_H2O(i)=y_H2O(i-1);
        y_N2(i)=y_N2(i-1);
        y_O2(i)=y_O2(i-1);
        y_CH4(i)=y_CH4(i-1);
        
        T(i) = T_find_h(h_total(i),y_CO2(i),y_H2O(i),y_N2(i),y_O2(i),y_CH4(i));

        h_CO2(i) = 0;
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = h_O2_func(R_uni,T(i));
        h_CH4(i) = 0;
        
        s_CO2(i) = 0;
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = s_O2_func(R_uni,P(i),T(i),y_O2(i));
        s_CH4(i) = 0;

        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);

        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = 0;
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = ef_O2_func(h_O2(i),s_O2(i),y_O2(i));
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end

    %% STATE 2L
    if i==4
        
        T(i) = 280; % K
        
        cp_H2O = 4.179*M_H2O; % kj/kmol*K
        
        P(i)=P(i-1);
        P_CO2(i) = 0;
        P_H2O(i) = P(i);
        P_N2(i) = 0;
        P_O2(i) = 0;
        P_CH4(i) = 0;
        
        y_CO2(i)=0;
        y_H2O(i)=1;
        y_N2(i)=0;
        y_O2(i)=0;
        y_CH4(i)=0;
 
        
        h_formation_liquid = -285830; % kJ/kmol
        
        h_CO2(i) = 0;
        h_H2O(i) = cp_H2O*(T(i)-T_ref) + h_formation_liquid; 
        h_N2(i) = 0;
        h_O2(i) = 0;
        h_CH4(i) = 0;
        
        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);
        
        s_ref_liquid = 69.95; % kJ/kmol*K

        s_CO2(i) = 0;
        s_H2O(i) = cp_H2O*log(T(i)/T_ref)-R_uni*log(y_H2O(i)*P_H2O(i)/P_ref)+s_ref_liquid;
        s_N2(i) = 0;
        s_O2(i) = 0;
        s_CH4(i) = 0;

        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);
        
        h0_liquid = cp_H2O*(T_0-T_ref) + h_formation_liquid;
        s0_liquid = cp_H2O*log(T_0/T_ref)-R_uni*log(y_H2O(i)*P_0/P_ref)+s_ref_liquid;
        
        ef_CO2(i) = 0;
        ef_H2O(i) = h_H2O(i)-h0_liquid-T_0*(s_H2O(i)-s0_liquid)+900;
        ef_N2(i) = 0;
        ef_O2(i) = 0;
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end

    %% STATE 3
    if i==5
        
        P(i) = P(i-1);
        rel_hum(i) = 0.9;

        T_max = T(i-2);
        T_min = T(i-1);
        T_avg = (T_min+T_max)/2;

        M_air = (3.76/4.76)*M_N2+(1/4.76)*M_O2;
        w1 = (M_H2O/M_air)*(P_H2O(i-2)/(P(i-2)-P_H2O(i-2)));

        E_enter = (3.76/4.76)*h_N2(i-2)+(1/4.76)*h_O2(i-2)+w1*(M_air/M_H2O)*h_H2O(i-2)-w1*(M_air/M_H2O)*h_total(i-1);
        
        y_CO2(i)=0;
        y_H2O(i)=1;
        y_N2(i)=1;
        y_O2(i)=1;
        y_CH4(i)=0;

        error = 10^-9;

        while true
            h_H2O(i) = h_H2O_func(R_uni,T_avg);
            h_N2(i) = h_N2_func(R_uni,T_avg);
            h_O2(i) = h_O2_func(R_uni,T_avg);
            
            h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);

            P_sat(i) = P_sat_func(T_avg);
            E_exit = (3.76/4.76)*h_N2(i)+(1/4.76)*h_O2(i)-(rel_hum(i)*P_sat(i))/(P(i)-rel_hum(i)*P_sat(i))*h_total(i-1)+(rel_hum(i)*P_sat(i))/(P(i)-rel_hum(i)*P_sat(i))*h_H2O(i);


            if abs(E_exit-E_enter) < error
                T(i) = T_avg;
                break
            elseif E_exit > E_enter
                T_max = T_avg;
                T_avg = (T_min+T_max)/2;
            elseif E_exit < E_enter
                T_min = T_avg;
                T_avg = (T_min+T_max)/2;
            end
        end

        P_H2O(i) = rel_hum(i)*P_sat_func(T(i));
        P_N2(i) = (3.76/4.76)*(P(i)-P_H2O(i));
        P_O2(i) = (1/4.76)*(P(i)-P_H2O(i));

        y_CO2(i)=0;
        y_H2O(i)=P_H2O(i)/P(i);
        y_N2(i)=P_N2(i)/P(i);
        y_O2(i)=P_O2(i)/P(i);
        y_CH4(i)=0;

        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);

        s_CO2(i) = 0;
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = s_O2_func(R_uni,P(i),T(i),y_O2(i));
        s_CH4(i) = 0;

        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);
        
        ef_CO2(i) = 0;
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = ef_O2_func(h_O2(i),s_O2(i),y_O2(i));
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end

    %% STATE 4s
    if i==6
        
        comp_ratio = 4;
        P(i) = comp_ratio*P(i-1);
        
        P_CO2(i) = 0;
        P_H2O(i) =comp_ratio*P_H2O(i-1);
        P_N2(i) = comp_ratio*P_N2(i-1);
        P_O2(i) = comp_ratio*P_O2(i-1);
        P_CH4(i) = 0;
        
        y_CO2(i)=0;
        y_H2O(i)=P_H2O(i)/P(i);
        y_N2(i)=P_N2(i)/P(i);
        y_O2(i)=P_O2(i)/P(i);
        y_CH4(i)=0;
        
        s_total(i)=s_total(i-1);
        
        T(i) = T_find_s(s_total(i),P(i),y_CO2(i),y_H2O(i),y_N2(i),y_O2(i),y_CH4(i));
        
        h_CO2(i) = 0;
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = h_O2_func(R_uni,T(i));
        h_CH4(i) = 0;
        
        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);

        s_CO2(i) = 0;
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = s_O2_func(R_uni,P(i),T(i),y_O2(i));
        s_CH4(i) = 0;

        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = 0;
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = ef_O2_func(h_O2(i),s_O2(i),y_O2(i));
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end

    %% STATE 4
    if i==7
        
        h_total(i) = (h_total(i-1)-h_total(i-2)+h_total(i-2)*n_c)/n_c;
        
        P(i)=P(i-1);
        P_CO2(i) = 0;
        P_H2O(i) =P_H2O(i-1);
        P_N2(i) = P_N2(i-1);
        P_O2(i) = P_O2(i-1);
        P_CH4(i) = 0;
        
        y_CO2(i)=0;
        y_H2O(i)=y_H2O(i-1);
        y_N2(i)=y_N2(i-1);
        y_O2(i)=y_O2(i-1);
        y_CH4(i)=0;
        
        T(i) = T_find_h(h_total(i),y_CO2(i),y_H2O(i),y_N2(i),y_O2(i),y_CH4(i));
        
        s_CO2(i) = 0;
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = s_O2_func(R_uni,P(i),T(i),y_O2(i));
        s_CH4(i) = 0;

        h_CO2(i) = 0;
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = h_O2_func(R_uni,T(i));
        h_CH4(i) = 0;

        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);

        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = 0;
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = ef_O2_func(h_O2(i),s_O2(i),y_O2(i));
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end


    %% STATE 4F
    if i==8
        
        T(i) = T_0; % K;
        P(i)=P(i-1);

        P_CO2(i) = 0;
        P_H2O(i) = 0;
        P_N2(i) = 0;
        P_O2(i) = 0;
        P_CH4(i) = P(i);
        
        y_CO2(i)=0;
        y_H2O(i)=0;
        y_N2(i)=0;
        y_O2(i)=0;
        y_CH4(i)=1;
        
        h_CO2(i) = 0;
        h_H2O(i) = 0;
        h_N2(i) = 0;
        h_O2(i) = 0;
        h_CH4(i) = h_CH4_func(R_uni,T(i));
        
        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);
        
        s_CO2(i) = 0;
        s_H2O(i) = 0;
        s_N2(i) = 0;
        s_O2(i) = 0;
        s_CH4(i) = s_CH4_func(R_uni,P(i),T(i),y_CH4(i));

        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);

        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = 0;
        ef_H2O(i) = 0;
        ef_N2(i) = 0;
        ef_O2(i) = 0;
        ef_CH4(i) = ef_CH4_func(h_CH4(i),s_CH4(i),y_CH4(i))+831650;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end

    %% STATE 5
    if i==9
        
        T(i) = 1000; % K;
        P(i) = P(i-2)*0.95; % kPa

        x = P_H2O(i-2)/(P(i-2)-P_H2O(i-2));

        y_CO2(i)=1/(7.52+1+2+2*4.76*x);
        y_H2O(i)=(2+2*4.76*x)/(7.52+1+2+2*4.76*x);
        y_N2(i)=7.52/(7.52+1+2+2*4.76*x);
        y_O2(i)=0;
        y_CH4(i)=0;
        
        P_CO2(i) = P(i)*y_CO2(i);
        P_H2O(i) = P(i)*y_H2O(i);
        P_N2(i) = P(i)*y_N2(i);
        P_O2(i) = 0;
        P_CH4(i) = 0;

        h_CO2(i) = h_CO2_func(R_uni,T(i));
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = 0;
        h_CH4(i) = 0;

        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);

        s_CO2(i) = s_CO2_func(R_uni,P(i),T(i),y_CO2(i));
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = 0;
        s_CH4(i) = 0;

        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);

        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = ef_CO2_func(h_CO2(i),s_CO2(i),y_CO2(i));
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = 0;
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end

    %% STATE 6s
    if i==10
       
        P(i) = P(i-1)/(comp_ratio^2*0.95);
        
        P_CO2(i) = P_CO2(i-1)/(comp_ratio^2*0.95);
        P_H2O(i) =P_H2O(i-1)/(comp_ratio^2*0.95);
        P_N2(i) = P_N2(i-1)/(comp_ratio^2*0.95);
        P_O2(i) = 0;
        P_CH4(i) = 0;
        
        y_CO2(i)=P_CO2(i)/P(i);
        y_H2O(i)=P_H2O(i)/P(i);
        y_N2(i)=P_N2(i)/P(i);
        y_O2(i)=0;
        y_CH4(i)=0;
        
        s_total(i)=s_total(i-1);
        
        T(i) = T_find_s(s_total(i),P(i),y_CO2(i),y_H2O(i),y_N2(i),y_O2(i),y_CH4(i));
        
        h_CO2(i) = h_CO2_func(R_uni,T(i));
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = 0;
        h_CH4(i) = 0;
        
        h_total(i) = y_CO2(i)*h_CO2(i)+y_H2O(i)*h_H2O(i)+y_N2(i)*h_N2(i)+y_O2(i)*h_O2(i)+y_CH4(i)*h_CH4(i);

        s_CO2(i) = s_CO2_func(R_uni,P(i),T(i),y_CO2(i));
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = 0;
        s_CH4(i) = 0;

        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = ef_CO2_func(h_CO2(i),s_CO2(i),y_CO2(i));
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = 0;
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end
    %% STATE 6
    if i==11
        
        h_total(i) = (h_total(i-2)-h_total(i-2)*n_t+h_total(i-1)*n_t);
        
        P(i)=P(i-1);

        P_CO2(i) = P_CO2(i-1);
        P_H2O(i) =P_H2O(i-1);
        P_N2(i) = P_N2(i-1);
        P_O2(i) = 0;
        P_CH4(i) = 0;
        
        y_CO2(i)=y_CO2(i-1);
        y_H2O(i)=y_H2O(i-1);
        y_N2(i)=y_N2(i-1);
        y_O2(i)=0;
        y_CH4(i)=0;
        
        T(i) = T_find_h(h_total(i),y_CO2(i),y_H2O(i),y_N2(i),y_O2(i),y_CH4(i));
        
        s_CO2(i) = s_CO2_func(R_uni,P(i),T(i),y_CO2(i));
        s_H2O(i) = s_H2O_func(R_uni,P(i),T(i),y_H2O(i));
        s_N2(i) = s_N2_func(R_uni,P(i),T(i),y_N2(i));
        s_O2(i) = 0;
        s_CH4(i) = 0;

        h_CO2(i) = h_CO2_func(R_uni,T(i));
        h_H2O(i) = h_H2O_func(R_uni,T(i));
        h_N2(i) = h_N2_func(R_uni,T(i));
        h_O2(i) = 0;
        h_CH4(i) = 0;

        s_total(i) = y_CO2(i)*s_CO2(i)+y_H2O(i)*s_H2O(i)+y_N2(i)*s_N2(i)+y_O2(i)*s_O2(i)+y_CH4(i)*s_CH4(i);

        P_sat(i)=P_sat_func(T(i));
        rel_hum(i)=P_H2O(i)/P_sat(i);

        ef_CO2(i) = ef_CO2_func(h_CO2(i),s_CO2(i),y_CO2(i));
        ef_H2O(i) = ef_H2O_func(h_H2O(i),s_H2O(i),y_H2O(i));
        ef_N2(i) = ef_N2_func(h_N2(i),s_N2(i),y_N2(i));
        ef_O2(i) = 0;
        ef_CH4(i) = 0;

        ef_total(i) = ef_CO2(i)*y_CO2(i)+ef_H2O(i)*y_H2O(i)+ef_N2(i)*y_N2(i)+ef_O2(i)*y_O2(i)+ef_CH4(i)*y_CH4(i);
    end
end
rel_hum(:) = rel_hum(:)*100;

ndot=zeros(1,num_state);
W_turbine = 175*10^3; % kW

ndot(11) = W_turbine/(h_total(9)-h_total(11));
ndot(10) = ndot(11);
ndot(9) = ndot(10);
ndot(8) = ndot(9)*y_CO2(9);
ndot(7) = ndot(11)-ndot(8);
ndot(6) = ndot(7);
ndot(5) = ndot(6);
ndot(4) = ndot(5)*(y_H2O(5)-y_H2O(3));
ndot(3) = ndot(5)-ndot(4);
ndot(2) = ndot(3);
ndot(1) = ndot(2);

%% PROCESS

num_process=5;

%% Qdot
Qdot=zeros(1,num_process);
Qcycle = 0;
for i=1:num_process
    if i==1
        Qdot(i)=0;
    elseif i==2
        Qdot(i)=0;
    elseif i==3
        Qdot(i)=0;
    elseif i==4
        Qdot(i)=ndot(i+5)*h_total(i+5)-ndot(i+3)*h_total(i+3)-ndot(i+4)*h_total(i+4);
    elseif i==5
        Qdot(i)=0;
    end
    Qcycle = Qcycle+Qdot(i);
end

%% Wdot
Wdot=zeros(1,num_process);
Wcycle = 0;
for i=1:num_process
    if i==1
        Wdot(i)=ndot(i)*(h_total(i)-h_total(i+2));
    elseif i==2
        Wdot(i)=0;
    elseif i==3
        Wdot(i)=ndot(i+2)*(h_total(i+2)-h_total(i+4));
    elseif i==4
        Wdot(i)=0;
    elseif i==5
        Wdot(i)=175*10^3;
    end
    Wcycle = Wcycle+Wdot(i);
end

%% Sdot_gen
Sdot_gen=zeros(1,num_process);
Sdot_cycle = 0;
for i=1:num_process
    if i==1
        Sdot_gen(i)=ndot(i)*s_total(i+2)-ndot(i)*s_total(i);
    elseif i==2
        Sdot_gen(i)=ndot(i+3)*s_total(i+3)-ndot(i+1)*s_total(i+1)-ndot(i+2)*s_total(i+2);
    elseif i==3
        Sdot_gen(i)=ndot(i+4)*s_total(i+4)-ndot(i+2)*s_total(i+2);
    elseif i==4
        Sdot_gen(i)=ndot(i+5)*s_total(i+5)-ndot(i+4)*s_total(i+4)-ndot(i+3)*s_total(i+3)-Qdot(i)/T(i+5);
    elseif i==5
        Sdot_gen(i)=ndot(i+6)*s_total(i+6)-ndot(i+4)*s_total(i+4);
    end
    Sdot_cycle = Sdot_cycle + Sdot_gen(i);
end

%% ex_flow_diff
ex_flow_diff=zeros(1,num_process);
ex_flow_cycle = 0;
for i=1:num_process
    if i==1
        ex_flow_diff(i) = ndot(i+2)*ef_total(i+2)-ndot(i)*ef_total(i);
    elseif i==2
        ex_flow_diff(i) = ndot(i+3)*ef_total(i+3)-ndot(i+1)*ef_total(i+1)-ndot(i+2)*ef_total(i+2);
    elseif i==3
        ex_flow_diff(i) = ndot(i+4)*ef_total(i+4)-ndot(i+2)*ef_total(i+2);
    elseif i==4
        ex_flow_diff(i) = ndot(i+5)*ef_total(i+5)-ndot(i+4)*ef_total(i+4)-ndot(i+3)*ef_total(i+3);
    elseif i==5
        ex_flow_diff(i) = ndot(i+6)*ef_total(i+6)-ndot(i+4)*ef_total(i+4);
    end
    ex_flow_cycle = ex_flow_cycle + ex_flow_diff(i);
end

%% ex_dest
ex_dest=zeros(1,num_process);
ex_dest_cycle = 0;
for i=1:num_process
    if i==1
        ex_dest(i) = -ex_flow_diff(i)-Wdot(i)+(1-T_0/T(i))*Qdot(i);
    elseif i==2
        ex_dest(i) = -ex_flow_diff(i)-Wdot(i)+(1-T_0/T(i))*Qdot(i);
    elseif i==3
        ex_dest(i) = -ex_flow_diff(i)-Wdot(i)+(1-T_0/T(i))*Qdot(i);
    elseif i==4
        ex_dest(i) = -ex_flow_diff(i)-Wdot(i)+(1-T_0/T(i))*Qdot(i);
    elseif i==5
        ex_dest(i) = -ex_flow_diff(i)-Wdot(i)+(1-T_0/T(i+4))*Qdot(i);
    end
    ex_dest_cycle = ex_dest_cycle + ex_dest(i);
end


%% eta
eta_array=zeros(1,num_process);
for i=1:num_process
    if i==1
        eta_array(i) = ndot(i)*(ef_total(i+1)-ef_total(i)) / -Wdot(i);
    elseif i==2
        eta_array(i) = 0;
    elseif i==3
        eta_array(i) = ndot(i+2)*(ef_total(i+4)-ef_total(i+2)) / -Wdot(i);
    elseif i==4
        eta_array(i) = ndot(i+5)*ef_total(i+5)/(ndot(i+3)*ef_total(i+3)+ndot(i+4)*ef_total(i+4));
    elseif i==5
        eta_array(i) = Wdot(i) / (ndot(i+5)*(ef_total(i+4)-ef_total(i+6))) ;
    end
end
eta_cycle = Wcycle/(ndot(i+3)*ef_total(i+3));

%% STATE TABLE

State = {'1';'2';'2L';'3';'4';'4F';'5';'6'};
Description = {'Low Pressure Compressor Inlet';'Intercooler Air Inlet';'Intercooler Water Inlet';'High Pressure Compressor Inlet';'Combustor Air Inlet';'Combustor Fuel Inlet';'Turbine Inlet';'Turbine Outlet'};
T_table = [T(1);T(3);T(4);T(5);T(7);T(8);T(9);T(11)];
P_table = [P(1);P(3);P(4);P(5);P(7);P(8);P(9);P(11)];
n_dot_table = [ndot(1);ndot(3);ndot(4);ndot(5);ndot(7);ndot(8);ndot(9);ndot(11)];
rel_hum_table = [rel_hum(1);rel_hum(3);rel_hum(4);rel_hum(5);rel_hum(7);rel_hum(8);rel_hum(9);rel_hum(11)];
h_table = [h_total(1);h_total(3);h_total(4);h_total(5);h_total(7);h_total(8);h_total(9);h_total(11)];
s_table = [s_total(1);s_total(3);s_total(4);s_total(5);s_total(7);s_total(8);s_total(9);s_total(11)];
ef_table = [ef_total(1);ef_total(3);ef_total(4);ef_total(5);ef_total(7);ef_total(8);ef_total(9);ef_total(11)];
TA1 = table(State,Description,T_table,P_table,n_dot_table,rel_hum_table,h_table,s_table,ef_table)


%% PROCESS TABLE
eta = 1;

Process = {'1->2';'2 + 2L ->3';'3->4';'4 + 4F->5';'5->6';'Cycle'};
Description = {'Low Pressure Compressor';'Intercooler';'High Pressure Compressor';'Combustor';'Turbine';'Entire Cycle'};
Q_dot = [Qdot(1);Qdot(2);Qdot(3);Qdot(4);Qdot(5);Qcycle];
W_dot = [Wdot(1);Wdot(2);Wdot(3);Wdot(4);Wdot(5);Wcycle];
sigma_gen_dot = [Sdot_gen(1);Sdot_gen(2);Sdot_gen(3);Sdot_gen(4);Sdot_gen(5);Sdot_cycle];
del_ef = [ex_flow_diff(1);ex_flow_diff(2);ex_flow_diff(3);ex_flow_diff(4);ex_flow_diff(5);ex_flow_cycle];
ex_dest_table = [ex_dest(1);ex_dest(2);ex_dest(3);ex_dest(4);ex_dest(5);ex_dest_cycle];
eta_table = {eta_array(1);'-';eta_array(3);eta_array(4);eta_array(5);eta_cycle};
TA2 = table(Process,Description,Q_dot,W_dot,sigma_gen_dot,del_ef,ex_dest_table,eta_table)

T_min = 700;
T_max = 3000;
T_avg = (T_min+T_max)/2;
        
    while true
       
        

        h_CO2_e = h_CO2_func(R_uni,T_avg);
        h_H2O_e = h_H2O_func(R_uni,T_avg);
        h_N2_e = h_N2_func(R_uni,T_avg);
        h_O2_e = 0;
        h_CH4_e = 0;
        
        h_total_e = y_CO2(9)*h_CO2_e+y_H2O(9)*h_H2O_e+y_N2(9)*h_N2_e+y_O2(9)*h_O2_e+y_CH4(9)*h_CH4_e;
        h_total_e = ndot(9)*h_total_e;

        h_total_i = ndot(7)*h_total(7)+ndot(8)*h_total(8);

        error = 10^-3;

        if abs(h_total_e-h_total_i) < error
            T_adiabatic = T_avg;
            break
        elseif h_total_e < h_total_i
            T_min = T_avg;
            T_avg = (T_min+T_max)/2;
        elseif h_total_e > h_total_i
            T_max = T_avg;
            T_avg = (T_min+T_max)/2;
        end
    end
T_adiabatic %K


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

