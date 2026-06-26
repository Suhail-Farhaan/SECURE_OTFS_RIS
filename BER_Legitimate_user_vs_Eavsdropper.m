tic
clc;
% clear variables;
close all;
%OTFS Frame Parameters
M=32;%Number of SubCarriers
N=4;%Number of Time Symbols
N_frames=1e5;
Pd=2;%Number of DD Paths in LOS Link
P1=3;%Number of DD Paths in the first Link
P2=2;%Number of DD Paths in the second Link
fc=(4*10^9);%carrier frequency
delta_f=15e3;%Sub_Carrier Spacing
T=1/delta_f;
U=8;%Number of RIS Elements
%% li and ki random (synthetic channel) 
%Delay_Doppler Corresponding to Direct Link
% l_d=randi([0,Pd-1],1,Pd);%l_d=[0,1];
% k_d=randi([0,Pd-1],1,Pd);%k_d=[0,1];
% tau=l_d./(M*delta_f);
% neu=k_d./(N*T);
% %Delay_Doppler Corresponding to First Link
% l_g=randi([0,P1-1],1,P1);%l_g=[0,1];
% k_g=randi([0,P1-1],1,P1);%k_g=[1,0];
% tau_1=l_g./(M*delta_f);
% neu_1=k_g./(N*T);
% %Delay_Doppler Corresponding to Second Link
% l_f=randi([0,P2-1],1,P2);%l_f=[1,1];
% k_f=randi([0,P2-1],1,P2);%k_f=[0,0];
% tau_2=l_f./(M*delta_f);
% neu_2=k_f./(N*T);
% 
% tau_eff=[tau_1(1)+tau_2(1),tau_1(1)+tau_2(2),tau_1(2)+tau_2(1),tau_1(2)+tau_2(2),tau_1(3)+tau_2(1),tau_1(3)+tau_2(2)];
% neu_eff=[neu_1(1)+neu_2(1),neu_1(1)+neu_2(2),neu_1(2)+neu_2(1),neu_1(2)+neu_2(2),neu_1(3)+neu_2(1),neu_1(3)+neu_2(2)];
% l_cascaded=floor(tau_eff.*(M*delta_f));
% k_cascaded=floor(neu_eff.*(N*T));

% l_cascaded=[2,3,1,2,0,1];%M=4,N=4
% k_cascaded=[1,1,3,3,3,3];
% l_cascaded=[3,3,3,3,1,1];%M=8,N=4
% k_cascaded=[3,2,1,0,1,0];
% l_cascaded=[0,1,2,3,2,3];%M=16,N=4
% k_cascaded=[1,2,3,3,1,2];
% l_cascaded=[0,0,1,1,1,1];
% k_cascaded=[1,2,1,2,1,2];%M=32,N=8
% l_cascaded=[1,0,5,3,0,1];
% k_cascaded=[1,1,3,2,4,3];%M=16,N=4
% l_cascaded=randi([0,(P1*P2)-1],1,(P1*P2));
% k_cascaded=randi([0,(P1*P2)-1],1,(P1*P2));
%% Excess tap delay model,jakes formula for li and ki
tau1=[(2.08*10^-6),(5.20*10^-6),(8.328*10^-6)];  %Excess Tap Delay of Prem sir model
tau2=[(11.46*10^-6),(20.8*10^-6)];

l_g=[];
for i=1:P1
    l_g(i)=ceil(tau1(i)*M*delta_f);  
end
tau_1=l_g./(M*delta_f);

l_f=[];
for j=1:P2
    l_f(i)=ceil(tau2(j)*M*delta_f);  
end
tau_2=l_f./(M*delta_f);
tau_eff=[tau_1(1)+tau_2(1),tau_1(1)+tau_2(2),tau_1(2)+tau_2(1),tau_1(2)+tau_2(2),tau_1(3)+tau_2(1),tau_1(3)+tau_2(2)];
l_cascaded=ceil(tau_eff.*(M*delta_f));

l_d=[];
for i=1:Pd
    l_d(i)=ceil(tau2(i)*M*delta_f);  
end
%%
%Finding Doppler indices(Integer) Using Jakes
%Formula(v_i=v_max*(cos(theta_i)), where theta_i is uniformly distributed
%from (-pi to pi)

c=(3*10^8);%velocity of light
s=500;%speed in kmph(120kmph<-->33.33m/s)                            
v_max=(fc*s*(5/18))/c;     
alpha1=(rand(1,P1));          
alpha1=(2*pi*alpha1-pi);%uniformly distributed theta from -pi to pi
%histogram(alpha1);
k1=[];
for i=1:P1
    k1(i)=ceil(v_max*N*T*abs(cos(alpha1(i))));
end
neu_1=k1./(N*T);

alpha2=(rand(1,P2));          
alpha2=(2*pi*alpha2-pi);%uniformly distributed theta from -pi to pi
%histogram(alpha2);
k2=[];
for j=1:P2
    k2(j)=ceil(v_max*N*T*abs(cos(alpha2(j))));
end
neu_2=k2./(N*T);
neu_eff=[neu_1(1)+neu_2(1),neu_1(1)+neu_2(2),neu_1(2)+neu_2(1),neu_1(2)+neu_2(2),neu_1(3)+neu_2(1),neu_1(3)+neu_2(2)];
doppler_spread=2*v_max;
k_cascaded=ceil(neu_eff.*(N*T));

alpha3=(rand(1,Pd));
k_d=[];
for j=1:P2
    k_d(j)=ceil(v_max*N*T*abs(cos(alpha3(j))));
end
%%
ze=zeros(1,M*N);
pie=circshift(eye(M*N),1,1);%Pie Matrix
z1=1j*2*pi; z2=M*N;
z=exp(z1/z2);
for i=1:M*N
    ze(i)=z.^(i-1);
end
Delta=diag(ze);%Delta Matrix
% m=randi([0 1],1,2*M*N); %Message(QPSK)
G_rx=eye(M);
G_tx=eye(M);
F_M=(1/sqrt(M))*dftmtx(M);
F_N=(1/sqrt(N))*dftmtx(N);
g=sqrt(1/(2*U*P1))*(randn((U*P1),1)+sqrt(-1)*randn((U*P1),1));%Fading Coefficient of the First Link i.e Source to RIS g~CN(0,1/Lg)
fc=sqrt(1/(2*P2))*(randn(P2,1)+sqrt(-1)*randn(P2,1));%Fading Coefficient of the Second Link i.e RIS to User f~CN(0,1/Lf)
hd=sqrt(1/2*Pd)*(randn(Pd,1)+sqrt(-1)*randn(Pd,1));%Fading Coefficient of the Direct Link i.e Source to User h~CN(0,1/Ld)
% h=[h1;zeros((Q*Lg*Lf)-length(h1),1)];
h_casc=kron(g,fc);%Cascaded RIS Channel with Q*Lg*Lf paths
h_q=reshape(h_casc,U,P1*P2);


H0=zeros(M*N);%Initialising Channel Matrix for non RIS
for i=1:Pd
    H0=H0+hd(i)*pie^(l_d(i))*Delta^(k_d(i));
end
% H0=kron(dftmtx(N),G_rx)*H0*kron(transpose(conj(dftmtx(N))),G_tx);
H1=zeros(M*N);
for i=1:(P1*P2)
    for q=1:6
        H1=H1+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
% H1=kron(F_N,G_rx)*H1*kron(transpose(conj(F_N)),G_tx);
H2=zeros(M*N);
for i=1:(P1*P2)
    for q=7:12
        H2=H2+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
% H2=kron(F_N,G_rx)*H2*kron(transpose(conj(F_N)),G_tx);
H3=zeros(M*N);
for i=1:(P1*P2)
    for q=13:18
        H3=H3+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
% H3=kron(F_N,G_rx)*H3*kron(transpose(conj(F_N)),G_tx);
H4=zeros(M*N);
for i=1:(P1*P2)
    for q=19:24
        H4=H4+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
% H4=kron(F_N,G_rx)*H4*kron(transpose(conj(F_N)),G_tx);
H5=zeros(M*N);
for i=1:(P1*P2)
    for q=25:30
        H5=H5+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end

H6=zeros(M*N);
for i=1:(P1*P2)
    for q=31:36
        H6=H6+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end

H7=zeros(M*N);
for i=1:(P1*P2)
    for q=36:42
        H7=H7+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end

H8=zeros(M*N);
for i=1:(P1*P2)
    for q=43:48
        H8=H8+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
%% Strongest DDCR Phase Optimization
for p=1:P1*P2
    for d=1:Pd
        L_ang_1(p)=(sum(abs(h_q(1,p)))+abs(hd(d)))^2;
        L_ang_2(p)=(sum(abs(h_q(2,p)))+abs(hd(d)))^2;
        L_ang_3(p)=(sum(abs(h_q(3,p)))+abs(hd(d)))^2;
        L_ang_4(p)=(sum(abs(h_q(4,p)))+abs(hd(d)))^2;
        L_ang_5(p)=(sum(abs(h_q(5,p)))+abs(hd(d)))^2;
        L_ang_6(p)=(sum(abs(h_q(6,p)))+abs(hd(d)))^2;
        L_ang_7(p)=(sum(abs(h_q(7,p)))+abs(hd(d)))^2;
        L_ang_8(p)=(sum(abs(h_q(8,p)))+abs(hd(d)))^2;
    end
end
[~,k1]=max(L_ang_1);
[~,k2]=max(L_ang_2);
[~,k3]=max(L_ang_3);
[~,k4]=max(L_ang_4);
[~,k5]=max(L_ang_5);
[~,k6]=max(L_ang_6);
[~,k7]=max(L_ang_7);
[~,k8]=max(L_ang_8);

t_star_1=angle(h_q(1,k1))-angle(max(abs(hd)));
t_star_2=angle(h_q(2,k2))-angle(max(abs(hd)));
t_star_3=angle(h_q(3,k3))-angle(max(abs(hd)));
t_star_4=angle(h_q(4,k4))-angle(max(abs(hd)));
t_star_5=angle(h_q(5,k5))-angle(max(abs(hd)));
t_star_6=angle(h_q(6,k6))-angle(max(abs(hd)));
t_star_7=angle(h_q(7,k7))-angle(max(abs(hd)));
t_star_8=angle(h_q(8,k8))-angle(max(abs(hd)));
t_star=[t_star_1,t_star_2,t_star_3,t_star_4,t_star_5,t_star_6,t_star_7,t_star_8];

t_set=0:(2*pi)/8:(14*pi)/8;
for i=1:length(t_set)
    val1(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(1)));
    val2(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(2)));
    val3(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(3)));
    val4(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(4)));
    val5(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(5)));
    val6(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(6)));
    val7(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(7)));
    val8(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(8)));
end
[~,v1]=min(val1);
[~,v2]=min(val2);
[~,v3]=min(val3);
[~,v4]=min(val4);
[~,v5]=min(val5);
[~,v6]=min(val6);
[~,v7]=min(val7);
[~,v8]=min(val8);

t1=[t_set(v1),t_set(v2),t_set(v3),t_set(v4)];
THETA1=[exp(1j*t1(1))*eye(M*N),exp(1j*t1(2))*eye(M*N),exp(1j*t1(3))*eye(M*N),exp(1j*t1(4))*eye(M*N)]';%THETA USING DDCR PHASE OPTIMIZATION
H_eff_1=kron(dftmtx(N),G_rx)*(H0+[H1,H2,H3,H4]*THETA1)*kron(transpose(conj(dftmtx(N))),G_tx);
[Ub1,Eb1,Vb1_H]=svd(H_eff_1);
Vb1=transpose(conj(Vb1_H));
Ub1_H=transpose(conj(Ub1));

t2=[t_set(v1),t_set(v2),t_set(v3),t_set(v4),t_set(v5),t_set(v6),t_set(v7),t_set(v8)];
THETA2=[exp(1j*t2(1))*eye(M*N),exp(1j*t2(2))*eye(M*N),exp(1j*t2(3))*eye(M*N),exp(1j*t2(4))*eye(M*N),exp(1j*t2(5))*eye(M*N),exp(1j*t2(6))*eye(M*N),exp(1j*t2(7))*eye(M*N),exp(1j*t2(8))*eye(M*N)]';%THETA USING DDCR PHASE OPTIMIZATION
H_eff_2=kron(dftmtx(N),G_rx)*(H0+[H1,H2,H3,H4,H5,H6,H7,H8]*THETA2)*kron(transpose(conj(dftmtx(N))),G_tx);
[Ub2,Eb2,Vb2_H]=svd(H_eff_2);
Vb2=transpose(conj(Vb2_H));
Ub2_H=transpose(conj(Ub2));
%% Projection Matrices
n=(1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N,1));% Auxilary Artificial Noise
Q1=H_eff_1*inv(transpose(conj(H_eff_1))*H_eff_1)*transpose(conj(H_eff_1));
T_1=(eye(M*N)-Q1)/norm((eye(M*N)-Q1));
o1=H_eff_1*T_1*n;
Q2=H_eff_2*inv(transpose(conj(H_eff_2))*H_eff_2)*transpose(conj(H_eff_2));
T_2=(eye(M*N)-Q2)/norm((eye(M*N)-Q2));
o2=H_eff_2*T_2*n;
%modulation size
mod_size=16;
% number of information symbols in one frame
N_syms_per_frame=N*M;
% number of information bits in one frame
N_bits_per_frame=N_syms_per_frame*log2(mod_size);
ber_IRS_legit_1 = [];    % Empty array to store BER with IRS.
ber_IRS_eavs_1 = [];    % Empty array to store BER with IRS.
ber_IRS_legit_2 = [];    % Empty array to store BER with IRS.
ber_IRS_eavs_2 = []; 
for SNRDB = -25:5:10             % Assigning SNR_db values
    SNR = 10^(SNRDB/10);         % Conversion from dB to linear scale
    ber_avg_IRS_legit_1=[];
    ber_avg_IRS_legit_2=[];
    ber_avg_IRS_eavs_1=[];
    ber_avg_IRS_eavs_2=[];
    for b=1:N_frames
        % generate random bits
        tx_info_bits=randi([0,1],N_bits_per_frame,1);
        % QAM modulation
        d=qammod(tx_info_bits,mod_size,'gray','InputType','bit');
        x1=Vb1*d;
        x2=Vb2*d;
        sigma_IRS_1 = sqrt(1./(((U/2)^2)*SNR));
        sigma_IRS_2 = sqrt(1./((U^2)*SNR));
        n_i_1=kron(transpose(conj(F_N)),G_tx)*sigma_IRS_1*(1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N,1)); % Noise at Rxr for IRS with optimised phase
        n_i_2=kron(transpose(conj(F_N)),G_tx)*sigma_IRS_2*(1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N,1)); % Noise at Rxr for IRS with optimised phase
        
        yl1=(Ub1*Eb1*Vb1_H)*x1+n_i_1;
        yl_tilde1=Ub1_H*yl1;
        db1=inv(Eb1)*yl_tilde1;
        detect1=qamdemod(db1,mod_size,'gray','OutputType','bit');

        ber1=0; %Calculating BER for RIS Legitimate User
        for i=1:length(tx_info_bits)
            if(tx_info_bits(i)~=detect1(i))
                ber1=ber1+1;
            end
        end
        ber_avg_IRS_legit_1=[ber_avg_IRS_legit_1 ber1];

        yl2=(Ub2*Eb2*Vb2_H)*x2+n_i_2;
        yl_tilde2=Ub2_H*yl2;
        db2=inv(Eb2)*yl_tilde2;
        detect2=qamdemod(db2,mod_size,'gray','OutputType','bit');

        ber2=0; %Calculating BER for RIS Legitimate User
        for i=1:length(tx_info_bits)
            if(tx_info_bits(i)~=detect2(i))
                ber2=ber2+1;
            end
        end
        ber_avg_IRS_legit_2=[ber_avg_IRS_legit_2 ber2];

        ye1=(Ue1*Ee1*Ve1_H)*x1+n_i_1;
        ye1_tilde=Ue1_H*ye1;
        de1=Ve1*inv(Ee1)*ye1_tilde;
        detected1=qamdemod(de1,mod_size,'gray','OutputType','bit');

        ber3=0; %Calculating BER for RIS Eavsdropper
        for i=1:length(tx_info_bits)
            if(tx_info_bits(i)~=detected1(i))
                ber3=ber3+1;
            end
        end
        ber_avg_IRS_eavs_1=[ber_avg_IRS_eavs_1 ber3];

        ye2=(Ue2*Ee2*Ve2_H)*x2+n_i_2;
        ye2_tilde=Ue2_H*ye2;
        de2=Ve2*inv(Ee2)*ye2_tilde;
        detected2=qamdemod(de2,mod_size,'gray','OutputType','bit');

        ber4=0; %Calculating BER for RIS Eavsdropper
        for i=1:length(tx_info_bits)
            if(tx_info_bits(i)~=detected2(i))
                ber4=ber4+1;
            end
        end
        ber_avg_IRS_eavs_2=[ber_avg_IRS_eavs_2 ber4];

    end
    BER_RIS_legit = mean(ber_avg_IRS_legit_1);
    No = 4*M*N;
    ber_IRS_legit_1 = [ber_IRS_legit_1 BER_RIS_legit/No];

    BER_RIS_legit_2 = mean(ber_avg_IRS_legit_2);
    No = 4*M*N;
    ber_IRS_legit_2 = [ber_IRS_legit_2 BER_RIS_legit_2/No];

    BER_RIS_eavs = mean(ber_avg_IRS_eavs_1);
    No = 4*M*N;
    ber_IRS_eavs_1 = [ber_IRS_eavs_1 BER_RIS_eavs/No];

    BER_RIS_eavs_2 = mean(ber_avg_IRS_eavs_2);
    No = 4*M*N;
    ber_IRS_eavs_2 = [ber_IRS_eavs_2 BER_RIS_eavs_2/No];
end
%%
SNRDB = -25:5:10;
% figure(1)
semilogy(SNRDB, ber_IRS_legit_1, 'b','Linewidth',2,LineStyle='-.');% Plotting BER vs SNRDB with RIS Optimised phase
hold on;
semilogy(SNRDB, ber_IRS_legit_2, 'bo-','Linewidth',2);% Plotting BER vs SNRDB with RIS Optimised phase
hold on;
semilogy(SNRDB, ber_IRS_eavs_1, 'r','Linewidth',2,LineStyle='-.');% Plotting BER vs SNRDB with RIS Optimised phase
hold on;
semilogy(SNRDB, ber_IRS_eavs_2, 'r*-','Linewidth',2);
xlabel('SNR(dB)');
ylabel('Bit Error Rate');
title('BER curve for 16QAM with Single Tap Equalization (M=32, N=8, P1=3,P2=2)');
legend('Legitimate User(U=4)','Legitimate User(U=8)','Eavesdropper(U=4)','Eavesdropper(U=8)',Location='southwest');
grid on;
% figure(2)
% scatterplot(db);
% title('Constellation Received by Legitimate User');
% figure(3)
% scatterplot(de);
% title('Constellation Received by Eavsdropper');
toc
