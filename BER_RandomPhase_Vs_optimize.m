tic
clc;
clear variables;
close all;
%OTFS Frame Parameters
M=32;%Number of SubCarriers
N=8;%Number of Time Symbols
N_frames=1e4;%Number of OTFS Frames
Pd=2;%Number of DD Paths in LOS Link
P1=3;%Number of DD Paths in the first Link
P2=2;%Number of DD Paths in the second Link
delta_f=15e3;%Sub_Carrier Spacing
T=1/delta_f;
U=4;%Number of RIS Elements
%Delay_Doppler Corresponding to Direct Link
l_d=randi([0,Pd-1],1,Pd);%l_d=[0,1];
k_d=randi([0,Pd-1],1,Pd);%k_d=[0,1];
tau=l_d./(M*delta_f);
neu=k_d./(N*T);
%Delay_Doppler Corresponding to First Link
l_g=randi([0,P1-1],1,P1);%l_g=[0,1];
k_g=randi([0,P1-1],1,P1);%k_g=[1,0];
tau_1=l_g./(M*delta_f);
neu_1=k_g./(N*T);
%Delay_Doppler Corresponding to Second Link
l_f=randi([0,P2-1],1,P2);%l_f=[1,1];
k_f=randi([0,P2-1],1,P2);%k_f=[0,0];
tau_2=l_f./(M*delta_f);
neu_2=k_f./(N*T);
%Delay_Doppler Corresponding to Cascaded Link
% l_cascaded=[l_g(1)+l_f(1),l_g(1)+l_f(2),l_g(2)+l_f(1),l_g(2)+l_f(2)];
% k_cascaded=[k_g(1)+k_f(1),k_g(1)+k_f(2),k_g(2)+k_f(1),k_g(2)+k_f(2)];
tau_eff=[tau_1(1)+tau_2(1),tau_1(1)+tau_2(2),tau_1(2)+tau_2(1),tau_1(2)+tau_2(2),tau_1(3)+tau_2(1),tau_1(3)+tau_2(2)];
neu_eff=[neu_1(1)+neu_2(1),neu_1(1)+neu_2(2),neu_1(2)+neu_2(1),neu_1(2)+neu_2(2),neu_1(3)+neu_2(1),neu_1(3)+neu_2(2)];
l_cascaded=floor(tau_eff.*(M*delta_f));
k_cascaded=floor(neu_eff.*(N*T));
% l_cascaded=[2,2,2,2,1,1];
% k_cascaded=[1,0,1,0,1,0];
ze=zeros(1,M*N);
pie=circshift(eye(M*N),1,1);%Pie Matrix
z1=1j*2*pi; z2=M*N;
z=exp(z1/z2);
for i=1:M*N
    ze(i)=z.^(i-1);
end
Delta=diag(ze);%Delta Matrix
G_rx=eye(M);
G_tx=eye(M);
F_M=(1/sqrt(M))*dftmtx(M);
F_N=(1/sqrt(N))*dftmtx(N);
g=sqrt(1/(2*U*P1))*(randn((U*P1),1)+sqrt(-1)*randn((U*P1),1));%Fading Coefficient of the First Link i.e Source to RIS g~CN(0,1/Lg)
f=sqrt(1/(2*P2))*(randn(P2,1)+sqrt(-1)*randn(P2,1));%Fading Coefficient of the Second Link i.e RIS to User f~CN(0,1/Lf)
hd=sqrt(1/2*Pd)*(randn(Pd,1)+sqrt(-1)*randn(Pd,1));%Fading Coefficient of the Direct Link i.e Source to User h~CN(0,1/Ld)
% h=[h1;zeros((Q*Lg*Lf)-length(h1),1)];
h_casc=kron(g,f);%Cascaded RIS Channel with Q*Lg*Lf paths
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
    for q=6:12
        H2=H2+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
% H2=kron(F_N,G_rx)*H2*kron(transpose(conj(F_N)),G_tx);
H3=zeros(M*N);
for i=1:(P1*P2)
    for q=12:18
        H3=H3+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
% H3=kron(F_N,G_rx)*H3*kron(transpose(conj(F_N)),G_tx);
H4=zeros(M*N);
for i=1:(P1*P2)
    for q=18:24
        H4=H4+h_casc(q)*pie^(l_cascaded(i))*Delta^(k_cascaded(i));
    end
end
% H4=kron(F_N,G_rx)*H4*kron(transpose(conj(F_N)),G_tx);

%% Strongest DDCR Phase Optimization
for p=1:P1*P2
    for d=1:Pd
        L_ang_1(p)=(sum(abs(h_q(1,p)))+abs(hd(d)))^2;
        L_ang_2(p)=(sum(abs(h_q(2,p)))+abs(hd(d)))^2;
        L_ang_3(p)=(sum(abs(h_q(3,p)))+abs(hd(d)))^2;
        L_ang_4(p)=(sum(abs(h_q(4,p)))+abs(hd(d)))^2;
    end
end
[~,k1]=max(L_ang_1);
[~,k2]=max(L_ang_2);
[~,k3]=max(L_ang_3);
[~,k4]=max(L_ang_4);

t_star_1=angle(h_q(1,k1))-angle(max(abs(hd)));
t_star_2=angle(h_q(2,k2))-angle(max(abs(hd)));
t_star_3=angle(h_q(3,k3))-angle(max(abs(hd)));
t_star_4=angle(h_q(4,k4))-angle(max(abs(hd)));
t_star=[t_star_1,t_star_2,t_star_3,t_star_4];

t_set=0:(2*pi)/8:(14*pi)/8;
for i=1:length(t_set)
    val1(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(1)));
    val2(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(2)));
    val3(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(3)));
    val4(i)=abs(exp(1j*t_set(i))-exp(1j*t_star(4)));
end
[~,v1]=min(val1);
[~,v2]=min(val2);
[~,v3]=min(val3);
[~,v4]=min(val4);
t=[t_set(v1),t_set(v2),t_set(v3),t_set(v4)];
THETA=[exp(1j*t(1))*eye(M*N),exp(1j*t(2))*eye(M*N),exp(1j*t(3))*eye(M*N),exp(1j*t(4))*eye(M*N)]';%THETA USING DDCR PHASE OPTIMIZATION
H_eff=kron(dftmtx(N),G_rx)*(H0+[H1,H2,H3,H4]*THETA)*kron(transpose(conj(dftmtx(N))),G_tx);
[UL,EL,VL_H]=svd(H_eff);
VL=transpose(conj(VL_H));
UL_H=transpose(conj(UL));
% H_eff=kron(F_N,G_rx)*H_eff*kron(transpose(conj(F_N)),G_tx);
% Generate random phase values in the range [0,2*pi]
phase = (2 * pi) * rand(U,1);
% Calculate the real and imaginary components of the vectors
realComponent = cos(phase);
imaginaryComponent = sin(phase);
theta = complex(realComponent, imaginaryComponent);%e power j theta
THETA_R=[theta(1)*eye(M*N),theta(2)*eye(M*N),theta(3)*eye(M*N),theta(4)*eye(M*N)]';%RIS Matrix for Random Phase Shift
H_ran=kron(F_N,G_rx)*(H0+[H1,H2,H3,H4]*THETA_R)*kron(transpose(conj(F_N)),G_tx);
[Ur,Er,Vr_H]=svd(H_ran);
Vr=transpose(conj(Vr_H));
Ur_H=transpose(conj(Ur));

[U0,E0,V0_H]=svd(H0);
V0=transpose(conj(V0_H));
U0_H=transpose(conj(U0));
%% Projection Matrices
w=(1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N, 1)); % Artificial Noise
Q1=H_eff*inv(transpose(conj(H_eff))*H_eff)*transpose(conj(H_eff));
T_1=(eye(M*N)-Q1)/norm((eye(M*N)-Q1));
o1=H_eff*T_1*w;

Q2=H_ran*inv(transpose(conj(H_ran))*H_ran)*transpose(conj(H_ran));
T_2=(eye(M*N)-Q2)/norm((eye(M*N)-Q2));
o2=H_ran*T_2*w;

Q3=H0*inv(transpose(conj(H0))*H0)*transpose(conj(H0));
T_3=(eye(M*N)-Q3)/norm((eye(M*N)-Q3));
o3=H0*T_3*w;
% H_ran=kron(F_N,G_rx)*H_ran*kron(transpose(conj(F_N)),G_tx);
%%
%modulation size
mod_size=16;
% number of information symbols in one frame
N_syms_per_frame=N*M;
% number of information bits in one frame
N_bits_per_frame=N_syms_per_frame*log2(mod_size);
ber_IRS_Ran = [];    % Empty array to store BER with IRS.
ber_IRS_DDCR = [];    % Empty array to store BER with IRS.
ber_NIRS = [];% Empty array to store BER W/o IRS.
for SNRDB = -25:5:10             % Assigning SNR_db values
    SNR = 10^(SNRDB/10);         % Conversion from dB to linear scale
    ber_avg_IRS_legit=[];
    ber_avg_IRS_legit_ran=[];
    ber_avg_NIRS=[];
    for b=1:N_frames
        % generate random bits
        tx_info_bits=randi([0,1],N_bits_per_frame,1);
        % QAM modulation
        d=qammod(tx_info_bits,mod_size,'gray','InputType','bit');
        x1=VL*d;
        x2=Vr*d;
        x3=V0*d;
        sigma = sqrt(1/SNR);   % Calculating variance
        sigma_IRS_1 = sqrt(1./((U^2)*SNR));
        sigma_IRS_2 = sqrt(1./((U^2)*SNR));
        n = kron(transpose(conj(F_N)),G_tx)*sigma * (1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N, 1)); % Noise
        n_i_1=kron(transpose(conj(F_N)),G_tx)*sigma_IRS_1*(1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N, 1)); % Noise at Rxr for IRS with random phase
        n_i_2=kron(transpose(conj(F_N)),G_tx)*sigma_IRS_2*(1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N,1)); % Noise at Rxr for IRS with optimised phase
        
        yl=(UL*EL*VL_H)*x1+o1+n_i_2;
        yl_tilde=UL_H*yl;
        dl=inv(EL)*yl_tilde;
        detect=qamdemod(dl,mod_size,'gray','OutputType','bit');

        ber=0; %Calculating BER for RIS Legitimate User(Optimized phase)
        for i=1:length(tx_info_bits)
            if(tx_info_bits(i)~=detect(i))
                ber=ber+1;
            end
        end
        ber_avg_IRS_legit=[ber_avg_IRS_legit ber];

        yr=(Ur*Er*Vr_H)*x2+o2+n_i_1;
        yr_tilde=Ur_H*yr;
        dr=inv(Er)*yr_tilde;
        detect_r=qamdemod(dr,mod_size,'gray','OutputType','bit');

        ber1=0; %%Calculating BER for RIS Legitimate User(Random phase)
        for i=1:length(tx_info_bits)
            if(tx_info_bits(i)~=detect_r(i))
                ber1=ber1+1;
            end
        end
        ber_avg_IRS_legit_ran=[ber_avg_IRS_legit_ran ber1];

        y0=(U0*E0*V0_H)*x3+o3+n;
        y0_tilde=U0_H*y0;
        d0=inv(E0)*y0_tilde;
        detected=qamdemod(d0,mod_size,'gray','OutputType','bit');

        ber2=0; %Calculating BER for NRIS Legitimate User
        for i=1:length(tx_info_bits)
            if(tx_info_bits(i)~=detected(i))
                ber2=ber2+1;
            end
        end
        ber_avg_NIRS=[ber_avg_NIRS ber2];

    end
    BER_mmse_DDCR = mean(ber_avg_IRS_legit);
    No = M*N;
    ber_IRS_DDCR = [ber_IRS_DDCR BER_mmse_DDCR/No];

    BER_mmse_Ran = mean(ber_avg_IRS_legit_ran);
    No = M*N;
    ber_IRS_Ran = [ber_IRS_Ran BER_mmse_Ran/No];

    BER_mmse_Without = mean(ber_avg_NIRS);
    No = M*N;
    ber_NIRS = [ber_NIRS BER_mmse_Without/No];
end
%%
SNRDB = -25:5:10;
semilogy(SNRDB, ber_NIRS, 'black*-','Linewidth',2);% Plotting BER vs SNRDB without RIS
hold on;
semilogy(SNRDB, ber_IRS_Ran, 'r*-','Linewidth',2);% Plotting BER vs SNRDB with RIS Random phase
hold on;
semilogy(SNRDB, ber_IRS_DDCR, 'b*-','Linewidth',2);% Plotting BER vs SNRDB with RIS Optimised phase
xlabel('SNR(dB)');
ylabel('Bit Error Rate');
title('BER curve for 16-QAM with Single Tap Equalization (M=4, N=4, P1=3,P2=2)');
legend('OTFS W/o RIS','RIS Aided OTFS(Random Phase)','RIS Aided OTFS(Optimised Phase(SDDCR))',Location='southwest');
grid on;
% figure(1);
% spy(H_eff);
% title('Non zero element positions in H_{eff} Matrix without phase optimisation');
% xlabel('Column Index');
% ylabel('Row Index');
toc
