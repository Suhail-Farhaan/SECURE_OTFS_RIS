tic
M=4;
N=4;
Q=4;
G_rx=eye(M);
G_tx=eye(M);
F_M=(1/sqrt(M))*dftmtx(M);
F_N=(1/sqrt(N))*dftmtx(N);
N_frames=1e4;
%modulation size
mod_size=16;
% number of information symbols in one frame
N_syms_per_frame=N*M;
% number of information bits in one frame
N_bits_per_frame=N_syms_per_frame*log2(mod_size);
% d_eavs=[];
% d_legit=[];
SNRDB=20;
SNR = 10^(SNRDB/10);         % Conversion from dB to linear scale
for b=1:N_frames
    % generate random bits
    tx_info_bits=randi([0,1],N_bits_per_frame,1);
    % QAM modulation
    d=qammod(tx_info_bits,mod_size,'gray','InputType','bit');
    x=Vb*d;
    sigma_IRS_2 = sqrt(1/((Q^2)*SNR));
    n_i_2=kron(transpose(conj(F_N)),G_tx)*sigma_IRS_2*(1/sqrt(2)) * (randn(M*N, 1) + 1j*randn(M*N,1)); % Noise at Rxr for IRS with optimised pha
    ye=(Ue*Ee*Ve_H)*x+n_i_2;
    ye_tilde=Ue_H*ye;
    de=inv(Ee)*ye_tilde;
    
    yl=(Ub*Eb*Vb_H)*x+n_i_2;
    yl_tilde=Ub_H*yl;
    db=inv(Eb)*yl_tilde;
    % Store estimated symbols for each frame
    de_all_frames(:, b) = de;
    db_all_frames(:, b) = db;

end
% Reshape the matrices into vectors
de_combined = reshape(de_all_frames, [], 1);
db_combined = reshape(db_all_frames, [], 1);
figure(1)
scatterplot(de_combined);
title('Constellation Received by Eavsdropper');
figure(2)
scatterplot(db_combined);
title('Constellation Received by Legitimate User');
toc