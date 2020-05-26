% MIMO sphere decoder
% works for 2-PAM only (it can be very easily modified)

clear all
close all
clc
Nt = 6; % number of transmit antennas
Nr = 8; % number of receive antennas
fade_var = 1; % fade variance of the MIMO channel
SNR_dB = 28; % SNR per bit (dB)
epsilon = 0.000001; % for choosing radius
num_frames = 10^3; % simulation runs

% SNR parameters
SNR = 10^(0.1*SNR_dB);
noise_var = 1*fade_var*Nt*Nr/(2*SNR*1*Nt); %awgn variance

% radius (square) of the sphere
d1_2 = chi2inv(1-epsilon,Nr)*noise_var;

C_Ser=0;
tic()
for i1=1:num_frames
% source
a = randi([0 1],1,Nt);

% BPSK (2PAM)
seq = 1-2*a;

% (Nr x Nt) MIMO flat channel
fade_chan = normrnd(0,sqrt(fade_var),Nr,Nt);

% AWGN
noise = normrnd(0,sqrt(noise_var),Nr,1);

% channel output
chan_op = fade_chan*seq.' + noise;

%--------------- sphere decoder-----------------------------------
[Q,R1] = qr(fade_chan); % QR decomposition
dec_seq = zeros(1,Nt); 
R = R1(1:Nt,:);% R matrix
Q1 = Q(:,1:Nt);
Q2 = Q(:,Nt+1:end);
y = Q1'*chan_op;
dec_seq1 = zeros(1,Nt); % decoded sequence
up_bound = zeros(1,Nt);
low_bound = zeros(1,Nt);
d2 = zeros(1,Nt);
y1 = zeros(1,Nt);
caseno = 1;
while (caseno~=0)
switch (caseno)
    case 1
          k = Nt;
          d2(k) = d1_2 - norm(Q2'*chan_op)^2;
          y1(k) = y(Nt);
          caseno=2;
    case 2
        temp = [(-sqrt(d2(k))+y1(k))/R(k,k) (sqrt(d2(k))+y1(k))/R(k,k)];
         up_bound(k) = 1-2*(max(temp)<1) ;
         low_bound(k) = 1-2*(min(temp)<-1);
         dec_seq(k) = low_bound(k)-1; % decoded sequence
         caseno=3;  
    case 3
        dec_seq(k) = dec_seq(k)+1;
              
         if dec_seq(k)<=up_bound(k)
             caseno = 5;
         else
             caseno = 4;
         end

    case 4
         k = k+1;
         if k==Nt+1
             break;
         else
             caseno = 3;
         end
    case 5
        if k==1
            caseno = 6;
        else
            k=k-1;
            d2(k) = d2(k+1)-(y1(k+1) - R(k+1,k+1)*dec_seq(k+1))^2;
            y1(k) = y(k)-R(k,k+1:Nt)*dec_seq(k+1:end).';
            caseno = 2;
        end
    case 6
        dec_seq1=dec_seq; % dec_seq1 is decoded sequence
        caseno =3;
        
end % for switch
end % for while
C_Ser = C_Ser +nnz(dec_seq1-seq);
end
toc()
% symbol error rate
SER = C_Ser/(num_frames*Nt)