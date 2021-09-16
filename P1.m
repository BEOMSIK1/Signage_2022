clear, clc, close all


FFT_Size = 256;
GI_Size = FFT_Size / 4;
D = (FFT_Size == 128 || FFT_Size == 256) * 14   +   (FFT_Size == 512) * 28;
Modulation_Order = 6;
Symbol_Size = (FFT_Size == 128) * 108   +   (FFT_Size == 256) * 234   +  (FFT_Size == 512) * 468;
Nt = 8;
Nr_Relay = 16;
Nu = Nt;
N_Relay = 2;
Multi_Path = 7;

Code_Rate = 2/3;
Z = 27; % Subblock Size : 27, 54, 81
N = 24 * Z; % 648, 1296, 1944
M = round(N * (1 - Code_Rate)); % ceil : 올림할려는 것이 아닌 소수 .000을 제거하기 위함
Data_Size = N - M;
LDPC_Iteration = 10;

SNR = 0 : 3 : 30;
Iteration = 10;

c = qammod(0 : 2^Modulation_Order-1, 2^Modulation_Order, 'UnitAveragePower', true);
L = length(c);

BER = zeros(1,length(SNR));

Data_Rate = Symbol_Size * Modulation_Order * Code_Rate / 3.6 * Nt / (N_Relay+1);


for SNR_Index = 1 : length(SNR)
    
    tic
    
    for Iteration_Index = 1 : Iteration
        
        Data = randi([0 1] , [Nt Data_Size]);
        
        Decoding = Data;
        
        for Relay_Index = 1 : N_Relay+1
            
            if Relay_Index < N_Relay+1 % Source -> Relay : SU-MIMO,  Relay -> Relay : SU-MIMO
                
                Decoding = SM_SCM(Decoding, Multi_Path, Nt, Nr_Relay, N_Relay, Data_Size, FFT_Size, Code_Rate, Modulation_Order, D, GI_Size, L, c, M, N, Z, LDPC_Iteration, SNR(SNR_Index));
                
            else % Relay -> Destination : MU-MIMO
                
                Decoding = ZF_Precoding_SCM(Decoding, Multi_Path, Nu, Nr_Relay, N_Relay, Data_Size, FFT_Size, Code_Rate, GI_Size, Modulation_Order, D, M, N, Z, LDPC_Iteration, SNR(SNR_Index));
                
            end
            
        end
        
        BER(SNR_Index) = BER(SNR_Index) + mean(mean(Data ~= Decoding)) / Iteration;
        
    end
    
    toc
    
end

Throughput = (1 - BER).^Data_Size * Data_Rate;


figure(1)
semilogy(SNR,BER,'-o','LineWidth',2)
title('BER Performance'), xlabel('SNR (dB)'), ylabel('BER'), grid on, set(gca,'FontSize',12)

figure(2)
plot(SNR,Throughput,'-o','LineWidth',2)
title('Throughput Performance'), xlabel('SNR (dB)'), ylabel('Throughput [Mbps]'), grid on, set(gca,'FontSize',12)