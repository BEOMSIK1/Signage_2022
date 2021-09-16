function Decoding = ZF_Precoding_SCM_ZF(Decoding, Multi_Path, Nu, Nt, N_Relay, Data_Size, FFT_Size, Code_Rate, GI_Size, Modulation_Order, D, M, N, Z, LDPC_Iteration, SNR)

%% SCM parameters
model = SCM();
model.n_mray = 20;  % scatter
model.n_path = Multi_Path+3;
model.asd = 25;
model.zsd = 10;
model.asa = 25;
model.zsa = 10;

rx_ant = 1;
model.ant(rx_ant,Nt);

%% channel generate
for d = 1:Nu
    [temp,rx_angle] = model.FD_channel(FFT_Size + GI_Size);
    h(:,:,1+(d-1)*rx_ant:d*rx_ant,:) = temp;
%     [~, idx] = max( abs( temp(:,1,1,1) ) );
%     sel_angle(:,d) = rx_angle(:,idx);
end
h_(:,:,:) = h(:,1,:,:);
H_ = fft(h_,FFT_Size,1);


%% simulate

H = Matrix_Generate(Z,Code_Rate);

Encoding = LDPC_Encoding(Decoding, H, Data_Size);

Interleaving = interleaver(Encoding,1);

Modulate = qammod(Interleaving.', 2^Modulation_Order, 'InputType', 'bit', 'UnitAveragePower', true).';

T = size(Modulate,2);
U = floor(T/FFT_Size);
V = rem(T,FFT_Size);



H_Reshape = reshape(reshape(H_,FFT_Size, Nt*Nu).',Nu,Nt*FFT_Size);


n = 1;

for t = 1 : U+1
    
    F = (t == U+1) * V  +  (t ~= U+1) * FFT_Size;
    
    
    for k = 1 : F
        
        H_k = H_Reshape(: , Nt*(k-1)+1 : Nt*k);
        
        G = H_k' * inv(H_k*H_k');

        NF(n) = trace(G'*G);
        
        Precoding(:,k) = G * (Modulate(:,n) / sqrt(NF(n)));
        
        n = n + 1;
        
    end
    
    IFFT = ifft(Precoding,FFT_Size,2) * (FFT_Size/sqrt(FFT_Size-D));
    
    CP_Add(: , (FFT_Size+GI_Size)*(t-1)+1 : (FFT_Size+GI_Size)*t) = [IFFT(: , FFT_Size-GI_Size+1 : FFT_Size) IFFT];
    
end

for u = 1 : Nu
    
    for Tx = 1 : Nt
        
        Receive(Tx,:) = conv(CP_Add(Tx,:) , h_(:,u,Tx).');
        
    end
    
    hx(u,:) = sum(Receive,1);
    
end

[y , N0] = awgn_noise(hx,SNR);

for t = 1 : U+1
    
    Index = FFT_Size*(t-1) + GI_Size*t + 1  :  (FFT_Size + GI_Size) * t;
    
    CP_Remove = y(:,Index);
    
    X(: , FFT_Size*(t-1)+1 : FFT_Size*t) = fft(CP_Remove,FFT_Size,2) / (FFT_Size/sqrt(FFT_Size-D));
    
end

X = X(: , 1:T) .* sqrt(NF);

Demodulate = qamdemod(X.', 2^Modulation_Order, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', N0).'; % approxllr

Deinterleaving = interleaver(Demodulate,2);

for Tx = 1 : Nu
    
    Decoding(Tx,:) = LDPC_Decoding(Deinterleaving(Tx,:), H, M, N, LDPC_Iteration);
    
end