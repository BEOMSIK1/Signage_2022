function Decoding = SM_SCM(Decoding, Multi_Path, Nt, Nr, N_Relay, Data_Size, FFT_Size, Code_Rate, Modulation_Order, D, GI_Size, L, c, M, N, Z, LDPC_Iteration, SNR)
%% SCM Parameters
model = SCM();
model.n_mray = 20 ;
model.n_path = Multi_Path+3;
model.los = 0;
model.asd = 25;
model.zsd = 10;
model.asa = 25;
model.zsa = 10;
model.ant(Nr,Nt);
d = 1;

%% Channel generate
[temp,rx_angle] = model.FD_channel(FFT_Size + GI_Size);
h(:,:,1+(d-1)*Nr:d*Nr,:) = temp;
%     [~, idx] = max( abs( temp(:,1,1,1) ) );

h_(:,:,:) = h(:,1,:,:);
H_ = fft(h_,FFT_Size,1);


% h = Rayleigh_channel([Multi_Path Nt*Nr], 3/(N_Relay+1), 1); % [h11 ; h12; h21 ; h22]
%
% H_ = fft(h,FFT_Size,2);

H = Matrix_Generate(Z,Code_Rate);

Encoding = LDPC_Encoding(Decoding, H, Data_Size);

Interleaving = interleaver(Encoding,1);

Modulate = qammod(Interleaving.', 2^Modulation_Order, 'InputType', 'bit', 'UnitAveragePower', true).';

T = size(Modulate,2);
U = floor(T/FFT_Size);
V = rem(T,FFT_Size);

for t = 1 : U+1
    
    F = (t == U+1) * V  +  (t ~= U+1) * FFT_Size;
    
    IFFT = ifft(Modulate(: , FFT_Size*(t-1)+1 : FFT_Size*(t-1)+F), FFT_Size, 2) * (FFT_Size/sqrt(FFT_Size-D));
    
    CP_Add(: , (FFT_Size+GI_Size)*(t-1)+1 : (FFT_Size+GI_Size)*t) = [IFFT(: , FFT_Size-GI_Size+1 : FFT_Size) IFFT];
    
end

for Rx = 1 : Nr
    
    for Tx = 1 : Nt
        
        Receive(Tx,:) = conv(CP_Add(Tx,:) / sqrt(Nt) , h_(:,Rx,Tx).');
        %Receive(Tx,:) = conv(CP_Add(Tx,:) , h_(:,Rx,Tx).');
        
    end
    
    hx(Rx,:) = sum(Receive,1);
    
end

[y , N0] = awgn_noise(hx,SNR);

m = 1;

for t = 1 : U+1
    
    Index = FFT_Size*(t-1) + GI_Size*t + 1  :  FFT_Size*t + GI_Size*t;
    
    CP_Remove = y(:,Index);
    
    Y = [fft(CP_Remove,FFT_Size,2) / (FFT_Size/sqrt(FFT_Size-D)); zeros(Nt,FFT_Size)] * sqrt(Nt);
    
    %Y = [fft(CP_Remove,FFT_Size,2) / (FFT_Size/sqrt(FFT_Size-D))] * sqrt(Nt);
    
    F = (t == U+1) * V  +  (t ~= U+1) * FFT_Size;
    
    H_Reshape = reshape(reshape(H_,FFT_Size, Nt*Nr).',Nr,Nt*FFT_Size);
    
    for k = 1 : F
        
        H_k = H_Reshape(: , Nt*(k-1)+1 : Nt*k);
        
        H_k = [H_k; sqrt(N0)*eye(Nt)];
        
        G = inv(H_k'*H_k) * H_k';
        
        G_Square = sum(abs(G).^2 , 2);
        
        [~ , Order] = sort(G_Square , 'descend');
        
        H_Sorting = H_k(:,Order);
        
        [Q , R] = qr(H_Sorting);
        
        Z_k = Q' * Y(:,k);
        
        % Threshold
        Temp_Candidate = zeros(Nt,1);
        
        Past = 0; Error_Candidate = 0;
        
        for Tx = Nt : -1 : 1
            
            Present = R(Tx,Tx) * c;
            
            if Tx ~= Nt
                Past = R(Tx,Tx+1:Nt) * Temp_Candidate(Tx+1:Nt,1);
            end
            
            Cumulative_Error = real(Z_k(Tx) - (Past + Present)).^2 + imag(Z_k(Tx) - (Past + Present)).^2 + Error_Candidate;
            
            [Error_Candidate , Min_Index] = min(Cumulative_Error);
            Temp_Candidate(Tx) = c(Min_Index);
            
        end
        
        Threshold = min(Error_Candidate);
        
        % Path Eliminations
        X_Candidate = [];
        
        Past = 0; Error_Candidate = 0;
        
        for Tx = Nt : -1 : 1
            
            Present = R(Tx,Tx) * c;
            
            if Tx ~= Nt
                
                Past = R(Tx,Tx+1:Nt) * X_Candidate(Tx+1:Nt,:);
                
                Error_Candidate = repelem(Error_Candidate,L);
                
            end
            
            Present = repmat(Present, 1, length(Past));
            
            Past = repelem(Past,L);
            
            Cumulative_Error = real(Z_k(Tx) - (Past + Present)).^2 + imag(Z_k(Tx) - (Past + Present)).^2 + Error_Candidate;
            
            [Error_Sorting , Error_Order] = sort(Cumulative_Error);
            
            Survival_Length = sum(Error_Sorting(1:L) <= Threshold);
            
            if Survival_Length == 0
                Survival_Length = L;
            end
            
            Error_Candidate = Error_Sorting(1:Survival_Length);
            
            Error_Index = Error_Order(1:Survival_Length);
            
            if Tx == Nt
                
                X_Candidate(Nt,:) = c(Error_Order(1:Survival_Length));
                
            else
                
                S = size(X_Candidate,2);
                
                Candidate_Tree = repelem(X_Candidate,1,L);
                Candidate_Tree(Tx,:) = repmat(c,1,S);
                
                X_Candidate = Candidate_Tree(:,Error_Index);
                
            end
            
        end
        
        X(Order,m) = X_Candidate(:,1);
        
        m = m + 1;
        
    end
    
end

Demodulate = qamdemod(X.', 2^Modulation_Order, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', N0).'; % approxllr

Deinterleaving = interleaver(Demodulate,2);

for Tx = 1 : Nt
    
    Decoding(Tx,:) = LDPC_Decoding(Deinterleaving(Tx,:), H, M, N, LDPC_Iteration);
    
end