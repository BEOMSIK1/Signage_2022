function Decoding = LDPC_Decoding(Lc, H, M, N, LDPC_Iteration)


%% 1. Initialization
Lq = ones(M,1) * Lc;
Lr = zeros(M,N);

for Iteration_Index = 1 : LDPC_Iteration
    
    %% 2. First Half Round Iteration (Horizontal Step)
    A = sign(Lq);
    B = abs(Lq);
    
    for i = 1 : M
        
        R = find(H(i,:) == 1);
        
        for k = 1 : length(R)
            
            Index = R([1 : k-1  k+1 : end]);
            
            Min = min( B(i,Index) );
            
            Lr(i,R(k)) = prod( A(i,Index) ) * Min;
            
        end
        
    end
    
    
    %% 3. Second Half Round Iteration (Vertical Step)
    for j = 1 : N
        
        Index = find(H(:,j) == 1);
        
        Lq(Index,j) = Lc(j) + sum(Lr(Index,j)) - Lr(Index,j);
        
    end
    
    
    %% 4. Soft Decision
    LQ = Lc + sum(Lr,1);
    
    
    %% 5. Hard Decision
    Decoding = LQ < 0;
    
    
    %% 6. Parity Check
    if rem(H*Decoding.' , 2) == 0
        break
    end
    
end

Decoding = Decoding(:,1:N-M);