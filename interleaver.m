function Output = interleaver(Encoding,Mode)


N = 10;
L = length(Encoding);

if Mode == 1
    
    for k = 1 : L/(10*N)
        
        for m = 1 : 10*N
            
            Output(: , (k-1)*10*N+m) = Encoding(: , floor((m-1)/N) + 10*rem(m-1,N) + 1 + (k-1)*10*N);
            
        end
        
        for m = floor(L/(10*N)) * 10*N + 1 : L
            
            Output(:,m) = Encoding(:,m);
            
        end
        
    end
    
else
    
    for k = 1 : L / (10*N)
        
        for m = 1 : 10*N
            
            Output(: , floor((m-1)/N) + 10*rem(m-1,N) + 1 + (k-1)*10*N ) = Encoding(: , (k-1)*10*N+m);
            
        end
        
        for m = floor(L/(10*N)) * 10*N + 1 : L
            
            Output(:,m) = Encoding(:,m);
            
        end
        
    end
    
end
