function Modulate = base_mod(Data,Modulation_Order)


if Modulation_Order == 1
    
    Modulate = 2*Data - 1;
    
else
    
    L = size(Data,2);
    
    Remainder = rem(L,Modulation_Order);
    
    if Remainder ~= 0
        Data(: , L+1 : L+Remainder) = 0;
        L = L + Remainder;
    end
    
    if Modulation_Order == 2
        
        Data = 2*Data - 1;
        
        I = Data(: , 2:Modulation_Order:L);
        Q = Data(: , 1:Modulation_Order:L);
        
        Modulate = (I + Q*1j) / sqrt(2);
        
    elseif Modulation_Order == 4
        
        I1_ = Data(: , 1:Modulation_Order:L);
        I1 = 4*I1_ - 2;
        
        I2_ = Data(: , 2:Modulation_Order:L);
        I2 = 2*(I1_ ~= I2_) - 1;
        
        Q1_ = Data(: , 3:Modulation_Order:L);
        Q1 = 4*Q1_ - 2;
        
        Q2_ = Data(: , 4:Modulation_Order:L);
        Q2 = 2*(Q1_ ~= Q2_) - 1;
        
        Modulate = ((I1+I2) + (Q1+Q2)*1j) / sqrt(10);
        
    elseif Modulation_Order == 6
        
        I1_ = Data(:,1:Modulation_Order:L);
        I1 = 8*I1_ - 4;
        
        I2_ = Data(:,2:Modulation_Order:L);
        I2 = 4*(I1_ ~= I2_) - 2;
        
        I3_ = Data(:,3:Modulation_Order:L);
        I3 = 2*xor(xor(I1_,I2_),I3_) - 1;
        
        Q1_ = Data(:,4:Modulation_Order:L);
        Q1 = 8*Q1_ - 4;
        
        Q2_ = Data(:,5:Modulation_Order:L);
        Q2 = 4*(Q1_ ~= Q2_) - 2;
        
        Q3_ = Data(:,6:Modulation_Order:L);
        Q3 = 2*xor(xor(Q1_,Q2_),Q3_) - 1;
        
        Modulate = ((I1+I2+I3) + ((Q1+Q2+Q3)*1j)) / sqrt(42);
        
    elseif Modulation_Order == 8
        
        I1_ = Data(:,1:Modulation_Order:L);
        I1 = 16*I1_ - 8;
        
        I2_ = Data(:,2:Modulation_Order:L-6);
        I2 = 8*(I1_ ~= I2_) - 4;
        
        I3_ = Data(:,3:Modulation_Order:L-5);
        I3 = 4*xor(xor(I1_,I2_),I3_) - 2;
        
        I4_ = Data(:,4:Modulation_Order:L-4);
        I4 = 2*xor(xor(xor(I1_,I2_),I3_),I4_) - 1;
        
        Q1_ = Data(:,5:Modulation_Order:L-3);
        Q1 = 16*Q1_ - 8;
        
        Q2_ = Data(:,6:Modulation_Order:L-2);
        Q2 = 8*(Q1_ ~= Q2_) - 4;
        
        Q3_ = Data(:,7:Modulation_Order:L-1);
        Q3 = 4*xor(xor(Q1_,Q2_),Q3_) - 2;
        
        Q4_ = Data(:,8:Modulation_Order:L);
        Q4 = 2*xor(xor(xor(Q1_,Q2_),Q3_),Q4_) - 1;
        
        Modulate = ((I1+I2+I3+I4) + ((Q1+Q2+Q3+Q4)*1j)) / sqrt(170);
        
    end
    
end