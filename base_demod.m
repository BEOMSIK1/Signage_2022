% QPSK성상도
%  10 | 11
%-----------
%  00 | 01
%
% 16QAM 성상도
% 0010 0110 | 1110 1010
% 0011 0111 | 1111 1011
%-----------------------
% 0001 0101 | 1101 1001
% 0000 0100 | 1100 1000

function Demodulate = base_demod(Modulate,Modulation_Order)


if Modulation_Order == 1
    
    Demodulate = Modulate > 0;
    
else
    
    Real = real(Modulate);
    Imag = imag(Modulate);
    
    [p , q] = size(Modulate);
    
    if Modulation_Order == 2
        
        I = Real > 0;
        Q = Imag > 0;
        
        Temp = [Q; I];
        
    elseif Modulation_Order == 4
        
        I1 = Real > 0;
        I2 = abs(Real) < 2/sqrt(10);
        
        Q1 = Imag > 0;
        Q2 = abs(Imag) < 2/sqrt(10);
        
        Temp = [I1; I2; Q1; Q2];
        
    elseif Modulation_Order == 6
        
        I1 = Real > 0;
        I2 = abs(Real) < 4/sqrt(42);
        I3 = (abs(Real) > 2/sqrt(42)) & (abs(Real) < 6/sqrt(42));
        
        Q1 = Imag > 0;
        Q2 = abs(Imag) < 4/sqrt(42);
        Q3 = (abs(Imag) > 2/sqrt(42)) & (abs(Imag) < 6/sqrt(42));
        
        Temp = [I1; I2; I3; Q1; Q2; Q3];
        
    elseif Modulation_Order == 8
        
        I1 = Real > 0;
        I2 = abs(Real) < 8/sqrt(170);
        I3 = (abs(Real) > 4/sqrt(170)) & (abs(Real) < 12/sqrt(170));
        I4 = ((abs(Real) > 2/sqrt(170)) & (abs(Real) < 6/sqrt(170)))  |  ((abs(Real) > 10/sqrt(170)) & (abs(Real) < 14/sqrt(170)));
        
        Q1 = Imag > 0;
        Q2 = abs(Imag) < 8/sqrt(170);
        Q3 = (abs(Imag) > 4/sqrt(170)) & (abs(Imag) < 12/sqrt(170));
        Q4 = ((abs(Imag) > 2/sqrt(170)) & (abs(Imag) < 6/sqrt(170)))  |  ((abs(Imag) > 10/sqrt(170)) & (abs(Imag) < 14/sqrt(170)));
        
        Temp = [I1; I2; I3; I4; Q1; Q2; Q3; Q4];
        
    end
    
    Demodulate = reshape(Temp,p,Modulation_Order*q);
    
end