function Encoding = LDPC_Encoding(Data, H, Data_Size)


A = H(: , 1:Data_Size);

T = H(: , Data_Size+1:end);

p2 = mod( mod(-inv(T)*A,2) * Data' , 2)';

Encoding = [Data p2];