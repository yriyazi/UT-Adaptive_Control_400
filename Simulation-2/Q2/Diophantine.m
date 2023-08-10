function [R , S ] = Diophantine(A , B , Ac)    
   %%Logical Definitions
    ord_A = numel(A)-1 ;    ord_B = numel(B)-1 ;    ord_Ac = numel(Ac)-1 ;

    B  = [zeros(1,ord_A-ord_B)        , B]      ;
    Ac = [zeros(1,2*ord_A-(ord_Ac+1)) , Ac]     ;
    E  =  zeros(2*ord_A               , 2*ord_A);
    %%       Silvester Matrix 
    for i = 1:ord_A
       E(1:2*ord_A,i  ) = [zeros(i-1,1); A';zeros(ord_A-i,1)];   
       E(1:2*ord_A,i+ord_A) = [zeros(i-1,1); B';zeros(ord_A-i,1)];   
    end
    %%       Result
    RS = E\ Ac'; 
    R = RS(1:ord_A)'; 
    S = RS(ord_A+1:end)';
end