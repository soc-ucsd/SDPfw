

function CDD_new = cdd_shuffle(L1, L2)
    %this one is much trickier. I'm ditching 
    i1 = [1,2,4];
    i2 = [2,3,5];
    
    %start with  one clique at a time    
    Bi = [1 0 1  1;
          0 1 1  1;
          0 0 1 -1];
    
    CDD1 =  zeros(6, 4);
    CDD2 =  zeros(6, 4);
    %CDD1 = zeros(size(DD_ext)); 
    %CDD2 = zeros(size(DD_ext));
    for i = 1:size(Bi, 2)
       bi = Bi(:, i) ; 
       Mi =  [bi(1) bi(3); bi(3) bi(2)];       
       M1_out =  L1 * Mi * L1';
       M2_out =  L2 * Mi * L2';
              
       CDD1(i1, i) =  M1_out([1,4,2]); 
       CDD2(i2, i) =  M2_out([1,4,2]);
    end
    
    %now for the fiber product
    Xp = [1 1; 1 1];
    Xn = [1 -1 ; -1 1];
    
    Xp1 = L1 * Xp * L1'; 
    Xp2 = L2 * Xp * L2'; 
    
    Xn1 = L1 * Xn * L1'; 
    Xn2 = L2 * Xn * L2'; 
    
    CDD_fiber  = zeros(6, 4);
    
    for i = 0:3
        p1 = bitand(i, 1);
        p2 = bitand(i, 2);
        
        if p1
            X1 = Xp1;
        else
            X1 = Xn1;
        end
        
         if p2
             X2 =  Xp2;
         else
             X2 =  Xn2;
         end
         
         scale = X1(2, 2)/X2(1, 1);
         
         X2_scale = scale*X2;
         
         v_out =  [X1(1), X1(4), X2_scale(4), X1(2), X2_scale(2), 0]';
         CDD_fiber(:, i+1) = v_out;                  
    end
    
    CDD_new = [CDD1 CDD2 CDD_fiber];
    
end

