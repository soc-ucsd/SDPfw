
function v = m2v(m)
    v = zeros(6, 1);
    v(1) = m(1, 1); 
    v(2) = m(2, 2);
    v(3) = m(3 ,3);
    v(4) = m(1, 2); 
    v(5) = m(2, 3);    
    v(6) = m(1, 3);
end