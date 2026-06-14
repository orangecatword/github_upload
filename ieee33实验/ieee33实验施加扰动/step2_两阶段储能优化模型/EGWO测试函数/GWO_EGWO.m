f_GWO = zeros(100,1);
f_EGWO = zeros(100,1); 
f_EGWO2 = zeros(100,1); 
f_EGWO3 = zeros(100,1); 
for i = 1:100
    % f_GWO(i) = GWO_test(); 
    % f_EGWO(i) = EGWO_test(); 
    f_EGWO2(i) = EGWO2_test(); 
    % f_EGWO3(i) = EGWO3_test(); 
end