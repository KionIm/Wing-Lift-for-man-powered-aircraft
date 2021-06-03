function[chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition)


chord= zeros(1,npartition);

for i = 1:npartition
    if i <= q1 
        chord(1,i) = (chord0-chord1)/(0-q1)*(i-0) + chord0;
    
    elseif i <= q2 
        chord(1,i) = (chord1-chord2)/(q1-q2)*(i-q1) + chord1;
        
    elseif i <= q3 
        chord(1,i) = (chord2-chord3)/(q2-q3)*(i-q2) + chord2;    

    elseif i <= q4 
        chord(1,i) = (chord3-chord4)/(q3-q4)*(i-q3) + chord3;
        
    elseif i <= q5 
        chord(1,i) = (chord4-chord5)/(q4-q5)*(i-q4) + chord4;    

    elseif i <= q6
        chord(1,i) = (chord5-chord6)/(q5-q6)*(i-q5) + chord5;
    
 
        
    end
    
end      