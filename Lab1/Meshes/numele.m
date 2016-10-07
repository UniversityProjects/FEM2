% numeriamo gli elementi

nele = size(vertices,1);

for iele=1:nele

    % trasformo iele in stringa
    s = num2str(iele);
    
    % (xb,yb) = baricentro di iele
    xb = sum(xv(vertices(iele,:)))/3;
    yb = sum(yv(vertices(iele,:)))/3;
    
    text(xb,yb,0,s,'Color','b','HorizontalAlignment','center');
end

    
    
    
    
    
    
    
    
    
    
    
    












