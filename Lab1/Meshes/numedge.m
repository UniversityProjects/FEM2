% numeriamo gli edge

nedge = size(endpoints,1);

for iedge=1:nedge
    %
    v1 = endpoints(iedge,1);
    v2 = endpoints(iedge,2);
    %
    x1 = xv(v1);
    y1 = yv(v1);
    %
    x2 = xv(v2);
    y2 = yv(v2);
    %
    siedge = num2str(iedge);
    %
    s = 0.35;
    %
    % if edgemarker(iedge)==1
    %    text(s*x1+(1-s)*x2,s*y1+(1-s)*y2,0,siedge,...
    %        'Color','k','BackGroundColor','w','FontWeight','bold','HorizontalAlignment','center');
    %else
        text(s*x1+(1-s)*x2,s*y1+(1-s)*y2,0,siedge,...
            'Color','k','BackGroundColor','w','HorizontalAlignment','center');
    %end
end



    
    
    
    
    
    
    
    
    
    
    
    












