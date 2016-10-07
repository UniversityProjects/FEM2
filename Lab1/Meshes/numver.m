% numeriamo i nodi

nver=size(xv,1);

for iv=1:nver
    s = num2str(iv);
    % if vertexmarker(iv)==1
    %    text(xv(iv),yv(iv),0,s,...
    %        'Color','r','BackGroundColor','w','FontWeight','bold');
    %else
        text(xv(iv),yv(iv),0,s,...
            'Color','r','BackGroundColor','w');
    %end
end

