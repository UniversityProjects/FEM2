function Z=uexactGmine(x,y,k1, k2)

% G = [...
%         -(pi/2)*sin((pi/2)*x)*cos((pi/2)*y), ... % du1/dx
%         -(pi/2)*sin((pi/2)*x)*cos((pi/2)*y); ... % du2/dx
%         -(pi/2)*cos((pi/2)*x)*sin((pi/2)*y), ... % du1/dy
%         -(pi/2)*cos((pi/2)*x)*sin((pi/2)*y)  ... % du2/dy
%     ];

G = [...
        x*y*(x^2 - 1)*(y^2 - 1), ...                          % du1/dx
        -((x^2 - 1)*(y^2 - 1)^2)/4 - (x^2*(y^2 - 1)^2)/2; ... % du2/dx
        ((x^2 - 1)^2*(y^2 - 1))/4 + (y^2*(x^2 - 1)^2)/2, ...  % du1/dy
        -x*y*(x^2 - 1)*(y^2 - 1) ...                          % du2/dy
    ];

Z = G(k1, k2);

end

