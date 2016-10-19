function [u]= initial(x,y);
% date le coordinate globali X=(x,y) calcola valore iniziale
% ... in quel punto

u = sin(pi*x)*sin(pi*y);

