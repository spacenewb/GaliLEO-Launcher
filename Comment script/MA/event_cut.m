function [position,isterminal,direction] = event_cut(~,Y,~, ~,~,para)

position = Y(3) - para.z_cut ; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zer