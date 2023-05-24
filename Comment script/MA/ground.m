function [position,isterminal,direction] = ground(~,Y, varargin)

position = Y(2); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached decreasing 