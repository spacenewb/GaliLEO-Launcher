function [value, isterminal, direction] = event_opening(~, Y, varargin)


value = 85549 - Y(2);

isterminal = 1;
direction = 1;