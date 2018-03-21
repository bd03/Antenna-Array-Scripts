function [ out ] = generate_Chebyshev( range, n )
%GENERATE_CHEBYSHEV returns length(range)xn matrix of 
%n Chebyshev polynomials over defined range
% dayi
% 11/21/2017

if floor(n) ~= n
    error('n should be an integer')
end

pol(1,:)=ones(1,length(range));
pol(2,:)=range;

for ii=3:n+1
    pol(ii,:)=2*range.*pol(ii-1,:)-pol(ii-2,:);
end

out=pol(1:(n+1),:);