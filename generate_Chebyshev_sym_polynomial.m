function [ eq ] = generate_Chebyshev_sym_polynomial( n )
%GENERATE_CHEBYSHEV_SYM_POLYNOMIAL returns the symbolic equation of n'th
%Chebyshev polynomial
% dayi
% 11/23/2017

if (floor(n) ~= n) || (abs(n) ~= n)
    error('n should be a non-negative integer')
end

syms t;

T=cell(1,n+1);
T{1}=1;
T{2}=t;

for ii=3:(n+1)
    T{ii} = 2*t*T{ii-1}-T{ii-2};
end

eq = simplify(T{n+1});

end

