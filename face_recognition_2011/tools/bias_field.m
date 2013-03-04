function [ Mb ] = bias_field(M)
%BIAS_FIELD Summary of this function goes here
%   Detailed explanation goes here

[m_h, m_w, m_d] = size(M);

a = flipud(hankel(linspace(1, 0.5, m_w)));
b = fliplr(hankel(linspace(0.001, 0.5, m_w)));
a(a == 0) = b(logical(triu(b, 1)));
a = flipud(a);
a = repmat(a(1:m_h, 1:m_w), [1 1 m_d]);

Mb = M .* a;

