function res = median_adjusted(b)
bc    = b(:);
bc(bc == 0) = NaN;
m_b   = nanmedian(bc);

res = b + m_b;