function y = ModifiedCompetitionRankings(x)
%            Author: Liber Eleutherios                                             %
%            E-Mail: libereleutherios@gmail.com                             %
%            Date: 8 April 2008                                                       %
% Prepare data
ctrl = isvector(x) & isnumeric(x);
if ctrl
  x = x(:);
  x = x(~isnan(x) & ~isinf(x));
else
  error('x is not a vector of numbers! The Modified Competition Rankings could not be calculated')
end
% Find the Frequency Distribution
[y, ind] = sort(x);
FreqTab(:, 1) = y([find(diff(y)); end]);
N1 = length(x);
N2 = length(FreqTab(:, 1));
if N1 == N2
  y(ind) = 1:N1;
  return
end
FreqTab(:, 2) = histc(y, FreqTab(:, 1));
% Find the rankings
y = zeros(N1, 1);
k = 1;
for i = 1:N2
  y(k:(k + FreqTab(i, 2) - 1)) = k + FreqTab(i, 2) - 1;
  k = k + FreqTab(i, 2);
end
y = sortrows([y, ind], 2);
y(:, 2) = [];
end