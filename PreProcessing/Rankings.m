function y = Rankings(x, meth)
%            Author: Liber Eleutherios                                             %
%            E-Mail: libereleutherios@gmail.com                             %
%            Date: 8 April 2008                                                       %
if nargin == 1;
  y = FractionalRankings(x);
  return;
end;
switch meth
    case 'fractional'
        y = FractionalRankings(x);
    case 'dense'
        y = DenseRankings(x);
    case 'competition1'
        y = StandardCompetitionRankings(x);
    case 'competition2'
        y = ModifiedCompetitionRankings(x);
    case 'ordinal'
        y = OrdinalRankings(x);
end