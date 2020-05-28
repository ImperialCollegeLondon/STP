function [ G ] = NHPO( G , WAR , obsOCC )
%NHPO Non-homogeneous probability of precipitation cccurrence
%   Inputs:
%   G - Gaussian rain field [-]
%   WAR - wet area ration [-]
%   obsOCC - normalized rainfall ocurrence field [-]
%   Output:
%   G - Gaussian field for rain interpolation. 0 indicate no rain.
%% Find WAR pixels
A=(WAR.*obsOCC);
if max(A)>1
   s=sum(A(A>1))-sum(A>1);
   A(A>1)=1;
   while s>0
       if isempty(max(A(A<1)-1))
           break
       end
       ind = max(A(A<1)-1)==(A-1);
       if 1-A(ind)<s
           s=s-(1-A(ind));
           s = unique(s);
           A(ind)=1;
       else
           A(ind)=A(ind)+s;
           break
       end
   end
end
C=1-A;
C=reshape(C,size(G,1),size(G,2));
G(G<C)=0;
%% Normalize the precentiles
qMin=min(G(G>0));
if ~isempty(qMin)
    G(G>0)=(G(G>0)-qMin)./(1-qMin);
end
G(G==1)=0.9999;
end