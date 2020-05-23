function [ beta_A ] = AODexp( beta_A , Zbas , zAOD , corrAOD )
beta_A=beta_A.*exp(-((Zbas-zAOD)./1000)./corrAOD);
end