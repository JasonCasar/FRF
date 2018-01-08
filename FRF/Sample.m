function [n, obs_index, Z_obs, St, M] = Sample(Z_t,S,ML,MW)

N = length(Z_t);
obs_index=ones([N 1]);

last=ML+MW-1;
if last<=N
    obs_index(ML:last,1)=0;
end

Z_obs=Z_t(obs_index>0);
St = S(obs_index(:,1)>0,:);
n = sum(obs_index(:,1));

M=zeros([N n]);
if ML>0
    for ijk=1:N
        if ijk<ML
              M(ijk,ijk)=1;
        elseif ijk>last
              M(ijk,ijk-MW)=1;
        end
    end
else
    for ijk=(last+1):N
        M(ijk,ijk-MW)=1;
    end 
end

end

