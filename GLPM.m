function [GLPM] = GLPM(n,s)
%好格子点法生成超立方体上均匀设计
h = 1:n;
ind = find(gcd(h,n)==1);
hm = h(ind);
m = length(hm);
udt = mod(h'*hm,n);
ind0 = find(udt==0);
udt(ind0)=n;
if s > m
    disp("s必须小于或等于m");
    return;
else
    mcs = nchoosek(m,s);
    if mcs < 1e5
        tind = nchoosek(1:m,s);
        [p,q] = size(tind);
        cd2 = zeros(p,1);
        for k=1:p
            UT = udt(1:n,tind(k,:));
            cd2(k,1) = CD2(UT);
        end
        tc = tind(find(abs(cd2-min(cd2))<1e-5),:);
        for r=1:size(tc,1)
            GLPM(:,:,r) = udt(:,tc(r,:));
        end
    else
        for k=1:n
            a = k;
            UT = mod(h'*a.^(0:s-1),n);
            cd2(k,1)=CD2(UT);
        end
        tc=find(abs(cd2-min(cd2))<1e-5);
        for r=1:size(tc,1)
            GLPM(:,:,r)=mod(h'*tc(r).^(0:s-1),n);
        end
        ind0=find(GLPM==0);
        GLPM(ind0)=n;
    end
 end

