function v = vfhcl_desk_col(rhf,rhcl)

%% definitions from desk. pot
%        hf         clf           hcl
rmat = [rhf(:)  rhf(:)+rhcl(:) rhcl(:)]; % in bond coordinates and assuming collinear geometry


%% two body terms

data_2body = importdata('DESK_2BODY.txt');

for j=1:3
    b2(j,:) = data_2body(1+(j-1)*12,:);
    c2(:,j) = data_2body((3:12)+12*(j-1),2);
end

vhf = v2(b2(1,:),c2(:,1),rmat(:,1));
vclf = v2(b2(2,:),c2(:,2),rmat(:,2));
vhcl = v2(b2(3,:),c2(:,3),rmat(:,3));


%% three body term

data_3body = importdata('DESK_3BODY.txt');

b3 = data_3body(1,1:3);
ijk = data_3body(2:end,1:3);
dijk = data_3body(2:end,4);

R31 = bsxfun(@power,rmat(:,1).*exp(-b3(1)*rmat(:,1)),ijk(:,1)');
R32 = bsxfun(@power,rmat(:,2).*exp(-b3(2)*rmat(:,2)),ijk(:,2)');
R33 = bsxfun(@power,rmat(:,3).*exp(-b3(3)*rmat(:,3)),ijk(:,3)');

v3body = sum(bsxfun(@times,dijk',R31.*R32.*R33),2);

%% add two and three body terms

v = reshape(vhf+vclf+vhcl+v3body,size(rhf));

%% zero the potential to the minimum of the hcl potential

if size(rhf,2)>1    
    v = v- min(v(:,end));
end

end

% two body potential
function v = v2(b,c,r)
R21 = exp(-b(1)*r)./r;
R22 = exp(-b(2)*r).*r;
v = c(1)*R21+sum(bsxfun(@times,c',bsxfun(@power,R22,0:9)),2);
end