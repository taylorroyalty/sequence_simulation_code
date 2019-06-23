%Applying equation 6 to get analytical solution for coupon collector problem (Figure 2)
g=round(linspace(1,100,10));%number of unique clades
k=1e6;%average genome size
l=100;%read length

G_mg=g*(k-l+1);%calculate size of metagenome
outfile=nan(numel(G_mg)*4,4);
count=1;

%Coupon Collector General Formulation; Flajolet et al. 1992; Eq. 13
fun=@(p,t) 1-prod(1-exp(-p.*t));

%Lognormal Abundance Distribution; Magurran 1988; Eq 2.11
log_fun=@(S,a,R)S*exp(-1*a^2*R.^2);
a=0.02; S=1; R=linspace(0,99,10);%define parameters for lognormal distribution
log_prob=log_fun(S,a,R)/sum(log_fun(S,a,R));%calculate proportion of distribution

for i=G_mg%iterate through     
%%Define Community Structure Vectors
p1=ones(1,i)/i;
indx1=round(i*0.5);
p2=[ones(1,indx1)*0.9/indx1 ones(1,i-indx1)*0.1/(i-indx1)];
indx2=round(i*0.1);
p3=[ones(1,indx2)*0.9/indx2 ones(1,i-indx2)*0.1/(i-indx2)];
indx3=round(diff(linspace(1,i,11)));
p4=[];
for j=1:1:10
    if j<10
        p4=[p4 ones(1,indx3(1))*log_prob(j)/indx3(1)];
    else
        p4=[p4 ones(1,i-length(p4))*log_prob(j)/(i-length(p4))];
    end
end

%%Calculate Shannon Diversity
h_index1=sum(-p1.*log(p1));%high diversity
h_index2=sum(-p2.*log(p2));%moderate diversity
h_index3=sum(-p3.*log(p3));%low diversity
h_index4=sum(-p4.*log(p4));%lognormal 

expected1=integral(@(t)fun(p1,t),0,inf,'ArrayValued',true);
expected2=integral(@(t)fun(p2,t),0,inf,'ArrayValued',true);
expected3=integral(@(t)fun(p3,t),0,inf,'ArrayValued',true);
expected4=integral(@(t)fun(p4,t),0,inf,'ArrayValued',true);

%%Write to output matrix
outfile(count,:) = [i,h_index1,double(expected1),1];
outfile(count+1,:) = [i,h_index2,double(expected2),2];
outfile(count+2,:) = [i,h_index3,double(expected3),3];
outfile(count+3,:) = [i,h_index4,double(expected4),4];
count=count+4;
end

filename = fullfile('.', 'modeled_sequencing_analytical.txt');
fid=fopen(filename, 'wt');
fprintf(fid, '%s,%s,%s,%s\n', 's','shannon_index','sequences','diversity_category');
fclose(fid);
dlmwrite('modeled_sequencing_analytical.txt',outfile,'precision',15,'-append')
