target_coverage=[0.5 0.7 0.9 0.95 0.99 1];%target coverages to log for
even_flag=4;%arbitrary flag reflecting weight matrix used
l=100;%Read length; only important for determining coverate
num_seq=1e7;%number of sequences sampled at a time.
rep_num=3;%replicates
s=[30 40 50 70 90 100];%number of operationally defined genomes for sample environment; each index a different environment
Mbp=1e6;%size of genomes; could be made into a vector in future
k=Mbp-l+1;%unique sample points for a given genome
outfile=nan(numel(s)*rep_num*numel(target_coverage),6); %predefined logging matrix
count=1;%index for writing to log matrix

%For log normal weight distribution
log_fun=@(S,a,R)S*exp(-1*a^2*R.^2);
a=0.02; S=1; R=linspace(1,100,10);

%normalize values
log_prob=log_fun(S,a,R)/sum(log_fun(S,a,R));

for rep=1:1:rep_num%iterate through triplicates
    for i=s%iterate through
        
%uncomment one weight vector and leave other three commented
%         w=[ones(1,k*i)*(1/k*i)];%weight vector; note probabilities are vectorized for the randsample function below
%         w=[ones(1,k*i/2)*0.9./(k*i/2) ones(1,k*i/2)*0.1./(k*i/2)]; %50:50%
%         w=[ones(1,k*i/10)*0.9./(k*i/10) ones(1,k*i/10*9)*0.1./(k*i/10*9)]; %1:9
        w=[ones(1,k*i/10)*log_prob(1)./(k*i/10) ones(1,k*i/10)*log_prob(2)./(k*i/10) ones(1,k*i/10)*log_prob(3)./(k*i/10) ones(1,k*i/10)*log_prob(4)./(k*i/10) ones(1,k*i/10)*log_prob(5)./(k*i/10) ones(1,k*i/10)*log_prob(6)./(k*i/10) ones(1,k*i/10)*log_prob(7)./(k*i/10) ones(1,k*i/10)*log_prob(8)./(k*i/10) ones(1,k*i/10)*log_prob(9)./(k*i/10) ones(1,k*i/10)*log_prob(10)./(k*i/10)];
        population=round(linspace(1,k*i,k*i));%indices sampled for population;vectorized
        
        max_coverage=i*k;%number of unique sized DNA fragments
        G_mg=zeros(k,i); %columns correspond to unique genome; row are l-sized DNA fragments
        act_coverage=0;%set sampled coverage to 0
        sequences=0;%set sequences sampled to 0
        for target=target_coverage %iterate through target coverages
            tic %time iterations
            while act_coverage<target%while the sampled coverage is less than the target coverage
                s_indx=randsample(population,num_seq,true,w);%sample num_seq number of indicies from population (replacement=true) with weight matrix
                G_mg(s_indx)=G_mg(s_indx)+1;%sampled indices go from 0 to 1
                act_coverage=numel(find(G_mg~=0))/max_coverage;%determine number of actual indicies not sampled
                sequences=sequences+num_seq;%add number of sequences to total sequenced
            end
            toc%output time it took to get target fraction
            fprintf('\nTarget Fraction: %d \nReplicate: %d \nGenome Number: %d \n',target,rep,i);
            outfile(count,:)=[i target even_flag sequences mean(mean(G_mg*l,1)) rep]; %write data to outfile; the double mean is to calculate mean coverage
            count=count+1;
        end
    end
end


filename = fullfile('.', 'modeled_sequencing_simulation_lognormal.txt');
fid=fopen(filename, 'wt');
fprintf(fid, '%s,%s,%s,%s,%s,%s\n', 'N_genomes','Fraction','Evenness','Sequences','Mean_Coverage','Replicate');  % header
fclose(fid);
dlmwrite("modeled_sequencing_simulation_lognormal.txt",outfile,'precision',15,'-append')

