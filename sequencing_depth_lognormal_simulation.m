%Matlab code generating data in Figure 3
target_coverage=[1];%target coverages to log for
even_flag=0;%arbitrary flag reflecting weight matrix used
l=100;%Read length; only important for determining coverate
num_seq=1e7;%number of sequences sampled at a time.
rep_num=1;%replicates
s=1;%number of operationally defined genomes for sample environment; each index a different environment
Mbp=1e6;%size of genomes; could be made into a vector in future
k=Mbp;%-l+1;%unique sample points for a given genome
outfile=nan(numel(s)*rep_num*numel(target_coverage),7); %predefined logging matrix
count=1;%index for writing to log matrix

%For log normal weight distribution
log_fun=@(S,a,R)S*exp(-1*a.^2*R.^2);
a=[0.5]; S=1; R=0:20:100-20;%octaves grouping four genomes together %linspace(0,100,20);


actual_timeseries=[];
groups=length(R);
for rep=1:1:rep_num%iterate through triplicates
    for i=s%iterate through
        for a_value=a
            log_prob=log_fun(S,a_value,R)/sum(log_fun(S,a_value,R));
            w=[];
            for it=1:1:groups
                w=[w ones(1,k*i/groups)*log_prob(it)./(k*i/groups)];
            end
%             w=[ones(1,k*i/10)*log_prob(1)./(k*i/10) ones(1,k*i/10)*log_prob(2)./(k*i/10) ones(1,k*i/10)*log_prob(3)./(k*i/10) ones(1,k*i/10)*log_prob(4)./(k*i/10) ones(1,k*i/10)*log_prob(5)./(k*i/10) ones(1,k*i/10)*log_prob(6)./(k*i/10) ones(1,k*i/10)*log_prob(7)./(k*i/10) ones(1,k*i/10)*log_prob(8)./(k*i/10) ones(1,k*i/10)*log_prob(9)./(k*i/10) ones(1,k*i/10)*log_prob(10)./(k*i/10)];
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
                    actual_timeseries=[actual_timeseries;act_coverage a_value];
                    sequences=sequences+num_seq;%add number of sequences to total sequenced
                end
                fprintf('\nTarget Fraction: %d \nReplicate: %d \nGenome Number: %d \na_value Number: %d \nSequences Needed: %d \n',target,rep,i,a_value,sequences);
                toc%output time it took to get target fraction
                shannon=sum(-w.*log(w));
                outfile(count,:)=[i target a_value sequences mean(mean(G_mg*l,1)) rep shannon]; %write data to outfile; the double mean is to calculate mean coverage
                count=count+1;
            end
        end
    end
end


filename = fullfile('.', 'modeled_sequencing_simulation_lognormal.txt');
fid=fopen(filename, 'wt');
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s\n', 'N_genomes','Fraction','a_value','Sequences','Mean_Coverage','Replicate','shannon');  % header
fclose(fid);
dlmwrite("modeled_sequencing_simulation_lognormal.txt",outfile,'precision',15,'-append')

