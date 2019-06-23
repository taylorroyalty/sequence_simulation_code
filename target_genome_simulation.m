%Matlab code generating data for Figure 3
%Code inputs
target_coverage=linspace(0.5,1,100);%[0.5 0.7 0.9 0.95 0.99 1];%fraction of exhaustion (genome coverage)
fraction_range=logspace(-2,0,30);%fraction that target genome represents of entire metagenome
num_seq=5e6;%number of sequences sampled at a time. impacts computational time/resolution
rep_num=1;%replicates
l=100;%Read length
filename_str='modeled_sequencing_target_simulation.csv';
Mbp=[15e6 20e6];%size of genomes;

%Simulation
outfile=[];
fraction_range=fliplr(fraction_range);
for g_size=Mbp
    k=g_size-l+1;%unique sample points for a given genome
    for rep=1:1:rep_num%iterate through triplicates
        for p=fraction_range
            fraction=p/k;
            G_mg=zeros(1,k+1);
            w=[ones(1,k)*fraction 1-p]; %this is a weight distribution (i.e., rareness/relative abundance)
            act_coverage=0;
            sequences=0;
            population=linspace(1,k+1,k+1);%indices sampled for population;vectorized
            tic
            for target=target_coverage
                while act_coverage<target
                    s_indx=randsample(population,num_seq,true,w);%sample num_seq number of indicies from population (replacement=true) with weight matrix
                    G_mg(s_indx)=G_mg(s_indx)+1;%sampled indices go from 0 to 1
                    sequences=sequences+num_seq;%add number of sequences to total sequenced
                    act_coverage=numel(find(G_mg(1:end-1)~=0))/k;%determine number of actual indicies not sampled
                end
                outfile=[outfile;[rep sequences p target g_size]];
            end
            toc
        end
    end
end


filename = fullfile('.', filename_str);
fid=fopen(filename, 'wt');
fprintf(fid, '%s,%s,%s,%s,%s\n', 'replicate','sequences','probability','target','genome_size');  % header
fclose(fid);
dlmwrite(filename_str,outfile,'precision',15,'-append')
