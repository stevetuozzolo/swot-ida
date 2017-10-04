function Truth = ReadTruthNetCDF (fname,DAll)

if strcmp(fname(end-6:end-3),'sity')
    group_name='/Reach_Timeseries_Simulator/';
    true_name='/Reach_Timeseries_True/';
else
    group_name='/Reach_Timeseries/';
end

gdrch=ncread(fname,'River_Info/gdrch');
Truth.A=ncread(fname,[true_name 'A']);
Truth.A=Truth.A(gdrch,:);
Truth.A0=Truth.A(:,1)';
Truth.q=0; %filler
Truth.n=NaN;%mean(mean(ncread(fname,'XS_Timeseries/n'))); %filler
Truth.Q=ncread(fname,[group_name 'Q']);
Truth.Q=Truth.Q(gdrch,:);
Truth.W=ncread(fname,[group_name 'W']);
Truth.W=Truth.W(gdrch,:);
Truth.h=ncread(fname,[group_name 'H']);
Truth.h=Truth.h(gdrch,:);

for i=1:size(Truth.A,2)
    Truth.dA(:,i)=Truth.A(:,i)-Truth.A0';
end

Truth.dAv=reshape(Truth.dA',DAll.nR*DAll.nt,1);

Truth.hv=reshape(Truth.h',DAll.nR*DAll.nt,1);

Truth.Wv=reshape(Truth.W',DAll.nR*DAll.nt,1);

return