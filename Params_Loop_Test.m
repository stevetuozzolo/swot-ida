Discharge=[50:50:1000]
for i=1:length(Discharge)
    fid=fopen('RunFile.txt');
    RunName=fgetl(fid);
    pfile=[RunName '\params.txt']
    A=regexp(fileread(pfile),'\n','split');
    A{10}=num2str(Discharge(i));
    fid=fopen(pfile,'w');
    fprintf(fid,'%s\n',A{:});
    fclose(fid)  
    
    RunExpsNetCDF;
    load([RunName '\RunData.mat'])
    rRMSE(i)=Err.Stats.rRMSE;
    bias(i)=Err.Stats.bias;
    RMSE(i)=Err.Stats.RMSE;
end