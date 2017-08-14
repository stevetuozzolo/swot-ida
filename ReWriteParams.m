%% Function to rewrite params file
function ReWriteParams(Qrng)
    fid=fopen('Olentangy\Params.txt','w');
    fprintf(fid,'%s\n','Chain length');
    fprintf(fid,'%f\n',10000);
    fprintf(fid,'%s\n','Burn in');
    fprintf(fid,'%f\n',2000);
    fprintf(fid,'%s\n','Random Number Seed');
    fprintf(fid,'%f\n',7);
    fprintf(fid,'%s\n','Number of times to use in estimate');
    fprintf(fid,'%i\t%i\n',[20 50]);
    fprintf(fid,'%s\n','Prior mean on mean discharge [m3/s]');
    fprintf(fid,'%f\n',Qrng);
    fprintf(fid,'%s\n','Prior coefficient of variation on mean flow [fraction]');
    fprintf(fid,'%f\n',.5);
end