function write_warning_message(FatalErrorFile, err_message)

fprintf('%s - Warning: %s\n',SESSIONS(jSession).name,err_message);
fp=fopen(FatalErrorFile,'at');
fprintf(fp,'%s - Warning: %s\n',SESSIONS(jSession).name,err_message);
fclose(fp);