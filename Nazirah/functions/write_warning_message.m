function write_warning_message(FatalErrorFile, session_name, err_message)

fprintf('%s - Warning: %s\n',session_name,err_message);
fp=fopen(FatalErrorFile,'at');
fprintf(fp,'%s - Warning: %s\n',session_name,err_message);
fclose(fp);