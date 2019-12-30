@echo off
title Converting STRUCTURE files into GENEPOP files using PGDSpider in a for loop
:: Put this batch file in the folder where PDGSPider is installed for it to run properly

FOR /L %%A IN (1,1,1000) DO (
  ECHO %%A
  PGDSpider2-cli.exe -inputfile bottleneck_files\raw\bottleneck_in%%A.stru -inputformat STRUCTURE -outputfile bottleneck_files\bottleneck_in%%A.txt -outputformat GENEPOP -spid SPID_struct_genepop_39msats.spid
 
)
