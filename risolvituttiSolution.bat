cd input
for %%f in (*.txt) do (..\bin\Debug\Assignement.exe -i %%f -o ..\output\summary.csv -s ..\output\solution.csv)