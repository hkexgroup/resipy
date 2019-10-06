#read in 2D survey files and compile them into one file 
#to do this we need to correct the electrode numbering for each survey 
#import modules 
import os, sys
import numpy as np
sys.path.append('..')
from resipy.Survey import Survey

#setup environmental variables 
data_dir = r'C:\Users\jamyd91\Documents\2PhD_projects\Extra\Willington\Data\STG' # Directory where the survey files are stored 
ext = '.stg' # file extension associated with the files 
ftype = 'Sting' # file format of the survey files 
#%% filter out relevant files 
all_files = os.listdir(data_dir)
files = []
for f in all_files:
    if f.endswith(ext):
        files.append(f)
files = sorted(files)

#%% go through and import each survey 
s = Survey(os.path.join(data_dir,files[0]),ftype=ftype) # import first survey 
print('read ' + os.path.join(data_dir,files[0]))
df = s.df
ABMN = np.array((df['a'],df['b'],df['m'],df['n']))
df_master = df.copy()

count = 0
for i in range(1,len(files)):
    s_temp = Survey(os.path.join(data_dir,files[i]),ftype=ftype)
    print('read ' + os.path.join(data_dir,files[i]))
    df_temp = s_temp.df
    #abmn = np.array(df_temp['a'],df_temp['b'],df_temp['m'],df_temp['n'])
    ABMN = np.array((df_master['a'],df_master['b'],df_master['m'],df_master['n']))
    max_idx = np.max(ABMN)
    df_temp['a'] = df_temp['a'] + max_idx
    df_temp['b'] = df_temp['b'] + max_idx
    df_temp['m'] = df_temp['m'] + max_idx
    df_temp['n'] = df_temp['n'] + max_idx
    df_master = df_master.append(df_temp)
    count += 1

#%% write to protocol like file (which can be parsed by ResIPy)
s.df = df_master
fh = open(os.path.join(data_dir,'compiled_protocol.dat'),'w')
fh.write("%i\n"%len(df_master))
df_master = df_master.reset_index()

for i in range(len(df_master)):
    line = "{:4d}\t{:4d}\t{:4d}\t{:4d}\t{:4d}\t{:16f}\n".format(i+1,
            df_master['a'][i],
            df_master['b'][i],
            df_master['m'][i],
            df_master['n'][i],
            df_master['resist'][i])
    fh.write(line)
fh.close()
#s.write2protocol(os.path.join(data_dir,'compiled_protocol.dat'))

#%% test 
s_test = Survey(os.path.join(data_dir,'compiled_protocol.dat'),ftype='Protocol')

    