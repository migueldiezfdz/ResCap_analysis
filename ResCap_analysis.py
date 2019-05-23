#!/usr/bin/python3

import subprocess
import glob
import sys
import os
import pandas
import argparse

#Parsing the script, description and selection of initila file and number of cores

parser = argparse.ArgumentParser(description=""
                                             "ResCap_analysis_1.0\n"
                                             "Hello world!\t"
                                             "This script allows you to analyze the results from the ResCap capture platform \t"
                                             "It performs an on-target selection, then a subsampling to normalize the samples and "
                                             "erase outliers, and then a selection of the resistance genes with KMA\t"
                                             "\t\tLet's hope it works!  ¯\_(ツ)_/¯\t\t")

parser.add_argument("-c","--clusters",default=1, help= 'Number of clusters available to use')
parser.add_argument("-t","--table", help='CSV files (comma separated) with two columns, called R1 and R2')
args = parser.parse_args()

cores = args.clusters
cores=str(cores)
print('number of cores =', cores)

scriptpath = os.path.dirname(sys.argv[0])

bowtie_db = '/dbs/bowtie_db/plataforma'
bowtie_db_path= scriptpath+bowtie_db
current_directory1 = os.getcwd()
extractR1R2 = '/*R1.fastq.gz'
loop_extract_path = current_directory1+extractR1R2


if args.table is not None:
    table_name = args.table
    df = pandas.read_csv(table_name, sep=',',
                         dtype=str)
    dict_of_samples = (df.set_index('R1')['R2'].to_dict())
    for key, value in dict_of_samples.items():
        salida_t = key+'.bow.gz'
        subprocess.call(['bowtie2', '-x', bowtie_db_path, '-p', cores, '-1', key, '-2', value, '-q', '--al-conc-gz', salida_t])
else:
    for R1 in glob.iglob(loop_extract_path):
        R2 = R1.replace('R1', 'R2')
        salida = R1.replace('.R1.fastq.gz', '.bow.gz')
        subprocess.call(['bowtie2','-x',bowtie_db_path,'-p',cores,'-1',R1, '-2', R2,'-q', '--al-conc-gz', salida])


dict_of_reads={}
list_of_reads=[]



seqkit='/bin/seqkit'
seqkit_path = scriptpath+seqkit

cmd= "{0} stats *b* -T -j {1} -o Tabla.txt".format(seqkit_path, cores)
subprocess.call(cmd, shell=True)


df = pandas.read_csv('Tabla.txt', sep='\t', dtype=str)
df1 = df[['file','num_seqs']]
df1.num_seqs = df1.num_seqs.astype(float)

dict_of_reads = (df1.set_index('file')['num_seqs'].to_dict())

for key, value in dict_of_reads.items():
    list_of_reads.append(value)

lowest_num = min(list_of_reads)
lowest_num = str(int(lowest_num))

sum=0
for num in list_of_reads:
    sum = sum +num
average= sum/len(list_of_reads)
outliers_value_norm=(average*0.1)
outliers_value=int(average*0.1)
outliers_value=str(outliers_value)


outliers_list=[]
outliers_value_list=[]

for key, value in dict_of_reads.items():
    if value <= outliers_value_norm:
        outliers_list.append(key)
        outliers_value_list.append(value)

os.makedirs('./outliers')
for f in outliers_list:
    subprocess.call(['mv', f,'./outliers'])

current_directory1 = os.getcwd()
extract_b = '/*.bow*'
loop_extract_b_path = current_directory1+extract_b

seqkit='/bin/seqkit'
seqkit_path = scriptpath+seqkit

for original in glob.iglob(loop_extract_b_path):
    salida1 = original.replace('.bow.', '.kma.')
    subprocess.call([seqkit_path, 'sample', original, '-n', lowest_num, '-j', cores, '-o', salida1])


current_directory1 = os.getcwd()
extract_c = '/*kma.1.gz'
loop_extract_c_path = current_directory1+extract_c

kma = '/bin/kma/kma'
kma_path= scriptpath+kma

kma_db = '/dbs/kma_db/data-base-kma'
kma_db_path= scriptpath+kma_db
print(kma_db_path)


for c1 in glob.iglob(loop_extract_c_path):
    c2 = c1.replace('kma.1.gz', 'kma.2.gz')
    salida1 = c1.replace('.kma.1.gz', '')
    print(c1,c2,salida1)
    subprocess.call([kma_path,'-ipe', c1, c2, '-o', salida1, '-t_db', kma_db_path, '-t', cores])