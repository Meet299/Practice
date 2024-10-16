import sys
import pandas as pd

run_df = pd.read_csv('data/run_metadata.tsv',sep='\t',header='infer')
sample_df = pd.read_csv('data/samples.tsv',sep='\t',header='infer')
#print(run_df)

run_df.index = run_df['run']
sample_df.index = sample_df['sample']

#print(run_df)
#print(sample_df)
chrom_selection = [l.strip() for l in open ("metadata/selection.tsv")]

def get_all_runs_for_a_sample(wildcards):
    #print(wildcards)
    #print(sample_df)
    all_runs = sample_df.loc[wildcards.sample,"runs"].split(",")
    
    run_path_list = []
    for r in all_runs:
        p= run_df.loc[r,'file_path']
        run_path_list.append(p)
    return run_path_list

rule split:
   input:
       all_runs = lambda wildcards: get_all_runs_for_a_sample(wildcards)
   params:
       selected_chromosomes = ",".join (chrom_selection)
   output:
       all_chrom_files = temp(expand("splitted/{{sample}}_chrom_{chrom}.bed", chrom = chrom_selection))
   shell:
       "zcat {input.all_runs} | python scripts/split_by_chromosome.py {params.selected_chromosomes} {wildcards.sample}"
rule sort:
   input:
       chrom_file = "splitted/{sample}_chrom_{chrom}.bed"
   params:
   output:
       sorted_chrom = temp("sorted_by_chrom/{sample}_chrom_{chrom}.bed")
   shell:
       "sort -S4G --parallel=4 -k2,2n -k3,3n {input.chrom_file} > {output.sorted_chrom}" 


rule merge:    
   input:
       all_sorted_files_in_order = lambda wildcards:  ["sorted_by_chrom/{sample}_chrom_{chrom}.bed".format (sample = wildcards.sample, chrom = c) for c in chrom_selection]
   params:
   output:
       sample_based_sorted_file = "sample_based_sorted/{sample}.bed.gz"
   shell:
       "cat {input.all_sorted_files_in_order} | gzip - > {output.sample_based_sorted_file}"

