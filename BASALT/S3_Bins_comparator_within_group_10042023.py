#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S3 (CheckM2 branch): Within-group bin comparison and selection.

This module compares bins from the same assembly group based on
coverage and PE connections, and selects the best representative
bins for each group.
"""

from tempfile import TemporaryFile
from Bio import SeqIO
from S2_BinsAbundance_PE_connections_multiple_processes_pool_10032023 import *
import copy
import os

FASTA_SUFFIXES = ('.fa', '.fasta', '.fna', '.fas', '.fsa')


def _is_fasta_file(filename):
    lower = filename.lower()
    return any(lower.endswith(suffix) for suffix in FASTA_SUFFIXES)


def _strip_fasta_suffix(filename):
    lower = filename.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return filename[:-len(suffix)]
    return os.path.splitext(filename)[0]


def _resolve_bin_filename(folder, bin_id):
    if not os.path.isdir(folder):
        return None

    candidate_path = os.path.join(folder, bin_id)
    if os.path.isfile(candidate_path):
        return bin_id

    target_base = _strip_fasta_suffix(bin_id)
    for file in os.listdir(folder):
        if _is_fasta_file(file) and _strip_fasta_suffix(file) == target_base:
            return file
    return None


def _is_legacy_singlecontig_zero_bin(bin_id):
    """
    Match the historical placeholder-style ``*_genomes.0`` bins only.

    Older BASALT code used the broad substring check ``'_genomes.0' in id``,
    which also matches normal zero-padded IDs like ``*_genomes.001`` and
    silently drops most real bins. We only want to exclude the exact
    ``..._genomes.0`` placeholder case here.
    """
    base = _strip_fasta_suffix(bin_id)
    return base.endswith('_genomes.0') and '_semibin_genomes.0' not in base


def contig_id_recorder(genome_folder):
    """
    Record contig IDs and lengths for all bins in a genome folder list.

    Parameters
    ----------
    genome_folder : list of str
        List of binset folder name prefixes.

    Returns
    -------
    tuple
        (relation, best_hit_genome, bins_extract)
        relation and best_hit_genome are initially empty and populated
        later; bins_extract is the list of bins not selected as best hits.
    """
    genomes_sum={}
    pwd=os.getcwd()
    n=0
    for item in genome_folder:
        n+=1
        print('Parsing '+item+' group of bin-set')
        print('--------------------------')
        if n == 1:
            genomes_sum[str(item)]={}
            relation, bin_len, bins1={}, {}, []
            frist=item
            for root, dirs, files in os.walk(item+'_genomes'):
                os.chdir(item+'_genomes')
                for file in files:
                    if _is_fasta_file(file):
                        bins1.append(file)
                        genomes_sum[item][file]={}
                        genomes_sum[item][file]['contig']={}
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            genomes_sum[item][file]['contig'][record.id]=len(record.seq)
                    else:
                        continue
            os.chdir(pwd)
        else:
            for root, dirs, files in os.walk(item+'_genomes'):
                os.chdir(item+'_genomes')
                bins2=[]
                for file in files:
                    if _is_fasta_file(file):
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue
                    else:
                        continue

                os.chdir(pwd)

    a, b, best_hit_genome, bins1_selected, bins2_selected, bins_extract={}, {}, {}, {}, {}, []
    for item in relation.keys():
        bin_set1=item.split('---')[0]
        bin_set2=item.split('---')[1]
        # print(bin_set1
        if bin_set1 not in a.keys():
            a[bin_set1]=relation[item]
            b[bin_set1]=item+'\t'+str(relation[item])+'\t'+str(round(float(relation[item])*float(bin_len[bin_set1])/float(bin_len[bin_set2]), 2))
        elif float(relation[item]) > float(a[bin_set1]):
            # print(bin_set1, item, float(relation[item]), float(a[bin_set1])
            a[bin_set1]=relation[item]
            b[bin_set1]=item+'\t'+str(relation[item])+'\t'+str(round(float(relation[item])*float(bin_len[bin_set1])/float(bin_len[bin_set2]), 2))
        else:
            continue

    for item in b.keys():
        bin_set1=b[item].split('\t')[0].split('---')[0]
        bin_set2=b[item].split('\t')[0].split('---')[1]
        bin1_score=float(b[item].split('\t')[1])
        bin2_score=float(b[item].split('\t')[2])
        total_score=bin1_score+bin2_score
        # if total_score >= 120 or bin1_score >= 80 or bin2_score >= 80:
        if total_score >= 100 or bin1_score >= 50 or bin2_score >= 50:
            best_hit_genome[item]=str(b[item])
            bins1_selected[bin_set1]=0
            bins2_selected[bin_set2]=0
    
    for item in bins1:
        if item not in bins1_selected.keys():
            bins_extract.append(item)

    for item in bins2:
        if item not in bins2_selected.keys():
            bins_extract.append(item)

    return relation, best_hit_genome, bins_extract

def checkm(genome_folder):
    """
    Read CheckM-like quality statistics for all bins in the given folders.

    Parameters
    ----------
    genome_folder : list of str
        List of binset folder name prefixes.

    Returns
    -------
    dict
        Mapping bin_id -> quality metrics (completeness, contamination,
        genome size, N50, optional marker lineage).
    """
    pwd=os.getcwd()
    bins_checkm={}
    print(pwd)
    for item in genome_folder:
        print('Reading checkm results of '+item)
        print('-------------------------')
        try:
            f=open(item+'_genomes/'+item+'_quality_report.tsv', 'r')
            print('Reading '+item+'_quality_report.tsv')
            n=0
            for line in f:
                n+=1
                if n >= 2:
                    genome_ids=str(line).strip().split('\t')[0]
                    if not _is_legacy_singlecontig_zero_bin(genome_ids):
                        bins_checkm[genome_ids]={}
                        bins_checkm[genome_ids]['Completeness']=float(str(line).strip().split('\t')[2].strip())
                        bins_checkm[genome_ids]['Genome size']=int(str(line).strip().split('\t')[1].strip())
                        bins_checkm[genome_ids]['N50']=float(str(line).strip().split('\t')[4].strip())
                        bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('\t')[3].strip())
                    elif '_semibin_genomes.0' in _strip_fasta_suffix(genome_ids):
                        bins_checkm[genome_ids]={}
                        bins_checkm[genome_ids]['Completeness']=float(str(line).strip().split('\t')[2].strip())
                        bins_checkm[genome_ids]['Genome size']=int(str(line).strip().split('\t')[1].strip())
                        bins_checkm[genome_ids]['N50']=float(str(line).strip().split('\t')[4].strip())
                        bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('\t')[3].strip())
        except:
            print('CheckM output file reading error!')

        # print('Reading dataframe-connections of', item
        # print('-------------------------'
        
        try:
            f=open(item+'_genomes/Bins_total_connections_'+str(item)+'.txt', 'r')
            print('Reading dataframe-connections of Bins_total_connections_'+item+'.txt')
            n=0
            for line in f:
                n+=1
                if n >= 2 and not _is_legacy_singlecontig_zero_bin(str(line).strip().split('\t')[0]):
                    bins_checkm[str(line).strip().split('\t')[0]]['Connections']=int(str(line).strip().split('\t')[1])
        except:
            print('Please make sure Bins_total_connections_'+str(item)+'.txt under the folder.')

        for item2 in bins_checkm.keys():
            if 'Connections' not in bins_checkm[item2].keys():
                bins_checkm[item2]['Connections']=0

    return bins_checkm

def genome_selector(best_hit_genome, bin_set_checkm):
    """
    Select best representative bins between candidate bin pairs.

    Parameters
    ----------
    best_hit_genome : dict
        Mapping representing best hit pairs between bins.
    bin_set_checkm : dict
        Mapping bin_id -> CheckM-derived metrics.

    Returns
    -------
    dict
        Mapping selected_bin_id -> summary string of comparison.
    """
    print('Selecting bin-set')
    print('------------------')
    bin_selected={}
    # marker_score={'root':0, 'k':1, 'p':1.5, 'c':2.3, 'o':3.4, 'f':5.1, 'g':7.6, 's':11.4}
    for item in best_hit_genome.keys():
        set1=str(best_hit_genome[item]).split('\t')[0].split('---')[0]
        set2=str(best_hit_genome[item]).split('\t')[0].split('---')[1]
        set1=_strip_fasta_suffix(set1)
        set2=_strip_fasta_suffix(set2)

        set1_cpn=bin_set_checkm[set1]['Completeness']
        set2_cpn=bin_set_checkm[set2]['Completeness']
        set1_ctn=bin_set_checkm[set1]['Contamination']
        set2_ctn=bin_set_checkm[set2]['Contamination']

        set1_cpn_ctn=float(set1_cpn)-float(set1_ctn)
        set2_cpn_ctn=float(set2_cpn)-float(set2_ctn)
        if float(set1_cpn_ctn) == float(set2_cpn_ctn):
            if float(bin_set_checkm[set1]['Genome size']) > float(bin_set_checkm[set2]['Genome size']):
                bin_selected[set1]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
            else:
                bin_selected[set2]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])

        elif float(set1_cpn_ctn) > float(set2_cpn_ctn):
            bin_selected[set1]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
        else:
            bin_selected[set2]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])

    print('bin-selecting accomplished!')
    print('---------------------------')
    return bin_selected

def two_groups_comparator(assembly, binset1, binset2, num):
    """
    Compare two binsets belonging to the same assembly group.

    Parameters
    ----------
    assembly : str
        Assembly name.
    binset1 : str
        First binset folder prefix.
    binset2 : str
        Second binset folder prefix.
    num : int
        Iteration index used in filenames.

    Returns
    -------
    None
        Writes comparison outputs and best-bin selections to disk.
    """
    pwd=os.getcwd()
    try:
        fx=open('Basalt_log.txt','a')
    except:
        ttttx=0
    print('Iteration '+str(num))
    print('Comparing '+binset1+' and '+binset2)
    fx.write('Iteration '+str(num)+'\n')
    fx.write('Comparing '+binset1+' and '+binset2+'\n')
    f=open('All_possible_bin_sets_iteration_'+str(num)+'.txt','w')
    f2=open('Best_hit_bin_sets_iteration_'+str(num)+'.txt','w')

    genome_folder=[binset1, binset2]

    a=contig_id_recorder(genome_folder)
    all_hit_genome=a[0]
    best_hit_genome=a[1]
    bins_extract=a[2]

    for item in all_hit_genome.keys():
        f.write(item+'\t'+str(all_hit_genome[item])+'\n')
    f.close()

    for item in best_hit_genome.keys():
        f2.write(str(best_hit_genome[item])+'\n')
    f2.close()

    f=open('Bins_marker_lineage_completeness_contamination_iteration_'+str(num)+'.txt', 'w')
    
    bin_set_checkm=checkm(genome_folder)

    for item in bin_set_checkm.keys():
        try:
            f.write(item+'\t'+str('Completeness:'+str(bin_set_checkm[item]['Completeness']))+'\t'+'Contaomination:'+str(bin_set_checkm[item]['Contamination'])+'\t'+'Genome size:'+str(bin_set_checkm[item]['Genome size'])+'\t'+'N50:'+str(bin_set_checkm[item]['N50'])+'\n')
        except:
            print(item+' error')
    f.close()

    bin_selected=genome_selector(best_hit_genome, bin_set_checkm)
    if len(bin_selected) == 0 and len(bin_set_checkm) != 0:
        print(
            'No bins were selected in iteration '
            + str(num)
            + '. Restoring all candidate bins for downstream selection.'
        )
        fx.write(
            'No bins were selected in iteration '
            + str(num)
            + '. Restoring all candidate bins for downstream selection.\n'
        )
        for item in bin_set_checkm.keys():
            bin_selected[item]='fallback candidate from '+binset1+' and '+binset2
    
    f=open('Extract_bins_in_iteration_'+str(num)+'.txt', 'w')
    for item in bins_extract:
        f.write(item+'\n')
    f.close()

    f=open('Best_bin_set_iteration_'+str(num)+'.txt','w')
    for item in bin_selected.keys():
        f.write(item+'\t'+str(bin_selected[item])+'\n')

    if len(bins_extract) >= 1:
        for item in bins_extract:
            name=_strip_fasta_suffix(item)

            f.write(item+'\t'+str(bin_set_checkm[name])+'\n')
            bin_selected[name]='unique genome in', binset2
            print(item+' unique genome in '+binset2)
            print('----------------')
            fx.write(item+' unique genome in '+binset2+'\n')
            fx.write('----------------'+'\n')
    f.close()

    try:
        os.mkdir('Iteration_'+str(num)+'_genomes')
    except:
        os.system('rm -rf Iteration_'+str(num)+'_genomes')
        os.mkdir('Iteration_'+str(num)+'_genomes')
        print('Iteration_'+str(num)+'_genomes exist')
        print('Re-created folder of Iteration_'+str(num)+'_genomes')
        fx.write('Iteration_'+str(num)+'_genomes exist'+'\n')
        fx.write('Re-created folder of Iteration_'+str(num)+'_genomes'+'\n')
    
    f=open('Iteration_'+str(num)+'_genomes/Iteration_'+str(num)+'_quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')

    f2=open('Iteration_'+str(num)+'_genomes/Bins_total_connections_Iteration_'+str(num)+'.txt','w')
    f2.write('Bin'+'\t'+'Total_connections'+'\n')

    copied_bins=0
    for item in bin_selected.keys():
        f.write(item+'\t'+str(bin_set_checkm[item]['Genome size'])+'\t'+str(bin_set_checkm[item]['Completeness'])+'\t'+str(bin_set_checkm[item]['Contamination'])+'\t'+str(bin_set_checkm[item]['N50'])+'\n')
        try:
            folder=item.split('_genomes.')[0]
            source_folder = pwd+'/'+folder+'_genomes'
            os.chdir(source_folder)
            bin_file = _resolve_bin_filename(source_folder, item)
            if bin_file is None:
                print('Copy bin-set error! Missing FASTA for '+str(item)+' in '+str(source_folder))
                continue
            os.system('cp '+bin_file+' '+pwd+'/Iteration_'+str(num)+'_genomes')
            copied_bins += 1
        except:
            print('Copy bin-set error!')
    f.close()  
    f2.close()
    print(
        'Copied '
        + str(copied_bins)
        + ' bin FASTA file(s) into Iteration_'
        + str(num)
        + '_genomes'
    )

    os.chdir(pwd)
    fx.write(binset1+' and '+binset2+' comparison done'+'\n')
    try:
        f_s3=open('S3_checkpoint.txt','a')
    except:
        xyzzzz=0
    f_s3.write(str(assembly)+'\t'+str(num)+'\t'+binset1+' and '+binset2+' comparison done'+'\n')

def bin_within_a_group_comparitor(binset, assembly, num):
    """
    Compare bins within a single binset and select best representatives.

    Parameters
    ----------
    binset : str
        Binset folder prefix.
    assembly : str
        Assembly name.
    num : int
        Iteration index used in filenames.

    Returns
    -------
    None
        Writes intermediate and best-bin files to disk.
    """
    pwd=os.getcwd()
    print('Comparing bins in final iteration')
    print('Parsing bins')
    print('---------------------------------')
    bins_sum={}
    for root, dirs, files in os.walk(binset):
        os.chdir(binset)
        for file in files:
            hz=file.split('.')[-1].lower()
            if hz in ('fa', 'fasta', 'fna', 'fas', 'fsa'):
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)

    relation, relation2, processed={}, {}, {}
    for item in bins_sum.keys():
        processed[item]=''
        for contig in bins_sum[item]['contig'].keys():
            for item2 in bins_sum.keys():
                if item2 != item and item2 not in processed.keys():
                     if contig in bins_sum[item2]['contig'].keys():
                        if str(item+'---'+item2) not in relation.keys():
                            relation[str(item+'---'+item2)]=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item]['totallen'])),2)
                            relation2[str(item+'---'+item2)]=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item2]['totallen'])),2)
                        else:
                            relation[str(item+'---'+item2)]+=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item]['totallen'])),2)
                            relation2[str(item+'---'+item2)]+=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item2]['totallen'])),2)
                else:
                    continue
        #del bins_sum[item]

    final_iteration_checkm={}
    f=open('Iteration_'+str(num)+'_quality_report.tsv', 'r')
    n=0
    for line in f:
        n+=1
        if n >= 2:
            ids=str(line).strip().split('\t')[0].strip()
            completeness=float(str(line).strip().split('\t')[2].strip())
            contamination=float(str(line).strip().split('\t')[3].strip())
            N50=float(str(line).strip().split('\t')[4].strip())
            genome_size=int(str(line).strip().split('\t')[1].strip())
            bin_id = _resolve_bin_filename(os.path.join(pwd, binset), ids)
            if bin_id is None:
                bin_id = ids+'.fa'
            final_iteration_checkm[bin_id]={'Completeness': completeness, 'Genome size': genome_size, 'N50': N50, 'Contamination': contamination}

    os.chdir(pwd)
    try:
        os.system('mkdir '+str(assembly)+'_BestBinsSet')
    except:
        os.system('rm -rf '+str(assembly)+'_BestBinsSet')
        os.system('mkdir '+str(assembly)+'_BestBinsSet')
        print(str(assembly)+'_BestBinsSet exist')
        print('Re-created folder of '+str(assembly)+'_BestBinsSet')

    os.chdir(pwd+'/'+str(assembly)+'_BestBinsSet')
    pos_bins={}
    f=open('Possible_same_bin.txt','w')
    f2=open('Highly_possible_same_bin.txt','w')
    for item in relation.keys():
        # if float(relation[item]) >= 50:
        f.write(item+'\t'+str(relation[item])+'\t'+str(relation2[item])+'\n')
        sim=float(relation[item]) + float(relation2[item])
        if sim >= 100 or float(relation[item]) >= 50 or float(relation2[item]) >= 50:
        # if sim >= 120 or float(relation[item]) >= 80 or float(relation2[item]) >= 80:
            f2.write(item+'\t'+str(relation[item])+'\t'+str(relation2[item])+'\t'+str(final_iteration_checkm[str(item).split('---')[0]])+'\t'+str(final_iteration_checkm[str(item).split('---')[1]])+'\n')
            pos_bins[item]=str(relation[item])+'\t'+str(relation2[item])

    f.close()
    f2.close()

    original_final_iteration_checkm=copy.deepcopy(final_iteration_checkm)
    remain_bin, del_bin={}, {}
    for item in pos_bins.keys():
        set1=str(item).split('---')[0]
        set2=str(item).split('---')[1]

        set1_cpn=final_iteration_checkm[set1]['Completeness']
        set2_cpn=final_iteration_checkm[set2]['Completeness']
        set1_ctn=final_iteration_checkm[set1]['Contamination']
        set2_ctn=final_iteration_checkm[set2]['Contamination']

        set1_cpn_ctn=float(set1_cpn)-float(set1_ctn)
        set2_cpn_ctn=float(set2_cpn)-float(set2_ctn)

        if float(set1_cpn_ctn) == float(set2_cpn_ctn):
            if float(final_iteration_checkm[set1]['Genome size']) > float(final_iteration_checkm[set2]['Genome size']):
                remain_bin[set1]=0
                del_bin[set2]=0
            else:
                remain_bin[set2]=0
                del_bin[set1]=0
        elif float(set1_cpn_ctn) > float(set2_cpn_ctn):
            remain_bin[set1]=0
            del_bin[set2]=0
        else:
            remain_bin[set2]=0
            del_bin[set1]=0

    for item in del_bin.keys():
        if item in final_iteration_checkm.keys():
            del final_iteration_checkm[item]
            print('Deleted '+item+' from the selected bins set')
        else:
            continue

    if len(final_iteration_checkm) == 0 and len(original_final_iteration_checkm) != 0:
        final_iteration_checkm=copy.deepcopy(original_final_iteration_checkm)
        print('Final iteration selection removed every bin for '+str(assembly)+'. Restoring the pre-filter candidate set.')

    bin_selected={}
    f=open(assembly+'_BestBinSet_quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')

    for item in final_iteration_checkm.keys():
        checkm_id_list=item.split('.')
        checkm_id_list.remove(checkm_id_list[-1])
        checkm_id='.'.join(checkm_id_list)
        # f.write(checkm_id+'\t'+str(final_iteration_checkm[item])+'\n')
        f.write(checkm_id+'\t'+str(final_iteration_checkm[item]['Genome size'])+'\t'+str(final_iteration_checkm[item]['Completeness'])+'\t'+str(final_iteration_checkm[item]['Contamination'])+'\t'+str(final_iteration_checkm[item]['N50'])+'\n')
        bin_selected[checkm_id]=0
    f.close()

    selected_bin_files={}
    expected_selected={_strip_fasta_suffix(item): item for item in final_iteration_checkm.keys()}
    os.chdir(pwd+'/'+'Iteration_'+str(num)+'_genomes')
    for root, dirs, files in os.walk(pwd+'/'+'Iteration_'+str(num)+'_genomes'):
        for file in files:
            if _is_fasta_file(file):
                file_base=_strip_fasta_suffix(file)
                if file_base in expected_selected.keys():
                    selected_bin_files[file]=expected_selected[file_base]
                    os.system('cp '+file+' '+pwd+'/'+str(assembly)+'_BestBinsSet')

    os.chdir(pwd+'/'+str(assembly)+'_BestBinsSet')
    f3=open(str(assembly)+'_BestBinsSet.depth.txt','w')
    f4=open('prebinned_genomes_output_for_dataframe_'+str(assembly)+'_BestBinsSet.txt','w')
    f5=open('Genome_group_all_list_'+str(assembly)+'_BestBinsSet.txt','w')

    if len(selected_bin_files) == 0:
        print('No selected FASTA bins were copied into '+str(assembly)+'_BestBinsSet')
        fallback_selected={_strip_fasta_suffix(item): item for item in original_final_iteration_checkm.keys()}
        os.chdir(pwd+'/'+'Iteration_'+str(num)+'_genomes')
        for root, dirs, files in os.walk(pwd+'/'+'Iteration_'+str(num)+'_genomes'):
            for file in files:
                if _is_fasta_file(file):
                    file_base=_strip_fasta_suffix(file)
                    if file_base in fallback_selected.keys():
                        selected_bin_files[file]=fallback_selected[file_base]
                        os.system('cp '+file+' '+pwd+'/'+str(assembly)+'_BestBinsSet')
        os.chdir(pwd+'/'+str(assembly)+'_BestBinsSet')
        if len(selected_bin_files) != 0:
            final_iteration_checkm=copy.deepcopy(original_final_iteration_checkm)
            print('Fallback copied '+str(len(selected_bin_files))+' FASTA bin(s) into '+str(assembly)+'_BestBinsSet')

    binset_folders={}
    for item in selected_bin_files.keys():
        if '_genomes.' in item:
            binset_folders[item.split('_genomes.')[0]]=0

    print('----------------------------')
    print('Parsing bins in best bin-set')
    contig_bin={}
    for root, dirs, files in os.walk(pwd+'/'+str(assembly)+'_BestBinsSet'):
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if _is_fasta_file(hz):
                    for record in SeqIO.parse(file, 'fasta'):
                        contig_bin[record.id]=str(file)

    title1, title2, title3=[], [], []
    for item in binset_folders.keys():
        print('Parsing files in folder '+item)
        for root, dirs, files in os.walk(pwd+'/'+str(item)+'_genomes'):
            os.chdir(pwd+'/'+str(item)+'_genomes')
            for file in files:
                if '.depth.txt' in file:
                    n=0
                    for line in open(file, 'r'):
                        n+=1
                        if n == 1 and len(title1) == 0:
                            title1.append(str(line))
                            f3.write(str(line))
                        else:
                            if str(line).strip().split('\t')[0] in contig_bin.keys():
                                f3.write(str(line))

                if 'prebinned_genomes_output_for_dataframe_' in file:
                    n=0
                    for line in open(file, 'r'):
                        n+=1
                        if n == 1 and len(title2) == 0:
                            title2.append(str(line))
                            f4.write(str(line))
                        else:
                            if str(line).strip().split('\t')[1] in bin_selected.keys():
                                f4.write(str(line))
                
                if 'Genome_group_all_list_' in file:
                    n=0
                    for line in open(file, 'r'): 
                        n+=1
                        if n == 1 and len(title3) == 0:
                            f5.write(str(line))
                            title3.append(str(line))
                        else:
                            if str(line).strip().split('\t')[0] in bin_selected.keys():
                                f5.write(str(line))

    f3.close()
    f4.close()
    f5.close()
    os.chdir(pwd)
    try:
        f_s3=open('S3_checkpoint.txt','a')
    except:
        xyzzzz=0
    num2=num+1
    f_s3.write(str(assembly)+'\t'+str(num2)+'\t'+str(assembly)+'_BestBinsSet done'+'\n')
    return str(assembly)+'_BestBinsSet'

def binset_filtration(binset):
    """
    Filter bins based on connections, coverage and quality metrics.

    Parameters
    ----------
    binset : str
        Binset folder prefix.

    Returns
    -------
    None
        Modifies binset contents in-place (e.g. removing low-quality bins).
    """
    print('Parsing '+binset)
    pwd=os.getcwd()
    os.chdir(pwd+'/'+binset)
    del_bin, bin_checkm=[], []
    for root, dirs, files in os.walk(pwd+'/'+binset):
        for file in files:
            if '_quality_report.tsv' in file:
                n=0
                for line in open(file, 'r'):
                    n+=1
                    if n >= 2:
                        bin_id=str(line).strip().split('\t')[0]
                        bin_id_f=_resolve_bin_filename(pwd+'/'+binset, bin_id)
                        if bin_id_f is None:
                            bin_id_f=bin_id+'.fa'
                        bin_checkm.append(bin_id_f)

                        try:
                            Completeness=float(str(line).strip().split('\t')[2].strip())
                            genome_size=int(str(line).strip().split('\t')[1].strip())
                        except:
                            Completeness=0
                            genome_size=0
    
                        if Completeness <= 5:
                            del_bin.append(bin_id_f)
                    
                        if genome_size <= 200000:
                            del_bin.append(bin_id_f)
            
        for file in files:
            if _is_fasta_file(file):
                if file not in bin_checkm:
                    del_bin.append(file)
            else:
                continue

    # fx=open('test.txt','w')
    # for item in bin_checkm:
    #     fx.write(str(item)+'\n')
    # fx.close()

    os.mkdir('Remove_bins')
    for root, dirs, files in os.walk(pwd+'/'+binset):
        for file in files:
            if file in del_bin:
                # print(file)
                os.system('mv '+file+' Remove_bins')
    
    os.system('tar zcvf Remove_bins.tar.gz Remove_bins')
    os.system('rm -rf Remove_bins')

    os.chdir(pwd)

def bins_comparator_multiple_groups(genome_folder, assembly):
    """
    Top-level entry for S3: compare bins across multiple groups.

    Parameters
    ----------
    genome_folder : list of str
        List of binset folder prefixes for all groups.
    assembly : str
        Assembly name.

    Returns
    -------
    str
        Name of the final best-binset folder for the assembly.
    """
    try:
        f_s3=open('S3_checkpoint.txt','a')
    except:
        f_s3=open('S3_checkpoint.txt','w')
    
    finished_step=[]
    for line in open('S3_checkpoint.txt','r'):
        ass=str(line).strip().split('\t')[0].strip()
        if ass == assembly:
            finished_step.append(str(line).strip().split('\t')[1].strip())

    for item in genome_folder:
        binset_filtration(item)

    num_binset=len(genome_folder)
    num=0
    for num in range(0,num_binset-1):
        num+=1
        if num == 1:
            if str(num) not in finished_step:
                two_groups_comparator(assembly, str(genome_folder[0]).split('_genomes')[0], str(genome_folder[-1]).split('_genomes')[0], num)
        else:
            genome_folder.remove(genome_folder[0])
            genome_folder.remove(genome_folder[-1])
            genome_folder.append('Iteration_'+str(num-1))
            if str(num) not in finished_step:
                two_groups_comparator(assembly, str(genome_folder[-1]).split('_genomes')[0], str(genome_folder[0]).split('_genomes')[0], num)
    i=0
    while i < num_binset-2:
        i+=1
        os.system('rm -rf Iteration_'+str(i)+'_genomes')

    bestbinset = str(assembly)+'_BestBinsSet'
    num2=num+1
    if str(num2) not in finished_step:
        bestbinset=bin_within_a_group_comparitor('Iteration_'+str(num)+'_genomes', str(assembly), num)    
    elif not os.path.isdir(bestbinset):
        raise RuntimeError(
            'Expected existing best-binset folder is missing after S3 checkpoint resume: '
            + str(bestbinset)
        )
    
    os.system('rm -rf Iteration_'+str(num)+'_genomes')
    os.system('mkdir '+str(assembly)+'_comparison_files')
    os.system('mv Bins_marker_lineage_completeness_contamination_iteration_* All_possible_bin_sets_iteration_* Best_hit_bin_sets_iteration_* Extract_bins_in_iteration_* Best_bin_set_iteration_* '+str(assembly)+'_comparison_files')

    # try:
    #     f_s3=open('S3_checkpoint.txt','a')
    # except:
    #     xyztt=0
    # f_s3.write('3rd bin selection within multiple groups done!')
    print('Done!')
    return bestbinset

if __name__ == '__main__': 
    # bins_folders_name_list=['1_adh_dn1_contigs.fasta_0.3_maxbin2_genomes', '1_adh_dn1_contigs.fasta_0.5_maxbin2_genomes', '1_adh_dn1_contigs.fasta_0.7_maxbin2_genomes', '1_adh_dn1_contigs.fasta_0.9_maxbin2_genomes', '1_adh_dn1_contigs.fasta_200_metabat_genomes', '1_adh_dn1_contigs.fasta_300_metabat_genomes', '1_adh_dn1_contigs.fasta_400_metabat_genomes', '1_adh_dn1_contigs.fasta_500_metabat_genomes', '1_adh_dn1_contigs.fasta_300_concoct_genomes', '1_adh_dn1_contigs.fasta_400_concoct_genomes']
    bins_folders_name_list=['1_assembly_sample1.fa_400_metabat_genomes','1_assembly_sample1.fa_500_metabat_genomes','1_assembly_sample1.fa_300_metabat_genomes','1_assembly_sample1.fa_200_metabat_genomes',
                            '1_assembly_sample1.fa_100_semibin_genomes','1_assembly_sample1.fa_1_SingleContig_genomes','1_assembly_sample1.fa_100_concoct_genomes',
                            '1_assembly_sample1.fa_0.3_maxbin2_genomes','1_assembly_sample1.fa_0.5_maxbin2_genomes','1_assembly_sample1.fa_0.7_maxbin2_genomes','1_assembly_sample1.fa_0.9_maxbin2_genomes']
    # depth_file='assembly.depth.txt'
    assembly='1_assembly_sample1.fa'
    # PE_connections_file='condensed.cytoscape.connections.tab'

    # binsabundance_pe_connections(bins_folders_name_list, depth_file, PE_connections_file)
    bins_comparator_multiple_groups(bins_folders_name_list, assembly)
