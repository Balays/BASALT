#!/usr/bin/env python

"""
Step S9 (CheckM2 branch): Short-read-based reassembly of bins.

This module performs reassembly of bins using short-read data to
improve contig structure after previous refinement steps.
"""

from Bio import SeqIO
import os, copy
from multiprocessing import Pool
import shutil
import subprocess

from qc_utils import run_checkm2_predict
from runtime_utils import (
    bounded_parallel_resources,
    build_bowtie2_index,
    ensure_min_free_gb,
    split_paired_sam_by_bin,
)


def mod_bin(binset_folder):
    """
    Create a re-indexed version of a binset with sequential bin IDs
    (CheckM2 branch).

    Parameters
    ----------
    binset_folder : str
        Folder containing the original per-bin FASTA and CheckM2
        ``quality_report.tsv`` output.

    Returns
    -------
    str
        Name of the new binset folder (``<binset_folder>_mod``).
    str
        Path to the combined ``Total_bins.fa`` in the working directory.
    list of str
        List of new bin IDs (``bin1``, ``bin2``, ...).
    dict
        Mapping ``bin_id -> quality_metrics`` parsed from the quality report.
    """
    pwd=os.getcwd()
    bins_checkm={}
    try:
        os.mkdir(str(binset_folder)+'_mod')
    except:
        print(str(binset_folder)+'_mod is exist. Re-create the folder')
        os.system('rm -rf '+str(binset_folder)+'_mod')
        os.mkdir(str(binset_folder)+'_mod')
    
    os.chdir(pwd+'/'+binset_folder)
    f=open('Mod_contig_id.txt','w')
    f1=open('Bin_name_mod.txt','w')
    f2=open('Bin_name_mod_quality_report.tsv','w')
    f2.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    n, record_seq, mod_bin_list, mod_bin_dict, bin_contig_num = 0, {}, [], {}, {}
    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            # print(hz
            if 'fa' in hz:
                n+=1
                m=0
                record_seq['bin'+str(n)]={}
                mod_bin_list.append('bin'+str(n))
                checkm_name_list=str(file).split('.')
                checkm_name_list.remove(checkm_name_list[-1])
                checkm_name='.'.join(checkm_name_list)
                mod_bin_dict[checkm_name]='bin'+str(n)
                f1.write(str(file)+'\t'+'bin'+str(n)+'.fa'+'\n')
                for record in SeqIO.parse(file, 'fasta'):
                    m+=1
                    f.write('bin'+str(n)+'_'+str(m)+'\t'+str(record.id)+'\n')
                    record_seq['bin'+str(n)]['bin'+str(n)+'_'+str(m)]=str(record.seq)
                bin_contig_num['bin'+str(n)]=m

    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            if 'quality_report.tsv' in str(file) and str(file) != 'Bin_name_mod_quality_report.tsv':
                n=0
                for line in open(str(file),'r'):
                    n+=1
                    if n >= 2:
                        checkm_name=str(line).strip().split('\t')[0]
                        mod_bin_checkm_name=str(mod_bin_dict[checkm_name])
                        bins_checkm[mod_bin_checkm_name]={}

                        try:
                            genome_size=str(line).strip().split('\t')[8].strip()
                            completeness=str(line).strip().split('\t')[1].strip()
                            contamination=str(line).strip().split('\t')[2].strip()
                            N50=str(line).strip().split('\t')[6].strip()
                        except:
                            genome_size=str(line).strip().split('\t')[1].strip()
                            completeness=str(line).strip().split('\t')[2].strip()
                            contamination=str(line).strip().split('\t')[3].strip()
                            N50=str(line).strip().split('\t')[4].strip()

                        bins_checkm[mod_bin_checkm_name]['N50']=int(N50)
                        bins_checkm[mod_bin_checkm_name]['Completeness']=float(completeness)
                        bins_checkm[mod_bin_checkm_name]['Genome size']=int(eval(genome_size))
                        bins_checkm[mod_bin_checkm_name]['Contamination']=float(contamination)
                        f2.write(mod_bin_checkm_name+'\t'+str(genome_size)+'\t'+str(completeness)+'\t'+str(contamination)+'\t'+str(N50)+'\n')
    f.close()
    f1.close()
    f2.close()

    os.system('mv Mod_contig_id.txt Bin_name_mod.txt Bin_name_mod_quality_report.tsv '+pwd+'/'+str(binset_folder)+'_mod')
    
    os.chdir(pwd+'/'+str(binset_folder)+'_mod')
    f1=open('Total_bins.fa','w')
    for bin_id in record_seq.keys():
        f=open(str(bin_id)+'.fa','w')
        for contigs in record_seq[bin_id].keys():
            f.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
            f1.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
        f.close()
    f1.close()
    os.system('mv Total_bins.fa '+pwd)
    os.chdir(pwd)
    return str(binset_folder)+'_mod', 'Total_bins.fa', mod_bin_list, bins_checkm


def parse_sam(sam_file, fq, pair, n):
    """
    Parse Bowtie2 SAM output and split paired reads by bin (CheckM2 branch).

    Parameters
    ----------
    sam_file : str
        Path to the SAM file produced by Bowtie2.
    fq : dict
        Mapping ``bin_id -> dict(read_id -> int)`` tracking how many
        mates of a read pair have been written.
    pair : dict
        Temporary mapping ``bin_id -> dict(read_id -> int)`` used to
        count paired hits before writing FASTQ records.
    n : int
        Dataset index, used to prefix read IDs.

    Returns
    -------
    None
        Writes bin-specific paired-end FASTQ files and a summary file
        ``Bin_reads_summary.txt`` in the working directory.
    """
    print('Parsing', sam_file)
    pair_counts, unmatched = split_paired_sam_by_bin(sam_file, n)
    with open('Bin_reads_summary.txt', 'w') as f_summary:
        for bin_id in sorted(pair_counts):
            f_summary.write(str(bin_id)+' SEQ number:'+str(pair_counts[bin_id])+'\n')
    if unmatched:
        print('Ignored', unmatched, 'read/bin assignments without both mates')


def parse_lr_sam(sam_file, long_read, sn):
    """
    Group long reads by bin based on SAM alignments (CheckM2 branch).

    Parameters
    ----------
    sam_file : str
        Path to the SAM file produced by mapping long reads to bins.
    long_read : str
        Path to the original long-read FASTQ file.
    sn : int
        Index of the long-read file, used in output naming.

    Returns
    -------
    dict
        Mapping ``bin_id -> dict(long_read_id -> None)`` indicating which
        long reads map to each bin.
    """
    print('Reading long reads id '+str(long_read))
    # f_not_mapped_reads=open('Not_mapped_reads.txt','w')
    bin_lr, bin_lr2, lr_bin, lr_bin2 = {}, {}, {}, {}
    m, m1, m2 = 0, 0, 0
    for line in open(sam_file,'r'):
        m1+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            read_id=flist[0]
            bin_id=flist[2].split('_')[0]
            bin_lr[bin_id]={}
            lr_bin[read_id]={}

        if m1 % 1000000 == 0:
            print('Read', m1,'lines')

    print('Collecting long reads '+str(long_read))
    for line in open(sam_file,'r'):
        m2+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id=flist[2].split('_')[0]
            seqs=flist[9]
            if 'bin' in bin_id and len(seqs) >= 100:
                read_id=flist[0]
                bin_lr[bin_id][read_id]=''
                lr_bin[read_id][bin_id]=''

        if m2 % 1000000 == 0:
            print('Read', m2,'lines')

    bin_lr2=copy.deepcopy(bin_lr)
    lr_bin2=copy.deepcopy(lr_bin)

    for item in bin_lr2.keys():
        if len(bin_lr2[item]) == 0:
            del bin_lr[item]

    for item in lr_bin2.keys():
        if len(lr_bin2[item]) == 0:
            del lr_bin[item]

    f=open('Bin_long_read'+str(sn)+'.txt','w')
    for bin_id in bin_lr.keys():
        f.write(str(bin_id)+'\t'+str(bin_lr[bin_id])+'\n')
        try:
            f1=open(str(bin_id)+'_lr.fq','a')
        except:
            f1=open(str(bin_id)+'_lr.fq','w')
        f1.close()
    f.close()

    f=open('Long_read_bin'+str(sn)+'.txt','w')
    for lr in lr_bin.keys():
        f.write(str(lr)+'\t'+str(lr_bin[lr])+'\n')
    f.close()

    n, record_bin_line, record_bin_line2 = 0, {}, {}
    for line in open(long_read,'r'):
        n+=1
        record_bin_line[n]=[]

    print('Splitting long reads to different bins '+str(long_read))
    n = 0
    for line in open(long_read,'r'):
        n+=1
        m=n-1
        if m % 4 == 0:
            seq_id=str(line).strip().split(' ')[0].split('@')[1]
            if seq_id in lr_bin.keys():
                for bin_id in lr_bin[seq_id].keys():
                    record_bin_line[n].append(bin_id)
                    record_bin_line[n+1].append(bin_id)
                    record_bin_line[n+2].append(bin_id)
                    record_bin_line[n+3].append(bin_id)

    record_bin_line2=copy.deepcopy(record_bin_line)
    for line in record_bin_line2.keys():
        if len(record_bin_line2[line]) == 0:
            del record_bin_line[line]
    seq_num=int(len(record_bin_line)/4)
    print(str(seq_num)+' reads from '+str(long_read)+' will be splitted into different bins')

    n1=0
    for line in open(long_read,'r'):
        n1+=1
        if n1 in record_bin_line.keys():
            for bin_id in record_bin_line[n1]:
                f1=open(str(bin_id)+'_lr.fq','a')
                f1.write(line)
                f1.close()
    if n1 % 1000000 == 0:
        print('Parse', n1,'lines')
    print('Long reads '+str(long_read)+' splitting done!')
    return bin_lr


def mapping_sr(total_fa, datasets_list, fq, pair, num_threads):
    """
    Map short reads to bins and prepare bin-specific FASTQ files
    for reassembly (CheckM2 branch).

    Parameters
    ----------
    total_fa : str
        Path to the concatenated bin FASTA file used as mapping reference.
    datasets_list : dict
        Mapping ``sample_id -> [r1_fastq, r2_fastq]`` with paired-end reads.
    fq : dict
        Mapping ``bin_id -> dict(read_id -> int)`` used to track mates.
    pair : dict
        Temporary mapping used when counting paired hits.
    num_threads : int
        Number of Bowtie2 threads to use.

    Returns
    -------
    None
        Writes bin-specific FASTQ files and summary statistics to disk.
    """
    build_bowtie2_index(total_fa, num_threads)
    n = 0
    for item in datasets_list.keys():
        n+=1
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(datasets_list[item][0])+' -2 '+str(datasets_list[item][1])+' -S '+str(item)+'.sam -q --no-unal')
        parse_sam(str(item)+'.sam', fq, pair, n)
        os.system('rm '+str(item)+'.sam')

def parse_checkm(checkm_containing_folder,pwd):
    #pwd=os.getcwd()
    bins_checkm={}

    os.chdir(pwd+'/'+checkm_containing_folder)
    for root, dirs, files in os.walk(pwd+'/'+checkm_containing_folder):
        for file in files:        
            if 'quality_report.tsv' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2:
                        genome_ids=str(line).strip().split('\t')[0]
                        if '_genomes.0' not in genome_ids:
                            bins_checkm[genome_ids]={}
                            try:
                                genome_size=str(line).strip().split('\t')[8].strip()
                                completeness=str(line).strip().split('\t')[1].strip()
                                contamination=str(line).strip().split('\t')[2].strip()
                                N50=str(line).strip().split('\t')[6].strip()
                            except:
                                genome_size=str(line).strip().split('\t')[1].strip()
                                completeness=str(line).strip().split('\t')[2].strip()
                                contamination=str(line).strip().split('\t')[3].strip()
                                N50=str(line).strip().split('\t')[4].strip()

                            bins_checkm[genome_ids]['N50']=int(N50)
                            bins_checkm[genome_ids]['Completeness']=float(completeness)
                            bins_checkm[genome_ids]['Genome size']=int(eval(genome_size))
                            bins_checkm[genome_ids]['Contamination']=float(contamination)
    return bins_checkm


def _run_external(command):
    """Run a per-bin assembler without hiding an unavailable executable."""
    try:
        return subprocess.run(command, check=False).returncode
    except OSError as error:
        print('Could not run {}: {}'.format(command[0], error))
        return 127


def _filter_contigs(source, destination, minimum_length=1000):
    """Write non-empty filtered assembly output and return its contig count."""
    count = 0
    with open(destination, 'w') as output:
        for record in SeqIO.parse(source, 'fasta'):
            if len(record.seq) >= minimum_length:
                count += 1
                output.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    if count == 0:
        os.remove(destination)
    return count


def _short_reassembly_worker(job):
    item, reads, output_folder, sequence_folder, threads, ram, pwd = job
    r1 = reads[0] if os.path.isabs(reads[0]) else os.path.join(pwd, reads[0])
    r2 = reads[1] if os.path.isabs(reads[1]) else os.path.join(pwd, reads[1])
    output_folder = os.path.abspath(os.path.join(pwd, output_folder))
    sequence_folder = os.path.abspath(os.path.join(pwd, sequence_folder))
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(sequence_folder, exist_ok=True)
    ensure_min_free_gb(pwd, os.environ.get('BASALT_MIN_FREE_GB', '0'),
                       context='short-read reassembly')

    spades_work = os.path.join(pwd, str(item)+'_spades_reassembly')
    spades_final = os.path.join(
        output_folder, str(item)+'_SPAdes_re-assembly_contigs.fa'
    )
    if not os.path.isfile(spades_final) or os.path.getsize(spades_final) == 0:
        command = [
            'spades.py', '-1', r1, '-2', r2, '-o', spades_work,
            '--careful', '-t', str(threads), '-m', str(ram),
        ]
        _run_external(command)
        contigs = os.path.join(spades_work, 'contigs.fasta')
        if not os.path.isfile(contigs):
            shutil.rmtree(spades_work, ignore_errors=True)
            command.remove('--careful')
            _run_external(command)
        if os.path.isfile(contigs):
            temporary = spades_final+'.tmp'
            if _filter_contigs(contigs, temporary) >= 2:
                os.replace(temporary, spades_final)
            elif os.path.exists(temporary):
                os.remove(temporary)
    if os.path.isfile(spades_final) and os.path.getsize(spades_final) > 0:
        shutil.rmtree(spades_work, ignore_errors=True)

    idba_work = os.path.join(pwd, str(item)+'_idba_reassembly')
    idba_input = os.path.join(pwd, str(item)+'_idba.fa')
    idba_final = os.path.join(output_folder, str(item)+'_IDBA_re-assembly_contigs.fa')
    if not os.path.isfile(idba_final) or os.path.getsize(idba_final) == 0:
        if _run_external(['fq2fa', '--merge', '--filter', r1, r2, idba_input]) == 0:
            first = next(SeqIO.parse(idba_input, 'fasta'), None)
            mode = '-r' if first is not None and len(first.seq) <= 125 else '-l'
            _run_external([
                'idba_ud', mode, idba_input, '-o', idba_work,
                '--num_threads', str(threads), '--min_contig', '1000',
            ])
            idba_contigs = os.path.join(idba_work, 'contig.fa')
            if os.path.isfile(idba_contigs) and os.path.getsize(idba_contigs) > 0:
                os.replace(idba_contigs, idba_final)
    if os.path.isfile(idba_final) and os.path.getsize(idba_final) > 0:
        shutil.rmtree(idba_work, ignore_errors=True)
        if os.path.exists(idba_input):
            os.remove(idba_input)

    succeeded = any(
        os.path.isfile(path) and os.path.getsize(path) > 0
        for path in (spades_final, idba_final)
    )
    if succeeded:
        for read_path in (r1, r2):
            if os.path.isfile(read_path):
                shutil.move(read_path, os.path.join(sequence_folder, os.path.basename(read_path)))
    return item, succeeded


def reassembly(bin_seq, reassembly_bin_folder, num_threads, bins_seq_folder,
               long_read, ram, pwd):
    """Run bounded parallel SPAdes/IDBA reassembly with divided resources."""
    jobs, threads_per_job, ram_per_job = bounded_parallel_resources(
        num_threads, ram, len(bin_seq)
    )
    if jobs == 0:
        return
    work = [
        (item, bin_seq[item], reassembly_bin_folder, bins_seq_folder,
         threads_per_job, ram_per_job, pwd)
        for item in bin_seq
    ]
    print('Short-read reassembly jobs:', jobs,
          'threads/job:', threads_per_job, 'RAM/job:', ram_per_job, 'GB')
    if jobs == 1:
        for job in work:
            _short_reassembly_worker(job)
    else:
        with Pool(processes=jobs) as pool:
            pool.map(_short_reassembly_worker, work)

def unicycler_mul(item, pwd, sr_folder, bin_seq, bin_lr, threads,
                  reassembly_bin_folder, bins_seq_folder):
    """Run one restartable Unicycler job and remove only its scratch folder."""
    output_folder = os.path.join(pwd, reassembly_bin_folder)
    final_output = os.path.join(
        output_folder, str(item)+'_UNICYCLER_re-assembly_contigs.fa'
    )
    long_read = os.path.join(pwd, str(bin_lr[item]))
    long_read_archive = os.path.join(pwd, bins_seq_folder, os.path.basename(long_read))
    if os.path.isfile(final_output) and os.path.getsize(final_output) > 0:
        return item, True

    ensure_min_free_gb(
        pwd,
        float(os.environ.get('BASALT_MIN_FREE_GB', '0')),
        context='Unicycler reassembly for '+str(item),
    )
    work_folder = os.path.join(pwd, str(item)+'_unicycler_reassembly')
    shutil.rmtree(work_folder, ignore_errors=True)
    command = [
        'unicycler',
        '-1', os.path.join(pwd, sr_folder, str(bin_seq[item][0])),
        '-2', os.path.join(pwd, sr_folder, str(bin_seq[item][1])),
        '-l', long_read,
        '-o', work_folder,
        '-t', str(threads),
        '--mode', 'conservative',
        '--no_pilon',
    ]
    result = subprocess.run(command, check=False)
    assembly = os.path.join(work_folder, 'assembly.fasta')
    succeeded = (
        result.returncode == 0 and os.path.isfile(assembly)
        and os.path.getsize(assembly) > 0
    )
    if succeeded:
        os.replace(assembly, final_output)
    else:
        print('Unicycler did not produce an assembly for '+str(item))
    shutil.rmtree(work_folder, ignore_errors=True)
    if os.path.isfile(long_read):
        os.makedirs(os.path.dirname(long_read_archive), exist_ok=True)
        shutil.move(long_read, long_read_archive)
    return item, succeeded

def reassembly_lr(bin_seq, bin_lr, reassembly_bin_folder, num_threads,
                  bins_seq_folder, sr_folder, ram, pwd):
    """
    Perform hybrid reassembly of bins using both short and long reads
    (CheckM2 branch).

    Parameters
    ----------
    bin_seq : dict
        Mapping ``bin_id -> [r1_fastq, r2_fastq]`` with bin-specific
        short reads.
    bin_lr : dict
        Mapping ``bin_id -> long_read_fastq`` with bin-specific long reads.
    reassembly_bin_folder : str
        Folder where hybrid reassembled bin FASTA files will be written.
    num_threads : int
        Number of threads to use for the hybrid assembler.
    bins_seq_folder : str
        Folder holding bin-specific short-read FASTQ files.
    sr_folder : str
        Folder used for intermediate short-read assemblies.
    ram : int
        Maximum RAM (in GB) available to assemblers.
    pwd : str
        Working directory path.

    Returns
    -------
    None
        Writes hybrid reassembled contig FASTA files into
        ``reassembly_bin_folder``.
    """
    items = [item for item in bin_lr if item in bin_seq]
    jobs, threads_per_job, _ = bounded_parallel_resources(
        num_threads, ram, len(items), target_threads=24, target_ram=48,
    )
    if jobs == 0:
        return
    work = [
        (
            item, pwd, sr_folder, bin_seq, bin_lr, threads_per_job,
            reassembly_bin_folder, bins_seq_folder,
        )
        for item in items
    ]
    print(
        'Unicycler jobs:', str(jobs),
        'threads/job:', str(threads_per_job),
    )
    if jobs == 1:
        for job in work:
            unicycler_mul(*job)
    else:
        with Pool(processes=jobs) as pool:
            pool.starmap(unicycler_mul, work)

def bin_comparison(paired_bins, bin_checkm):
    pwd=os.getcwd()
    f=open('Reassembled_bins_comparison.txt','w')
    best_bin, best_bin_checkm={}, {}
    for item in paired_bins.keys():
        best_bin_checkm_name_list=item.split('.')
        best_bin_checkm_name_list.remove(best_bin_checkm_name_list[-1])
        best_bin_checkm_name='.'.join(best_bin_checkm_name_list)
        f.write(str(item)+'\t'+str(bin_checkm[best_bin_checkm_name])+'\n')
        for item2 in paired_bins[item]:
            reass_bin_checkm_name_list=item2.split('.')
            reass_bin_checkm_name_list.remove(reass_bin_checkm_name_list[-1])
            reass_bin_checkm_name='.'.join(reass_bin_checkm_name_list)
            f.write(str(item2)+'\t'+str(bin_checkm[reass_bin_checkm_name])+'\n')
            best_bin_cpn=bin_checkm[best_bin_checkm_name]['Completeness']
            best_bin_ctn=bin_checkm[best_bin_checkm_name]['Contamination']
            best_bin_ml=bin_checkm[best_bin_checkm_name]['N50']
            reass_bin_cpn=bin_checkm[reass_bin_checkm_name]['Completeness']
            reass_bin_ctn=bin_checkm[reass_bin_checkm_name]['Contamination']
            reass_bin_ml=bin_checkm[reass_bin_checkm_name]['N50']
        
            delta_cpn_ctn_bestbin=float(best_bin_cpn)-5*float(best_bin_ctn)
            delta_cpn_ctn_reass_bin=float(reass_bin_cpn)-float(5*reass_bin_ctn)

            if '_SPAdes_' in reass_bin_checkm_name or '_UNICYCLER_' in reass_bin_checkm_name:
                if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
                    best_bin_checkm_name=best_bin_checkm_name
                elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
                    best_bin_checkm_name=reass_bin_checkm_name
                elif delta_cpn_ctn_bestbin == delta_cpn_ctn_reass_bin:
                    if reass_bin_ml > best_bin_ml:
                        best_bin_checkm_name=reass_bin_checkm_name
                    else:
                        best_bin_checkm_name=best_bin_checkm_name
                else:
                    continue
            elif '_IDBA_' in reass_bin_checkm_name:
                delta_idba = delta_cpn_ctn_reass_bin-delta_cpn_ctn_bestbin
                if delta_idba < 3:
                    best_bin_checkm_name=best_bin_checkm_name
                elif delta_idba >= 3:
                    best_bin_checkm_name=reass_bin_checkm_name
                else:
                    continue
                # if best_bin_ml > reass_bin_ml:
                #    best_bin_checkm_name=best_bin_checkm_name
                # elif best_bin_ml < reass_bin_ml:
                #    best_bin_checkm_name=reass_bin_checkm_name
                # else:
                #    best_bin_checkm_name=reass_bin_checkm_name
            # elif '_UNICYCLER_' in reass_bin_checkm_name:
            #     if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=best_bin_checkm_name
            #     elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=reass_bin_checkm_name
            #     elif delta_cpn_ctn_bestbin == delta_cpn_ctn_reass_bin:
            #         if reass_bin_ml > best_bin_ml:
            #             best_bin_checkm_name=reass_bin_checkm_name
            #         else:
            #             best_bin_checkm_name=best_bin_checkm_name
            #     else:
            #         continue
            # elif '_FLYe_' in reass_bin_checkm_name:
            #     if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=best_bin_checkm_name
            #     elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=reass_bin_checkm_name
            #     else:
            #         continue

        best_bin[best_bin_checkm_name+'.fa']=best_bin_checkm_name
        best_bin_checkm[best_bin_checkm_name]=bin_checkm[best_bin_checkm_name].copy()
    f.close()
    return best_bin, best_bin_checkm

def re_assembly_main(binset_folder, datasets_list, long_read,
                     hybri_reassembly, ram, num_threads):
    """
    Entry point for S9 short-read reassembly (CheckM2 branch).

    Parameters
    ----------
    binset_folder : str
        Folder containing the starting binset.
    datasets_list : dict
        Mapping ``sample_id -> [r1_fastq, r2_fastq]`` with paired-end reads.
    long_read : list of str
        List of long-read FASTQ files used for optional hybrid reassembly.
    hybri_reassembly : {'y', 'n'}
        Flag controlling whether hybrid reassembly should be performed.
    ram : int
        Maximum RAM (in GB) available to assemblers.
    num_threads : int
        Number of threads to use for mapping and assembly.

    Returns
    -------
    None
        Coordinates mapping, short-read reassembly, optional long-read
        reassembly, and CheckM-based bin selection; results are written
        to ``<binset_folder>_re-assembly_binset`` and related folders.
    """
    pwd=os.getcwd()
    try:
        f=open('Assembly_status.txt','a')
    except:
        f=open('Assembly_status.txt','w')
    f.close()

    try:
        os.mkdir(binset_folder+'_re-assembly_binset')
    except:
        print(binset_folder+'_re-assembly_binset exists')
    
    assembled_bins={}
    try:
        os.mkdir(binset_folder+'_sr_bins_seq')
    except:
        print(binset_folder+'_sr_bins_seq exists')
        for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_sr_bins_seq'):
            for file in files:        
                if '.fq' in file and '_seq_' in file and 'bin' in file:
                    assembled_bins[str(file).split('_seq_')[0].strip()]=0

    A=mod_bin(binset_folder)
    mod_bin_folder=A[0]
    total_fa=A[1]
    mod_bin_list=A[2]
    original_bins_checkm=A[3]
    fq, pair, bin_seq, bin_checkm={}, {}, {}, {}
    bin_checkm=original_bins_checkm.copy()
    for bin_id in mod_bin_list:
        fq[str(bin_id)]={}
        pair[str(bin_id)]={}
        bin_seq[str(bin_id)]=[]
        bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R1.fq')
        bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R2.fq')
        if len(assembled_bins) == 0:
            f1=open(str(bin_id)+'_seq_R1.fq','w')
            f2=open(str(bin_id)+'_seq_R2.fq','w')
            f1.close()
            f2.close()

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Short-read mapping done!' in line:
            x=1

    if x == 0:
        f_not_mapped_reads=open('Not_mapped_reads.txt','w')
        f_not_mapped_reads.close()

        mapping_sr(total_fa, datasets_list, fq, pair, num_threads)
        f=open('Assembly_status.txt','a')
        f.write('Short-read mapping done!'+'\n')
        f.close()

    bin_seq2={}
    for item in bin_seq:
        if item not in assembled_bins.keys():
            bin_seq2[item]=bin_seq[item]

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Short-read assembly done!' in line:
            x=1

    if x == 0:
        reassembly(bin_seq2, binset_folder+'_re-assembly_binset', num_threads, binset_folder+'_sr_bins_seq', long_read, ram, pwd)
        f=open('Assembly_status.txt','a')
        f.write('Short-read assembly done!'+'\n')
        f.close()

    if len(long_read) != 0 and hybri_reassembly == 'y':
        print('Reassemblying long reads')
        assembled_bins={}
        try:
            os.mkdir(binset_folder+'_lr_bins_seq')
        except:
            print(binset_folder+'_lr_bins_seq exists')
            for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_lr_bins_seq'):
                for file in files:        
                    if '_lr.fq' in file and 'bin' in file:
                        assembled_bins[str(file).split('_lr')[0].strip()]=0
        os.chdir(pwd)

        if len(assembled_bins) == 0:
            print('Mapping long reads')
            x=0
            for line in open('Assembly_status.txt','r'):
                if 'Long-read mapping done!' in line:
                    x=1

            if x == 0:
                n, bin_lr=0, {}
                for lrs in long_read:
                    n+=1
                    os.system('minimap2 -t '+str(num_threads)+' -ax map-ont Total_bins.fa '+str(lrs)+' > lr'+str(n)+'.sam')
                    print('Splitting long reads '+str(n))
                    bin_lr.update(parse_lr_sam('lr'+str(n)+'.sam', lrs, n))

                f=open('Assembly_status.txt','a')
                f.write('Long-read mapping done!'+'\n')
                f.close()

                for i in range(1, n+1):
                    os.system('rm lr'+str(n)+'.sam')

        x=0
        for line in open('Assembly_status.txt','r'):
            if 'Long-read assembly done!' in line:
                x=1

        if x == 0:
            n, bin_lr=0, {}
            for lrs in long_read:
                n+=1
                for line in open('Bin_long_read'+str(n)+'.txt','r'):
                    bin_id=str(line).strip().split('\t')[0]
                    bin_lr[bin_id]=0
 
            bin_lr2={}
            for item in bin_lr.keys():
                if item not in assembled_bins.keys():
                    bin_lr2[item]=str(item)+'_lr.fq'
            
            reassembly_lr(bin_seq, bin_lr2, binset_folder+'_re-assembly_binset', num_threads, binset_folder+'_lr_bins_seq', binset_folder+'_sr_bins_seq', ram, pwd)
            f=open('Assembly_status.txt','a')
            f.write('Long-read assembly done!'+'\n')
            f.close()

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Assembly quality evaluation done!' in line:
            x=1
    
    if x == 0:
        print('Checking reassembled bins')
        os.chdir(str(binset_folder)+'_re-assembly_binset')
        for root, dirs, files in os.walk(pwd+'/'+str(binset_folder)+'_re-assembly_binset'):
            for file in files:
                try:
                    ff=open(file+'_filtrated.fa','w')
                    for record in SeqIO.parse(file, 'fasta'):
                        if len(record.seq) != 0:
                            ff.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    ff.close()
                    os.system('mv '+str(file)+'_filtrated.fa '+str(file))
                except:
                    xyzzzz=0
        os.chdir(pwd)
        
        # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset_folder)+'_re-assembly_binset '+str(binset_folder)+'_re-assembly_binset_checkm')
        run_checkm2_predict(
            str(binset_folder)+'_re-assembly_binset',
            'fa',
            str(binset_folder)+'_re-assembly_binset_checkm',
            num_threads,
        )
        f=open('Assembly_status.txt','a')
        f.write('Assembly quality evaluation done!'+'\n')
        f.close()

    os.chdir(pwd)
    reassembly_bin_checkm=parse_checkm(binset_folder+'_re-assembly_binset_checkm', pwd)
    bin_checkm.update(reassembly_bin_checkm)

    # checkm_t=open('bin_checkm_t.txt','w')
    # for item in bin_checkm.keys():
    #     checkm_t.write(str(item)+'\t'+str(bin_checkm[item])+'\n')
    # checkm_t.close()

    paired_bins, best_bin, best_bin_checkm={}, {}, {}
    os.chdir(pwd+'/'+binset_folder+'_re-assembly_binset')
    for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_re-assembly_binset'):
        for file in files:
            if '_SPAdes_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_SPAdes_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
            elif '_IDBA_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_IDBA_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
            elif '_UNICYCLER_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_UNICYCLER_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
            elif '_FLYe_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_FLYe_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
    os.chdir(pwd)

    # paired_bin_t=open('Paired_bins_t.txt','w')
    # for item in paired_bins.keys():
    #     paired_bin_t.write(str(item)+'\t'+str(paired_bins[item])+'\n')
    # paired_bin_t.close()

    best_bin_after_c=bin_comparison(paired_bins, bin_checkm)
    best_bin=best_bin_after_c[0].copy()
    best_bin_checkm=best_bin_after_c[1].copy()

    try:
        os.mkdir(binset_folder+'_re-assembly')
    except:
        print(binset_folder+'_re-assembly exists')
    
    selected_bins={}
    for item in best_bin.keys():
        if '_re-assembly_contigs.fa' in item:
            s_name=item.split('_')[0]+'.fa'
            selected_bins[s_name]=1
            os.system('cp '+pwd+'/'+binset_folder+'_re-assembly_binset/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')
        else:
            selected_bins[item]=1
            os.system('cp '+pwd+'/'+str(binset_folder)+'_mod/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')

    os.chdir(pwd+'/'+str(binset_folder)+'_mod')
    for root, dirs, files in os.walk(pwd+'/'+str(binset_folder)+'_mod'):
        for file in files:
            if '.fa' in file:
                if file not in selected_bins.keys():
                    os.system('cp '+pwd+'/'+str(binset_folder)+'_mod/'+file+' '+pwd+'/'+binset_folder+'_re-assembly')
                    item_checkm_name=file.split('.fa')[0]
                    best_bin_checkm[item_checkm_name]=bin_checkm[item_checkm_name]
    os.chdir(pwd)

    os.chdir(pwd+'/'+binset_folder+'_re-assembly')
    f=open('Best_binset_after_re-assembly_quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    for item in best_bin_checkm.keys():
        # f.write(str(item)+'\t'+str(best_bin_checkm[item])+'\n')
        f.write(str(item)+'\t'+str(best_bin_checkm[item]['Genome size'])+'\t'+str(best_bin_checkm[item]['Completeness'])+'\t'+str(best_bin_checkm[item]['Contamination'])+'\t'+str(best_bin_checkm[item]['N50'])+'\n')
    f.close()
    os.chdir(pwd)

    # try:
    #     os.system('rm -rf '+binset_folder+'_sr_bins_seq')
    # except:
    #     x=1

    # try:
    #     os.system('rm -rf '+binset_folder+'_lr_bins_seq')
    # except:
    #     x=1
    print('Re-assembly done!')

if __name__ == '__main__': 
    binset_folder='1_Opera_unpolished_cat_contigs.fasta_BestBinsSet_outlier_refined_filtrated_retrieved'
    datasets_list={'1':['PE_r1_RH_S001_insert_270_mate1.fq','PE_r2_RH_S001_insert_270_mate2.fq'], '2':['PE_r1_RH_S002_insert_270_mate1.fq','PE_r2_RH_S002_insert_270_mate2.fq'], '3':['PE_r1_RH_S003_insert_270_mate1.fq','PE_r2_RH_S003_insert_270_mate2.fq'],'4':['PE_r1_RH_S004_insert_270_mate1.fq','PE_r2_RH_S004_insert_270_mate2.fq'],'5':['PE_r1_RH_S005_insert_270_mate1.fq','PE_r2_RH_S005_insert_270_mate2.fq']}
    long_read=['anonymous_reads1.fq','anonymous_reads2.fq','anonymous_reads3.fq','anonymous_reads4.fq','anonymous_reads5.fq'] ### Write the name of long read here, if there is not, just let it to be blank
    hybri_reassembly='n' ### Use Unicycler to re-assembly. e.g. --hybrid y / --hybrid n; defalt no
    num_threads=20
    ram=250
    re_assembly_main(binset_folder, datasets_list, long_read, hybri_reassembly, ram, num_threads)
