import os
import shutil
import pandas as pd

class foldseek:
    """
    Class to hold methods to work with foldseek executable.

    Methods
    -------
    runEasySearch()
    parseFoldSeekOutput()

    ## TODO ##
    More options to manage databases (create custom db...)
    More search parameters
    More options to manage output of search (customized dataframe)
    Foldseek search option (db against db?)
    Use TM score (check github tutorial)
    ##

    """
    def runEasySearch(target_pdb,output_folder,db,format_output=None,sensitivity=9.5,exhaustive_search=False,max_seqs=1000,
                    e_value_threshold=0.001,alignment_type='Smith–Waterman',coverage_threshold=0,coverage_mode='query',
                    output_file='fs_output.txt',local=False,overwrite=False):
        """
        Run structural search of target against specified database.

        Parameters
        ----------
        target_pdb : str
            Path to the target pdb or pdb folder to perform the structural search.
        output_folder : str
            Folder to run foldseek.
        db: str
            Absolute path to the search database. Available databases at bubbles:
                /home/bubbles/foldseek_databases/af50 (AlphaFold UniProt Protein Structure Database clustered with MMseqs2 at 50% sequence identity and 90% bidrectional coverage.)
                /home/bubbles/foldseek_databases/afsp (AlphaFold Swissprot Protein Structure Database.)
                /home/bubbles/foldseek_databases/pdb (Protein Data Bank.)
        format_output : list
            List with the data to write in the columns of the output. The possible arguments are as follows:
            query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen
            tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,mismatch,qcov,tcov
            qset,qsetid,tset,tsetid,taxid,taxname,taxlineage,
            lddt,lddtfull,qca,tca,t,u,qtmscore,ttmscore,alntmscore,rmsd,prob
        sensitivity : float
            Adjust sensitivity to speed trade-off; lower is faster, higher more sensitive.(fast: 7.5, default: 9.5)
        exhaustive_search : bool
            Skips prefilter and performs an all-vs-all alignment (more sensitive but much slower)
        max_seqs : int
            Adjust the amount of prefilter handed to alignment; increasing it can lead to more hits
        e_value_threshold : float
            List matches below this E-value (range 0.0-inf, default: 0.001); increasing it reports more distant structures
        alignment_type : str
            0: 3Di Gotoh-Smith-Waterman (local, not recommended), 1: TMalign (global, slow), 2: 3Di+AA Gotoh-Smith-Waterman (local, default)
        coverage_threshold : float
            List matches above this fraction of aligned (covered) residues (see --cov-mode) (default: 0.0); higher coverage = more global alignment
        coverage_mode : str
            0: coverage of query and target, 1: coverage of target, 2: coverage of query
        output_file: str
            File to write foldseek output.
        local: bool
            Run foldseek locally.
        overwrite: bool
            Overwrite results of previous run.
        """

        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        if os.path.exists(output_folder+'/'+output_file) and not overwrite:
            print('Output file already exists. Use overwrite argument to run FoldSeek again.')
        else:
            if format_output == None:
                format_output = ['query','target','fident','alnlen','mismatch','gapopen','evalue','qseq','tseq']

            if local:
                print('WARNING: This is not currently supported. Proceed with caution.')
                if not os.path.exists(db):
                    raise ValueError('The database does not exist in the specified path. Download the desired database with foldseek databases command or use one of the existing ones.')
                os.chdir(output_folder)
                os.system('foldseek easy-search ../'+target_pdb+' ../'+db+' '+output_file+' tmp')
                os.chdir('..')
            else:
                shutil.copyfile(target_pdb,output_folder+'/'+target_pdb.split('/')[-1])
                line = 'foldseek easy-search '+target_pdb.split('/')[-1]+' '+db+' '+output_file+' '
                line += '--format-output '+','.join(format_output)+' '
                line += '-s '+str(sensitivity)+' '
                if exhaustive_search:
                    line += '--exhaustive-search '
                line += '--max-seqs '+str(max_seqs)+' '
                line += '-e '+str(e_value_threshold)+' '
                if alignment_type == 'Smith–Waterman':
                    line += '--alignment-type '+str(2)+' '
                elif alignment_type == 'TMalign':
                    line += '--alignment-type '+str(1)+' '
                else:
                    raise ValueError('Alignment type not recognised.')
                line += '-c '+str(coverage_threshold)+' '
                if coverage_mode == 'query':
                    line += '--cov-mode '+str(2)+' '
                elif coverage_mode == 'target':
                    line += '--cov-mode '+str(1)+' '
                elif coverage_mode == 'both':
                    line += '--cov-mode '+str(0)+' '
                else:
                    raise ValueError('Coverage mode not recognised.')

                line += 'tmp \n'

                with open(output_folder+'/'+'command.sh','w') as f:
                    f.write(line)

                print('Output folder has been generated. Execute command.sh to run foldseek.')


    def parseFoldSeekOutput(output_folder,output_file='fs_output.txt',format_output=None):
        """
        Parse foldseek easysearch output.

        Parameters
        ----------
        output_folder : str
            Path to the folder with easysearch output.
        output_file : str
            File with the easysearch output.

        Returns
        -------
        df : pandas.DataFrame
            Dataframe with the results of FoldSeek easysearch run.
        """

        if format_output == None:
            format_output = ['query','target','fident','alnlen','mismatch','gapopen','evalue','qseq','tseq']

        lines = []
        f = open(output_folder+'/'+output_file,'r')
        lines = f.readlines()
        spl_lines = [x.split() for x in lines]

        df = pd.DataFrame(spl_lines,columns=format_output)
        df.set_index(format_output[0],inplace=True)

        return df

    def downloadDatabase(db_folder,db,db_output_name):
        """
        Download database in foldseek format.

        Parameters
        ----------
        db_folder : str
            Name of folder to contain db files.
        db : str
            Name of the database to download.
        db_output_name : str
            Name for the db files.
        """

        available_db = ['Alphafold/UniProt','Alphafold/UniProt50','Alphafold/Proteome','Alphafold/Swiss-Prot','ESMAtlas30','PDB']
        if db not in available_db:
            raise ValueError(db +'database is not currently supported please select one of the following: '+','.join(available_db))

        if not os.path.exists(db_folder):
            os.mkdir(db_folder)

        if not os.path.exists(db_folder+'/'+db_output_name):
            os.chdir(db_folder)
            os.system('foldseek databases '+db+' '+db_output_name+' tmp')
            os.chdir('..')
        else:
            print('A file with name '+db_output_name+' already exists in '+db_folder+'. Assuming database is correctly downloaded.')
