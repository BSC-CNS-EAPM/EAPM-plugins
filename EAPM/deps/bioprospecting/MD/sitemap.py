import os
import mdtraj as md
from pkg_resources import resource_stream, Requirement, resource_listdir
import io
import pandas as pd

def calculateSiteMapForTrajectory(job_folder, trajectory_files, topology_file, target_residues,
                                  site_box=10, resolution='fine', reportsize=100, frames=None,
                                  overwrite=False, pH=7.0):
    """
    Calculate SiteMap for a full trajectory. The function
    creates separate PDBs for each frame to parallelise the
    calculation

    Parameters
    ==========
    job_folder : str
        Name of the job folder
    trajectory_files : (list, dict)
        Paths to the trajectory files
    topology_file : str
        Path to the topology file
    """

    if frames != None:
        raise ValueError('Not yet implemented!')

    # Check input's format
    if isinstance(trajectory_files, str):
        trajectory_files = [trajectory_files]

    if isinstance(trajectory_files, list):
        trajectory_files_tmp = {}
        for t in trajectory_files:
            fn = t.split('/')[-1].split('.')[0]
            trajectory_files_tmp[fn] = t
        trajectory_files = trajectory_files_tmp

    if isinstance(target_residues, int):
        target_residues = [target_residues]
    if not isinstance(target_residues, (list, tuple)):
        raise ValueError('target_residues should be an integer, or a list or a tuple of integers')

    if not isinstance(trajectory_files, dict):
        raise ValueError('trajectory_files should be a list or dictionary of trajectory paths.')

    # Create job folder
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Copy script file
    _copyScriptFile(job_folder, 'prepareForSiteMap.py')

    # Iterate trajectories
    jobs = []
    for t in trajectory_files:
        tdir = job_folder+'/'+t

        # Create folder for trajectory
        if not os.path.exists(tdir):
            os.mkdir(tdir)

        # Create folder for input pdbs
        if not os.path.exists(tdir+'/input_models'):
            os.mkdir(tdir+'/input_models')

        if not os.path.exists(tdir+'/output_models'):
            os.mkdir(tdir+'/output_models')

        traj = md.load(trajectory_files[t], top=topology_file)

        zfill = len(str(traj.n_frames))

        # Write input PDBs
        for i in range(1, traj.n_frames+1):

            pdb_name = t+'_'+str(i).zfill(zfill)
            input_pdb = tdir+'/input_models/'+pdb_name+'.pdb'

            # Save PDB models if they do not exists
            if not os.path.exists(input_pdb):
                traj[i-1].save(input_pdb)

            # Set up Prepwizard
            output_folder = tdir+'/output_models/'+pdb_name
            command = 'cd '+output_folder+'\n'
            command += '"${SCHRODINGER}/utilities/prepwizard" '
            command += '../../input_models/'+pdb_name+'.pdb '
            command += pdb_name+'.pdb '
            command += '-fillsidechains '
            command += '-disulfides '
            command += '-rehtreat '
            command += '-propka_pH '+str(pH)+' '
            command += '-f 2005 '
            command += '-rmsd 0.3 '
            command += '-JOBNAME '+pdb_name+' '
            command += '-HOST localhost:1 '
            command += '-WAIT\n'

            input_mae = output_folder+'/'+pdb_name+'_protein.mae'
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)

            # Convert PDB into MAE for sitemap
            command += '"${SCHRODINGER}/run" ../../../._prepareForSiteMap.py '
            command += pdb_name+'.pdb '
            command += '. '
            command += '--protein_only\n'

            for r in target_residues:

                rdir = output_folder+'/'+str(r)
                if not os.path.exists(rdir):
                    os.mkdir(rdir)

                # Add site map command
                command += 'cd '+str(r)+'\n'
                command += '"${SCHRODINGER}/sitemap" '
                command += '-j '+pdb_name+' '
                command += '-prot ../'+pdb_name+'_protein.mae'+' '
                command += '-sitebox '+str(site_box)+' '
                command += '-resolution '+str(resolution)+' '
                command += '-keepvolpts yes '
                command += '-keeplogs yes '
                command += '-reportsize '+str(reportsize)+' '
                command += '-siteasl \"res.num {'+str(r)+'}\" '
                command += '-HOST localhost:1 '
                command += '-TMPLAUNCHDIR '
                command += '-WAIT\n'
                command += 'cd ../\n'

            command += 'cd ../../../../\n'
            jobs.append(command)

    return jobs

def readSiteMapForTrajectory(job_folder, verbose=True):
    """
    Read the results from a SiteMap calculation run over trajectories.
    See:
        bioprospecting.MD.sitemap.calculateSiteMapForTrajectory()

    Parameters
    ==========
    job_folder : str
        Path to the output job folder from calculateSiteMapForTrajectory().
    """

    def checkIfCompleted(log_file):
        """
        Check log file for calculation completition.

        Parameters
        ==========
        log_file : str
            Path to the standard sitemap log file.

        Returns
        =======
        completed : bool
            Did the simulation end correctly?
        """
        with open(log_file) as lf:
            for l in lf:
                if 'Version' in l:
                    version = l.split()[2]
        completed = False
        with open(log_file) as lf:
            for l in lf:
                if '2021' in version:
                    if 'Site points found' in l:
                        completed =  True
        return completed

    def checkIfFound(log_file):
        """
        Check log file for found sites.

        Parameters
        ==========
        log_file : str
            Path to the standard sitemap log file.

        Returns
        =======
        found : bool
            Did the simulation end correctly?
        """
        found = True
        with open(log_file) as lf:
            for l in lf:
                if 'Number of sites' in l:
                    found = int(l.split()[-1])
        return found

    def parseVolumeInfo(eval_log):
        """
        Parse eval log file for site scores.

        Parameters
        ==========
        eval_log : str
            Eval log file from sitemap output

        Returns
        =======
        pocket_data : dict
            Scores for the given pocket
        """
        with open(eval_log) as lf:
            c = False
            for l in lf:
                if 'SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc' in l:
                    c = True
                    labels = l.split()
                    continue
                if c:
                    values = [float(x) for x in l.split()]
                    pocket_data = {x:y for x,y in zip(labels, values)}
                    c = False
        return pocket_data

    for m in sorted(os.listdir(job_folder)):
        mdir = job_folder+'/'+m
        if not os.path.isdir(mdir):
            continue

        sm_data = {}
        sm_data['Trajectory'] = []
        sm_data['Residue'] = []
        sm_data['Frame'] = []

        for t in sorted(os.listdir(mdir+'/output_models')):
            tpath = mdir+'/output_models/'+t
            ti = int(t.split('_')[-1])

            for r in sorted(os.listdir(tpath)):
                try: int(r)
                except: continue
                rdir = tpath+'/'+r

                log_file = rdir+'/'+t+'.log'
                if os.path.exists(log_file):
                    completed = checkIfCompleted(log_file)
                else:
                    if verbose:
                        message = 'Log file for model %s, residue %s, and trajectory index %s was not found!\n' % (m, r, ti)
                        message += 'It seems the calculation has not run yet...'
                        print(message)
                    continue

                if not completed:
                    if verbose:
                        print('There was a problem with model %s, residue %s, and trajectory index %s' % (m, r, ti))
                    continue
                else:
                    found = checkIfFound(log_file)
                    if not found:
                        if verbose:
                            print('No sites were found for model %s, residue %s, and trajectory index %s' % (m, r, ti))
                        continue

                pocket_data = parseVolumeInfo(log_file)

                sm_data['Trajectory'].append(m)
                sm_data['Residue'].append(r)
                sm_data['Frame'].append(ti)

                for k in pocket_data:
                    sm_data.setdefault(k, [])
                    sm_data[k].append(pocket_data[k])


    sm_data = pd.DataFrame(sm_data)
    sm_data.set_index(['Trajectory', 'Frame', 'Residue'], inplace=True)

    return sm_data

def _copyScriptFile(output_folder, script_name, no_py=False, subfolder=None, hidden=True):
    """
    Copy a script file from the prepare_proteins package.

    Parameters
    ==========

    """
    # Get script
    path = "bioprospecting/MD/scripts"
    if subfolder != None:
        path = path+'/'+subfolder

    script_file = resource_stream(Requirement.parse("bioprospecting"),
                                     path+'/'+script_name)
    script_file = io.TextIOWrapper(script_file)

    # Write control script to output folder
    if no_py == True:
        script_name = script_name.replace('.py', '')

    if hidden:
        output_path = output_folder+'/._'+script_name
    else:
        output_path = output_folder+'/'+script_name

    with open(output_path, 'w') as sof:
        for l in script_file:
            sof.write(l)
