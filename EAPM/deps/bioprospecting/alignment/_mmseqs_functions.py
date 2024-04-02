import subprocess
import shutil
import os

class mmseqs2:

    def clusterSequences(sequences, min_seq_id=0.5, c=0.8, cov_mode=1, stdout=True, stderr=True):
        """
        Cluster sequences using MMseqs2 algorithm.

        Paramters
        =========
        sequences : dict
            The sequences to cluster.
        """

        # Manage stdout and stderr
        if stdout:
            stdout = None
        else:
            stdout = subprocess.DEVNULL

        if stderr:
            stderr = None
        else:
            stderr = subprocess.DEVNULL

        # Write input fasta
        with open('mmseqs_tmp.fasta', 'w') as of:
            for name in sequences:
                of.write('>'+name+'\n')
                of.write(sequences[name]+'\n')

        # Execute mmseqs2
        command =  'mmseqs easy-cluster '
        command += 'mmseqs_tmp.fasta '
        command += 'clusterRes '
        command += 'tmp '
        command += '--min-seq-id '+str(min_seq_id)+' '
        command += '-c '+str(c)+' '
        command += '--cov-mode '+str(cov_mode)+' '
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        # Read clustering results
        clusters = {}
        with open('clusterRes_cluster.tsv') as cr:
            for l in cr:
                cn, m = l.split()
                clusters.setdefault(cn, [])
                clusters[cn].append(m)

        os.remove('mmseqs_tmp.fasta')
        os.remove('clusterRes_cluster.tsv')
        os.remove('clusterRes_all_seqs.fasta')
        os.remove('clusterRes_rep_seq.fasta')
        shutil.rmtree('tmp')

        return clusters
