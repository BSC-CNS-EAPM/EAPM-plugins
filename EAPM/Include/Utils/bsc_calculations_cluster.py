localIPs = {"cactus": "84.88.51.217", "alma": "84.88.51.250", "bubbles": "84.88.51.219"}


def setup_bsc_calculations_based_on_horus_remote(
    remote_name, remote_host: str, jobs, partition, scriptName, cpus, job_name, program
):
    import bsc_calculations

    cluster = "local"

    if remote_name != "local":
        cluster = remote_host

    if remote_host in localIPs.values():
        cluster = "powerpuff"

    ## Define cluster
    # cte_power
    if cluster == "plogin1.bsc.es":
        bsc_calculations.cte_power.jobArrays(
            jobs,
            job_name=job_name,
            partition=partition,
            program=program,
            script_name=scriptName,
            gpus=cpus,
        )
    # marenostrum
    elif cluster == "mn1.bsc.es":
        bsc_calculations.marenostrum.jobArrays(
            jobs,
            job_name=job_name,
            partition=partition,
            program=program,
            script_name=scriptName,
            cpus=cpus,
        )
    # minotauro
    elif cluster == "mt1.bsc.es":
        bsc_calculations.minotauro.jobArrays(
            jobs,
            job_name=job_name,
            partition=partition,
            program=program,
            script_name=scriptName,
            gpus=cpus,
        )
    # powerpuff
    elif cluster == "powerpuff":
        print("Generating powerpuff girls jobs...")
        bsc_calculations.local.parallel(
            jobs,
            cpus=min(40, len(jobs)),
            script_name=scriptName,
        )
    # local
    elif cluster == "local":
        print("Generating local jobs...")
        bsc_calculations.local.parallel(
            jobs,
            cpus=min(40, len(jobs)),
            script_name=scriptName,
        )
    else:
        raise Exception("Cluster not supported.")

    return cluster
