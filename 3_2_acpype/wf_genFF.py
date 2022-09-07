from dflow import config, s3_config

config["host"] = ""
s3_config["endpoint"] = ""

from dflow.utils import copy_artifact, copy_file

from dflow import (
    Inputs,
    InputParameter,
    Outputs,
    OutputArtifact,
    upload_artifact,
    download_artifact,
    Workflow,
    Step,
    Steps,
    argo_range
)

from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices,
    upload_packages
)

if "__file__" in locals():
    upload_packages.append(__file__)

import os
from pathlib import Path


def os_shellcmd(cmdtxt, name='shellcmd.sh', extraparam=None):
    import os
    with open(name, 'w') as f:
        f.write(cmdtxt)
    os.system('bash %s %s' % (name, extraparam))


def shell_source(cmd):
    """Sometime you want to emulate the action of "source" in bash,
    settings some environment variables. Here is a way to do it."""
    import subprocess, os
    pipe = subprocess.Popen("%s ; env" % cmd, stdout=subprocess.PIPE, shell=True, executable="bash", encoding="utf8")
    output = pipe.communicate()[0]
    env = dict((line.split("=", 1) for line in output.splitlines() if "=" in line))
    os.environ.update(env)


class genGAFF2(OP):
    """
    genGAFF2 using Ambertools and acpype
    """

    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "mol2": Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "top": Artifact(Path)
        })

    @OP.exec_sign_check
    def execute(self, op_in):
        import os
        shell_source("source /opt/Miniconda/bin/activate && conda activate md")
        print(os.environ["PATH"])
        #resname=op_in["resname"]
        begindir = os.getcwd()
        os.chdir(op_in["mol2"].parent)
        #mol2name=op_in["mol2"].name
        mol2name="out.mol2"
        cmd='''\
#!/bin/bash
antechamber -i \$1 -fi mol2 -o Lig.mol2 -fo mol2 -pf y
parmchk2 -i Lig.mol2 -f mol2 -o Lig.frcmod
cat > leap.in << OPEOF
source leaprc.gaff2
loadamberparams Lig.frcmod
lig=loadmol2 Lig.mol2
check lig
saveamberparm lig  Lig.prmtop Lig.inpcrd
quit
OPEOF

tleap -f leap.in
#python acpype.py -p Lig.prmtop -x Lig.inpcrd -d
acpype -p Lig.prmtop -x Lig.inpcrd -d
'''
        os_shellcmd(cmd,name="gen_FF.sh",extraparam=" {}".format(mol2name))
        #os_shellcmd(cmd, name="gen_FF.sh")
        #topfile=glob.glob("*.gro")[0]
        topfile=Path(".").absolute()
        print("#####topfile",topfile)
        op_out=OPIO(
                {
            "top": topfile
             })
        os.chdir(begindir)
        return op_out

if __name__ == "__main__":
    input_mol2_artifact = upload_artifact("filesIO/out.mol2")
    input_dir_artifact = upload_artifact(Path("filesIO"))
    from dflow import RemoteExecutor

    remote_executor = RemoteExecutor(
        host="", port=22, username="", password="",
        remote_command="pyhton3"
    )

    genff = Step(name='s-genff',
                 template=PythonOPTemplate(genGAFF2,
                                           image='python:3.8',
                                           image_pull_policy='IfNotPresent',
                                           command="python3"),
                 artifacts={"mol2": input_mol2_artifact},
                 executor=remote_executor
                 )


    wf = Workflow(name="w-genff")  # step name and workflow name cannot using uppercase characters
    wf.add(genff)
    wf.submit()

    import time

    while wf.query_status() in ["Pending", "Running"]:
        time.sleep(5)
    step_todownload = wf.query_step(name="s-genff")[0]
    download_artifact(step_todownload.outputs.artifacts["top"])
