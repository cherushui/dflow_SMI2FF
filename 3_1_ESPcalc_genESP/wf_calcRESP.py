from dflow import config, s3_config

from dflow.utils import copy_artifact,copy_file

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

def os_shellcmd(cmdtxt,name='shellcmd.sh',extraparam=None):
    with open(name,'w') as f:
        f.write(cmdtxt)
    os.system('bash %s %s'%(name,extraparam))

def shell_source(cmd):
    """Sometime you want to emulate the action of "source" in bash,
    settings some environment variables. Here is a way to do it."""
    import subprocess, os
    pipe = subprocess.Popen("%s ; env" % cmd, stdout=subprocess.PIPE, shell=True,executable="bash", encoding="utf8")
    output = pipe.communicate()[0]
    env = dict((line.split("=", 1) for line in output.splitlines() if "=" in line))
    os.environ.update(env)

class calcRESP(OP):
    """
    using Multiwfn to generate RESP charge.
    """
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "fchk": Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "chg": Artifact(Path)
        })
            
    @OP.exec_sign_check
    def execute(self,op_in):
        begindir=os.getcwd()
        fchkname=op_in["fchk"]
        os.chdir(fchkname.parent)
        multiwfncmdtxt="""\
7
18
5
11
a
y
q
1
eqvcons_PG.txt
1
y
0
0
q
"""
        txtfilename="multiwfncmd.txt"
        with open(txtfilename,"w") as f:
            f.write(multiwfncmdtxt)
        os.environ["PATH"]="/root/software/Multiwfn_3.8_dev_bin_Linux_noGUI::"+os.environ["PATH"]
        #shell_source("export PATH=/root/software/Multiwfn_3.8_dev_bin_Linux_noGUI:$PATH")
        os_shellcmd("Multiwfn {} < {}".format(fchkname,txtfilename))
        chgout=Path(fchkname.stem+".chg").absolute()
        #print("####chgout",chgout)
        op_out=OPIO(
                {
            "chg": chgout
             })
        os.chdir(begindir)
        return op_out


class genmol2(OP):
    """
    using openbabel to generate mol2 file, containing RESP charge.
    """

    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "fchk": Artifact(Path),
            "chg":Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "mol2": Artifact(Path)
        })

    @OP.exec_sign_check
    def execute(self, op_in):
        begindir = os.getcwd()
        fchkname = op_in["fchk"]
        chgname = op_in["chg"]
        os.chdir(fchkname.parent)
        #shell_source("source /opt/Miniconda/bin/activate && conda activate md")
        from openbabel import pybel
        mol2string=self.genmol2(str(fchkname.absolute()))
        mol2stirng_resp=self.dumpchg(mol2string,str(chgname.absolute()))
        mol2name=Path(fchkname.stem+".mol2").absolute()
        with open(str(mol2name),"w") as f:
            f.write(mol2stirng_resp)
        op_out=OPIO(
                {
            "mol2": mol2name
             })
        os.chdir(begindir)
        return op_out
    def genmol2(self,fchkname):
        from openbabel import pybel
        mol=next(pybel.readfile("fchk",fchkname))
        mol2string=mol.write("mol2")
        return mol2string
    def dumpchg(self,mol2string,chgfile):
        #readresp
        with open(chgfile, "r") as f:
            chargetxt = f.readlines()
        respcharges = []
        for linetxt in chargetxt:
            respcharges.append(linetxt.strip().split()[-1])
        #dump charge
        tag_change_atom = False
        atom_fmt = "{} {} {} {} {} {} {} {} {}"
        count = 0
        mol2string_resp = ""
        for idx, linetxt in enumerate(mol2string.splitlines()):
            if idx == 4:
                linetxt = 'USER_CHARGES'
            if linetxt == "@<TRIPOS>BOND":
                tag_change_atom = False
            if tag_change_atom == True:
                (atom_id, atom_name, x, y, z, atom_type, subst_id, subst_name, charge) = linetxt.strip().split()
                charge = respcharges[count]
                count += 1
                # subst_name=?
                txtstring = atom_fmt.format(atom_id, atom_name, x, y, z, atom_type, subst_id, subst_name, charge)
                txtstring += '\n'
            else:
                txtstring = linetxt + '\n'
            if linetxt == "@<TRIPOS>ATOM":
                tag_change_atom = True
            mol2string_resp += txtstring
        return mol2string_resp

if __name__ == "__main__":
    input_fchk_artifact=upload_artifact("filesIO/out.fchk")
    from dflow import RemoteExecutor
    remote_executor = RemoteExecutor(
        host="", port=22, username="", password="",
        remote_command="python3"
    )

    cresp=Step(name='s-resp',
            template=PythonOPTemplate(calcRESP,
                                      image='python:3.8',
                                      image_pull_policy='IfNotPresent'),
                artifacts={"fchk": input_fchk_artifact},
                executor=remote_executor
    )
    gmol2=Step(name='s-mol2',
            template=PythonOPTemplate(genmol2,
                                      image='python:3.8',
                                      image_pull_policy='IfNotPresent'),
                artifacts={"fchk": input_fchk_artifact,
                           "chg":cresp.outputs.artifacts["chg"]},
                executor=remote_executor
    )
    wf = Workflow(name="w-resp")  #step name and workflow name cannot using uppercase characters
    wf.add(cresp)
    wf.add(gmol2)
    wf.submit()

    import time
    while wf.query_status() in ["Pending", "Running"]:
        time.sleep(5)
    step_todownload = wf.query_step(name="s-mol2")[0]
    download_artifact(step_todownload.outputs.artifacts["mol2"])
