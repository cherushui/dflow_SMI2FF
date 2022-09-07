from dflow import config, s3_config
config["host"] = ""
s3_config["endpoint"] = ""
from dflow.plugins.lebesgue import LebesgueContext
lebesgue_context = LebesgueContext(
    executor="lebesgue_v2",
    extra={"scass_type":"c4_m4_cpu","program_id":},
    username="",
    password="",
    app_name='Default',
    org_id='123',
    user_id='456',
    tag='',
)


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
    Slices
)
from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
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

class g16calc(OP):
    def __init__(self):
        pass
    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "input": Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "log": Artifact(Path),
            "fchk": Artifact(Path)
        })
    @OP.exec_sign_check
    def execute(self,op_in):
        cwd=os.getcwd()
        os.chdir(op_in["input"])

        cmd="source /root/g16.sh ; g16 < DME_1.gjf > out.log ; formchk DME_1.chk ; mv DME_1.fchk out.fchk"
        os_shellcmd(cmd)
        op_out=OPIO(
                {
            "log": op_in["input"]/"out.log",
            "fchk": op_in["input"]/"out.fchk"
             })
        os.chdir(cwd)
        return op_out 

if __name__ == "__main__":
    input_artifact=upload_artifact("filesIO")
    qmcalc=Step(name='g16calc',
            template=PythonOPTemplate(g16calc,image='LBG_gau16_V1', command=["/opt/Miniconda/bin/python"],image_pull_policy='IfNotPresent'),
                artifacts={"input": input_artifact},
                parameters={}
    )
    wf = Workflow(name="wf-ele",context=lebesgue_context)
    wf.add(qmcalc)
    wf.submit()

    import time
    while wf.query_status() in ["Pending", "Running"]:
        time.sleep(10)
    step_todownload = wf.query_step(name="g16calc")[0]
    download_artifact(step_todownload.outputs.artifacts["log"])
    download_artifact(step_todownload.outputs.artifacts["fchk"])
