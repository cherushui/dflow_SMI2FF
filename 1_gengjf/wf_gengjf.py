from dflow import config, s3_config
config["host"] = ""
s3_config["endpoint"] = ""

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

import os,subprocess
from pathlib import Path

def os_shellcmd(cmdtxt,name='shellcmd.sh',extraparam=None):
    with open(name,'w') as f:
        f.write(cmdtxt)
    os.system('bash %s %s'%(name,extraparam))

import json
class workflowset:
    def __init__(self,filename="electrolyte.json"):
        with open(filename,"r") as f:
            self.infodict=json.load(f)
        #self.smi=self.infodict["input"]["smi"]
        #self.gausskeyword=self.json["input"]["gausskeyword"]

class gengjf(OP):
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "name":str,
            "smi": str,
            "gauss_keyword":str,
            "charge":int,
            "multiplicity":int,
            "workflowset": workflowset
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "gengjf": Artifact(Path),
            "workflowset": workflowset
        })
            
    @OP.exec_sign_check
    def execute(self,op_in):
        name = op_in["name"]
        smi = op_in["smi"]
        gauss_keyword = op_in["gauss_keyword"]
        charge = op_in["charge"]
        multiplicity = op_in["multiplicity"]

        workflowset=op_in["workflowset"]
        #convert smi to xyz
        from openbabel import pybel
        mol=pybel.readstring('smi',smi)
        mol.make3D()
        xyzstring=mol.write("xyz")
        coor="\n".join(xyzstring.splitlines()[2:])  #delete first two lines
        #name="test" #change to a input?

        gjfhead="%%chk=%s\n%%mem=4GB\n%%nproc=4\n\n" %name
        gjfkeyword="{keyword}\n\n{name}\n\n{charge} {multiplicity}\n".format(
            keyword=gauss_keyword,
            name=name,
            charge=charge,
            multiplicity=multiplicity)
        gjftail="\n\n"
        with open("%s.gjf"%name,"w") as f:
            txt=gjfhead+gjfkeyword+coor+gjftail
            f.write(txt)
        op_out=OPIO(
                {
            "gengjf": Path("%s.gjf" %name),
            "workflowset": workflowset
             })
        return op_out 

    # @classmethod
    #def formatinput(cls,*wargs,**kwargs):
    # def formatinput(cls,dictinput)
    #     length = len(keylist)
    #     js_validation = None
    #     for i in range(length):
    #         if i == 0 and keylist[i] in self.scanjs:
    #             js_validation = copy.deepcopy(self.scanjs[keylist[i]])
    #         elif i > 0 and keylist[i] in js_validation:
    #             js_validation = js_validation[keylist[i]]
    #         else:
    #             raise RuntimeError('Key {} is missing, plz check the json file.'.format(keylist))
    #         if js_validation == '' or js_validation == None:
    #             raise RuntimeError('Key {} is empty, plz check the json file.'.format(keylist))
    #     return
    # @classmethod
    # def check_necessary(cls,keylist):
    #     length=len(keylist)
    #     for i in range(length):

if __name__ == "__main__":
    #inputartifact=upload_artifact("ready_files")
    workflowset_instance=workflowset("electrolyte.json")

    infodict = workflowset_instance.infodict
    smi=infodict["input"]["smi"][0]
    name=infodict["input"]["name"][0]
    formatted_input={
        "name":name,
        "smi": smi,
        "gauss_keyword": infodict["input"]["OPT"]["gauss_keyword"],
        "charge": int(infodict["input"]["charge"][0]),
        "multiplicity": int(infodict["input"]["multiplicity"][0]),
        "workflowset": workflowset_instance
                     }
    print("formatted_input",formatted_input)
    from dflow import RemoteExecutor
    remote_executor = RemoteExecutor(
        host="", port=22, username="", password="",
        remote_command="python3"
    )

    qmcalc=Step(name='step-gengjf',
            template=PythonOPTemplate(gengjf, image='python:3.8',image_pull_policy='IfNotPresent'),
                artifacts={},
                parameters=formatted_input,
                executor=remote_executor
    )


    wf = Workflow(name="test-gengjf")
    wf.add(qmcalc)
    wf.submit()

    #print('wf.id:',wf.id)
    import time
    while wf.query_status() in ["Pending", "Running"]:
        time.sleep(5)
    step_todownload = wf.query_step(name="step-gengjf")[0]
    download_artifact(step_todownload.outputs.artifacts["gengjf"])
