"""Check prescribed-path evidence; never implies production acceptance."""
import argparse
import json
from pathlib import Path

import numpy as np

NAME="RIVASAND02BranchReversalResearch"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results",type=Path)
    HERE=parser.parse_args().results.resolve()
    checks={}
    def load_csv(path):
        values=np.loadtxt(path,delimiter=",",skiprows=1)
        assert np.isfinite(values).all(), f"Nonfinite output in {path}"
        return values
    for n in [128,512,2048,8192]:
        load=lambda name,v: load_csv(HERE/f"rotating_results/{name}_v{v}_n{n}.csv")
        original=load("RIVASAND02","0.0");research=load(NAME,"0.0")
        np.testing.assert_array_equal(original,research)
        checks[f"pure_shear_identical_{n}"]=True
        mixed=load(NAME,"0.0001")
        assert mixed[-1,-1]==4
        checks[f"mixed_path_four_reversals_{n}"]=True
        old=load("RIVASAND02","0.0001")
        assert old[-1,-1]==(4 if n==128 else 0)
    for layer in ["bottom","loose"]:
        a=load_csv(HERE/f"rotation_results/{layer}_0.0.csv")
        b=load_csv(HERE/f"rotation_results/{layer}_0.61.csv")
        c,s=np.cos(.61),np.sin(.61)
        q=np.array([[c,-s,0],[s,c,0],[0,0,1]])
        tensors=np.zeros((len(a),3,3))
        for k,(i,j) in enumerate([(0,0),(1,1),(2,2),(0,1),(1,2),(0,2)]):
            tensors[:,i,j]=tensors[:,j,i]=a[:,k+1]
        rotated=q@tensors@q.T
        expected=np.column_stack([rotated[:,i,j] for i,j in [(0,0),(1,1),(2,2),(0,1),(1,2),(0,2)]])
        error=float(np.max(abs(expected-b[:,1:7])))
        np.testing.assert_allclose(expected,b[:,1:7],atol=1e-8,rtol=1e-10)
        np.testing.assert_array_equal(a[:,-1],b[:,-1])
        checks[f"{layer}_rotation_max_stress_error_kPa"]=error
    assert "PASS: research latch rejection, trial-path independence, checkpoint rejection" in (HERE/"guard_checks.log").read_text()
    checks["latch_and_checkpoint_rejection_and_trial_path_independence"]=True
    result={"status":"passed", "scope":"Prescribed strain paths and adapter guards; no production acceptance", "checks":checks}
    (HERE/"verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__": main()
