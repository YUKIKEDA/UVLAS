"""Finalize a separate fine runner only after corrected loads/regressions are audited."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    output=F/'avl-span1280-fine-ks-preparation.json'
    if output.exists():raise FileExistsError(output)
    audit=F/'avl-span1280-capacity-ks-artifact-audit.json';assert json.loads(audit.read_text())['all_criteria_pass']
    fine=json.loads((F/'avl-span1280-fine-ks-pending.json').read_text());cap=json.loads((F/'avl-span1280-capacity-ks-preparation.json').read_text())
    fine.update(capacity_audit_sha256=digest(audit),parent_runs_sha256=digest(F/'avl-span1280-capacity-ks-runs.json'),
        executable_sha256=digest(ROOT/cap['target']/'bin/avl'),source_sha256=digest(__file__),load_sha256=digest(F/'avl-span1280-load-result.json'))
    output.write_text(json.dumps(fine,indent=2)+'\n')
    text=rd_paths.sibling('run_avl_span1280_fine.py').read_text()
    for a,b in [('avl-span1280-fine-runs','avl-span1280-fine-ks-runs'),('avl-span1280-fine-preparation','avl-span1280-fine-ks-preparation'),
        ('avl-span1280-capacity-artifact-audit','avl-span1280-capacity-ks-artifact-audit'),('avl-span1280-capacity-preparation','avl-span1280-capacity-ks-preparation')]:text=text.replace(a,b)
    text=text.replace("'avl-span1280-resource-guard-result.json']","'avl-span1280-resource-guard-result.json','avl-span1280-load-result.json']")
    text=text.replace("protocol=ROOT/prep['protocol'];assert digest(protocol)==prep['protocol_sha256']", "protocol=ROOT/prep['protocol'];assert digest(protocol)==prep['protocol_sha256']\n    assert prep['correction_protocol_sha256']==digest(FOLDER/'avl-span1280-ksmax-addendum.md')")
    rd_paths.sibling('run_avl_span1280_fine_ks.py').write_text(text)
    print('Fine runner finalized after corrected gate')
if __name__=='__main__':main()
