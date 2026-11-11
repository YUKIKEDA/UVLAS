"""Fixed candidate interface; only independently screened cells may be used."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib
import json
import numpy as np
import mq3i1_moments_v2 as candidate

base=candidate.family.base
raw=candidate.family.raw
data=candidate.family.data
volume=candidate.family.volume
moments=candidate.family.moments
coefficients=candidate.coefficients
system=candidate.system
evaluate=candidate.evaluate
DIM=117
WORK=True
METRIC='covariant'
ROOT=rd_paths.ROOT
KEYS=['q','T','rows','U','Qs','hard_U','length','geometry_projection','geometry_vectors','A','b']
STATES={}


def configure(metric,work):
    assert metric=='covariant' and work
    candidate.family.configure(metric,work)


def register():
    for name in ['mq3i1-local-v2-result','mq3i1-additional-local-result']:
        result=json.loads((rd_paths.RECORDS/f'{name}.json').read_text())
        assert result['all_checks_pass'], name
        for source,sha in result['sources'].items():
            assert hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==sha
        for entry in result['states']:
            path=ROOT/entry['artifact']
            assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['artifact_sha256']
            saved=np.load(path)
            state={key:saved[key] for key in KEYS}
            state.update(space_ranks=entry['preparation']['space_ranks'],
                constraint_ranks=entry['preparation']['base_constraint_ranks'])
            STATES[hashlib.sha256(state['q'].tobytes()).hexdigest()]=state


def choose(q):
    key=hashlib.sha256(q.tobytes()).hexdigest()
    if key not in STATES:
        raise ValueError('reference_cell_not_screened:'+key)
    return STATES[key]
