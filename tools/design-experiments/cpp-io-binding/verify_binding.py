from pathlib import Path
import sys,json,gc
sys.path.insert(0,str(Path(sys.argv[1]).resolve()))
import _uvlas_design_binding as b
import numpy as np
r=b.compute(1000);view=r.view();assert not view.flags.writeable;del r;gc.collect();assert np.array_equal(view,np.arange(1000)*.25)
try:view[0]=3.;raise AssertionError('writable')
except ValueError:pass
assert b.compute(0).view().shape==(0,)
good='{"schema":"uvlas.design-io/1","id":"18446744073709551615","scale":1.0}'
assert json.loads(b.validate(good))['id']=='18446744073709551615'
bad=[good.replace('1.0','0.0'),good.replace('1.0','NaN'),good.replace('"id":','"extra":0,"id":'),good[:-1]+',"scale":2}',good.replace('18446744073709551615','18446744073709551616'),good.replace('"18446744073709551615"','9007199254740993'),good.replace('18446744073709551615','01')]
for value in bad:
    try:b.validate(value);raise AssertionError('accepted invalid input '+value)
    except (ValueError,RuntimeError,OverflowError,IndexError):pass
try:b.compute(1000001);raise AssertionError('budget')
except ValueError:pass
print(json.dumps(dict(passed=True,read_only_owner_lifetime=True,empty_result=True,max_uint64_string=True,negative_controls=len(bad)+2,scope='binding and metadata spike; not product parser or solver')))
