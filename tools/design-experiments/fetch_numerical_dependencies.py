"""Acquire pinned official header sources into the workspace, never install globally."""
from pathlib import Path
import urllib.request,hashlib,json,tarfile,shutil,time
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/numerical-deps'
PACKAGES=[('eigen','5.0.1','https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz'),('spectra','1.2.0','https://github.com/yixuan/spectra/archive/refs/tags/v1.2.0.tar.gz'),('boost','1.92.0','https://archives.boost.io/release/1.92.0/source/boost_1_92_0.tar.bz2')]
def main():
    OUT.mkdir(exist_ok=True);records=[]
    for name,version,url in PACKAGES:
        target=OUT/url.rsplit('/',1)[-1];print('FETCH',name,url,flush=True);start=time.monotonic()
        if not target.exists():
            request=urllib.request.Request(url,headers={'User-Agent':'UVLAS-numerical-dependency-spike/1.0'})
            with urllib.request.urlopen(request,timeout=60) as response,target.with_suffix(target.suffix+'.partial').open('wb') as stream:shutil.copyfileobj(response,stream)
            target.with_suffix(target.suffix+'.partial').replace(target)
        digest=hashlib.sha256(target.read_bytes()).hexdigest();count=0;roots=set()
        with tarfile.open(target,'r:*') as archive:
            for m in archive:
                parts=Path(m.name).parts
                if not parts:continue
                roots.add(parts[0])
                if not m.isfile():continue
                if name=='boost' and not ('boost' in parts[1:2] or parts[-1] in ['LICENSE_1_0.txt','boost.css']):continue
                p=(OUT/m.name).resolve()
                if not p.is_relative_to(OUT.resolve()):raise ValueError('archive_path_escape')
                p.parent.mkdir(parents=True,exist_ok=True)
                with archive.extractfile(m) as src,p.open('wb') as dst:shutil.copyfileobj(src,dst)
                count+=1
        record=dict(name=name,version=version,url=url,archive=target.relative_to(ROOT).as_posix(),sha256=digest,roots=sorted(roots),extracted_files=count,seconds=time.monotonic()-start);records.append(record);print(json.dumps(record),flush=True)
    dest=OUT/'header-sources.lock.json';assert not dest.exists();dest.write_text(json.dumps(dict(packages=records),indent=2),encoding='utf-8')
if __name__=='__main__':main()
