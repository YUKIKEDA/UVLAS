"""Pinned official sources for the preimplementation I/O/binding spike only."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import urllib.request,hashlib,json,tarfile,shutil
ROOT=rd_paths.ROOT;OUT=ROOT/'build/numerical-deps'
def main():
    records=[];lock=OUT/'io-sources.lock.json';assert not lock.exists()
    for name,version,url in [('pybind11','3.1.0','https://github.com/pybind/pybind11/archive/refs/tags/v3.1.0.tar.gz'),('json','3.12.0','https://github.com/nlohmann/json/archive/refs/tags/v3.12.0.tar.gz')]:
        archive=OUT/f'{name}-{version}.tar.gz'
        if not archive.exists():
            with urllib.request.urlopen(url,timeout=60) as src,archive.open('wb') as dst:shutil.copyfileobj(src,dst)
        files=[]
        with tarfile.open(archive,'r:gz') as tar:
            for item in tar:
                if not item.isfile():continue
                target=(OUT/item.name).resolve();assert target.is_relative_to(OUT.resolve())
                target.parent.mkdir(parents=True,exist_ok=True)
                with tar.extractfile(item) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
                files.append(target.relative_to(ROOT).as_posix())
        record=dict(name=name,version=version,url=url,archive=archive.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),extracted_files=len(files));records.append(record);print(json.dumps(record),flush=True)
    lock.write_text(json.dumps(dict(packages=records),indent=2),encoding='utf-8')
if __name__=='__main__':main()
