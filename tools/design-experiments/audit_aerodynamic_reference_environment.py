"""Record distribution metadata and every installed reference source hash."""
from pathlib import Path
import hashlib,json,sys,importlib.metadata as md
ROOT=Path(__file__).resolve().parents[2]
def main():
    name=sys.argv[1];dist=md.distribution(name);sources={}
    for f in dist.files:
        if f.suffix in ['.py','.txt','.md'] or f.name in ['METADATA','RECORD','WHEEL']:
            p=Path(dist.locate_file(f)).resolve();sources[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    dependencies={d.metadata['Name']:d.version for d in md.distributions()};out=ROOT/'build'/('reference-aerosandbox' if name.lower()=='aerosandbox' else 'reference-ptera')/'environment-manifest.json'
    out.write_text(json.dumps(dict(distribution=dist.metadata['Name'],version=dist.version,python=sys.version,dependencies=dependencies,sources=sources),indent=2),encoding='utf-8');print(name,dist.version,len(sources))
if __name__=='__main__':main()
