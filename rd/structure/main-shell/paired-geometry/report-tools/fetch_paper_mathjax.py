from pathlib import Path
import urllib.request,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
def main():
    directory=ROOT/'build/paper-render';directory.mkdir(parents=True,exist_ok=True)
    url='https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js';p=directory/'mathjax-3.2.2-tex-svg.js'
    if not p.exists():
        with urllib.request.urlopen(url,timeout=60) as response:p.write_bytes(response.read())
    assert p.stat().st_size>100000
    report=dict(version='3.2.2',url=url,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),license='Apache-2.0',purpose='local mathematical typesetting, no runtime CDN fetch')
    (directory/'mathjax-lock.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))
if __name__=='__main__':main()
