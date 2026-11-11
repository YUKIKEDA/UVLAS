"""Load all eight decks without solving; record the standard console I4 limitation."""
from pathlib import Path
import datetime,json,os,re,shutil,subprocess,time
import avl_reference_resources as resource
ROOT=Path(__file__).resolve().parents[2];FOLDER=ROOT/'docs/verification';digest=resource.digest
def main():
    output=FOLDER/'avl-span1280-load-result.json'
    if output.exists():raise FileExistsError(output)
    prep=json.loads((FOLDER/'avl-span1280-fine-preparation.json').read_text())
    exe=ROOT/'build/reference-avl/source-symmetry16000-ks1000-spacing-pi64/bin/avl'
    root=ROOT/'build/analysis'/('avl-span1280-load-'+datetime.datetime.now().strftime('%Y%m%d-%H%M%S'));root.mkdir()
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(ROOT/'build/reference-avl/toolchain/lib');begin=time.monotonic();records=[]
    for deck in prep['decks']:
        folder=root/f'{deck["case"]}-M{deck["m"]}';folder.mkdir();shutil.copyfile(ROOT/deck['path'],folder/'wing.avl')
        commands='load wing.avl\nquit\n';(folder/'session.txt').write_text(commands)
        result=subprocess.run(['/usr/bin/time','-v','-o','resource.txt',str(exe)],input=commands.encode(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,cwd=folder,env=env,timeout=10)
        (folder/'console.txt').write_bytes(result.stdout);text=result.stdout.decode();timing=(folder/'resource.txt').read_text()
        peak=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',timing)[1])*1024
        assert result.returncode==0 and peak<256*1024**2 and time.monotonic()-begin<60
        assert not re.search(r'overflow|increase|REJECTED|NaN|Infinity',text,re.I)
        assert re.search(r'\b640 Strips',text) and text.count('AVL   c>')>=2
        expected=str(deck['native_unknowns']) if deck['native_unknowns']<=9999 else r'\*{4}'
        assert re.search(expected+r' Vortices',text)
        records.append(dict(case=deck['case'],m=deck['m'],directory=folder.relative_to(ROOT).as_posix(),strips=640,
            vortices_console_numeric=deck['native_unknowns'] if deck['native_unknowns']<=9999 else None,
            i4_display_overflow=deck['native_unknowns']>9999,peak_bytes=peak,
            files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in folder.iterdir()}))
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file());assert size<8*1024**2
    output.write_text(json.dumps(dict(source_sha256=digest(__file__),executable_sha256=digest(exe),
        protocol_sha256=digest(FOLDER/'avl-span1280-ksmax-addendum.md'),display_protocol_sha256=digest(FOLDER/'avl-span1280-load-diagnostics-addendum.md'),
        cases=records,seconds=time.monotonic()-begin,total_bytes=size,no_solve=True,all_criteria_pass=True),indent=2)+'\n')
    print(json.dumps(dict(loads=len(records),numeric_vortex_counts=6,i4_limit_recorded=2,all_criteria_pass=True)))
if __name__=='__main__':main()
