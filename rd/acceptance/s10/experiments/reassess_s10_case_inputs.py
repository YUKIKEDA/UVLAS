"""Inspect representative public source files before defining experimental gates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import sys,urllib.request,json,zipfile,hashlib,io,xml.etree.ElementTree as ET
import numpy as np
from scipy.io import loadmat
from inspect_s10_public_resources import RemoteZip
ROOT=rd_paths.ROOT;OUT=ROOT/'build/analysis/s10-reassessment'
class BoundedRaw(RemoteZip):
    def read(self,n=-1):
        if n<0:n=self.size-self.pos
        n=min(n,self.size-self.pos)
        if n==0:return b''
        if self.transferred+n>128*1024**2:raise ValueError('128 MiB transfer ceiling')
        chunks=[]
        while n:
            count=min(n,8*1024**2);start=self.pos;end=start+count-1
            request=urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}','User-Agent':'UVLAS-public-case-audit/1.0'})
            with urllib.request.urlopen(request,timeout=60) as response:
                if response.status!=206 or response.headers.get('Content-Range')!=f'bytes {start}-{end}/{self.size}':raise ValueError('range response')
                data=response.read(count+1)
            assert len(data)==count;chunks.append(data);self.transferred+=count;self.pos+=count;self.ranges.append([start,end]);n-=count
        return b''.join(chunks)
def describe(value):
    if isinstance(value,dict):return {k:describe(v) for k,v in value.items()}
    a=np.asarray(value)
    out=dict(shape=list(a.shape),dtype=str(a.dtype))
    if a.dtype.kind in 'fiu':
        out['finite']=bool(np.isfinite(a).all())
        if a.size:out.update(min=float(np.nanmin(a)),max=float(np.nanmax(a)))
        if a.size<=40:out['values']=a.tolist()
    elif a.dtype.kind in 'US' and a.size<=20:out['values']=a.tolist()
    return out
def spreadsheet():
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(OUT/'Run_Numbers.xlsx') as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for item in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns):strings.append(''.join(item.itertext()))
        sheets={}
        for name in z.namelist():
            if not(name.startswith('xl/worksheets/sheet') and name.endswith('.xml')):continue
            rows=[]
            for row in ET.fromstring(z.read(name)).findall('.//s:row',ns):
                cells={}
                for cell in row.findall('s:c',ns):
                    v=cell.find('s:v',ns)
                    if v is not None:cells[cell.attrib['r']]=strings[int(v.text)] if cell.attrib.get('t')=='s' else v.text
                    elif cell.attrib.get('t')=='inlineStr':cells[cell.attrib['r']]=''.join(cell.itertext())
                if cells:rows.append(cells)
            sheets[name]=rows
        return sheets
def main():
    dest=OUT/'case-input-audit.json';assert not dest.exists();meta=json.loads((ROOT/'build/reference-audit/s10-metadata.json').read_text());raw=next(f for f in meta['files'] if f['key']=='raw_data.zip');stream=BoundedRaw(raw['links']['self'],raw['size']);report=dict(run_spreadsheet=spreadsheet(),raw_files={})
    selectors=['Type_GVT_ID_4003_rAoA_0_vel_0_TipWeight_FPS_TE_II_ver_1.mat','Type_GVT_ID_4004_rAoA_0_vel_0_TipWeight_FPS_LE_II_ver_1.mat','Type_static_ID_4012_rAoA_3_vel_20_TipWeight_FPS_TE_II_ver_1.mat','Type_Flutter_ID_4030_rAoA_0_vel_20p_TipWeight_FPS_LE_II_ver_1.mat']
    with zipfile.ZipFile(stream) as z:
        for selection in selectors:
            names=[n for n in z.namelist() if n.endswith(selection)];assert len(names)==1,selection;name=names[0];p=OUT/Path(name).name
            if not p.exists():p.write_bytes(z.read(name))
            assert p.stat().st_size==z.getinfo(name).file_size
            data=loadmat(p,simplify_cells=True);fields={k:describe(v) for k,v in data.items() if not k.startswith('__')};report['raw_files'][name]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),crc=z.getinfo(name).CRC,fields=fields);print('INSPECTED',name,list(fields),flush=True)
    report['transferred_bytes']=stream.transferred;report['ranges']=stream.ranges
    with zipfile.ZipFile(ROOT/'build/reference-audit/s10-analyzed.zip') as z:
        report['analyzed_modes']={}
        for config in ['LE','TE']:
            data=loadmat(io.BytesIO(z.read(f'Analyzed_data/modes_{config}.mat')),simplify_cells=True)[f'modes_{config}'];report['analyzed_modes'][config]=describe(data)
    dest.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8');print('COMPLETE',stream.transferred,flush=True)
if __name__=='__main__':main()
