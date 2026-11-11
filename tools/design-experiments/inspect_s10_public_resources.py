"""Acquire small public run metadata and ZIP directory without downloading 1.9 GB."""
from pathlib import Path
import urllib.request,json,io,zipfile,hashlib
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'build/analysis/s10-reassessment';OUT.mkdir(exist_ok=True)
class RemoteZip(io.RawIOBase):
    def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.transferred=0;self.ranges=[]
    def seekable(self):return True
    def tell(self):return self.pos
    def seek(self,offset,whence=0):
        self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
        if not 0<=self.pos<=self.size:raise ValueError('seek range')
        return self.pos
    def read(self,n=-1):
        if n<0:n=self.size-self.pos
        n=min(n,self.size-self.pos)
        if n==0:return b''
        if n>8*1024**2 or self.transferred+n>32*1024**2:raise ValueError('bounded resource inspection exceeded')
        start=self.pos;end=start+n-1;request=urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}','User-Agent':'UVLAS-public-evidence-audit/1.0'})
        with urllib.request.urlopen(request,timeout=60) as response:
            if response.status!=206 or response.headers.get('Content-Range')!=f'bytes {start}-{end}/{self.size}':raise ValueError(f'range unsupported: {response.status} {response.headers.get("Content-Range")}')
            data=response.read(n+1)
        if len(data)!=n:raise ValueError('range length')
        self.transferred+=n;self.ranges.append([start,end]);self.pos+=n;return data
def main():
    meta=json.loads((ROOT/'build/reference-audit/s10-metadata.json').read_text());items={x['key']:x for x in meta['files']};report=dict(source='https://zenodo.org/records/19114860',version='v1',files={})
    item=items['Run_Numbers.xlsx'];p=OUT/'Run_Numbers.xlsx'
    if not p.exists():
        with urllib.request.urlopen(item['links']['self'],timeout=60) as response:p.write_bytes(response.read())
    assert len(p.read_bytes())==item['size'];assert 'md5:'+hashlib.md5(p.read_bytes()).hexdigest()==item['checksum'];report['files'][p.name]=dict(bytes=item['size'],sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    item=items['raw_data.zip'];stream=RemoteZip(item['links']['self'],item['size'])
    try:
        with zipfile.ZipFile(stream) as z:
            entries=[dict(name=x.filename,size=x.file_size,compressed=x.compress_size,crc=x.CRC) for x in z.infolist()];report['raw_inventory']=entries
            selected=[x for x in z.infolist() if x.filename.lower().endswith(('sensor_loc.m','readme.txt','readme.m')) and x.file_size<1024**2]
            for x in selected:
                data=z.read(x);target=OUT/Path(x.filename).name;assert not target.exists();target.write_bytes(data);report['files'][target.name]=dict(original_name=x.filename,bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
        report['range_inspection']='complete'
    except Exception as e:report['range_inspection']='unavailable';report['range_error']=repr(e)
    report['transferred_raw_bytes']=stream.transferred;report['ranges']=stream.ranges;dest=OUT/'public-inspection.json';assert not dest.exists();dest.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='raw_inventory'}));print('raw_inventory_entries',len(report.get('raw_inventory',[])))
if __name__=='__main__':main()
