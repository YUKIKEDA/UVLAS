"""Typeset the canonical Japanese manuscript without changing frozen evidence."""
from pathlib import Path
import sys,re,html,base64,json,hashlib
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'build/paper-tools'))
import markdown
def main():
    if (ROOT/'docs/papers/source-hashes.json').exists():
        raise SystemExit('Report v1 is frozen; prepare a separately named revision')
    folder=ROOT/'docs/papers';p=folder/'pg-ancf-technical-report-v1.md';source=p.read_text(encoding='utf-8');tables=(folder/'pg-ancf-assets/result-tables.md').read_text(encoding='utf-8')
    quality=tables.split('\n\n')[0];physical=tables.split('| Domain |',1)[1];physical='| Domain |'+physical
    for token,value in [('@@QUALITY_TABLE@@',quality),('@@PHYSICAL_TABLE@@',physical)]:
        if token in source:source=source.replace(token,value)
    assert '@@' not in source
    p.write_text(source,encoding='utf-8')
    equations=[]
    def replace_math(match):
        equations.append(match.group(0));return '\n\nMATHBLOCKTOKEN'+str(len(equations)-1)+'\n\n'
    protected=re.sub(r'\$\$[\s\S]*?\$\$',replace_math,source)
    body=markdown.markdown(protected,extensions=['tables','fenced_code'])
    for i,eq in enumerate(equations):body=body.replace('<p>MATHBLOCKTOKEN'+str(i)+'</p>','<div class="equation">'+html.escape(eq)+'</div>')
    assert 'MATHBLOCKTOKEN' not in body
    def image_replace(m):
        filename=m.group(1);target=(folder/filename).resolve();assert target.is_relative_to(folder.resolve()) and target.is_file();kind='image/svg+xml' if target.suffix=='.svg' else 'image/png';uri='data:'+kind+';base64,'+base64.b64encode(target.read_bytes()).decode();return 'src="'+uri+'"'
    body=re.sub(r'src="(pg-ancf-assets/[^\"]+)"',image_replace,body)
    # Markdown image paragraphs become indivisible figures with explicit captions.
    body=re.sub(r'<p><img alt="([^\"]*)" src="([^\"]*)"\s*/?></p>',r'<figure><img alt="\1" src="\2"><figcaption>\1</figcaption></figure>',body)
    js=(ROOT/'build/paper-render/mathjax-3.2.2-tex-svg.js').read_text(encoding='utf-8');lock=json.loads((ROOT/'build/paper-render/mathjax-lock.json').read_text());assert hashlib.sha256(js.encode()).hexdigest()==lock['sha256']
    css='''
    @page{size:A4;margin:18mm}*{box-sizing:border-box}html{background:#fff;color:#20252b}
    body{font-family:"Yu Mincho","Noto Serif JP",serif;font-size:10pt;line-height:1.63;margin:0}
    h1,h2,h3{font-family:"Yu Gothic","Noto Sans JP",sans-serif;line-height:1.4;color:#173e50;break-after:avoid}
    h1{font-size:20pt;margin:0 0 11pt}h2{font-size:13pt;margin:20pt 0 8pt;border-bottom:1px solid #c1d0d6;padding-bottom:4pt}
    h3{font-size:10.7pt;margin:14pt 0 5pt}p{margin:0 0 8pt;orphans:3;widows:3}strong{font-weight:700}
    a{color:#175a72;text-decoration:none;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;font-size:7.6pt;line-height:1.45;margin:10pt 0 7pt;break-inside:auto}
    thead{display:table-header-group}tr{break-inside:avoid}th{background:#eaf0f3;font-family:"Yu Gothic",sans-serif;font-weight:600;text-align:left}
    th,td{border-bottom:1px solid #ccd6da;padding:5pt 4pt;vertical-align:top;overflow-wrap:anywhere}
    .equation{break-inside:avoid;margin:10pt 0;font-size:10.2pt;text-align:center}
    mjx-container[jax="SVG"]{overflow:visible!important}mjx-container svg{max-width:100%;height:auto}
    figure{margin:12pt 0;break-inside:avoid}figure img{width:100%;height:auto;max-height:76mm;object-fit:contain}figcaption{font-size:8.3pt;line-height:1.5;margin-top:5pt;color:#43545e}
    pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f4f5;padding:8pt;font-size:8pt;line-height:1.5;break-inside:avoid}code{font-family:Consolas,monospace;font-size:.92em}ul,ol{padding-left:1.7em}li{margin:3pt 0}
    body>p:first-of-type{font-family:Cambria,serif;font-size:11.5pt;color:#38515d;line-height:1.4}
    '''
    config='''window.MathJax={tex:{displayMath:[['$$','$$']],inlineMath:[],tags:'none'},svg:{fontCache:'local'},options:{enableMenu:false},startup:{typeset:true}};'''
    document='<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>PG-ANCF Technical Report v1</title><style>'+css+'</style><script>'+config+'</script><script>'+js.replace('</script','<\\/script')+'</script></head><body>'+body+'</body></html>'
    out=folder/'pg-ancf-technical-report-v1.html';out.write_text(document,encoding='utf-8')
    print(json.dumps(dict(passed=True,equations=len(equations),html_bytes=out.stat().st_size,manuscript_sha256=hashlib.sha256(source.encode()).hexdigest())))
if __name__=='__main__':main()
