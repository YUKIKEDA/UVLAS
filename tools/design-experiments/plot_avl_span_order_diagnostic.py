"""Standalone scientific SVG charts using ReportLab's chart library."""
from pathlib import Path
import hashlib,json,math,xml.etree.ElementTree as ET
import reportlab
from reportlab.graphics.shapes import Drawing,String,Rect
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics import renderSVG
from reportlab.lib import colors
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    result=F/'avl-span-order-diagnostic-result.json';r=json.loads(result.read_text());assert r['all_criteria_pass']
    folder=F/'figures';folder.mkdir(exist_ok=True);files={};data={};palette=[colors.HexColor(x) for x in ['#0072B2','#999999','#D55E00']]
    draw=Drawing(550,340);draw.add(String(35,316,'FW-V01: spatial shape of consecutive span differences',fontName='Helvetica-Bold',fontSize=13))
    c=next(s for s in r['channels'] if s['case']=='FW-V01' and s['channel']=='tau0');centers=[(-3.8+.4*j)/4 for j in range(20)]
    series=[[(x,v/s['norm']) for x,v in zip(centers,s['difference'])] for s in c['stages']];data['normalized_tau0']=series
    chart=LinePlot();chart.x=65;chart.y=84;chart.width=445;chart.height=208;chart.data=series;chart.joinedLines=1
    chart.xValueAxis.valueMin=-1;chart.xValueAxis.valueMax=1;chart.xValueAxis.valueStep=.2
    chart.yValueAxis.valueMin=0;chart.yValueAxis.valueMax=.5;chart.yValueAxis.valueStep=.1
    for j,col in enumerate(palette):chart.lines[j].strokeColor=col;chart.lines[j].strokeWidth=1.7
    chart.lines[0].strokeDashArray=[5,3];chart.lines[1].strokeDashArray=[2,2];draw.add(chart)
    draw.add(String(220,61,'span coordinate y / 4 m',fontSize=10));draw.add(String(65,300,'weighted component / norm of its own difference',fontSize=9))
    for j,label in enumerate(['160 to 320','320 to 640','640 to 1280']):
        draw.add(Rect(60+j*155,38,12,6,fillColor=palette[j],strokeColor=None));draw.add(String(78+j*155,37,label,fontSize=10))
    draw.add(String(60,14,'Each curve normalized separately; this figure is not an error bound.',fontSize=9))
    p=folder/'avl-span-order-shape.svg';renderSVG.drawToFile(draw,str(p));ET.parse(p);files[p.name]=digest(p)
    draw=Drawing(640,370);draw.add(String(40,345,'Final span-difference energy: regional shares',fontName='Helvetica-Bold',fontSize=14))
    names=['tau0','tauP1','K0','Jline','Jpoint'];labels=[f'{case[-2:]} {name}' for case in ['FW-V01','FW-V02'] for name in names]
    values=[[100*c['stages'][-1]['regional_fractions'][reg] for c in r['channels']] for reg in ['root','interior','tip']];data['regional_percentages']=values
    bars=VerticalBarChart();bars.x=60;bars.y=105;bars.width=530;bars.height=215;bars.data=values;bars.categoryAxis.style='stacked';bars.categoryAxis.categoryNames=labels
    bars.categoryAxis.labels.angle=35;bars.categoryAxis.labels.boxAnchor='ne';bars.categoryAxis.labels.fontSize=8
    bars.valueAxis.valueMin=0;bars.valueAxis.valueMax=100;bars.valueAxis.valueStep=20
    for j,col in enumerate(palette):bars.bars[j].fillColor=col;bars.bars[j].strokeColor=colors.white
    draw.add(bars);draw.add(String(60,327,'Percent of squared weighted norm; FW-V01 followed by FW-V02',fontSize=9))
    for j,label in enumerate(['root: |y| <= 0.8','interior','tip: |y| >= 3.2']):
        draw.add(Rect(65+j*180,37,12,9,fillColor=palette[j],strokeColor=None));draw.add(String(83+j*180,37,label,fontSize=10))
    draw.add(String(60,15,'Energy concentration indicates sensitivity location, not a proven cause.',fontSize=9))
    p=folder/'avl-span-order-regions.svg';renderSVG.drawToFile(draw,str(p));ET.parse(p);files[p.name]=digest(p)
    (F/'avl-span-order-figures.json').write_text(json.dumps(dict(source_sha256=digest(__file__),parent_sha256=digest(result),
        reportlab_version=reportlab.Version,files=files,data=data),indent=2)+'\n')
    print(json.dumps(dict(svg_figures=2,xml_valid=True,files=files)))
if __name__=='__main__':main()
