const path = require('path');
const fs = require('fs');
const { pathToFileURL } = require('url');
const { chromium } = require('<home>/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
async function main() {
  const root = path.resolve(__dirname, '../..');
  if (fs.existsSync(path.join(root,'docs/papers/source-hashes.json'))) throw new Error('Report v1 is frozen; prepare a separately named revision');
  const browser = await chromium.launch({ executablePath: 'C:/Program Files/Google/Chrome/Application/chrome.exe', headless: true });
  try {
    const page = await browser.newPage({viewport:{width:960,height:1200}});
    const failures=[];page.on('pageerror',e=>failures.push(String(e)));
    // The HTML is standalone; no remote content is permitted during typesetting.
    await page.route(/^https?:\/\//, route => route.abort());
    await page.goto(pathToFileURL(path.join(root,'docs/papers/pg-ancf-technical-report-v1.html')).href);
    await page.evaluate(async()=>{await document.fonts.ready;await MathJax.startup.promise;});
    await page.emulateMedia({media:'print'});
    const math = await page.evaluate(()=>({equations:document.querySelectorAll('.equation').length,typeset:document.querySelectorAll('mjx-container').length,errors:Array.from(document.querySelectorAll('mjx-merror')).map(x=>x.textContent)}));
    if(failures.length||math.errors.length||math.equations!==32||math.typeset!==32)throw new Error(JSON.stringify({failures,math}));
    const out=path.join(root,'output/pdf/pg-ancf-technical-report-v1.pdf');fs.mkdirSync(path.dirname(out),{recursive:true});
    await page.pdf({path:out,format:'A4',preferCSSPageSize:true,printBackground:true,displayHeaderFooter:true,
      headerTemplate:'<div style="font-family:Arial;font-size:8px;width:100%;margin:0 18mm;color:#52616a">UVLAS TECHNICAL REPORT <span style="float:right">PG-ANCF / v1 / 2026-11-06</span></div>',
      footerTemplate:'<div style="font-family:Arial;font-size:8px;width:100%;text-align:center;color:#52616a"><span class="pageNumber"></span> / <span class="totalPages"></span></div>'});
    fs.writeFileSync(path.join(root,'build/paper-render/browser-check.json'),JSON.stringify({passed:true,math,remote_fetch_allowed:false},null,2));console.log(JSON.stringify({passed:true,output:out,math}));
  } finally {await browser.close();}
}
main().catch(e=>{console.error(e);process.exitCode=1;});
