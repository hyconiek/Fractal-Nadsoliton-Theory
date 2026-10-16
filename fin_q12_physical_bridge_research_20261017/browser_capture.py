"""Run inside browser-harness: capture only the current public research tab.

Body text remains in ignored .firecrawl cache; tracked log stores hashes and
URLs only. This is not a claim that a whole paper/proof has been read.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse

root=Path.cwd();label=CAPTURE_LABEL
assert label and all(c.isalnum() or c in '-_' for c in label)
info=page_info();host=urlparse(info['url']).hostname
assert host in {'arxiv.org','journals.plos.org','www.nature.com','pmc.ncbi.nlm.nih.gov'},host
body=js('document.body.innerText')
assert isinstance(body,str) and len(body)>200
cache=root/'.firecrawl/fin_pc_physics_20261017';cache.mkdir(parents=True,exist_ok=True)
path=cache/(label+'.txt');path.write_text(body)
entry={'label':label,'observed_utc':datetime.now(timezone.utc).isoformat(),
       'url':info['url'],'title':info['title'],'method':'local browser-harness DOM body; recording disabled',
       'text_sha256':hashlib.sha256(body.encode()).hexdigest(),'text_bytes':len(body.encode()),
       'cache_path':str(path.relative_to(root)),
       'scope':'browser access/capture, not automatic full proof reading or physical evidence'}
log=root/'fin_q12_physical_bridge_research_20261017/BROWSER_READ_LOG.jsonl'
with log.open('a') as f:f.write(json.dumps(entry,ensure_ascii=False)+'\n')
print(json.dumps(entry,ensure_ascii=False))
