"""Execute within browser-harness; capture public task-page text only."""
from datetime import datetime,timezone
import hashlib,json
from pathlib import Path
from urllib.parse import urlparse

root=Path.cwd();label=CAPTURE_LABEL;assert label and all(x.isalnum() or x in '-_' for x in label)
info=page_info();assert urlparse(info['url']).hostname in {'arxiv.org','www.cambridge.org','projecteuclid.org'}
body=js('document.body.innerText');assert isinstance(body,str) and len(body)>200
folder=root/'.firecrawl/fin_dh_trace_20261017';folder.mkdir(parents=True,exist_ok=True)
path=folder/(label+'.txt');path.write_text(body)
entry={'label':label,'observed_utc':datetime.now(timezone.utc).isoformat(),'url':info['url'],'title':info['title'],
    'method':'actual local browser-harness public DOM capture; recording off',
    'text_sha256':hashlib.sha256(body.encode()).hexdigest(),'text_bytes':len(body.encode()),
    'cache_path':str(path.relative_to(root)),
    'scope':'access and extracted-text provenance; not automatic full-paper reading or original PDF-byte hashing'}
with (root/'fin_intermediate_clock_bridge_research_20261017/BROWSER_READ_LOG.jsonl').open('a') as f:f.write(json.dumps(entry,ensure_ascii=False)+'\n')
print(json.dumps(entry,ensure_ascii=False))
